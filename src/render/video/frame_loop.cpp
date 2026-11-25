// The video frame loop (frame_loop.hpp): mdatk/render/video.py video_positions and render_video
// (MDATK's own code, ported), with the frames composed on a pool of threads.
#include "video/frame_loop.hpp"

#include <algorithm>
#include <condition_variable>
#include <deque>
#include <exception>
#include <map>
#include <mutex>
#include <string>
#include <thread>
#include <utility>

#include "mdbridge/core/log.hpp"
#include "mdbridge/core/paths.hpp"
#include "video/encoder.hpp"

namespace fs = std::filesystem;

namespace mdbridge::render::video {

FrameCompositor::~FrameCompositor() = default;

std::vector<std::size_t> videoPositions(std::size_t frames, std::int64_t stride) {
    const auto step = static_cast<std::size_t>(std::max<std::int64_t>(1, stride));
    std::vector<std::size_t> out;
    out.reserve(frames / step + 1);
    for (std::size_t pos = 0; pos < frames; pos += step) out.push_back(pos);
    return out;
}

unsigned frameLoopThreads(const FrameLoopOptions& options) {
    if (options.threads > 0) return options.threads;
    return std::clamp(std::thread::hardware_concurrency(), 1u, 8u);
}

namespace {

// Composes frames on worker threads; results are taken in any order by frame index.
class ComposePool {
public:
    ComposePool(const FrameCompositor& compositor, unsigned threads) : m_compositor(compositor) {
        m_workers.reserve(threads);
        for (unsigned t = 0; t < threads; ++t) m_workers.emplace_back([this] { work(); });
    }
    ~ComposePool() {
        {
            const std::lock_guard<std::mutex> lock(m_mutex);
            m_stop = true;
            m_jobs.clear();
        }
        m_jobReady.notify_all();
        for (std::thread& t : m_workers) t.join();
    }
    ComposePool(const ComposePool&) = delete;
    ComposePool& operator=(const ComposePool&) = delete;

    void submit(std::size_t k, DisplayList layer) {
        {
            const std::lock_guard<std::mutex> lock(m_mutex);
            m_jobs.emplace_back(k, std::move(layer));
        }
        m_jobReady.notify_one();
    }

    // Frame k once composed (rethrows the compositor's error for it).
    std::vector<std::uint8_t> take(std::size_t k) {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_doneReady.wait(lock, [&] { return m_done.count(k) != 0; });
        auto it = m_done.find(k);
        Done d = std::move(it->second);
        m_done.erase(it);
        lock.unlock();
        if (d.error) std::rethrow_exception(d.error);
        return std::move(d.rgba);
    }

private:
    struct Done {
        std::vector<std::uint8_t> rgba;
        std::exception_ptr error;
    };

    void work() {
        for (;;) {
            std::pair<std::size_t, DisplayList> job;
            {
                std::unique_lock<std::mutex> lock(m_mutex);
                m_jobReady.wait(lock, [&] { return m_stop || !m_jobs.empty(); });
                if (m_stop) return;
                job = std::move(m_jobs.front());
                m_jobs.pop_front();
            }
            Done d;
            try {
                d.rgba = m_compositor.compose(job.second);
            } catch (...) {
                d.error = std::current_exception();
            }
            {
                const std::lock_guard<std::mutex> lock(m_mutex);
                m_done.emplace(job.first, std::move(d));
            }
            m_doneReady.notify_all();
        }
    }

    const FrameCompositor& m_compositor;
    std::mutex m_mutex;
    std::condition_variable m_jobReady, m_doneReady;
    std::deque<std::pair<std::size_t, DisplayList>> m_jobs;
    std::map<std::size_t, Done> m_done;
    bool m_stop = false;
    std::vector<std::thread> m_workers;
};

}  // namespace

void runFrameLoop(panels::BuiltFigure& figure, const FontSet& fonts, std::span<const std::size_t> positions,
                  FrameCompositor& compositor, const FrameSink& sink, const core::Reporter& reporter,
                  const FrameLoopOptions& options) {
    // blitting: static content (axes, ticks, full traces, heatmaps) is drawn once; each frame
    // restores it and draws only the artists that move
    panels::markAnimated(figure.figure, figure.animatedArtists());
    figure.figure.draw(fonts);
    compositor.setBackground(buildDisplayList(figure.figure, fonts, DisplayLayer::Static));

    const std::size_t n = positions.size();
    const std::size_t every = std::max<std::size_t>(1, n / 200);
    const auto write = [&](std::size_t k, std::vector<std::uint8_t> rgba) {
        reporter.checkCancel();  // as MDATK: no frame is written once cancelled
        sink(k, rgba);
        if (k % every == 0 || k == n - 1) {
            reporter.progress(static_cast<double>(k + 1) / static_cast<double>(n),
                              "Video frame " + std::to_string(k + 1) + "/" + std::to_string(n));
        }
    };
    const unsigned threads = frameLoopThreads(options);
    if (threads <= 1 || n <= 1) {
        for (std::size_t k = 0; k < n; ++k) {
            reporter.checkCancel();
            figure.update(positions[k], fonts);
            write(k, compositor.compose(buildDisplayList(figure.figure, fonts, DisplayLayer::Animated)));
        }
        return;
    }
    // frames in flight (submitted, not yet written): enough to keep the pool busy while the calling
    // thread updates the panels and writes, few enough to bound the memory (a 1080p frame is 8 MB)
    const std::size_t inFlight = 2 * static_cast<std::size_t>(threads);
    ComposePool pool(compositor, threads);
    std::size_t next = 0;  // the next frame to write
    for (std::size_t k = 0; k < n; ++k) {
        reporter.checkCancel();
        figure.update(positions[k], fonts);
        pool.submit(k, buildDisplayList(figure.figure, fonts, DisplayLayer::Animated));
        while (k + 1 - next >= inFlight) {
            write(next, pool.take(next));
            ++next;
        }
    }
    for (; next < n; ++next) write(next, pool.take(next));
}

fs::path renderVideo(const core::Dataset& dataset, const core::Project& project, const fs::path& path,
                     const FontSet& fonts, FrameCompositor& compositor, const core::Reporter& reporter,
                     const VideoRenderOptions& options) {
    if (dataset.frameCount() == 0) {
        throw VideoError("no analysed frames to render");  // RN-02
    }
    const core::VideoParams& vp = project.output.video;
    panels::BuiltFigure figure = panels::buildFigure(dataset, project, fonts);
    // RN-06: the canvas is exactly figure_size(width, height), always even
    const int w = figure.figure.widthPx(), h = figure.figure.heightPx();
    const fs::path folder = fs::absolute(path).parent_path();
    std::error_code ec;
    fs::create_directories(folder, ec);
    if (ec && !fs::is_directory(folder)) {
        throw VideoError("video: cannot create the folder " + core::pathToUtf8(folder) + ": " + ec.message());
    }
    const std::vector<std::size_t> positions = videoPositions(dataset.frameCount(), vp.frame_stride);
    const std::string pathText = core::pathToUtf8(path);
    reporter.log("Rendering " + std::to_string(positions.size()) + " video frame(s) at " + std::to_string(w) + "x" +
                 std::to_string(h) + ", " + std::to_string(vp.fps) + " fps -> " + pathText);
    EncoderOptions eo;
    eo.width = w;
    eo.height = h;
    eo.fps = static_cast<double>(vp.fps);
    eo.codec = vp.codec;
    eo.crf = static_cast<int>(std::clamp<std::int64_t>(vp.crf, -1000, 1000));
    eo.threads = options.encoderThreads;
    VideoEncoder encoder(path, eo, [](const std::string& message) { core::logWarning(message); });
    runFrameLoop(
        figure, fonts, positions, compositor,
        [&](std::size_t, std::vector<std::uint8_t>& rgba) {
            if (rgba.size() != static_cast<std::size_t>(w) * static_cast<std::size_t>(h) * 4) {
                throw VideoError("video: a composed frame has " + std::to_string(rgba.size()) + " bytes, not " +
                                 std::to_string(w) + "x" + std::to_string(h) + "x4");
            }
            encoder.writeFrame(rgba.data());
        },
        reporter,
        options.loop);
    encoder.finish();
    reporter.log("Wrote " + pathText + " (" + std::to_string(encoder.framesWritten()) + " frames)");
    return path;
}

}  // namespace mdbridge::render::video
