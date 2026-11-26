#include "gui/display_path.hpp"

#include <algorithm>
#include <cstddef>
#include <mutex>
#include <string>
#include <vector>

#include <QDir>

#include "mdbridge/core/paths.hpp"

namespace mdbridge::app::gui {

namespace {

// Characters a shown path ends at: white space, quotes, characters Windows forbids in names, and
// ':' (after the drive's); control characters.
bool endsPath(QChar c) {
    return c.isSpace() || c.unicode() < 0x20 || c == QLatin1Char('"') || c == QLatin1Char('\'') ||
           c == QLatin1Char('<') || c == QLatin1Char('>') || c == QLatin1Char('|') || c == QLatin1Char('*') ||
           c == QLatin1Char('?') || c == QLatin1Char(':');
}

bool partOfWord(QChar c) { return c.isLetterOrNumber() || c == QLatin1Char('_'); }

// Text a path may follow in a message: the start, white space, a quote or an opening bracket, '='.
bool opensPath(const QString& text, qsizetype at) {
    if (at == 0) return true;
    const QChar c = text[at - 1];
    return c.isSpace() || c == QLatin1Char('"') || c == QLatin1Char('\'') || c == QLatin1Char('(') ||
           c == QLatin1Char('[') || c == QLatin1Char('=') || c == QLatin1Char('<');
}

struct Roots {
    std::mutex mutex;
    QStringList list;
};

Roots& roots() {
    static Roots r;
    return r;
}

}  // namespace

QString displayPathText(const QString& path, PathStyle style) {
    if (style == PathStyle::Posix) return path;
    QString out = path;
    out.replace(QLatin1Char('/'), QLatin1Char('\\'));
    return out;
}

QString displayPathsInText(const QString& text, PathStyle style, const QStringList& rootFolders) {
    if (style == PathStyle::Posix || text.isEmpty()) return text;
    QString generic = text;  // the same length: positions match
    generic.replace(QLatin1Char('\\'), QLatin1Char('/'));
    const qsizetype n = generic.size();
    std::vector<bool> convert(static_cast<std::size_t>(n), false);
    const auto mark = [&](qsizetype from, qsizetype to) {
        for (qsizetype i = from; i < to; ++i) convert[static_cast<std::size_t>(i)] = true;
    };
    // a path's rest: up to the next character that ends a path
    const auto rest = [&](qsizetype from) {
        qsizetype i = from;
        while (i < n && !endsPath(generic[i])) ++i;
        mark(from, i);
    };
    // the roots first (longest first: they may hold spaces)
    std::vector<QString> folders;
    for (const QString& root : rootFolders) {
        QString r = root;
        r.replace(QLatin1Char('\\'), QLatin1Char('/'));
        while (r.size() > 3 && r.endsWith(QLatin1Char('/'))) r.chop(1);
        if (r.size() >= 3) folders.push_back(r);
    }
    std::sort(folders.begin(), folders.end(), [](const QString& a, const QString& b) { return a.size() > b.size(); });
    for (const QString& r : folders) {
        for (qsizetype at = generic.indexOf(r, 0, Qt::CaseInsensitive); at >= 0;
             at = generic.indexOf(r, at + 1, Qt::CaseInsensitive)) {
            const qsizetype end = at + r.size();
            if (at > 0 && (partOfWord(generic[at - 1]) || generic[at - 1] == QLatin1Char('/'))) continue;
            if (end < n && partOfWord(generic[end])) continue;  // a longer name: "C:/data" in "C:/database"
            mark(at, end);
            rest(end);
        }
    }
    for (qsizetype i = 0; i + 2 < n; ++i) {
        const QChar c = generic[i];
        if (c.isLetter() && c.unicode() < 0x80 && generic[i + 1] == QLatin1Char(':') &&
            generic[i + 2] == QLatin1Char('/') && (i == 0 || !partOfWord(generic[i - 1]))) {
            rest(i + 2);  // "C:/...", also after a device prefix
            continue;
        }
        if (c == QLatin1Char('/') && generic[i + 1] == QLatin1Char('/') && opensPath(generic, i)) {
            if (i + 3 < n && (generic[i + 2] == QLatin1Char('?') || generic[i + 2] == QLatin1Char('.')) &&
                generic[i + 3] == QLatin1Char('/')) {
                mark(i, i + 4);  // a device path: "\\?\", "\\.\"
                rest(i + 4);
            } else if (!endsPath(generic[i + 2]) && generic[i + 2] != QLatin1Char('/')) {
                rest(i);  // "//server/share/..."
            }
        }
    }
    QString out = text;
    for (qsizetype i = 0; i < n; ++i) {
        if (convert[static_cast<std::size_t>(i)] && out[i] == QLatin1Char('/')) out[i] = QLatin1Char('\\');
    }
    return out;
}

QString displayPath(const std::filesystem::path& path) {
    const std::string text = core::pathToUtf8(path);
    return displayPathText(QString::fromUtf8(text.data(), static_cast<qsizetype>(text.size())), kNativePathStyle);
}

QString displayPath(const QString& path) { return displayPathText(path, kNativePathStyle); }

QString displayPathsIn(const QString& text) {
    if constexpr (kNativePathStyle == PathStyle::Posix) {
        return text;
    } else {
        QStringList all = displayRoots();
        all << QDir::homePath() << QDir::tempPath();
        return displayPathsInText(text, kNativePathStyle, all);
    }
}

void setDisplayRoots(const QStringList& list) {
    Roots& r = roots();
    const std::lock_guard<std::mutex> lock(r.mutex);
    r.list = list;
}

QStringList displayRoots() {
    Roots& r = roots();
    const std::lock_guard<std::mutex> lock(r.mutex);
    return r.list;
}

}  // namespace mdbridge::app::gui
