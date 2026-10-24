// Progress reporting and cancellation as MDATK's Reporter does them (data/topology/golden_model.json).

#include <cmath>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mdbridge/core/progress.hpp"
#include "test_support.hpp"

using namespace mdbridge::core;
using namespace mdbridge::test;

TEST(Progress, EventsMatchMdatk) {
    Array events;
    Reporter rep([&](const std::string& m) { events.push_back(Array{"log", m}); },
                 [&](double f, const std::string& m) { events.push_back(Array{"progress", f, m}); });
    rep.log("hello");
    rep.rawLog("raw");
    rep.progress(0.5, "half");
    rep.progress(-1);
    rep.progress(2.0, "over");
    rep.progress(std::nan(""), "nan");
    const Reporter s = rep.sub(0.2, 0.6);
    s.progress(0.5, "sub half");
    const Reporter ss = s.sub(0.5, 1.0);
    ss.progress(0.0, "subsub 0");
    ss.progress(1.0, "subsub 1");
    rep.sub(0.9, 1.5).progress(1.0, "beyond");

    const Value golden = *readJson(kData / "topology" / "golden_model.json").find("reporter_events");
    ASSERT_EQ(events.size(), golden.size());
    for (std::size_t i = 0; i < events.size(); ++i) {
        const Array& want = *(*golden.get<Array>())[i].get<Array>();
        const Array& got = *events[i].get<Array>();
        ASSERT_EQ(got.size(), want.size());
        EXPECT_EQ(got[0], want[0]);
        if (want.size() == 3) {
            EXPECT_NEAR(*got[1].get<double>(), *want[1].get<double>(), 1e-12) << i;
            EXPECT_EQ(got[2], want[2]);
        } else {
            EXPECT_EQ(got[1], want[1]);
        }
    }

    EXPECT_NO_THROW(rep.checkCancel());
    s.cancel();  // any reporter of the family cancels all
    try {
        rep.checkCancel();
        ADD_FAILURE() << "not cancelled";
    } catch (const Cancelled& error) {
        EXPECT_STREQ(error.what(), "operation cancelled");
    }
    EXPECT_TRUE(ss.cancelled());
    EXPECT_FALSE(Reporter().cancelled());
}
