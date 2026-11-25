// The plot panels (src/render/panels plots.cpp: timeseries, dual_axis, heatmap, contact_map) against
// MDATK's figures (data/render/figures, tests/parity/make_figure_cases.py).

#include <gtest/gtest.h>

#include "figure_records.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/render/panels/build.hpp"

using namespace mdbridge;
using namespace mdbridge::render;

TEST(PanelRegistry, UnknownTypeAsMdatk) {
    core::Dataset ds({0}, {0.0});
    core::Project p;
    panels::Context ctx(ds, p);
    Figure fig(100, 100, 72.0);
    core::PanelSpec spec;
    spec.type = "nope";
    try {
        (void)panels::createPanel(fig, {0, 1, 0, 1}, spec, ctx);
        FAIL();
    } catch (const core::ValueError& e) {
        EXPECT_STREQ(e.what(), "unknown panel type 'nope' (choose from structure, timeseries, dual_axis, heatmap, "
                               "contact_map, sequence, secondary)");
    }
}
