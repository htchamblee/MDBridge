#include "page_host.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exception>
#include <optional>
#include <set>
#include <string_view>
#include <utility>
#include <vector>

#include <QAbstractButton>
#include <QApplication>
#include <QLabel>
#include <QScrollArea>
#include <QTabWidget>
#include <QTest>
#include <QVBoxLayout>

#include "gui/forms/form_context.hpp"
#include "gui/pages/form_page.hpp"
#include "gui/pages/page_registry.hpp"
#include "mdbridge/core/json.hpp"
#include "mdbridge/core/paths.hpp"
#include "mdbridge/core/problems.hpp"

namespace fs = std::filesystem;

namespace mdbridge::guitest {

namespace {

QString qs(const std::string& text) { return QString::fromUtf8(text.data(), static_cast<qsizetype>(text.size())); }

std::string str(const core::Value& table, std::string_view key) {
    const core::Value* v = table.find(key);
    return v ? v->stringOr() : std::string();
}

// The page id an action key belongs to: the field's page, else its first part when that is a page
// id ("structure.reset_view"), else "".
std::string pageOfKey(const std::string& key) {
    const std::string_view page = core::problemPage(key);
    if (!page.empty()) {
        return std::string(page);
    }
    const std::string first = key.substr(0, key.find('.'));
    for (const core::ProblemPage& p : core::kProblemPages) {
        if (p.id == first) {
            return first;
        }
    }
    return {};
}

bool hostHasPage(const std::string& id) {
    const std::vector<QString> ids = app::pages::pageIds();
    return std::find(ids.begin(), ids.end(), qs(id)) != ids.end();
}

std::optional<double> numberOf(const core::Value& v) {
    if (const auto* d = v.get<double>()) return *d;
    if (const auto* i = v.get<std::int64_t>()) return static_cast<double>(*i);
    return std::nullopt;
}

}  // namespace

PageHost::PageHost() {
    setObjectName(QStringLiteral("pagehost"));
    m_context = std::make_unique<app::forms::FormContext>(m_document);
    auto* layout = new QVBoxLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);
    m_tabs = new QTabWidget(this);
    layout->addWidget(m_tabs);
    m_pages = app::pages::createPages(*m_context);
    for (QWidget* p : m_pages) {
        m_tabs->addTab(p, app::pages::pageTitleOf(p).replace(QLatin1Char('&'), QStringLiteral("&&")));  // no mnemonic
    }
}

PageHost::~PageHost() {
    // the pages' field widgets go before the context and the document
    delete m_tabs;
    m_tabs = nullptr;
    m_pages.clear();
}

QWidget* PageHost::page(const QString& id) const {
    for (QWidget* p : m_pages) {
        if (app::pages::pageIdOf(p) == id) {
            return p;
        }
    }
    return nullptr;
}

bool pageHostRuns(const GuiCase& c, std::string& why) {
    static const std::set<std::string> supported{"start", "lossless", "type", "table", "snapshot", "toggle"};
    for (const std::string& need : c.needs) {
        if (!supported.count(need)) {
            why = "needs what neither the harness nor the page host has yet: " + need;
            return false;
        }
    }
    if (c.start.empty()) {
        why = "starts without a project (the page host opens the case's start project)";
        return false;
    }
    for (const char* k : {"dialogs", "questions", "uncaught"}) {
        const core::Value* v = c.outcome.isTable() ? c.outcome.find(k) : nullptr;
        if (v && v->isArray() && !v->get<core::Array>()->empty()) {
            why = std::string("expects ") + k + " of the window";
            return false;
        }
    }
    bool onPage = false;
    for (const core::Value& action : c.actions) {
        const std::string what = str(action, "action");
        if (what == "snapshot") {
            // "slots" (the layout page's slot editors) is left out of the comparison until that page
            // exists (W2-A); the snapshot's values are compared
            for (const char* k : {"page", "modified"}) {
                if (action.find(k)) {
                    why = std::string("a snapshot of '") + k + "' needs the window";
                    return false;
                }
            }
            if (const core::Value* keys = action.find("keys")) {
                for (const core::Value& k : *keys->get<core::Array>()) {
                    if (!hostHasPage(pageOfKey(k.stringOr()))) {
                        why = "snapshots '" + k.stringOr() + "' of a page the host does not have yet";
                        return false;
                    }
                    onPage = true;
                }
            }
            continue;
        }
        if (what == "save" || what == "screenshot" || what == "open") {
            if (action.find("answers") || (what == "screenshot" && action.find("widget"))) {
                why = "a " + what + " with window parts";
                return false;
            }
            continue;
        }
        if (what != "set" && what != "type" && what != "table" && what != "click" && what != "toggle") {
            why = "action '" + what + "' needs the window";
            return false;
        }
        const std::string key = str(action, "key");
        if (!hostHasPage(pageOfKey(key))) {
            why = "key '" + key + "' is on a page the host does not have yet";
            return false;
        }
        onPage = true;
    }
    if (!onPage) {
        why = "no action on a settings page (a window case: open, wait, save)";
        return false;
    }
    return true;
}

std::string compareValues(const core::Value& actual, const core::Value& expected, const std::string& where) {
    const std::optional<double> a = numberOf(actual);
    const std::optional<double> e = numberOf(expected);
    if (a && e) {
        const bool same = *a == *e || (std::isnan(*a) && std::isnan(*e)) ||
                          std::abs(*a - *e) <= 1e-9 * std::max(std::abs(*a), std::abs(*e));
        return same ? std::string() : where + ": expected " + core::toJson(expected, -1) + ", got " + core::toJson(actual, -1);
    }
    if (const core::Table* et = expected.get<core::Table>()) {
        const core::Table* at = actual.get<core::Table>();
        if (!at) {
            return where + ": expected a table, got " + core::toJson(actual, -1);
        }
        for (const auto& [key, value] : *et) {
            const core::Value* got = actual.find(key);
            if (!got) {
                return where + "." + key + ": missing";
            }
            if (std::string d = compareValues(*got, value, where + "." + key); !d.empty()) {
                return d;
            }
        }
        for (const auto& [key, value] : *at) {
            if (!expected.find(key)) {
                return where + "." + key + ": not expected";
            }
        }
        return {};
    }
    if (const core::Array* ea = expected.get<core::Array>()) {
        const core::Array* aa = actual.get<core::Array>();
        if (!aa || aa->size() != ea->size()) {
            return where + ": expected " + core::toJson(expected, -1) + ", got " + core::toJson(actual, -1);
        }
        for (std::size_t k = 0; k < ea->size(); ++k) {
            if (std::string d = compareValues((*aa)[k], (*ea)[k], where + "." + std::to_string(k)); !d.empty()) {
                return d;
            }
        }
        return {};
    }
    return actual == expected ? std::string()
                              : where + ": expected " + core::toJson(expected, -1) + ", got " + core::toJson(actual, -1);
}

std::string runOnPageHost(const GuiCase& c, const fs::path& dataDir, const fs::path& caseDir, const fs::path& screenshots,
                          std::vector<fs::path>* shots) {
    std::string failures;
    std::string shotName = c.name;
    std::replace(shotName.begin(), shotName.end(), '/', '_');
    std::vector<core::Value> snapshots;
    std::set<std::string> pagesUsed;
    {
        PageHost host;
        host.resize(560, 900);
        host.show();
        if (!QTest::qWaitForWindowExposed(&host)) {
            return "the page host was not shown";
        }
        std::size_t k = 0;
        try {
            prepareCase(c, dataDir, caseDir);
            host.document().open(caseDir / core::pathFromUtf8(c.start));
            QApplication::processEvents();
            for (k = 0; k < c.actions.size(); ++k) {
                const core::Value& action = c.actions[k];
                const std::string what = str(action, "action");
                const std::string key = str(action, "key");
                QWidget* page = key.empty() ? nullptr : host.page(qs(pageOfKey(key)));
                if (page) {
                    host.tabs()->setCurrentWidget(page);
                    pagesUsed.insert(pageOfKey(key));
                }
                if (what == "set" || what == "toggle") {
                    const core::Value* value = action.find("value");
                    if (!value) throw CaseError("'value' is missing");
                    // a field that is not on screen is scrolled to first, as a user would
                    app::pages::revealFieldOn(page, qs(key));
                    setKey(page, qs(key), *value);
                } else if (what == "type") {
                    app::pages::revealFieldOn(page, qs(key));
                    typeKey(page, qs(key), qs(str(action, "text")));
                } else if (what == "table") {
                    const auto integer = [&](std::string_view name, int fallback) {
                        const core::Value* v = action.find(name);
                        return v && v->get<std::int64_t>() ? static_cast<int>(*v->get<std::int64_t>()) : fallback;
                    };
                    tableAction(page, qs(key), str(action, "do"), integer("row", -1), integer("by", 0));
                } else if (what == "click") {
                    auto* button = qobject_cast<QAbstractButton*>(findField(page, qs(key)));
                    if (!button) throw CaseError(key + " is no button");
                    if (!button->isEnabled()) throw CaseError(key + " is disabled");
                    QTest::mouseClick(button, Qt::LeftButton);
                } else if (what == "open") {
                    host.document().open(caseDir / core::pathFromUtf8(str(action, "path")));
                } else if (what == "save") {
                    host.context().commitPendingEdits();
                    const std::string path = str(action, "path");
                    if (path.empty()) {
                        host.document().save();
                    } else {
                        host.document().saveAs(caseDir / core::pathFromUtf8(path));
                    }
                } else if (what == "snapshot") {
                    host.context().commitPendingEdits();
                    core::Table values;
                    if (const core::Value* keys = action.find("keys")) {
                        for (const core::Value& kv : *keys->get<core::Array>()) {
                            values.emplace_back(kv.stringOr(), host.document().get(kv.stringOr()));
                        }
                    }
                    core::Table record{{"action", core::Value(static_cast<std::int64_t>(k))}};
                    if (!values.empty()) record.emplace_back("values", core::Value(std::move(values)));
                    if (action.find("density_status")) {
                        auto* label = host.findChild<QLabel*>(QStringLiteral("density.status"));
                        record.emplace_back("density_status",
                                            core::Value(label ? label->text().toStdString() : std::string("<no label>")));
                    }
                    snapshots.emplace_back(std::move(record));
                } else if (what == "screenshot") {
                    const fs::path shot = saveScreenshot(&host, screenshots, shotName + "-" + str(action, "name"));
                    if (shots) shots->push_back(shot);
                } else {
                    throw CaseError("the page host has no action '" + what + "'");
                }
                QApplication::processEvents();
            }
            // each page the case used, whole (the scroll area's content)
            for (const std::string& id : pagesUsed) {
                auto* page = qobject_cast<app::pages::FormPage*>(host.page(qs(id)));
                if (!page) continue;
                host.tabs()->setCurrentWidget(page);
                QApplication::processEvents();
                const fs::path shot = saveScreenshot(page->scrollArea()->widget(), screenshots, shotName + "-" + id);
                if (shots) shots->push_back(shot);
            }
        } catch (const std::exception& e) {
            failures += "action " + std::to_string(k + 1) + ": " + e.what() + "\n";
        }
    }
    if (!failures.empty()) {
        return failures;
    }
    if (const core::Value* expected = c.outcome.find("snapshots")) {
        const core::Array& want = *expected->get<core::Array>();
        if (want.size() != snapshots.size()) {
            failures += "expected " + std::to_string(want.size()) + " snapshots, took " + std::to_string(snapshots.size()) + "\n";
        }
        for (std::size_t k = 0; k < std::min(want.size(), snapshots.size()); ++k) {
            core::Value expectedSnapshot = want[k];
            if (core::Table* t = expectedSnapshot.get<core::Table>()) {
                std::erase_if(*t, [](const auto& entry) { return entry.first == "slots"; });  // needs the layout page
            }
            if (std::string d = compareValues(snapshots[k], expectedSnapshot, "snapshot " + std::to_string(k + 1));
                !d.empty()) {
                failures += d + "\n";
            }
        }
    }
    failures += checkExpected(c, dataDir, caseDir);
    return failures;
}

}  // namespace mdbridge::guitest
