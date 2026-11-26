#include "gui/pages/layout_page.hpp"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <QComboBox>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QSpinBox>
#include <QVBoxLayout>

#include "gui/forms/field_widgets.hpp"
#include "gui/forms/form_context.hpp"
#include "mdbridge/analysis/analysis.hpp"
#include "mdbridge/analysis/registry.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/json.hpp"
#pragma push_macro("slots")
#undef slots
#include "mdbridge/core/panels.hpp"
#pragma pop_macro("slots")
#include "mdbridge/core/problems.hpp"
#include "mdbridge/core/project_document.hpp"

namespace mdbridge::app::pages {

using forms::toQString;
using forms::toStdString;

namespace {

// MDATK's SLOT_DEFAULT_OPTIONS (gui/app.py:66): the options a new panel of these types gets.
core::Table slotDefaultOptions(const std::string& type) {
    if (type == "sequence") return core::Table{{"residues", ""}};
    if (type == "secondary") return core::Table{{"source", "dynamic"}, {"residues", ""}, {"dotbracket", ""}};
    return {};
}

std::string strip(const std::string& s) {
    const std::size_t a = s.find_first_not_of(" \t\r\n\f\v");
    if (a == std::string::npos) return {};
    const std::size_t b = s.find_last_not_of(" \t\r\n\f\v");
    return s.substr(a, b - a + 1);
}

std::vector<std::string> splitKeys(const std::string& text) {
    std::vector<std::string> out;
    std::size_t start = 0;
    for (;;) {
        const std::size_t end = text.find(',', start);
        std::string item = strip(text.substr(start, end == std::string::npos ? std::string::npos : end - start));
        if (!item.empty()) out.push_back(std::move(item));
        if (end == std::string::npos) return out;
        start = end + 1;
    }
}

std::string joined(const std::vector<std::string>& items, std::size_t from) {
    std::string out;
    for (std::size_t k = from; k < items.size(); ++k) {
        if (k > from) out += ", ";
        out += items[k];
    }
    return out;
}

QString optionsText(const core::Table& options) {
    return options.empty() ? QString() : toQString(core::toJson(core::Value(options), -1));
}

bool hasSeries(const std::string& type) {
    const core::PanelType* pt = core::findPanelType(type);
    return pt == nullptr || pt->maxSeries > 0;
}

bool isDefaultTitle(const std::string& title) {
    return std::find(core::kDefaultPanelTitles.begin(), core::kDefaultPanelTitles.end(), title) !=
           core::kDefaultPanelTitles.end();
}

std::string comboId(const QComboBox* box) {
    return box->currentIndex() < 0 ? std::string() : toStdString(box->currentData().toString());
}

void selectId(QComboBox* box, const QString& id) {
    const int at = box->findData(id);
    box->setCurrentIndex(at);
}

}  // namespace

// ---- automatic titles after a rename ----------------------------------------------------------

std::vector<std::string> retitleAutomaticPanels(core::Project& project, const core::Project& before,
                                                const std::vector<std::string>& changed) {
    std::vector<std::string> out;
    for (std::size_t k = 0; k < project.layout.panels.size(); ++k) {
        core::PanelSpec& spec = project.layout.panels[k];
        if (spec.series.empty() || spec.title.empty() || k >= before.layout.panels.size()) continue;
        const core::PanelSpec& old = before.layout.panels[k];
        const std::string oldLabel = old.series.empty() ? std::string() : analysis::keyLabel(old.series[0], before);
        // the editors' automatic title: the data's label, or a default title
        if (spec.title != oldLabel && !isDefaultTitle(spec.title)) continue;
        const std::string prefix = "layout.panels." + std::to_string(k) + ".series";
        const bool moved = std::any_of(changed.begin(), changed.end(),
                                       [&](const std::string& c) { return c.compare(0, prefix.size(), prefix) == 0; });
        // a default title ("RMSD") is not its data's label: it stays while the panel's keys do
        if (moved || spec.title == oldLabel) {
            const std::string title = analysis::keyLabel(spec.series[0], project);
            if (title != spec.title) {
                spec.title = title;
                out.push_back("layout.panels." + std::to_string(k) + ".title");
            }
        }
    }
    return out;
}

// ---- a slot's text fields ---------------------------------------------------------------------

// The More keys, Title and Options fields of a slot editor: a line edit stored when editing
// finishes (and by FormContext::commitPendingEdits), flagged with the form error of its text.
class SlotTextField : public forms::FieldWidget {
public:
    SlotTextField(forms::FormContext& context, int k, QString field, const QString& placeholder, const QString& hint,
                  LayoutPage* page)
        : forms::FieldWidget(context, "layout.panels." + std::to_string(k) + "." + toStdString(field), makeSpec(hint),
                             page),
          m_page(page), m_field(std::move(field)) {
        auto* layout = new QVBoxLayout(this);
        layout->setContentsMargins(0, 0, 0, 0);
        m_edit = new QLineEdit(this);
        m_edit->setPlaceholderText(placeholder);
        layout->addWidget(m_edit);
        setEditor(m_edit);
        connect(m_edit, &QLineEdit::textEdited, this, [this] { m_dirty = true; });
        connect(m_edit, &QLineEdit::editingFinished, this, [this] { commit(); });
    }

    QLineEdit* edit() const { return m_edit; }
    void refresh() override {}  // the page shows the values (LayoutPage::refreshSlots)
    void commit() override {
        if (!m_dirty) return;
        m_dirty = false;
        if (m_slot) m_page->commitText(*m_slot, m_field);
    }
    QString displayedText() const override { return m_edit->text(); }
    bool dirty() const { return m_dirty; }
    void showText(const QString& text) {
        const RefreshScope scope(*this);
        if (m_edit->text() != text) m_edit->setText(text);
        m_dirty = false;
    }
    void flagChanged() { updateFlag(); }
    LayoutPage::Slot* m_slot = nullptr;

protected:
    void onDocumentChanged(const core::DocumentChange&) override {}

private:
    static forms::FieldSpec makeSpec(const QString& hint) {
        forms::FieldSpec spec;
        spec.hint = hint;
        return spec;
    }
    LayoutPage* m_page;
    QString m_field;
    QLineEdit* m_edit = nullptr;
    bool m_dirty = false;
};

// A drop-down of the page as a field widget of `path` (problems flag it, the problems list goes to
// it); the drop-down itself is named `name` (the GUI cases' key: "layout.panels.1.data" edits the
// first of "layout.panels.1.series").
class SlotComboField : public forms::FieldWidget {
public:
    SlotComboField(forms::FormContext& context, const std::string& path, const QString& name, QWidget* parent)
        : forms::FieldWidget(context, path, forms::FieldSpec{}, parent) {
        auto* layout = new QVBoxLayout(this);
        layout->setContentsMargins(0, 0, 0, 0);
        m_combo = new QComboBox(this);
        m_combo->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
        m_combo->setMinimumContentsLength(8);
        forms::guardWheel(m_combo);
        layout->addWidget(m_combo);
        setEditor(m_combo);
        m_combo->setObjectName(name);
    }
    QComboBox* combo() const { return m_combo; }
    void refresh() override {}  // the page shows the values
    QString displayedText() const override { return m_combo->currentText(); }

protected:
    void onDocumentChanged(const core::DocumentChange& change) override {
        if (change.affects(path())) updateFlag();
    }

private:
    QComboBox* m_combo = nullptr;
};

struct LayoutPage::Slot {
    int k = 0;
    Section* section = nullptr;
    QComboBox* type = nullptr;
    QComboBox* data = nullptr;
    QComboBox* ga = nullptr;
    QComboBox* gb = nullptr;
    SlotTextField* more = nullptr;
    SlotTextField* title = nullptr;
    SlotTextField* options = nullptr;
    QString autoTitle;          // MDATK's row["auto"]: the title that follows the data
    QString optionsError;       // the located form error of the options text ("" none)
    std::string path(const char* field) const { return "layout.panels." + std::to_string(k) + "." + field; }
};

// ---- the page ---------------------------------------------------------------------------------

LayoutPage::LayoutPage(forms::FormContext& context, QWidget* parent)
    : FormPage(context, QStringLiteral("layout"), tr("Layout & panels"), parent) {
    Section* sec = addSection(QStringLiteral("figure"), tr("Figure layout"));
    FormGrid* f = sec->form();
    auto* presetField = new SlotComboField(context, "layout.preset", QStringLiteral("layout.preset"), this);
    m_preset = presetField->combo();
    m_preset->setAccessibleName(tr("Grid preset"));
    for (const core::LayoutPreset& p : core::layoutPresets()) {
        m_preset->addItem(toQString(std::string(p.label)), toQString(std::string(p.name)));
    }
    connect(m_preset, &QComboBox::activated, this, &LayoutPage::onPreset);
    f->add(tr("Grid preset"), presetField);
    f->addHint(tr("Panels fill the slots in order (slot 1 is the large one). Drag panels in the figure to swap "
                  "them. Choose the data and, for analyses that compare residue groups, which groups to use: each "
                  "group choice is computed and cached on its own. “More keys” adds series to a "
                  "time-series panel (comma separated data keys); options are JSON, e.g. {\"mode\": \"frequency\"} "
                  "for a contact map or {\"progressive\": true} for a time series."));
    auto* box = new QWidget;
    box->setObjectName(QStringLiteral("layout.slots"));
    m_slotBox = new QVBoxLayout(box);
    m_slotBox->setContentsMargins(0, 0, 0, 0);
    m_slotBox->setSpacing(6);
    f->addWide(box);
    m_stillFrame = new QSpinBox;
    m_stillFrame->setObjectName(QStringLiteral("layout.still_frame"));
    m_stillFrame->setRange(0, 1000000000);
    m_stillFrame->setAccessibleName(tr("Real-data still of frame"));
    forms::guardWheel(m_stillFrame);
    m_renderStill = makeButton(tr("Render still"), QStringLiteral("layout.render_still"),
                               tr("Analyse (or load cached results) and render one frame of the video"));
    connect(m_renderStill, &QPushButton::clicked, this, [this] { emit renderStillRequested(m_stillFrame->value()); });
    f->addRow(tr("Real-data still of frame"), {m_stillFrame, m_renderStill});

    connect(&context, &forms::FormContext::documentChanged, this, &LayoutPage::documentChanged);
    rebuildSlots();
}

LayoutPage::~LayoutPage() {
    for (const auto& slot : m_slots) {
        context().setFormError(slot->path("options"), nullptr);
    }
}

core::PanelSpec LayoutPage::defaultSpec() const { return core::PanelSpec{"timeseries", {}, "", {}, {}}; }

void LayoutPage::rebuildSlots() {
    const core::Project& project = context().document().project();
    const std::string preset = project.layout.preset;
    ++m_loading;
    selectId(m_preset, toQString(preset));
    if (m_preset->currentIndex() < 0) {
        m_preset->addItem(tr("%1 (unknown)").arg(toQString(preset)), toQString(preset));
        selectId(m_preset, toQString(preset));
    }
    for (const auto& slot : m_slots) {
        context().setFormError(slot->path("options"), nullptr);
        // at once: its field widgets belong to the FormContext, which may go before a later event
        // loop would delete them (rebuilds come from the preset, an undo or a new project, never
        // from a slot's own widgets)
        delete slot->section;
    }
    m_slots.clear();
    const core::LayoutPreset* layout = core::findLayoutPreset(preset);
    const int n = layout ? static_cast<int>(core::presetSlots(*layout).size()) : 0;
    const bool large = preset.rfind("large", 0) == 0;
    for (int k = 0; k < n; ++k) {
        auto slot = std::make_unique<Slot>();
        slot->k = k;
        const QString heading = tr("Slot %1").arg(k + 1) + (k == 0 && large ? tr(" (large)") : QString());
        slot->section = new Section(QStringLiteral("layout"), QStringLiteral("slot%1").arg(k + 1), heading);
        FormGrid* f = slot->section->form();
        const auto combo = [&](const std::string& path, const char* field, const QString& name) {
            auto* c = new SlotComboField(context(), path, toQString(slot->path(field)), this);
            c->combo()->setAccessibleName(name);
            return c;
        };
        SlotComboField* typeField = combo(slot->path("type"), "type", tr("Panel type"));
        slot->type = typeField->combo();
        for (const core::PanelType& t : core::panelTypes()) {
            slot->type->addItem(toQString(std::string(t.label)), toQString(std::string(t.name)));
        }
        f->add(tr("Panel type"), typeField);
        SlotComboField* dataField = combo(slot->path("series"), "data", tr("Data"));
        slot->data = dataField->combo();
        slot->data->setMinimumContentsLength(20);
        f->add(tr("Data"), dataField);
        slot->ga = combo(slot->path("group_a"), "group_a", tr("Group A"))->combo();
        slot->gb = combo(slot->path("group_b"), "group_b", tr("Group B"))->combo();
        f->addRow(tr("Groups"), {slot->ga->parentWidget(), new QLabel(tr("vs")), slot->gb->parentWidget()},
                  tr("Residue groups this panel's analysis compares; each choice is cached separately"));
        slot->more = new SlotTextField(context(), k, QStringLiteral("more"), tr("additional data keys (comma separated)"),
                                       QString(), this);
        f->add(tr("More keys"), slot->more);
        slot->title = new SlotTextField(context(), k, QStringLiteral("title"), tr("automatic"),
                                        tr("Follows the selected data until you type your own title"), this);
        f->add(tr("Title"), slot->title);
        slot->options = new SlotTextField(context(), k, QStringLiteral("options"), QStringLiteral("{}"), QString(), this);
        f->add(tr("Options (JSON)"), slot->options);
        Slot* raw = slot.get();
        for (SlotTextField* t : {raw->more, raw->title, raw->options}) t->m_slot = raw;
        connect(raw->type, &QComboBox::activated, this, [this, raw] { onType(*raw); });
        connect(raw->data, &QComboBox::activated, this, [this, raw] { onData(*raw, false); });
        connect(raw->ga, &QComboBox::activated, this, [this, raw] { onData(*raw, false); });
        connect(raw->gb, &QComboBox::activated, this, [this, raw] { onData(*raw, false); });
        m_slotBox->addWidget(slot->section);
        m_slots.push_back(std::move(slot));
    }
    fixTabOrder();
    --m_loading;
    refreshSlots(true);
}

// D10: widgets made after the page go to the end of the window's focus chain; Tab visits the slot
// editors between the preset and the still row, in their order, as the page shows them.
void LayoutPage::fixTabOrder() {
    QWidget* previous = m_preset;
    for (const auto& slot : m_slots) {
        for (QWidget* w : slot->section->findChildren<QWidget*>()) {
            // not a combo's pop-up (a window of its own), not a widget that hands the focus on
            if (w->window() != m_preset->window() || w->focusProxy() || !(w->focusPolicy() & Qt::TabFocus)) continue;
            QWidget::setTabOrder(previous, w);
            previous = w;
        }
    }
    QWidget::setTabOrder(previous, m_stillFrame);
    QWidget::setTabOrder(m_stillFrame, m_renderStill);
}

void LayoutPage::fillGroups(Slot& slot) {
    const core::Project& project = context().document().project();
    for (QComboBox* c : {slot.ga, slot.gb}) {
        const QString keep = c->currentIndex() >= 0 ? c->currentData().toString() : QString();
        c->clear();
        c->addItem(tr("(analysis default)"), QString());
        for (const core::ResidueGroup& g : project.groups) {
            if (c->findData(toQString(g.name)) < 0) c->addItem(toQString(g.name), toQString(g.name));
        }
        selectId(c, keep);
        if (c->currentIndex() < 0) c->setCurrentIndex(0);
    }
}

void LayoutPage::fillData(Slot& slot, const std::string& type, const std::string& select) {
    const core::Project& project = context().document().project();
    std::vector<std::string> keys;
    if (const core::PanelType* pt = core::findPanelType(type); pt != nullptr && pt->maxSeries > 0) {
        try {
            keys = analysis::expectedKeys(project).of(pt->dataKind);
        } catch (const std::exception&) {
            keys.clear();
        }
    }
    if (!select.empty() && std::find(keys.begin(), keys.end(), select) == keys.end()) {
        keys.insert(keys.begin(), select);
    }
    slot.data->clear();
    for (const std::string& key : keys) {
        slot.data->addItem(toQString(analysis::keyLabel(key, project) + "   [" + key + "]"), toQString(key));
    }
    if (!keys.empty()) selectId(slot.data, toQString(select.empty() ? keys.front() : select));
    slot.data->setEnabled(!keys.empty());
}

void LayoutPage::updateGroupWidgets(Slot& slot) {
    const std::string key = comboId(slot.data);
    int groupSlots = 0;
    if (key.find('/') != std::string::npos) {
        try {
            groupSlots = analysis::analysisFor(analysis::splitKey(key).first).groupSlots();
        } catch (const std::exception&) {
            groupSlots = 0;
        }
    }
    slot.ga->setEnabled(groupSlots >= 1);
    slot.gb->setEnabled(groupSlots >= 2);
}

std::string LayoutPage::primaryKey(const Slot& slot) const {
    const std::string key = comboId(slot.data);
    if (key.find('/') == std::string::npos) return key;
    const auto [inst, out] = analysis::splitKey(key);
    const std::string base = analysis::parseInstance(inst).base;
    const analysis::AnalysisDescriptor* cls = analysis::findAnalysis(base);
    if (cls == nullptr) return key;
    const auto n = static_cast<std::size_t>(std::max(0, cls->groupSlots()));
    std::vector<std::string> chosen{comboId(slot.ga), comboId(slot.gb)};
    chosen.resize(std::min(n, chosen.size()));
    if (std::all_of(chosen.begin(), chosen.end(), [](const std::string& c) { return c.empty(); })) return key;
    std::vector<std::string> defaults = cls->groups(context().document().project());
    if (defaults.empty()) defaults.assign(n, std::string());
    std::vector<std::string> groups;
    for (std::size_t k = 0; k < std::min(chosen.size(), defaults.size()); ++k) {
        groups.push_back(chosen[k].empty() ? defaults[k] : chosen[k]);
    }
    if (n == 1 && (groups.empty() || groups[0].empty())) return key;
    return analysis::instanceId(base, groups) + "/" + out;
}

std::string LayoutPage::autoTitle(const Slot& slot) const {
    if (comboId(slot.type) == "structure") return "Structure";
    const std::string key = primaryKey(slot);
    return key.empty() ? std::string() : analysis::keyLabel(key, context().document().project());
}

core::PanelSpec LayoutPage::editorSpec(const Slot& slot) const {
    core::PanelSpec spec;
    spec.type = comboId(slot.type);
    if (spec.type.empty()) spec.type = "timeseries";
    if (hasSeries(spec.type)) {
        const std::string primary = primaryKey(slot);
        if (!primary.empty()) spec.series.push_back(primary);
        for (std::string& key : splitKeys(toStdString(slot.more->edit()->text()))) spec.series.push_back(std::move(key));
    }
    spec.title = strip(toStdString(slot.title->edit()->text()));
    if (slot.optionsError.isEmpty()) {
        const std::string text = strip(toStdString(slot.options->edit()->text()));
        if (!text.empty()) {
            try {
                core::Value v = core::parseJson(text, "");
                if (const auto* t = v.get<core::Table>()) spec.options = *t;
            } catch (const std::exception&) {
            }
        }
    }
    return spec;
}

bool LayoutPage::editorMatches(const Slot& slot, const core::PanelSpec& spec) const {
    if (comboId(slot.type) != spec.type) return false;
    if (!hasSeries(spec.type)) return true;
    return primaryKey(slot) == (spec.series.empty() ? std::string() : spec.series.front());
}

void LayoutPage::loadSlot(Slot& slot, const core::PanelSpec& spec, bool keepText) {
    ++m_loading;
    const std::string primary = spec.series.empty() ? std::string() : spec.series.front();
    std::string baseKey = primary;
    std::vector<std::string> groups;
    if (primary.find('/') != std::string::npos) {
        const auto [inst, out] = analysis::splitKey(primary);
        analysis::Instance parsed = analysis::parseInstance(inst);
        groups = std::move(parsed.groups);
        baseKey = parsed.base + "/" + out;
    }
    selectId(slot.type, toQString(spec.type));
    if (slot.type->currentIndex() < 0) {
        slot.type->addItem(tr("%1 (unknown)").arg(toQString(spec.type)), toQString(spec.type));
        selectId(slot.type, toQString(spec.type));
    }
    fillData(slot, spec.type, baseKey);
    fillGroups(slot);
    selectId(slot.ga, toQString(groups.size() > 0 ? groups[0] : std::string()));
    if (slot.ga->currentIndex() < 0) slot.ga->setCurrentIndex(0);
    selectId(slot.gb, toQString(groups.size() > 1 ? groups[1] : std::string()));
    if (slot.gb->currentIndex() < 0) slot.gb->setCurrentIndex(0);
    if (!keepText) {
        slot.more->showText(toQString(joined(spec.series, 1)));
        if (slot.optionsError.isEmpty()) slot.options->showText(optionsText(spec.options));
        slot.title->showText(toQString(spec.title));
    }
    const std::string title = spec.title;
    slot.autoTitle = title == autoTitle(slot) || isDefaultTitle(title) ? toQString(title) : QString();
    updateGroupWidgets(slot);
    --m_loading;
}

void LayoutPage::refreshSlots(bool force) {
    const core::Project& project = context().document().project();
    const std::vector<core::PanelSpec>& panels = project.layout.panels;
    for (const auto& s : m_slots) {
        Slot& slot = *s;
        const auto k = static_cast<std::size_t>(slot.k);
        const core::PanelSpec spec = k < panels.size() ? panels[k] : defaultSpec();
        if (force || !editorMatches(slot, spec)) {
            loadSlot(slot, spec, false);
            continue;
        }
        // the choices' labels name the settings as they are now; the selections stay
        ++m_loading;
        const std::string chosen = comboId(slot.data);
        fillData(slot, spec.type, chosen);
        fillGroups(slot);
        updateGroupWidgets(slot);
        if (!slot.more->dirty()) slot.more->showText(toQString(joined(spec.series, 1)));
        if (!slot.title->dirty()) slot.title->showText(toQString(spec.title));
        if (!slot.options->dirty() && slot.optionsError.isEmpty()) slot.options->showText(optionsText(spec.options));
        --m_loading;
    }
}

void LayoutPage::documentChanged(const core::DocumentChange& change) {
    using Kind = core::DocumentChange::Kind;
    if (change.kind == Kind::Saved && change.fields.empty()) return;
    const core::Project& project = context().document().project();
    const core::LayoutPreset* layout = core::findLayoutPreset(project.layout.preset);
    const int n = layout ? static_cast<int>(core::presetSlots(*layout).size()) : 0;
    if (change.kind == Kind::Replaced || n != slotCount() || change.affects("layout.preset") ||
        comboId(m_preset) != project.layout.preset) {
        rebuildSlots();
        return;
    }
    refreshSlots(false);
}

// ---- the user's changes -----------------------------------------------------------------------

QString LayoutPage::optionsProblem() const {
    for (const auto& slot : m_slots) {
        if (!slot->optionsError.isEmpty()) return slot->optionsError;
    }
    return {};
}

void LayoutPage::storeSlot(int k, const QString& label, const std::function<void(core::PanelSpec&)>& change) {
    core::ProjectDocument& doc = context().document();
    const auto index = static_cast<std::size_t>(k);
    doc.edit(toStdString(label), [&](core::Project& p) {
        std::vector<core::PanelSpec>& panels = p.layout.panels;
        // AP-07: the empty slots up to this one hold what their editors show
        while (panels.size() <= index) panels.push_back(editorSpec(*m_slots[panels.size()]));
        change(panels[index]);
    });
    if (const QString problem = optionsProblem(); !problem.isEmpty()) {
        emit previewNotUpdated(problem);  // MDATK's refresh_canvas: the options still do not parse
    }
}

void LayoutPage::onPreset(int index) {
    if (m_loading) return;
    const core::Project& project = context().document().project();
    const std::string chosen = toStdString(m_preset->itemData(index).toString());
    if (chosen == project.layout.preset) return;
    if (const QString problem = optionsProblem(); !problem.isEmpty()) {
        // AP-09 (MDATK 0.1.0): refused, the preset put back
        emit logLine(tr("Grid preset not changed; fix this first: %1").arg(problem));
        ++m_loading;
        selectId(m_preset, toQString(project.layout.preset));
        --m_loading;
        return;
    }
    context().document().set("layout.preset", core::Value(chosen), "Set Grid preset");
}

void LayoutPage::onType(Slot& slot) {
    if (m_loading) return;
    const std::string type = comboId(slot.type);
    const core::Project& project = context().document().project();
    const auto k = static_cast<std::size_t>(slot.k);
    const std::string before = k < project.layout.panels.size() ? project.layout.panels[k].type : defaultSpec().type;
    if (type == before) return;  // GTK emits no "changed" for the active item
    ++m_loading;
    if (const core::Table defaults = slotDefaultOptions(type);
        !defaults.empty() && strip(toStdString(slot.options->edit()->text())).empty()) {
        slot.options->showText(optionsText(defaults));
    }
    fillData(slot, type, std::string());
    --m_loading;
    onData(slot, true);
}

void LayoutPage::onData(Slot& slot, bool newType) {
    if (m_loading && !newType) return;
    ++m_loading;
    updateGroupWidgets(slot);
    const QString automatic = toQString(autoTitle(slot));
    const QString current = slot.title->edit()->text().trimmed();
    // a new panel type always gets its automatic title; within a type, a typed title is kept
    if (newType || current.isEmpty() || current == slot.autoTitle) slot.title->showText(automatic);
    slot.autoTitle = automatic;
    --m_loading;
    const core::PanelSpec spec = editorSpec(slot);
    storeSlot(slot.k, tr("Edit panel %1").arg(slot.k + 1), [&](core::PanelSpec& p) {
        p.type = spec.type;
        p.series = spec.series;
        p.title = spec.title;
        if (slot.optionsError.isEmpty()) p.options = spec.options;
    });
}

void LayoutPage::commitText(Slot& slot, const QString& field) {
    const std::string text = toStdString(slot.options->edit()->text());
    const auto k = static_cast<std::size_t>(slot.k);
    if (field == QStringLiteral("options")) {
        const std::string stripped = strip(text);
        std::string message;
        core::Table options;
        if (!stripped.empty()) {
            try {
                core::Value v = core::parseJson(stripped, "");
                if (const auto* t = v.get<core::Table>()) options = *t;
                else message = core::pythonStringRepr(stripped) + " is not a JSON object";
            } catch (const std::exception& e) {
                std::string reason = e.what();
                const std::size_t at = reason.find("not valid JSON: ");
                if (at != std::string::npos) reason = reason.substr(at + 16);
                message = core::pythonStringRepr(stripped) + " is not valid JSON: " + reason;
            }
        }
        const std::string path = slot.path("options");
        const bool wasError = !slot.optionsError.isEmpty();
        if (!message.empty()) {
            // located as the slot's panel would be (AP-10)
            core::Project where = context().document().project();
            while (where.layout.panels.size() <= k) where.layout.panels.push_back(defaultSpec());
            where.layout.panels[k].type = comboId(slot.type);
            const core::Problem problem = core::locatedProblem(path, message, &where);
            slot.optionsError = toQString(problem.text());
            context().setFormError(path, &problem);
            slot.options->flagChanged();
            emit previewNotUpdated(slot.optionsError);
            return;
        }
        slot.optionsError.clear();
        context().setFormError(path, nullptr);
        slot.options->flagChanged();
        const core::Project& project = context().document().project();
        const bool same = k < project.layout.panels.size() ? project.layout.panels[k].options == options : options.empty();
        if (same) {
            if (wasError) emit previewRequested();
            if (const QString problem = optionsProblem(); !problem.isEmpty()) emit previewNotUpdated(problem);
            return;
        }
        storeSlot(slot.k, tr("Set Panel %1 options").arg(slot.k + 1), [&](core::PanelSpec& p) { p.options = options; });
        return;
    }
    if (field == QStringLiteral("title")) {
        const std::string title = strip(toStdString(slot.title->edit()->text()));
        const core::Project& project = context().document().project();
        if (k < project.layout.panels.size() && project.layout.panels[k].title == title) return;
        storeSlot(slot.k, tr("Set Panel %1 title").arg(slot.k + 1), [&](core::PanelSpec& p) { p.title = title; });
        return;
    }
    // more keys
    const core::PanelSpec spec = editorSpec(slot);
    const core::Project& project = context().document().project();
    if (k < project.layout.panels.size() && !hasSeries(project.layout.panels[k].type)) return;
    std::vector<std::string> series;
    if (k < project.layout.panels.size() && !project.layout.panels[k].series.empty()) {
        series.push_back(project.layout.panels[k].series.front());
        for (std::string& key : splitKeys(toStdString(slot.more->edit()->text()))) series.push_back(std::move(key));
    } else {
        series = spec.series;
    }
    if (k < project.layout.panels.size() && project.layout.panels[k].series == series) return;
    storeSlot(slot.k, tr("Set Panel %1 keys").arg(slot.k + 1), [&](core::PanelSpec& p) { p.series = series; });
}

bool LayoutPage::swapPanels(int from, int to) {
    if (const QString problem = optionsProblem(); !problem.isEmpty()) {
        emit logLine(tr("Panels not swapped; fix this first: %1").arg(problem));  // AP-09 (MDATK 0.1.0)
        return false;
    }
    const int n = slotCount();
    if (from < 0 || to < 0 || std::max(from, to) >= n || from == to) return false;
    context().document().edit(toStdString(tr("Swap panels %1 and %2").arg(from + 1).arg(to + 1)), [&](core::Project& p) {
        std::vector<core::PanelSpec>& panels = p.layout.panels;
        const auto slotCount = static_cast<std::size_t>(n);
        std::vector<core::PanelSpec> specs;
        for (std::size_t k = 0; k < slotCount; ++k) {
            specs.push_back(k < panels.size() ? panels[k] : editorSpec(*m_slots[k]));
        }
        std::swap(specs[static_cast<std::size_t>(from)], specs[static_cast<std::size_t>(to)]);
        std::size_t count = std::min(panels.size(), slotCount);
        // a panel moved into an empty slot: that slot (and those before it) hold panels now (AP-07)
        if (static_cast<std::size_t>(std::min(from, to)) < count) {
            count = std::max(count, static_cast<std::size_t>(std::max(from, to)) + 1);
        }
        std::vector<core::PanelSpec> next(specs.begin(), specs.begin() + static_cast<std::ptrdiff_t>(count));
        for (std::size_t k = slotCount; k < panels.size(); ++k) next.push_back(panels[k]);  // AP-08: kept
        panels = std::move(next);
    });
    emit logLine(tr("Swapped panels %1 and %2").arg(from + 1).arg(to + 1));
    return true;
}

LayoutPage::SlotState LayoutPage::slotState(int k) const {
    SlotState s;
    if (k < 0 || k >= slotCount()) return s;
    const Slot& slot = *m_slots[static_cast<std::size_t>(k)];
    s.frame = slot.section->title();
    s.type = slot.type->currentData().toString();
    if (slot.data->count() > 0 && slot.data->currentIndex() >= 0) s.data = slot.data->currentData().toString();
    for (int i = 0; i < slot.data->count(); ++i) s.choices.emplace_back(slot.data->itemData(i).toString(), slot.data->itemText(i));
    s.groupA = slot.ga->currentData().toString();
    s.groupB = slot.gb->currentData().toString();
    s.groupAEnabled = slot.ga->isEnabled();
    s.groupBEnabled = slot.gb->isEnabled();
    s.more = slot.more->edit()->text();
    s.title = slot.title->edit()->text();
    s.options = slot.options->edit()->text();
    return s;
}

}  // namespace mdbridge::app::pages
