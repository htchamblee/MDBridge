#include "gui/forms/list_editor.hpp"

#include <algorithm>
#include <exception>
#include <numeric>
#include <optional>
#include <utility>

#include <QBrush>
#include <QColor>
#include <QComboBox>
#include <QEvent>
#include <QHeaderView>
#include <QItemSelectionModel>
#include <QLineEdit>
#include <QPushButton>
#include <QStyle>
#include <QTableView>
#include <QToolButton>
#include <QTimer>
#include <QVBoxLayout>

#include "gui/forms/flow_layout.hpp"
#include "mdbridge/core/errors.hpp"
#include "mdbridge/core/python.hpp"

namespace mdbridge::app::forms {

namespace {

QColor tint(const QColor& base, const QString& level) {
    const QColor flag = level == QLatin1String("warning") ? QColor(0xe5, 0xa5, 0x0a) : QColor(0xe0, 0x1b, 0x24);
    return QColor::fromRgbF(static_cast<float>(base.redF() * 0.7 + flag.redF() * 0.3),
                            static_cast<float>(base.greenF() * 0.7 + flag.greenF() * 0.3),
                            static_cast<float>(base.blueF() * 0.7 + flag.blueF() * 0.3));
}

QString choiceText(const FieldSpec& spec, const QString& value) {
    for (const Choice& c : spec.choices) {
        if (c.value == value) {
            return c.label.isEmpty() ? c.value : c.label;
        }
    }
    return (value.isEmpty() ? QStringLiteral("\"\"") : value) + QStringLiteral(" (unknown)");
}

std::optional<core::FieldKind> kindAt(const core::Project& project, const std::string& path) {
    const auto info = core::fieldInfo(project, path);
    return info ? std::optional<core::FieldKind>(info->kind) : std::nullopt;
}

}  // namespace

// ---- ListModel -------------------------------------------------------------------------------

ListModel::ListModel(FormContext& context, std::string listPath, std::vector<ListColumn> columns, QObject* parent)
    : QAbstractTableModel(parent), m_context(context), m_listPath(std::move(listPath)), m_columns(std::move(columns)) {
    reload();
}

void ListModel::reload() {
    beginResetModel();
    m_rows = 0;
    if (core::fieldInfo(m_context.document().project(), m_listPath)) {
        const core::Value list = m_context.document().get(m_listPath);
        if (const core::Array* items = list.get<core::Array>()) {
            m_rows = static_cast<int>(items->size());
        }
    }
    endResetModel();
}

void ListModel::problemsChanged() {
    if (m_rows > 0 && !m_columns.empty()) {
        emit dataChanged(index(0, 0), index(m_rows - 1, static_cast<int>(m_columns.size()) - 1),
                         {Qt::BackgroundRole, Qt::ToolTipRole});
    }
}

int ListModel::rowCount(const QModelIndex& parent) const { return parent.isValid() ? 0 : m_rows; }

int ListModel::columnCount(const QModelIndex& parent) const {
    return parent.isValid() ? 0 : static_cast<int>(m_columns.size());
}

std::string ListModel::cellPath(int row, int column) const {
    std::string path = m_listPath + "." + std::to_string(row);
    const std::string& key = m_columns[static_cast<std::size_t>(column)].key;
    return key.empty() ? path : path + "." + key;
}

QVariant ListModel::headerData(int section, Qt::Orientation orientation, int role) const {
    if (orientation == Qt::Horizontal && section >= 0 && section < columnCount()) {
        const ListColumn& c = m_columns[static_cast<std::size_t>(section)];
        if (role == Qt::DisplayRole) {
            return c.title;
        }
        if (role == KeyRole) {
            return toQString(c.key);
        }
    }
    if (orientation == Qt::Vertical && role == Qt::DisplayRole) {
        return section + 1;
    }
    return {};
}

Qt::ItemFlags ListModel::flags(const QModelIndex& index) const {
    if (!index.isValid()) {
        return Qt::NoItemFlags;
    }
    const auto kind = kindAt(m_context.document().project(), cellPath(index.row(), index.column()));
    if (!kind) {
        return Qt::ItemIsSelectable;
    }
    Qt::ItemFlags out = Qt::ItemIsSelectable | Qt::ItemIsEnabled;
    if (*kind == core::FieldKind::Boolean) {
        return out | Qt::ItemIsUserCheckable;
    }
    return out | Qt::ItemIsEditable;
}

QVariant ListModel::data(const QModelIndex& index, int role) const {
    if (!index.isValid()) {
        return {};
    }
    const std::string path = cellPath(index.row(), index.column());
    const core::ProjectDocument& doc = m_context.document();
    const auto kind = kindAt(doc.project(), path);
    if (!kind) {
        return {};
    }
    const FieldSpec& spec = m_columns[static_cast<std::size_t>(index.column())].spec;
    QString text;
    try {
        text = toQString(doc.text(path));
    } catch (const core::FieldError&) {
        return {};
    }
    switch (role) {
    case PathRole:
        return toQString(path);
    case Qt::EditRole:
        return *kind == core::FieldKind::Boolean ? QVariant() : QVariant(text);
    case Qt::DisplayRole:
        if (*kind == core::FieldKind::Boolean) {
            return {};
        }
        return spec.editor == Editor::Combo ? choiceText(spec, text) : text;
    case Qt::CheckStateRole:
        if (*kind == core::FieldKind::Boolean) {
            return static_cast<int>(text == QLatin1String("true") ? Qt::Checked : Qt::Unchecked);
        }
        return {};
    case Qt::DecorationRole:
        if (spec.editor == Editor::Color) {
            return colorSwatch(matplotlibColor(text));
        }
        return {};
    case Qt::BackgroundRole:
    case Qt::ToolTipRole: {
        std::vector<core::Problem> problems = m_context.problemsAt(path, true);
        if (index.column() == 0) {  // problems of the whole row (item) on its first cell
            const std::vector<core::Problem> row = m_context.problemsAt(m_listPath + "." + std::to_string(index.row()));
            problems.insert(problems.end(), row.begin(), row.end());
        }
        if (problems.empty()) {
            return role == Qt::ToolTipRole && !spec.hint.isEmpty() ? QVariant(spec.hint) : QVariant();
        }
        if (role == Qt::BackgroundRole) {
            return QBrush(tint(QColor(255, 255, 255), toQString(problems.front().level)));
        }
        QStringList lines;
        for (const core::Problem& p : problems) {
            lines.append(toQString(p.message));
        }
        lines.removeDuplicates();
        return lines.join(QLatin1Char('\n'));
    }
    default:
        return {};
    }
}

bool ListModel::setData(const QModelIndex& index, const QVariant& value, int role) {
    if (!index.isValid() || !store) {
        return false;
    }
    const std::string path = cellPath(index.row(), index.column());
    if (role == Qt::CheckStateRole) {
        const bool checked = value.toInt() == Qt::Checked;
        return store(path, checked ? QStringLiteral("true") : QStringLiteral("false"));
    }
    if (role == Qt::EditRole) {
        return store(path, value.toString());
    }
    return false;
}

// ---- ListDelegate ----------------------------------------------------------------------------

ListDelegate::ListDelegate(const ListModel& model, QObject* parent) : QStyledItemDelegate(parent), m_model(model) {}

QWidget* ListDelegate::createEditor(QWidget* parent, const QStyleOptionViewItem& option, const QModelIndex& index) const {
    const FieldSpec& spec = m_model.columns()[static_cast<std::size_t>(index.column())].spec;
    QWidget* editor = nullptr;
    if (spec.editor == Editor::Combo || spec.editor == Editor::EditableCombo) {
        auto* combo = new QComboBox(parent);
        if (spec.editor == Editor::EditableCombo) {
            combo->setEditable(true);
            combo->setInsertPolicy(QComboBox::NoInsert);
            combo->setDuplicatesEnabled(true);
            combo->setCompleter(nullptr);
            for (const Choice& c : spec.choices) {
                combo->addItem(c.value, c.value);
            }
        } else {
            for (const Choice& c : spec.choices) {
                combo->addItem(c.label.isEmpty() ? c.value : c.label, c.value);
            }
        }
        guardWheel(combo);
        editor = combo;
    } else {
        auto* edit = new QLineEdit(parent);
        edit->setFrame(false);
        edit->setPlaceholderText(spec.placeholder);
        editor = edit;
    }
    (void)option;
    editor->setObjectName(index.data(ListModel::PathRole).toString());
    editor->setAccessibleName(m_model.headerData(index.column(), Qt::Horizontal, Qt::DisplayRole).toString());
    m_open.erase(std::remove_if(m_open.begin(), m_open.end(), [](const QPointer<QWidget>& w) { return w.isNull(); }),
                 m_open.end());
    m_open.emplace_back(editor);
    return editor;
}

void ListDelegate::setEditorData(QWidget* editor, const QModelIndex& index) const {
    const QString text = index.data(Qt::EditRole).toString();
    if (auto* combo = qobject_cast<QComboBox*>(editor)) {
        if (combo->isEditable()) {
            if (!combo->lineEdit()->isModified()) {
                combo->lineEdit()->setText(text);
            }
            return;
        }
        int i = combo->findData(text);
        if (i < 0) {  // an unknown value is kept (AP-05)
            combo->addItem(choiceText({}, text), text);
            i = combo->count() - 1;
        }
        combo->setCurrentIndex(i);
        return;
    }
    if (auto* edit = qobject_cast<QLineEdit*>(editor)) {
        if (!edit->isModified()) {
            edit->setText(text);
        }
    }
}

void ListDelegate::setModelData(QWidget* editor, QAbstractItemModel* model, const QModelIndex& index) const {
    if (auto* combo = qobject_cast<QComboBox*>(editor)) {
        const QString text = combo->isEditable() ? combo->lineEdit()->text() : combo->currentData().toString();
        if (text != index.data(Qt::EditRole).toString()) {
            model->setData(index, text, Qt::EditRole);
        }
        return;
    }
    if (auto* edit = qobject_cast<QLineEdit*>(editor)) {
        if (edit->isModified() && model->setData(index, edit->text(), Qt::EditRole)) {
            edit->setModified(false);
        }
    }
}

void ListDelegate::destroyEditor(QWidget* editor, const QModelIndex& index) const {
    m_open.erase(std::remove_if(m_open.begin(), m_open.end(),
                                [editor](const QPointer<QWidget>& w) { return w.isNull() || w.data() == editor; }),
                 m_open.end());
    QStyledItemDelegate::destroyEditor(editor, index);
}

void ListDelegate::commitOpenEditors() {
    const auto open = m_open;
    for (const QPointer<QWidget>& editor : open) {
        if (editor) {
            emit commitData(editor);
        }
    }
}

// ---- ListEditor ------------------------------------------------------------------------------

ListEditor::ListEditor(FormContext& context, std::string listPath, FieldSpec spec, QWidget* parent)
    : FieldWidget(context, std::move(listPath), std::move(spec), parent) {
    auto* layout = new QVBoxLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);
    layout->setSpacing(4);
    m_model = new ListModel(context, path(), columnsFor(this->spec()), this);
    m_model->store = [this](const std::string& cell, const QString& text) { return storeCell(cell, text); };
    m_view = new QTableView(this);
    m_view->setModel(m_model);
    m_delegate = new ListDelegate(*m_model, m_view);
    m_view->setItemDelegate(m_delegate);
    m_view->setSelectionBehavior(QAbstractItemView::SelectRows);
    m_view->setSelectionMode(QAbstractItemView::SingleSelection);
    m_view->setEditTriggers(QAbstractItemView::DoubleClicked | QAbstractItemView::EditKeyPressed |
                            QAbstractItemView::AnyKeyPressed | QAbstractItemView::SelectedClicked);
    m_view->setTabKeyNavigation(false);  // Tab leaves the table (keyboard users reach the buttons)
    m_view->horizontalHeader()->setStretchLastSection(false);  // fitColumns() sizes them all
    m_view->verticalHeader()->setVisible(false);
    m_view->setMinimumHeight(90);
    m_view->setProperty("mdbridgeList", true);
    layout->addWidget(m_view);
    setEditor(m_view);

    // the buttons wrap onto a further line when the field column is too narrow for all of them
    // (fonts wider than the one the side bar was laid out with), as the form rows do
    auto* buttonRow = new QWidget(this);
    auto* buttons = new FlowLayout(buttonRow, 4, 4);
    buttons->setContentsMargins(0, 0, 0, 0);
    m_buttonRow = buttons;
    const QString name = toQString(path());
    const auto pushButton = [&](const QString& text, const QString& id, QStyle::StandardPixmap icon, const QString& tip) {
        auto* b = new QPushButton(style()->standardIcon(icon), text, this);
        b->setObjectName(name + id);
        b->setToolTip(tip);
        b->setAccessibleName(tip);
        buttons->addWidget(b);
        return b;
    };
    const auto toolButton = [&](Qt::ArrowType arrow, const QString& id, const QString& tip) {
        auto* b = new QToolButton(this);
        b->setArrowType(arrow);
        b->setObjectName(name + id);
        b->setToolTip(tip);
        b->setAccessibleName(tip);
        b->setFocusPolicy(Qt::StrongFocus);
        buttons->addWidget(b);
        return b;
    };
    m_add = pushButton(QStringLiteral("Add"), QStringLiteral(":add"), QStyle::SP_FileDialogNewFolder,
                       QStringLiteral("Add a row"));
    m_remove = pushButton(QStringLiteral("Remove"), QStringLiteral(":remove"), QStyle::SP_DialogDiscardButton,
                          QStringLiteral("Remove the selected row"));
    m_up = toolButton(Qt::UpArrow, QStringLiteral(":up"), QStringLiteral("Move the selected row up"));
    m_down = toolButton(Qt::DownArrow, QStringLiteral(":down"), QStringLiteral("Move the selected row down"));
    if (this->spec().browse) {
        m_browse = pushButton(QStringLiteral("Browse…"), QStringLiteral(":browse"), QStyle::SP_DialogOpenButton,
                              QStringLiteral("Browse: ") +
                                  (this->spec().dialogTitle.isEmpty() ? label() : this->spec().dialogTitle));
    }
    layout->addWidget(buttonRow);
    connect(m_add, &QPushButton::clicked, this, &ListEditor::addRow);
    connect(m_remove, &QPushButton::clicked, this, &ListEditor::removeSelected);
    connect(m_up, &QToolButton::clicked, this, [this] { moveSelected(-1); });
    connect(m_down, &QToolButton::clicked, this, [this] { moveSelected(1); });
    if (m_browse) {
        connect(m_browse, &QPushButton::clicked, this, &ListEditor::browse);
    }
    connect(&context, &FormContext::problemsChanged, m_model, &ListModel::problemsChanged);
    m_view->viewport()->installEventFilter(this);
    connect(m_model, &QAbstractItemModel::modelReset, this, [this] {
        QTimer::singleShot(0, this, &ListEditor::fitColumns);  // after the header has its sections again
    });
    refresh();
}

bool ListEditor::eventFilter(QObject* watched, QEvent* event) {
    if (watched == m_view->viewport() && event->type() == QEvent::Resize) {
        fitColumns();
    }
    return FieldWidget::eventFilter(watched, event);
}

void ListEditor::fitColumns() {
    QHeaderView* header = m_view->horizontalHeader();
    const int columns = m_model->columnCount();
    const int available = m_view->viewport()->width();
    if (columns <= 0 || available <= 0) {
        return;
    }
    std::vector<int> natural(static_cast<std::size_t>(columns));
    const int rows = std::min(m_model->rowCount(), 200);
    const int grid = m_view->showGrid() ? 1 : 0;
    for (int c = 0; c < columns; ++c) {
        int fit = std::max(header->sectionSizeHint(c), header->minimumSectionSize());
        for (int r = 0; r < rows; ++r) {
            fit = std::max(fit, m_view->sizeHintForIndex(m_model->index(r, c)).width() + grid);
        }
        natural[static_cast<std::size_t>(c)] = fit;
    }
    int cap = *std::max_element(natural.begin(), natural.end());
    const auto total = [&natural](int limit) {
        return std::accumulate(natural.begin(), natural.end(), 0, [limit](int sum, int w) { return sum + std::min(w, limit); });
    };
    if (total(cap) > available) {
        // the largest cap at which everything fits (water filling), never below the header's minimum
        int low = header->minimumSectionSize();
        int high = cap;
        while (low < high) {
            const int mid = (low + high + 1) / 2;
            if (total(mid) <= available) {
                low = mid;
            } else {
                high = mid - 1;
            }
        }
        cap = low;
    }
    int used = 0;
    for (int c = 0; c + 1 < columns; ++c) {
        const int section = std::min(natural[static_cast<std::size_t>(c)], cap);
        header->resizeSection(c, section);
        used += section;
    }
    // the last column takes the rest (the header's own last-section stretch keeps it at least
    // its default size, which made narrow tables scroll)
    header->resizeSection(columns - 1, std::max(available - used, header->minimumSectionSize()));
}

void ListEditor::addButton(QWidget* button) {
    button->setParent(this);
    m_buttonRow->addWidget(button);
}

std::vector<ListColumn> ListEditor::columnsFor(const FieldSpec& spec) const {
    const core::Project& project = document().project();
    const auto info = core::fieldInfo(project, path());
    const std::string modelPath = info ? info->modelPath : path();
    const auto columnSpec = [&](const std::string& key) {
        return key.empty() ? FieldSpec{} : defaultFieldSpec(modelPath + "." + key);
    };
    std::vector<ListColumn> out;
    for (const Column& c : spec.columns) {
        out.push_back({c.key, c.title, columnSpec(c.key)});
    }
    if (!out.empty() || !info) {
        return out;
    }
    switch (info->kind) {
    case core::FieldKind::Records: {
        const core::Value item = newItem();
        if (const core::Table* table = item.get<core::Table>()) {
            for (const auto& [key, value] : *table) {
                QString title = toQString(key);
                title.replace(QLatin1Char('_'), QLatin1Char(' '));
                if (!title.isEmpty()) {
                    title[0] = title[0].toUpper();
                }
                out.push_back({key, title, columnSpec(key)});
            }
        }
        break;
    }
    case core::FieldKind::NumberRows: {
        std::size_t width = 3;
        const core::Value list = document().get(path());
        if (const core::Array* rows = list.get<core::Array>()) {
            for (const core::Value& row : *rows) {
                if (const core::Array* numbers = row.get<core::Array>()) {
                    width = std::max(width, numbers->size());
                }
            }
        }
        for (std::size_t k = 0; k < width; ++k) {
            out.push_back({std::to_string(k), QString::number(k + 1), FieldSpec{}});
        }
        break;
    }
    default:
        out.push_back({"", label(), FieldSpec{}});
        break;
    }
    return out;
}

core::Value ListEditor::newItem() const {
    if (!spec().newItem.isNull()) {
        return spec().newItem;
    }
    const auto info = core::fieldInfo(document().project(), path());
    if (!info) {
        return {};
    }
    switch (info->kind) {
    case core::FieldKind::TextList:
        return core::Value(std::string());
    case core::FieldKind::NumberList:
        return core::Value(0.0);
    case core::FieldKind::NumberRows:
        return core::Value(core::Array{core::Value(0.0), core::Value(0.0), core::Value(0.0)});
    case core::FieldKind::Records:
        try {  // the record's defaults: what an empty table reads as
            core::ProjectDocument scratch;
            scratch.replace(document().project());
            scratch.insertItem(path(), 0, core::Value(core::Table{}));
            return scratch.get(path() + ".0");
        } catch (const std::exception&) {
            return core::Value(core::Table{});
        }
    default:
        return {};
    }
}

QModelIndex ListEditor::cellIndex(const std::string& cellPath) const {
    for (int r = 0; r < m_model->rowCount(); ++r) {
        for (int c = 0; c < m_model->columnCount(); ++c) {
            if (m_model->cellPath(r, c) == cellPath) {
                return m_model->index(r, c);
            }
        }
    }
    return {};
}

int ListEditor::selectedRow() const {
    const QModelIndexList rows = m_view->selectionModel()->selectedRows();
    if (!rows.isEmpty()) {
        return rows.front().row();
    }
    const QModelIndex current = m_view->currentIndex();
    return current.isValid() ? current.row() : -1;
}

void ListEditor::selectRow(int row) {
    if (row < 0 || row >= m_model->rowCount()) {
        m_view->clearSelection();
        return;
    }
    m_view->setCurrentIndex(m_model->index(row, 0));
    m_view->selectRow(row);
}

bool ListEditor::storeCell(const std::string& cellPath, const QString& text) {
    // while a cell's editor commits, a change of the document only updates the cells (a reset of
    // the model would delete the editor under the view's feet)
    m_writing = true;
    bool stored = false;
    try {
        const QModelIndex index = cellIndex(cellPath);
        if (m_cellEdit && index.isValid() &&
            m_cellEdit(index.row(), m_model->columns()[static_cast<std::size_t>(index.column())].key, text)) {
            stored = true;
        } else {
            stored = commitText(text, cellPath);
        }
    } catch (...) {
        m_writing = false;
        throw;
    }
    m_writing = false;
    return stored;
}

void ListEditor::addRow() {
    commit();
    if (!fieldExists()) {
        return;
    }
    const int at = m_model->rowCount();
    try {
        document().insertItem(path(), static_cast<std::size_t>(at), newItem());
    } catch (const core::FieldError& error) {
        const core::Problem problem = error.problem(&document().project());
        context().setFormError(path(), &problem);
        return;
    }
    selectRow(at);
}

void ListEditor::removeSelected() {
    commit();
    const int row = selectedRow();
    if (row < 0 || row >= m_model->rowCount()) {
        return;
    }
    document().removeItem(path(), static_cast<std::size_t>(row));
    selectRow(std::min(row, m_model->rowCount() - 1));
}

void ListEditor::moveSelected(int delta) {
    commit();
    const int row = selectedRow();
    const int to = row + delta;
    if (row < 0 || to < 0 || to >= m_model->rowCount()) {
        return;
    }
    document().moveItem(path(), static_cast<std::size_t>(row), static_cast<std::size_t>(to));
    selectRow(to);
}

void ListEditor::browse() {
    commit();
    if (!fieldExists()) {
        return;
    }
    FileDialogRequest request;
    request.mode = spec().pathMode;
    request.title = spec().dialogTitle.isEmpty() ? label() : spec().dialogTitle;
    request.start = toQString(document().folder());
    request.filters = nameFilters(spec());
    const QStringList files = runFileDialog(window(), request);
    if (files.isEmpty()) {
        return;
    }
    core::ProjectDocument& doc = document();
    doc.beginGroup(toStdString(QStringLiteral("Add files to ") + label()));
    try {
        for (const QString& file : files) {
            const core::Value text(toStdString(pickedPathText(doc, file)));
            core::Value item = text;
            if (!spec().browseKey.empty()) {
                item = newItem();
                if (core::Table* table = item.get<core::Table>()) {
                    auto it = std::find_if(table->begin(), table->end(),
                                           [&](const auto& e) { return e.first == spec().browseKey; });
                    if (it != table->end()) {
                        it->second = text;
                    } else {
                        table->emplace_back(spec().browseKey, text);
                    }
                }
            }
            doc.insertItem(path(), static_cast<std::size_t>(m_model->rowCount()), item);
        }
    } catch (...) {
        doc.endGroup();
        throw;
    }
    doc.endGroup();
    selectRow(m_model->rowCount() - 1);
}

void ListEditor::refresh() {
    const int row = m_view ? selectedRow() : -1;
    m_model->reload();
    m_view->setEnabled(fieldExists());
    if (row >= 0) {
        selectRow(std::min(row, m_model->rowCount() - 1));
    }
}

void ListEditor::commit() { m_delegate->commitOpenEditors(); }

void ListEditor::onDocumentChanged(const core::DocumentChange& change) {
    if (!change.affects(path())) {
        return;
    }
    if (m_writing && change.kind == core::DocumentChange::Kind::Edit) {
        // our own cell edit: the rows stay, the open editor stays valid
        m_model->problemsChanged();
        if (m_model->rowCount() > 0 && m_model->columnCount() > 0) {
            emit m_model->dataChanged(m_model->index(0, 0),
                                      m_model->index(m_model->rowCount() - 1, m_model->columnCount() - 1));
        }
    } else {
        refresh();
    }
    updateFlag();
}

CellLocation findCell(QWidget* root, const std::string& cellPath) {
    if (!root) {
        return {};
    }
    for (ListEditor* editor : root->findChildren<ListEditor*>()) {
        const std::string& list = editor->path();
        if (cellPath.size() > list.size() && cellPath.compare(0, list.size(), list) == 0 && cellPath[list.size()] == '.') {
            const QModelIndex index = editor->cellIndex(cellPath);
            if (index.isValid()) {
                return {editor, index};
            }
        }
    }
    return {};
}

}  // namespace mdbridge::app::forms
