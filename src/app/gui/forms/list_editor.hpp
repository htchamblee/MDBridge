#pragma once
// The list editor (milestone 5 W0-C): a table of a list field's items (MDATK's ListEditor):
// lists of records (groups, inputs.external, analyses.rmsd.targets, visual.representations), of
// text (inputs.trajectories) and of number rows (density.anchor), with Add / Remove / Move up /
// Move down (and Browse… for files) as one undo step each (insertItem / removeItem / moveItem).
//
// A cell is the field "<list>.<row>.<key>" ("<list>.<row>" for plain items) and is stored with
// setText when its editor closes: Enter, and also when the cell loses the focus to a click
// elsewhere (AP-15, as MDATK 0.1.1's keep_edit_on_leave) or on commitPendingEdits(). Cells with a
// problem are tinted and carry the message as tooltip; problems of the list itself flag the table.
//
// Names: the QTableView is the list path ("groups"); the buttons "<list>:add", ":remove", ":up",
// ":down", ":browse"; the editor widget "<list>:field". findCell() maps a cell's field path to the
// view and index (the GUI harness edits cells through it).

#include <functional>
#include <string>
#include <vector>

#include <QAbstractTableModel>
#include <QPointer>
#include <QStyledItemDelegate>

#include "gui/forms/field_widgets.hpp"

class QPushButton;
class QTableView;
class QToolButton;

namespace mdbridge::app::forms {

struct ListColumn {
    std::string key;  // "" for plain items
    QString title;
    FieldSpec spec;   // the cell's editor: Combo, EditableCombo, Color, Check (Boolean cells), else text
};

class ListModel : public QAbstractTableModel {
    Q_OBJECT

public:
    // The role of header data that holds a column's key.
    static constexpr int KeyRole = Qt::UserRole + 1;
    // The role of cell data that holds the cell's field path.
    static constexpr int PathRole = Qt::UserRole + 2;

    ListModel(FormContext& context, std::string listPath, std::vector<ListColumn> columns, QObject* parent = nullptr);

    int rowCount(const QModelIndex& parent = {}) const override;
    int columnCount(const QModelIndex& parent = {}) const override;
    QVariant data(const QModelIndex& index, int role) const override;
    QVariant headerData(int section, Qt::Orientation orientation, int role) const override;
    Qt::ItemFlags flags(const QModelIndex& index) const override;
    bool setData(const QModelIndex& index, const QVariant& value, int role) override;

    const std::vector<ListColumn>& columns() const { return m_columns; }
    std::string cellPath(int row, int column) const;
    // Re-reads the list (after a change of the document).
    void reload();
    // Problems changed: tints and tooltips.
    void problemsChanged();
    // Stores a cell's text: (cell path, text) -> stored; set by ListEditor.
    std::function<bool(const std::string&, const QString&)> store;

private:
    FormContext& m_context;
    std::string m_listPath;
    std::vector<ListColumn> m_columns;
    int m_rows = 0;
};

// Opens text, combo and editable-combo editors for cells; commits when they lose the focus (Qt's
// delegates do, AP-15) and on commitOpenEditors().
class ListDelegate : public QStyledItemDelegate {
    Q_OBJECT

public:
    explicit ListDelegate(const ListModel& model, QObject* parent = nullptr);
    QWidget* createEditor(QWidget* parent, const QStyleOptionViewItem& option, const QModelIndex& index) const override;
    void setEditorData(QWidget* editor, const QModelIndex& index) const override;
    void setModelData(QWidget* editor, QAbstractItemModel* model, const QModelIndex& index) const override;
    void destroyEditor(QWidget* editor, const QModelIndex& index) const override;
    // Stores what the open editors hold (before a save, or a button of the editor).
    void commitOpenEditors();

private:
    const ListModel& m_model;
    mutable std::vector<QPointer<QWidget>> m_open;
};

class ListEditor : public FieldWidget {
    Q_OBJECT

public:
    ListEditor(FormContext& context, std::string listPath, FieldSpec spec, QWidget* parent = nullptr);

    QTableView* view() const { return m_view; }
    ListModel* model() const { return m_model; }
    int rowCount() const { return m_model->rowCount(); }
    // The index of the cell whose field path is `cellPath` (invalid when it shows none).
    QModelIndex cellIndex(const std::string& cellPath) const;
    int selectedRow() const;
    void selectRow(int row);

    // The buttons' actions (one undo step each; Browse… adds every file picked as one step).
    void addRow();
    void removeSelected();
    void moveSelected(int delta);
    void browse();
    QPushButton* addButton() const { return m_add; }
    QPushButton* removeButton() const { return m_remove; }
    QToolButton* upButton() const { return m_up; }
    QToolButton* downButton() const { return m_down; }
    QPushButton* browseButton() const { return m_browse; }

    // Called before a typed cell value is stored (row, column key, text); returning true means it
    // dealt with the edit itself (e.g. a group rename with its cascade, MDATK 0.1.1 on_edit).
    using CellEditHandler = std::function<bool(int row, const std::string& key, const QString& text)>;
    void setCellEditHandler(CellEditHandler handler) { m_cellEdit = std::move(handler); }

    // The item an Add inserts (spec.newItem, else the record's defaults, "" or zeros).
    core::Value newItem() const;

    void refresh() override;
    void commit() override;
    QString displayedText() const override { return {}; }

protected:
    void onDocumentChanged(const core::DocumentChange& change) override;

private:
    bool storeCell(const std::string& cellPath, const QString& text);
    std::vector<ListColumn> columnsFor(const FieldSpec& spec) const;

    QTableView* m_view = nullptr;
    ListModel* m_model = nullptr;
    ListDelegate* m_delegate = nullptr;
    QPushButton* m_add = nullptr;
    QPushButton* m_remove = nullptr;
    QToolButton* m_up = nullptr;
    QToolButton* m_down = nullptr;
    QPushButton* m_browse = nullptr;
    CellEditHandler m_cellEdit;
    bool m_writing = false;
};

// The list editor, view and index of the cell whose field path is `cellPath` ("groups.1.color",
// "inputs.trajectories.2") among the widgets below `root`; editor null when no list editor shows
// it.
struct CellLocation {
    ListEditor* editor = nullptr;
    QModelIndex index;
};
CellLocation findCell(QWidget* root, const std::string& cellPath);

}  // namespace mdbridge::app::forms
