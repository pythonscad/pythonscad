#include "gui/UnsavedChangesDialog.h"

#include <QApplication>
#include <QFileInfo>
#include <QFont>
#include <QHBoxLayout>
#include <QShortcut>
#include <QStyle>
#include <QToolButton>
#include <QVariant>
#include <QVBoxLayout>

#include "gui/Editor.h"
#include "gui/MainWindow.h"
#include "gui/TabManager.h"
#include "gui/parameter/ParameterWidget.h"

#include <genlang/genlang.h>

UnsavedChangesDialog::UnsavedChangesDialog(TabManager *tabManager, MainWindow *mainWindow,
                                           QWidget *parent)
  : QDialog(parent), tabManager(tabManager), mainWindow(mainWindow), dialogResult(Cancel)
{
  setWindowTitle(_("Unsaved Changes"));
  setModal(true);
  setMinimumWidth(400);
  setMinimumHeight(300);

  auto *layout = new QVBoxLayout(this);

  countLabel = new QLabel(this);
  layout->addWidget(countLabel);

  listWidget = new QListWidget(this);
  listWidget->setSelectionMode(QAbstractItemView::SingleSelection);
  listWidget->setSpacing(2);
  connect(listWidget, &QListWidget::itemClicked, this, &UnsavedChangesDialog::onItemClicked);
  layout->addWidget(listWidget);

  auto *warningLabel = new QLabel(_("If you quit PythonSCAD now, these changes will be lost."), this);
  layout->addWidget(warningLabel);

  auto *shortcutLabel = new QLabel(_("Press Ctrl+D to discard changes and quit."), this);
  QFont italicFont = shortcutLabel->font();
  italicFont.setItalic(true);
  italicFont.setPointSize(italicFont.pointSize() - 1);
  shortcutLabel->setFont(italicFont);
  layout->addWidget(shortcutLabel);

  layout->addSpacing(10);

  auto *buttonLayout = new QHBoxLayout();
  buttonLayout->addStretch();

  cancelButton = new QPushButton(_("Cancel"), this);
  connect(cancelButton, &QPushButton::clicked, this, &QDialog::reject);
  buttonLayout->addWidget(cancelButton);

  discardButton = new QPushButton(_("Discard Changes"), this);
  connect(discardButton, &QPushButton::clicked, this, &UnsavedChangesDialog::onDiscardClicked);
  buttonLayout->addWidget(discardButton);

  layout->addLayout(buttonLayout);

  auto *discardShortcut = new QShortcut(QKeySequence("Ctrl+D"), this);
  connect(discardShortcut, &QShortcut::activated, this, &UnsavedChangesDialog::onDiscardClicked);

  populateList();
  updateCountLabel();
}

void UnsavedChangesDialog::populateList()
{
  for (EditorInterface *editor : tabManager->editorList) {
    if (editor->isContentModified() || editor->parameterWidget->isModified()) {
      auto *item = new QListWidgetItem(listWidget);
      item->setSizeHint(QSize(0, 36));

      auto *widget = createItemWidget(editor);
      listWidget->setItemWidget(item, widget);

      itemToEditor[item] = editor;
      editorToItem[editor] = item;
    }
  }
}

QWidget *UnsavedChangesDialog::createItemWidget(EditorInterface *editor)
{
  auto *widget = new QWidget();
  auto *layout = new QHBoxLayout(widget);
  layout->setContentsMargins(8, 4, 8, 4);

  auto *nameLabel = new QLabel(getDisplayName(editor), widget);
  nameLabel->setToolTip(editor->filepath.isEmpty() ? _("Untitled") : editor->filepath);
  layout->addWidget(nameLabel, 1);

  auto *saveButton = new QToolButton(widget);
  saveButton->setIcon(QApplication::style()->standardIcon(QStyle::SP_DialogSaveButton));
  saveButton->setToolTip(_("Save this file"));
  saveButton->setProperty("editor", QVariant::fromValue(static_cast<void *>(editor)));
  connect(saveButton, &QToolButton::clicked, this, &UnsavedChangesDialog::onSaveButtonClicked);
  layout->addWidget(saveButton);

  return widget;
}

QString UnsavedChangesDialog::getDisplayName(EditorInterface *editor) const
{
  if (editor->filepath.isEmpty()) {
    return _("Untitled");
  }
  QFileInfo info(editor->filepath);
  return info.fileName();
}

void UnsavedChangesDialog::updateCountLabel()
{
  const int count = listWidget->count();
  if (count == 1) {
    countLabel->setText(_("There is 1 file with unsaved changes:"));
  } else {
    countLabel->setText(QString(_("There are %1 files with unsaved changes:")).arg(count));
  }
}

void UnsavedChangesDialog::onItemClicked(QListWidgetItem *item)
{
  if (!itemToEditor.contains(item)) return;

  EditorInterface *editor = itemToEditor[item];
  tabManager->switchToEditor(editor);
  mainWindow->actionRenderPreview();
}

void UnsavedChangesDialog::onSaveButtonClicked()
{
  auto *button = qobject_cast<QToolButton *>(sender());
  if (!button) return;

  auto *editor = static_cast<EditorInterface *>(button->property("editor").value<void *>());
  if (!editor) return;

  tabManager->switchToEditor(editor);

  const bool saved = tabManager->save(editor);
  if (saved) {
    removeListItem(editor);
    updateCountLabel();

    if (listWidget->count() == 0) {
      dialogResult = AllSaved;
      accept();
    }
  }
}

void UnsavedChangesDialog::onDiscardClicked()
{
  dialogResult = DiscardAll;
  accept();
}

void UnsavedChangesDialog::removeListItem(EditorInterface *editor)
{
  if (!editorToItem.contains(editor)) return;

  QListWidgetItem *item = editorToItem[editor];
  const int row = listWidget->row(item);

  QWidget *widget = listWidget->itemWidget(item);
  listWidget->removeItemWidget(item);
  delete widget;

  listWidget->takeItem(row);
  delete item;

  itemToEditor.remove(item);
  editorToItem.remove(editor);
}

void UnsavedChangesDialog::updateListItem(EditorInterface *editor)
{
  if (!editorToItem.contains(editor)) return;

  QListWidgetItem *item = editorToItem[editor];
  QWidget *oldWidget = listWidget->itemWidget(item);
  listWidget->removeItemWidget(item);
  delete oldWidget;

  QWidget *newWidget = createItemWidget(editor);
  listWidget->setItemWidget(item, newWidget);
}
