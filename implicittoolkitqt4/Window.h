#ifndef WINDOW_H
#define WINDOW_H

#include <QMainWindow>
#include <QMenu>
#include <QMenuBar>
#include <QApplication>

#include "ContentWidget.h"

class Window : public QMainWindow
{

  Q_OBJECT

public:

  Window( QWidget* parent = nullptr );

private slots:

  void loadMesh();
  void toggleGrid();
  void toggleNormals();
  void toggleSamplesVsConvex();
  void toggleSphereTree();
  void descendSphereTree();
  void ascendSphereTree();

private:

  std::string gridDisplayName() const;
  std::string normalDisplayName() const;
  std::string samplesVsConvexDisplayName() const;
  std::string toggleSphereTreeDisplayName() const;
  std::string descendSphereTreeDisplayName() const;
  std::string ascendSphereTreeDisplayName() const;

  ContentWidget* m_content_widget;
  QAction* m_toggle_grid_display;
  QAction* m_toggle_normal_display;
  QAction* m_toggle_samples_vs_convex_samples;
  QAction* m_toggle_sphere_tree_display;
  QAction* m_descend_sphere_tree_display;
  QAction* m_ascend_sphere_tree_display;

};

#endif
