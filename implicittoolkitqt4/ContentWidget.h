#ifndef CONTENT_WIDGET_H
#define CONTENT_WIDGET_H

#include <QWidget>
#include <QtGui>

#include "GLWidget.h"

class ContentWidget : public QWidget
{

  Q_OBJECT

public:

  ContentWidget( QWidget* parent = nullptr );
  
  GLWidget* glWidget();

  void toggleGridDisplay();
  bool gridDisplayStatus() const;

  void toggleNormalDisplay();
  bool normalDisplayStatus() const;

  void toggleSamplesVsConvex();
  bool samplesVsConvexDisplayStatus() const;

  void toggleSphereTreeDisplay();
  bool sphereTreeDisplayStatus() const;

  void descendSphereTree();
  unsigned nextLowestSphereTreeLevel() const;

  void ascendSphereTree();
  unsigned nextHighestSphereTreeLevel() const;

  void loadMesh();

public slots:

  void resetCamera();
  void computeSDF();
  void exportImage();

private slots:

  void setCellWidth( const double& cell_width );
  void setGridPadding( const double& grid_padding );
  void setLevelSetValue( const double& level_set_value );

  void renderInputMesh( const bool render );
  void renderLevelSet( const bool render );
  void renderSurfaceSampling( const bool render );

  void setEdgeSubsamples( const int num_subsamples );
  void setFaceSubsamples( const int num_subsamples );

  void saveProcessedMesh();

private:

  std::string getOpenFileNameFromUser();
  std::string getSaveFileNameFromUser( const std::string& display_text );

  GLWidget* m_glWidget;
  QDoubleSpinBox* m_cell_width_spinbox;
  QDoubleSpinBox* m_grid_padding_spinbox;

};

#endif
