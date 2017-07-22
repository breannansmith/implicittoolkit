#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QGLWidget>
#include <QtGui>
#include <QtOpenGL>
#include <cmath>
#include <iostream>
#include <fstream>

#include "implicitutilities/MathDefines.h"
#include "implicitutilities/ObjParser.h"
#include "implicitutilities/DistanceTools.h"
#include "implicitutilities/SignedDistanceField.h"
#include "implicitutilities/SurfaceSampling.h"
#include "implicitutilities/Moments.h"
#include "implicitutilities/SphereTree.h"

#include "marchingcubes/MarchingCubes.h"

#include "PerspectiveCameraController.h"
#include "OrthographicCameraController.h"

#include "OpenGL3DSphereRenderer.h"

class GLWidget : public QGLWidget
{

  Q_OBJECT

public:

  GLWidget( QWidget* parent );
  ~GLWidget() = default;

  QSize minimumSizeHint() const;
  QSize sizeHint() const;

  const scalar& cellWidth() const;
  const scalar& gridPadding() const;
  const scalar& isosurfaceValue() const;
  bool drawInputMesh() const;
  bool drawSDFMesh() const;
  bool drawSurfaceSampling() const;
  bool drawSampleNormals() const;
  void toggleDrawNormals();
  bool drawGrid() const;
  void toggleDrawGrid();
  bool samplesVsConvexDisplayStatus() const;
  void toggleSamplesVsConvex();
  bool sphereTreeDisplayStatus() const;
  void toggleSphereTreeDisplay();
  unsigned nextLowestSphereTreeLevel() const;
  void descendSphereTree();
  unsigned nextHighestSphereTreeLevel() const;
  void ascendSphereTree();
  void setCellWidth( const scalar& cell_width );
  void setGridPadding( const scalar& grid_padding );
  void setInputMeshRendering( const bool render );
  void setIsosurfaceRendering( const bool render );
  void setSurfaceSamplingRendering( const bool render );
  void setIsosurfaceValue( const scalar& value );
  void loadTriangleMesh( const std::string& file_name );
  void centerCamera();
  void computeSignedDistanceField();
  void saveScreenshot( const QString& file_name );
  unsigned numEdgeSubSamples() const;
  unsigned numFaceSubSamples() const;
  void setEdgeSubsamples( const unsigned num_samples );
  void setFaceSubsamples( const unsigned num_samples );

  void saveProcessedMesh( const std::string& file_name );

protected:

  void initializeGL();
  void resizeGL( int width, int height );
  void paintGL();
  void mousePressEvent( QMouseEvent* event );
  void mouseReleaseEvent( QMouseEvent* event );
  void mouseMoveEvent( QMouseEvent* event );
  void wheelEvent( QWheelEvent* event );

private slots:

  void enablePerspectiveCamera();
  void enableOrthographicXYCamera();
  void enableOrthographicZYCamera();
  void enableOrthographicZXCamera();

private:

  void computeMarchingCubesMesh();
  void clearSDF();
  void recomputeGrid();
  
  void displayAxes() const;
  void renderInputMesh() const;
  void renderTriangleMesh( const Matrix3Xuc& triangles, const Matrix3Xsc& vertices ) const;
  void renderSignedDistanceField() const;
  void renderGrid() const;
  void renderSurfaceSampling() const;

  bool m_use_perspective_camera;
  PerspectiveCameraController m_camera_controller;
  OrthographicCameraController m_orthographic_controller;

  OpenGL3DSphereRenderer m_sphere_renderer;

  QPoint m_lastPos;
  bool m_draw_axes;

  bool m_draw_sdf;
  bool m_draw_input_mesh;
  bool m_draw_sdf_mesh;
  bool m_draw_surface_sampling;
  bool m_draw_normals;
  bool m_draw_grid;
  bool m_samples_vs_convex;

  // Distance field
  SignedDistanceField m_sdf;
  scalar m_cell_width;
  scalar m_grid_padding;

  // Surface sampling
  SurfaceSampling m_surface_sampling;
  Matrix3Xsc m_sample_normals;

  // Moments of the input mesh
  Moments m_moments;

  // Input mesh
  Matrix3Xuc m_input_triangles;
  Matrix2Xuc m_input_edges;
  Matrix3Xsc m_input_vertices;

  // Vertices of the input's convex hull
  Matrix3Xsc m_convex_hull;

  // Sphere tree of surface samples
  SphereTree m_sphere_tree;
  // Whether or not to draw the sphere tree
  bool m_draw_tree;
  // Level of the sphere tree to draw
  unsigned m_level_to_draw;

  // Marching cubes visualization of the signed distance field
  Matrix3Xuc m_sdf_triangles;
  Matrix3Xsc m_sdf_vertices;
  // Value of the level set to visualize
  scalar m_level_set_value;

  // Max extent of the mesh's bounding box along any direction
  scalar m_mesh_max_bbox_width;

};

#endif
