#include "Window.h"

#include "implicitutilities/StringUtilities.h"

Window::Window( QWidget* parent )
: QMainWindow( parent )
, m_content_widget( new ContentWidget )
, m_toggle_grid_display( nullptr )
{
  QMenu* file;
  file = menuBar()->addMenu( tr( "File" ) );
  QMenu* view;
  view = menuBar()->addMenu( tr( "View" ) );
  QMenu* actions;
  actions = menuBar()->addMenu( tr( "Actions" ) );


  // File menu actions
  // Load the input mesh
  QAction* load_mesh = new QAction( tr( "Open Mesh..."), this );
  load_mesh->setShortcut( tr( "Ctrl+o" ) );
  file->addAction( load_mesh );
  connect( load_mesh, SIGNAL( triggered() ), this, SLOT( loadMesh() ) );

  // Save the all data
  QAction* save_processed_mesh = new QAction( tr( "Save Processed Mesh..." ), this );
  save_processed_mesh->setShortcut( tr( "Ctrl+s" ) );
  file->addAction( save_processed_mesh );
  connect( save_processed_mesh, SIGNAL( triggered() ), m_content_widget, SLOT( saveProcessedMesh() ) );

  {
    QAction* seperator = new QAction( this );
    seperator->setSeparator( true );
    file->addAction( seperator );
  }

  // Export an image of the scene
  QAction* export_image = new QAction( tr( "Save Screenshot..." ), this );
  file->addAction( export_image );
  connect( export_image, SIGNAL( triggered() ), m_content_widget, SLOT( exportImage() ) );


  // View menu actions
  QAction* reset_view = new QAction( tr( "Center Camera" ), this );
  connect(reset_view, SIGNAL(triggered()), m_content_widget, SLOT(resetCamera()));
  view->addAction(reset_view);

  // Perspective camera
  QAction* enable_perspective_camera = new QAction( tr( "Perspective Camera" ), this );
  view->addAction( enable_perspective_camera );
  connect( enable_perspective_camera, SIGNAL( triggered() ), m_content_widget->glWidget(), SLOT( enablePerspectiveCamera() ) );
  
  // Orthographic projections
  QAction* enable_xy_orthographic_camera = new QAction( tr( "XY Orthographic Camera" ), this );
  view->addAction( enable_xy_orthographic_camera );
  connect( enable_xy_orthographic_camera, SIGNAL( triggered() ), m_content_widget->glWidget(), SLOT( enableOrthographicXYCamera() ) );
  
  QAction* enable_zy_orthographic_camera = new QAction( tr( "ZY Orthographic Camera" ), this );
  view->addAction( enable_zy_orthographic_camera );
  connect( enable_zy_orthographic_camera, SIGNAL( triggered() ), m_content_widget->glWidget(), SLOT( enableOrthographicZYCamera() ) );
  
  QAction* enable_zx_orthographic_camera = new QAction( tr( "ZX Orthographic Camera" ), this );
  view->addAction( enable_zx_orthographic_camera );
  connect( enable_zx_orthographic_camera, SIGNAL( triggered() ), m_content_widget->glWidget(), SLOT( enableOrthographicZXCamera() ) );

  {
    QAction* seperator = new QAction( this );
    seperator->setSeparator( true );
    view->addAction( seperator );
  }

  // Render option toggles
  m_toggle_grid_display = new QAction( tr( gridDisplayName().c_str() ), this );
  view->addAction( m_toggle_grid_display );
  connect( m_toggle_grid_display, SIGNAL( triggered() ), this, SLOT( toggleGrid() ) );

  m_toggle_normal_display = new QAction( tr( normalDisplayName().c_str() ), this );
  view->addAction( m_toggle_normal_display );
  connect( m_toggle_normal_display, SIGNAL( triggered() ), this, SLOT( toggleNormals() ) );

  m_toggle_samples_vs_convex_samples = new QAction( tr( samplesVsConvexDisplayName().c_str() ), this );
  view->addAction( m_toggle_samples_vs_convex_samples );
  connect( m_toggle_samples_vs_convex_samples, SIGNAL( triggered() ), this, SLOT( toggleSamplesVsConvex() ) );

  m_toggle_sphere_tree_display = new QAction( tr( toggleSphereTreeDisplayName().c_str() ), this );
  view->addAction( m_toggle_sphere_tree_display );
  connect( m_toggle_sphere_tree_display, SIGNAL( triggered() ), this, SLOT( toggleSphereTree() ) );

  m_descend_sphere_tree_display = new QAction( tr( descendSphereTreeDisplayName().c_str() ), this );
  view->addAction( m_descend_sphere_tree_display );
  connect( m_descend_sphere_tree_display, SIGNAL( triggered() ), this, SLOT( descendSphereTree() ) );

  m_ascend_sphere_tree_display = new QAction( tr( ascendSphereTreeDisplayName().c_str() ), this );
  view->addAction( m_ascend_sphere_tree_display );
  connect( m_ascend_sphere_tree_display, SIGNAL( triggered() ), this, SLOT( ascendSphereTree() ) );

  // Action menu actions
  QAction* compute_sdf = new QAction( tr( "Compute SDF" ), this );
  connect( compute_sdf, SIGNAL( triggered() ), m_content_widget, SLOT( computeSDF() ) );
  actions->addAction( compute_sdf );

  setCentralWidget( m_content_widget );
}

void Window::loadMesh()
{
  m_content_widget->loadMesh();
  m_descend_sphere_tree_display->setText( tr( descendSphereTreeDisplayName().c_str() ) );
  m_ascend_sphere_tree_display->setText( tr( ascendSphereTreeDisplayName().c_str() ) );
}

std::string Window::gridDisplayName() const
{
  if( m_content_widget->gridDisplayStatus() )
  {
    return "Hide Grid";
  }
  return "Show Grid";
}

std::string Window::normalDisplayName() const
{
  if( m_content_widget->normalDisplayStatus() )
  {
    return "Hide Normals";
  }
  return "Show Normals";
}

std::string Window::samplesVsConvexDisplayName() const
{
  if( m_content_widget->samplesVsConvexDisplayStatus() )
  {
    return "Show Convex Hull Samples";
  }
  return "Show Surface Samples";
}

std::string Window::toggleSphereTreeDisplayName() const
{
  if( m_content_widget->sphereTreeDisplayStatus() )
  {
    return "Hide Sphere Tree";
  }
  return "Show Sphere Tree";
}

std::string Window::descendSphereTreeDisplayName() const
{
  return "Descend to Sphere Tree Level " + StringUtilities::convertToString( m_content_widget->nextLowestSphereTreeLevel() );;
}

std::string Window::ascendSphereTreeDisplayName() const
{
  return "Ascend to Sphere Tree Level " + StringUtilities::convertToString( m_content_widget->nextHighestSphereTreeLevel() );;
}

void Window::toggleGrid()
{
  m_content_widget->toggleGridDisplay();
  m_toggle_grid_display->setText( tr( gridDisplayName().c_str() ) );
}

void Window::toggleNormals()
{
  m_content_widget->toggleNormalDisplay();
  m_toggle_normal_display->setText( tr( normalDisplayName().c_str() ) );
}

void Window::toggleSamplesVsConvex()
{
  m_content_widget->toggleSamplesVsConvex();
  m_toggle_samples_vs_convex_samples->setText( tr( samplesVsConvexDisplayName().c_str() ) );
}

void Window::toggleSphereTree()
{
  m_content_widget->toggleSphereTreeDisplay();
  m_toggle_sphere_tree_display->setText( tr( toggleSphereTreeDisplayName().c_str() ) );
}

void Window::descendSphereTree()
{
  m_content_widget->descendSphereTree();
  m_descend_sphere_tree_display->setText( tr( descendSphereTreeDisplayName().c_str() ) );
  m_ascend_sphere_tree_display->setText( tr( ascendSphereTreeDisplayName().c_str() ) );
}

void Window::ascendSphereTree()
{
  m_content_widget->ascendSphereTree();
  m_descend_sphere_tree_display->setText( tr( descendSphereTreeDisplayName().c_str() ) );
  m_ascend_sphere_tree_display->setText( tr( ascendSphereTreeDisplayName().c_str() ) );
}
