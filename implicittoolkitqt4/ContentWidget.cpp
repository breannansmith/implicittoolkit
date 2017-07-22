#include "ContentWidget.h"

ContentWidget::ContentWidget( QWidget* parent )
: QWidget( parent )
, m_glWidget( new GLWidget( this ) )
, m_cell_width_spinbox( nullptr )
, m_grid_padding_spinbox( nullptr )
{
  QVBoxLayout* mainLayout = new QVBoxLayout;
  mainLayout->addWidget( m_glWidget );

  // Master grid layout below the OpenGL display
  QGridLayout* master_grid_layout  = new QGridLayout;

  // Leftmost spinboxes
  QGridLayout* left_spinbox_layout = new QGridLayout;

  // Top leftmost spinboxes
  QGridLayout* top_left_spinbox_layout = new QGridLayout;
  QLabel* cell_width_label = new QLabel( tr( "Cell width:" ) );
  m_cell_width_spinbox = new QDoubleSpinBox;
  m_cell_width_spinbox->setRange( 0.0001, 5.0 );
  m_cell_width_spinbox->setValue( m_glWidget->cellWidth() );
  m_cell_width_spinbox->setSingleStep( 0.0001 );
  m_cell_width_spinbox->setDecimals( 4 );
  connect( m_cell_width_spinbox, SIGNAL( valueChanged(double) ), this, SLOT( setCellWidth(double) ) );

  QLabel* grid_padding_label = new QLabel( tr( "Grid padding:" ) );
  m_grid_padding_spinbox = new QDoubleSpinBox;
  m_grid_padding_spinbox->setRange( 0.0, 5.0 );
  m_grid_padding_spinbox->setValue( m_glWidget->gridPadding() );
  m_grid_padding_spinbox->setSingleStep( 0.0001 );
  m_grid_padding_spinbox->setDecimals( 4 );
  connect( m_grid_padding_spinbox, SIGNAL( valueChanged(double) ), this, SLOT( setGridPadding(double) ) );

  top_left_spinbox_layout->setColumnStretch( 0, 0 );
  top_left_spinbox_layout->setColumnStretch( 1, 1 );
  top_left_spinbox_layout->setColumnStretch( 2, 0 );
  top_left_spinbox_layout->setColumnStretch( 3, 1 );
  top_left_spinbox_layout->addWidget( cell_width_label, 0, 0 );
  top_left_spinbox_layout->addWidget( m_cell_width_spinbox, 0, 1 );
  top_left_spinbox_layout->addWidget( grid_padding_label, 0, 2 );
  top_left_spinbox_layout->addWidget( m_grid_padding_spinbox, 0, 3 );

  // Bottom leftmost spinboxes
  QGridLayout* bottom_left_spinbox_layout = new QGridLayout;
  QLabel* level_set_label = new QLabel( tr( "Isosurface:" ) );

  QDoubleSpinBox* level_set_spinbox = new QDoubleSpinBox;
  level_set_spinbox->setRange( -10.0, 10.0 );
  level_set_spinbox->setValue( m_glWidget->isosurfaceValue() );
  level_set_spinbox->setSingleStep( 0.01 );
  connect( level_set_spinbox, SIGNAL( valueChanged(double) ), this, SLOT( setLevelSetValue(double) ) );

  bottom_left_spinbox_layout->setColumnStretch( 0, 0 );
  bottom_left_spinbox_layout->setColumnStretch( 1, 1 );
  bottom_left_spinbox_layout->addWidget( level_set_label, 0, 0 );
  bottom_left_spinbox_layout->addWidget( level_set_spinbox, 0, 1 );

  // Controls for surface sampling
  QGridLayout* sampling_layout = new QGridLayout;

  QLabel* edge_sampling_label = new QLabel( tr( "Edge Samps:" ) );
  QSpinBox* edge_sample_spinbox = new QSpinBox;
  edge_sample_spinbox->setRange( 0, 10 );
  edge_sample_spinbox->setValue( int(m_glWidget->numEdgeSubSamples()) );
  connect( edge_sample_spinbox, SIGNAL( valueChanged(int) ), this, SLOT( setEdgeSubsamples(int) ) );

  QLabel* face_sampling_label = new QLabel( tr( "Face Samps:" ) );
  QSpinBox* face_sample_spinbox = new QSpinBox;
  face_sample_spinbox->setRange( 0, 10 );
  face_sample_spinbox->setValue( int(m_glWidget->numFaceSubSamples()) );
  connect( face_sample_spinbox, SIGNAL( valueChanged(int) ), this, SLOT( setFaceSubsamples(int) ) );

  sampling_layout->setColumnStretch( 0, 0 );
  sampling_layout->setColumnStretch( 1, 1 );
  sampling_layout->setColumnStretch( 2, 0 );
  sampling_layout->setColumnStretch( 3, 1 );
  sampling_layout->addWidget( edge_sampling_label, 0, 0 );
  sampling_layout->addWidget( edge_sample_spinbox, 0, 1 );
  sampling_layout->addWidget( face_sampling_label, 0, 2 );
  sampling_layout->addWidget( face_sample_spinbox, 0, 3 );

  left_spinbox_layout->addLayout( top_left_spinbox_layout, 0, 0 );
  left_spinbox_layout->addLayout( bottom_left_spinbox_layout, 1, 0 );
  left_spinbox_layout->addLayout( sampling_layout, 2, 0 );


  // Rightmost checkboxes
  QGridLayout* right_checkbox_layout = new QGridLayout;
  
  // Controls for selecting/rendering the input mesh
  QCheckBox* input_mesh_checkbox = new QCheckBox( tr( "Display Input Mesh" ) );
  input_mesh_checkbox->setChecked( m_glWidget->drawInputMesh() );
  connect( input_mesh_checkbox, SIGNAL(toggled(bool)), this, SLOT( renderInputMesh(bool) ) );

  QCheckBox* level_set_checkbox = new QCheckBox( tr( "Display Isosurface" ) );
  level_set_checkbox->setChecked( m_glWidget->drawSDFMesh() );
  connect( level_set_checkbox, SIGNAL(toggled(bool)), this, SLOT( renderLevelSet(bool) ) );

  QCheckBox* surface_sampling_checkbox = new QCheckBox( tr( "Display Surface Sampling" ) );
  surface_sampling_checkbox->setChecked( m_glWidget->drawSurfaceSampling() );
  connect( surface_sampling_checkbox, SIGNAL(toggled(bool)), this, SLOT( renderSurfaceSampling(bool) ) );

  right_checkbox_layout->addWidget( input_mesh_checkbox, 0, 0 );
  right_checkbox_layout->addWidget( level_set_checkbox, 1, 0 );
  right_checkbox_layout->addWidget( surface_sampling_checkbox, 2, 0 );



  master_grid_layout->setColumnStretch( 0, 1 );
  master_grid_layout->setColumnStretch( 1, 0 );
  master_grid_layout->addLayout( left_spinbox_layout, 0, 0 );
  master_grid_layout->addLayout( right_checkbox_layout, 0, 1 );

  mainLayout->addLayout( master_grid_layout );

  setLayout(mainLayout);
}

GLWidget* ContentWidget::glWidget()
{
  return m_glWidget;
}

void ContentWidget::toggleGridDisplay()
{
  m_glWidget->toggleDrawGrid();
}

bool ContentWidget::gridDisplayStatus() const
{
  return m_glWidget->drawGrid();
}

void ContentWidget::toggleNormalDisplay()
{
  m_glWidget->toggleDrawNormals();
}

bool ContentWidget::normalDisplayStatus() const
{
  return m_glWidget->drawSampleNormals();
}

void ContentWidget::toggleSamplesVsConvex()
{
  m_glWidget->toggleSamplesVsConvex();
}

bool ContentWidget::samplesVsConvexDisplayStatus() const
{
  return m_glWidget->samplesVsConvexDisplayStatus();
}

void ContentWidget::toggleSphereTreeDisplay()
{
  m_glWidget->toggleSphereTreeDisplay();
}

bool ContentWidget::sphereTreeDisplayStatus() const
{
  return m_glWidget->sphereTreeDisplayStatus();
}

void ContentWidget::descendSphereTree()
{
  return m_glWidget->descendSphereTree();
}

unsigned ContentWidget::nextLowestSphereTreeLevel() const
{
  return m_glWidget->nextLowestSphereTreeLevel();
}

void ContentWidget::ascendSphereTree()
{
  return m_glWidget->ascendSphereTree();
}

unsigned ContentWidget::nextHighestSphereTreeLevel() const
{
  return m_glWidget->nextHighestSphereTreeLevel();
}

void ContentWidget::loadMesh()
{
  const std::string file_name = getOpenFileNameFromUser();
  if( !file_name.empty() )
  {
    m_glWidget->loadTriangleMesh( file_name );
  }
}

void ContentWidget::resetCamera()
{
  m_glWidget->centerCamera();
}

void ContentWidget::computeSDF()
{
  m_glWidget->computeSignedDistanceField();
}

void ContentWidget::exportImage()
{
  assert( m_glWidget != nullptr );
  const std::string file_name = getSaveFileNameFromUser( "Please Specify an Image Name" );
  if( file_name.empty() ) m_glWidget->saveScreenshot( QString::fromStdString( file_name ) );
}

void ContentWidget::setCellWidth( const double& cell_width )
{
  assert( cell_width > 0.0 );
  m_glWidget->setCellWidth( cell_width );
}

void ContentWidget::setGridPadding( const double& grid_padding )
{
  assert( grid_padding >= 0.0 );
  m_glWidget->setGridPadding( grid_padding );
}

void ContentWidget::setLevelSetValue( const double& level_set_value )
{
  m_glWidget->setIsosurfaceValue( level_set_value );
}

void ContentWidget::renderInputMesh( const bool render )
{
  m_glWidget->setInputMeshRendering( render );
}

void ContentWidget::renderLevelSet( const bool render )
{
  m_glWidget->setIsosurfaceRendering( render );
}

void ContentWidget::renderSurfaceSampling( const bool render )
{
  m_glWidget->setSurfaceSamplingRendering( render );
}

void ContentWidget::setEdgeSubsamples( const int num_subsamples )
{
  assert( num_subsamples >= 0 );
  m_glWidget->setEdgeSubsamples( unsigned( num_subsamples ) );
}

void ContentWidget::setFaceSubsamples( const int num_subsamples )
{
  assert( num_subsamples >= 0 );
  m_glWidget->setFaceSubsamples( unsigned(num_subsamples) );
}

void ContentWidget::saveProcessedMesh()
{
  assert( m_glWidget != nullptr );
  const std::string file_name = getSaveFileNameFromUser( "Please an Output File Name" );
  if( file_name.empty() ) return;
  try
  {
    m_glWidget->saveProcessedMesh( file_name );
  }
  catch( const std::string& error )
  {
    std::cerr << "Failed to save processed mesh to file " << file_name << ": " << error << std::endl;
  }
}

std::string ContentWidget::getOpenFileNameFromUser()
{
  QString file_name = QFileDialog::getOpenFileName( this, tr( "Open Triangle Mesh" ) );
  return file_name.toStdString();
}

std::string ContentWidget::getSaveFileNameFromUser( const std::string& display_text )
{
  QString file_name = QFileDialog::getSaveFileName( this, tr( display_text.c_str() ) );
  return file_name.toStdString();
}
