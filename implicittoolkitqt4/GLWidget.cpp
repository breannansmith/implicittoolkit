#include "GLWidget.h"

#include "implicitutilities/MeshTools.h"

#include "implicitutilities/ProcessedMeshWriter.h"

#include "implicitutilities/Qhull.h"

class Qt4PRogressWrapper final : public AlgorithmProgressTracker
{

public:

  Qt4PRogressWrapper( QProgressDialog& progress_dialog )
  : m_progress_dialog(progress_dialog)
  {}

  virtual void setProgressIndicator( int value )
  {
    assert( value >= 0 ); assert( value <= 100 );
    m_progress_dialog.setValue(value);
  }

  virtual bool algorithmCancelled()
  {
    return m_progress_dialog.wasCanceled();
  }

private:

  QProgressDialog& m_progress_dialog;

};

#ifndef NDEBUG
static std::string glErrorToString( const GLenum error_code )
{
  switch( error_code )
  {
    case GL_NO_ERROR:
      return "GL_NO_ERROR";
    case GL_INVALID_ENUM:
      return "GL_INVALID_ENUM";
    case GL_INVALID_VALUE:
      return "GL_INVALID_VALUE";
    case GL_INVALID_OPERATION:
      return "GL_INVALID_OPERATION";
    case GL_INVALID_FRAMEBUFFER_OPERATION:
      return "GL_INVALID_FRAMEBUFFER_OPERATION";
    case GL_OUT_OF_MEMORY:
      return "GL_OUT_OF_MEMORY";
    case GL_STACK_UNDERFLOW:
      return "GL_STACK_UNDERFLOW";
    case GL_STACK_OVERFLOW:
      return "GL_STACK_OVERFLOW";
    default:
      return "Unknown error. Please contact the maintainer of this code.";
  }
}

static bool checkGLErrors()
{
  const GLenum error_code = glGetError();
  if( error_code != GL_NO_ERROR )
  {
    std::cerr << "OpenGL error: " << glErrorToString( error_code ) << std::endl;
    return false;
  }
  return true;
}
#endif

GLWidget::GLWidget( QWidget* parent )
: QGLWidget( QGLFormat( QGL::SampleBuffers ), parent )
, m_use_perspective_camera( true )
, m_camera_controller()
, m_orthographic_controller()
, m_sphere_renderer( 3 )
, m_lastPos( 0, 0 )
, m_draw_axes( false )
, m_draw_sdf( false )
, m_draw_input_mesh( true )
, m_draw_sdf_mesh( true )
, m_draw_surface_sampling( false )
, m_draw_normals( false )
, m_draw_grid( true )
, m_samples_vs_convex( true )
, m_sdf()
, m_cell_width( 1.0 )
, m_grid_padding( 0.1 * m_cell_width )
, m_surface_sampling()
, m_moments()
, m_input_triangles()
, m_input_edges()
, m_input_vertices()
, m_convex_hull()
, m_sphere_tree()
, m_draw_tree( false )
, m_level_to_draw( 0 )
, m_sdf_triangles()
, m_sdf_vertices()
, m_level_set_value( 0.0 )
, m_mesh_max_bbox_width()
{
  assert( m_cell_width > 0.0 );
  assert( m_grid_padding >= 0.0 );
}

QSize GLWidget::minimumSizeHint() const
{
  return QSize( 50, 50 );
}

QSize GLWidget::sizeHint() const
{
  return QSize( 512, 512 );
}

const scalar& GLWidget::cellWidth() const
{
  return m_cell_width;
}

const scalar& GLWidget::gridPadding() const
{
  return m_grid_padding;
}

const scalar& GLWidget::isosurfaceValue() const
{
  return m_level_set_value;
}

bool GLWidget::drawInputMesh() const
{
  return m_draw_input_mesh;
}

bool GLWidget::drawSDFMesh() const
{
  return m_draw_sdf_mesh;
}

bool GLWidget::drawSurfaceSampling() const
{
  return m_draw_surface_sampling;
}

bool GLWidget::drawSampleNormals() const
{
  return m_draw_normals;
}

void GLWidget::toggleDrawNormals()
{
  m_draw_normals = !m_draw_normals;
  updateGL();
}

bool GLWidget::drawGrid() const
{
  return m_draw_grid;
}

void GLWidget::toggleDrawGrid()
{
  m_draw_grid = !m_draw_grid;
  updateGL();
}

bool GLWidget::samplesVsConvexDisplayStatus() const
{
  return m_samples_vs_convex;
}

void GLWidget::toggleSamplesVsConvex()
{
  m_samples_vs_convex = !m_samples_vs_convex;
  updateGL();
}

bool GLWidget::sphereTreeDisplayStatus() const
{
  return m_draw_tree;
}

void GLWidget::toggleSphereTreeDisplay()
{
  m_draw_tree = !m_draw_tree;
  updateGL();
}

unsigned GLWidget::nextLowestSphereTreeLevel() const
{
  return std::min( m_level_to_draw + 1, m_sphere_tree.depth() );
}

void GLWidget::descendSphereTree()
{
  if( m_level_to_draw + 1 <= m_sphere_tree.depth() )
  {
    m_level_to_draw += 1;
    updateGL();
  }
}

unsigned GLWidget::nextHighestSphereTreeLevel() const
{
  if( m_level_to_draw == 0 ) return 0;
  return m_level_to_draw - 1;
}

void GLWidget::ascendSphereTree()
{
  if( m_level_to_draw >= 1 )
  {
    m_level_to_draw -= 1;
    updateGL();
  }
}

void GLWidget::setCellWidth( const scalar& cell_width )
{
  assert( cell_width > 0.0 );
  m_cell_width = cell_width;
  clearSDF();
  recomputeGrid();
  updateGL();
}

void GLWidget::setGridPadding( const scalar& grid_padding )
{
  assert( grid_padding >= 0.0 );
  m_grid_padding = grid_padding;
  clearSDF();
  recomputeGrid();
  updateGL();
}

void GLWidget::clearSDF()
{
  m_sdf_triangles.resize( 3, 0 );
  m_sdf_vertices.resize( 3, 0 );
  m_sdf.clear();
  m_sample_normals.setZero();
}

void GLWidget::setInputMeshRendering( const bool render )
{
  m_draw_input_mesh = render;
  updateGL();
}

void GLWidget::setIsosurfaceRendering( const bool render )
{
  m_draw_sdf_mesh = render;
  updateGL();
}

void GLWidget::setSurfaceSamplingRendering( const bool render )
{
  m_draw_surface_sampling = render;
  updateGL();
}

void GLWidget::setIsosurfaceValue( const scalar& value )
{
  m_level_set_value = value;
  computeMarchingCubesMesh();
}

void GLWidget::loadTriangleMesh( const std::string& file_name )
{
  {
    // Load the user requested file
    Matrix3Xsc new_vertices;
    Matrix2Xuc new_edges;
    Matrix3Xuc new_triangles;
    try
    {
      ObjParser::loadMesh( file_name, new_vertices, new_edges, new_triangles );
    }
    catch( const std::string& error )
    {
      std::cerr << "Failed to load " << file_name << ": " << error << std::endl;
      return;
    }

    // Check the Euler characteristic to verify that we have a closed volumetric mesh
    //{
    //  const int euler_char = MeshTools::eulerCharacteristic( new_vertices, new_edges, new_triangles );
    //  if( euler_char != 2 )
    //  {
    //    std::cerr << "Euler characteristic incorect for " << file_name << " with a value of " << euler_char << std::endl;
    //    return;
    //  }
    //}

    // Compute moments
    Moments moments( new_vertices, new_triangles );
    // Translate and rotate the mesh to align the principal axes with the Cartesian axes at the origin
    new_vertices = moments.R().transpose() * ( new_vertices.colwise() - moments.x() );

    // If everything checked out, save the mesh
    m_input_vertices.swap( new_vertices );
    m_input_edges.swap( new_edges );
    m_input_triangles.swap( new_triangles );
    // ... and save the moments
    m_moments = std::move( moments );

    // Compute the vertices that compose the convex hull of the input vertices
    Qhull::computeConvexHull( m_input_vertices, m_convex_hull );

    const Eigen::Matrix<scalar,6,1> bbox = MeshTools::computeBoundingBox( m_convex_hull );
    m_mesh_max_bbox_width = ( bbox.segment<3>( 3 ) - bbox.segment<3>( 0 ) ).maxCoeff();
  }

  // Compute the grid for the implicit function
  clearSDF();
  recomputeGrid();

  // Compute the initial surface sampling
  m_surface_sampling.computeSurfaceSampling( m_input_vertices, m_input_edges, m_input_triangles );
  m_sphere_tree.buildTree( m_surface_sampling.samples(), m_sdf.gridStart(), m_sdf.gridEnd() );
  m_level_to_draw = 0;
  m_sdf.evaluateGradients( m_surface_sampling.samples(), m_sample_normals );

  // Center the camera at the object
  centerCamera();
}

void GLWidget::saveProcessedMesh( const std::string& file_name )
{
  // TODO: If SDF not computed, emit error
  ProcessedMeshWriter::write( file_name, m_sdf, m_surface_sampling, m_moments, m_input_vertices, m_input_edges, m_input_triangles, m_convex_hull, m_cell_width, m_grid_padding, m_surface_sampling.numEdgeSubSamples(), m_surface_sampling.numFaceSubSamples(), m_sphere_tree );
}

void GLWidget::centerCamera()
{
  const Vector3s start = m_sdf.gridStart();
  const Vector3s end = ( m_sdf.gridStart().array() + m_sdf.cellDelta().array() * ( m_sdf.gridDimensions().array() - 1 ).cast<scalar>() ).matrix();
  const Vector3s center = 0.5 * ( start + end );
  const scalar radius = 1.05 * ( end - center ).norm();

  if( m_use_perspective_camera )
  {
    m_camera_controller.centerCameraAtSphere( center, radius );
  }
  else
  {
    GLint viewport[4];
    glGetIntegerv( GL_VIEWPORT, viewport );
    const GLint width = viewport[2];
    const GLint height = viewport[3];

    const GLdouble ratio = GLdouble( height ) / GLdouble( width );
    const GLdouble scale = std::max( ratio * radius, radius );

    m_orthographic_controller.setCenterAndScale( center, scale );
  }
  
  updateGL();
}

void GLWidget::enablePerspectiveCamera()
{
  m_use_perspective_camera = true;
  GLint viewport[4];
  glGetIntegerv( GL_VIEWPORT, viewport );
  const GLint width = viewport[2];
  const GLint height = viewport[3];
  m_camera_controller.setPerspective( width, height );
  centerCamera();
}

void GLWidget::enableOrthographicXYCamera()
{
  m_use_perspective_camera = false;
  m_orthographic_controller.useXYView();
  GLint viewport[4];
  glGetIntegerv( GL_VIEWPORT, viewport );
  const GLint width = viewport[2];
  const GLint height = viewport[3];
  m_orthographic_controller.setPerspective( width, height );
  centerCamera();
}

void GLWidget::enableOrthographicZYCamera()
{
  m_use_perspective_camera = false;
  m_orthographic_controller.useZYView();
  GLint viewport[4];
  glGetIntegerv( GL_VIEWPORT, viewport );
  const GLint width = viewport[2];
  const GLint height = viewport[3];
  m_orthographic_controller.setPerspective( width, height );
  centerCamera();
}

void GLWidget::enableOrthographicZXCamera()
{
  m_use_perspective_camera = false;
  m_orthographic_controller.useZXView();
  GLint viewport[4];
  glGetIntegerv( GL_VIEWPORT, viewport );
  const GLint width = viewport[2];
  const GLint height = viewport[3];
  m_orthographic_controller.setPerspective( width, height );
  centerCamera();
}

void GLWidget::saveScreenshot( const QString& file_name )
{
  QImage frame_buffer = grabFrameBuffer();
  frame_buffer.save( file_name );
}

unsigned GLWidget::numEdgeSubSamples() const
{
  return m_surface_sampling.numEdgeSubSamples();
}

unsigned GLWidget::numFaceSubSamples() const
{
  return m_surface_sampling.numFaceSubSamples();
}

void GLWidget::setEdgeSubsamples( const unsigned num_samples )
{
  m_surface_sampling.setEdgeSubsamples( num_samples );
  m_surface_sampling.computeSurfaceSampling( m_input_vertices, m_input_edges, m_input_triangles );
  m_sphere_tree.buildTree( m_surface_sampling.samples(), m_sdf.gridStart(), m_sdf.gridEnd() );
  m_level_to_draw = 0;
  m_sdf.evaluateGradients( m_surface_sampling.samples(), m_sample_normals );
  updateGL();
}

void GLWidget::setFaceSubsamples( const unsigned num_samples )
{
  m_surface_sampling.setFaceSubsamples( num_samples );
  m_surface_sampling.computeSurfaceSampling( m_input_vertices, m_input_edges, m_input_triangles );
  m_sphere_tree.buildTree( m_surface_sampling.samples(), m_sdf.gridStart(), m_sdf.gridEnd() );
  m_level_to_draw = 0;
  m_sdf.evaluateGradients( m_surface_sampling.samples(), m_sample_normals );
  updateGL();
}

void GLWidget::computeSignedDistanceField()
{
  QProgressDialog progress( "Computing Distance...", "Abort Computation", 0, 100, this );
  Qt4PRogressWrapper progress_wrapper( progress );

  progress.setWindowModality( Qt::WindowModal );

  progress.setValue( 0 );

  DistanceTools::computeSignedDistanceField( m_sdf.gridStart(), m_sdf.gridDimensions(), m_sdf.cellDelta(), m_input_vertices, m_input_triangles, m_sdf.vals(), &progress_wrapper );

  progress.setValue( 100 );

  if( progress.wasCanceled() ) { m_sdf.clearValues(); }

  m_sdf.evaluateGradients( m_surface_sampling.samples(), m_sample_normals );
  assert( m_sample_normals.rows() == m_surface_sampling.samples().rows() );
  assert( m_sample_normals.cols() == m_surface_sampling.samples().cols() );

  computeMarchingCubesMesh();
}

void GLWidget::computeMarchingCubesMesh()
{
  // Extract a mesh from the signed distance field
  MarchingCubes mc( int(m_sdf.nx()), int(m_sdf.ny()), int(m_sdf.nz()) );
  // Use the original marching cubes algorithm
  mc.set_method( true );
  // Use the points we loaded
  mc.set_ext_data( m_sdf.vals().data() );
  // Initialize some storage before the solve
  mc.init_all();
  // Run marching cubes
  mc.run( m_level_set_value );

  // Allocate space for the mesh
  m_sdf_triangles.resize( 3, mc.ntrigs() );
  m_sdf_vertices.resize( 3, mc.nverts() );

  // Extract the vertices
  for( int i = 0; i < mc.nverts(); ++i )
  {
    const Vertex* cur_vert = mc.vert( i );
    assert( cur_vert != nullptr );
    m_sdf_vertices.col( i ) = Vector3s( cur_vert->x, cur_vert->y, cur_vert->z );
    // Rescale to be same size as mesh and re-center
    m_sdf_vertices.col( i ) = m_sdf_vertices.col( i ).array() * m_sdf.cellDelta().array() + m_sdf.gridStart().array();
  }

  // Extract the triangles
  for( int i = 0; i < mc.ntrigs(); ++i )
  {
    const Triangle* cur_tri = mc.trig(i);
    assert( cur_tri != nullptr );
    m_sdf_triangles.col( i ) << unsigned(cur_tri->v1), unsigned(cur_tri->v2), unsigned(cur_tri->v3);
  }

  // Clean up
  mc.clean_all();

  updateGL();
}

void GLWidget::recomputeGrid()
{
  assert( ( m_input_vertices.cols() == 0 && m_input_triangles.cols() == 0 ) || ( m_input_vertices.cols() != 0 && m_input_triangles.cols() != 0 ) );

  if( m_input_vertices.cols() == 0 ) return;

  const Eigen::Matrix<scalar,6,1> bounding_box = MeshTools::computeBoundingBox( m_input_vertices );

  m_sdf.initializeCenteredGrid( bounding_box, Vector3s::Constant( cellWidth() ), Vector3s::Constant( gridPadding() ) );
}

void GLWidget::initializeGL()
{
  glEnable( GL_DEPTH_TEST );
  qglClearColor( QColor(255,255,255,255) );
  glShadeModel( GL_SMOOTH );
  GLfloat global_ambient[] = { 0.5f, 0.5f, 0.5f, 1.0f };
  glLightModelfv( GL_LIGHT_MODEL_AMBIENT, global_ambient );
  GLfloat diffuse[] = { 1.0f, 1.0f, 1.0f , 1.0f };
  glLightfv( GL_LIGHT0, GL_DIFFUSE, diffuse );
  glEnable( GL_LIGHTING );
  glEnable( GL_LIGHT0 );
  glEnable( GL_NORMALIZE ); // Prpobably slow... but fine for now

  //m_sphere_renderer = new OpenGL3DSphereRenderer( 1.0, 4 );

  assert( checkGLErrors() );
}

void GLWidget::resizeGL( int width, int height )
{
  assert( width >= 0 ); assert( height >= 0 );

  if( m_use_perspective_camera )
  {
    m_camera_controller.setPerspective( width, height );
  }
  else
  {
    m_orthographic_controller.setPerspective( width, height );
  }

  assert( checkGLErrors() );
}

void GLWidget::paintGL()
{
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity();

  if( m_use_perspective_camera )
  {
    m_camera_controller.positionCamera();
  }
  else
  {
    m_orthographic_controller.positionCamera();
  }

  // Render the x, y, and z axis in the scene
  if( m_draw_axes ) { displayAxes(); }

  // Render the signed distance field
  if( m_draw_sdf ) { renderSignedDistanceField(); }

  // Render the surface sampling
  if( m_draw_surface_sampling ) { renderSurfaceSampling(); }

  // Render the input mesh
  if( drawInputMesh() )
  {
    GLfloat mcolorambient[] = { 0.2f, 0.2f, 0.2f, 1.0f };
    glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT, mcolorambient );
    GLfloat mcolordiffuse[] = { 0.0f, 0.8f, 0.0f, 1.0f };
    glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, mcolordiffuse );

    renderTriangleMesh( m_input_triangles, m_input_vertices );
  }

  // Render the sdf visualization mesh
  if( drawSDFMesh() )
  {
    GLfloat mcolorambient[] = { 0.2f, 0.2f, 0.2f, 1.0f };
    glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT, mcolorambient );
    GLfloat mcolordiffuse[] = { 0.8f, 0.0f, 0.0f, 1.0f };
    glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, mcolordiffuse );

    renderTriangleMesh( m_sdf_triangles, m_sdf_vertices );
  }

  if( m_draw_grid )
  {
    renderGrid();
  }

  //renderVolumeSampledVertices();

  assert( checkGLErrors() );
}

void GLWidget::mousePressEvent( QMouseEvent* event )
{
  if( event->buttons() & Qt::LeftButton )
  {
    m_draw_axes = true;
    updateGL();
  }
  m_lastPos = event->pos();
}

void GLWidget::mouseReleaseEvent( QMouseEvent* event )
{
  if( !( event->buttons() & Qt::LeftButton ) )
  {
    m_draw_axes = false;
    updateGL();
  }
}

void GLWidget::mouseMoveEvent( QMouseEvent* event )
{
  const int dx = event->x() - m_lastPos.x();
  const int dy = event->y() - m_lastPos.y();
  m_lastPos = event->pos();

  if( m_use_perspective_camera )
  {
    if( event->buttons() & Qt::LeftButton )
    {
      const GLdouble horizontal_amount = 0.01 * GLdouble( dx );
      const GLdouble vertical_amount =   0.01 * GLdouble( dy );
      m_camera_controller.addToZenithAngle( -vertical_amount );
      m_camera_controller.addToAzimuthAngle( -horizontal_amount );
      updateGL();
    }
    else if( event->buttons() & Qt::MidButton )
    {
      const GLdouble horizontal_amount = 0.002 * GLdouble( dx );
      const GLdouble vertical_amount =   0.002 * GLdouble( dy );
      m_camera_controller.trackCameraHorizontal( -m_camera_controller.getDistFromCenter()*horizontal_amount );
      m_camera_controller.trackCameraVertical( m_camera_controller.getDistFromCenter()*vertical_amount  );
      updateGL();
    }
    else if( event->buttons() & Qt::RightButton )
    {
      const GLdouble horizontal_amount = 0.01 * GLdouble( dx );
      const GLdouble vertical_amount =   0.01 * GLdouble( dy );
      m_camera_controller.addToDistFromCenter( m_camera_controller.getDistFromCenter()*(horizontal_amount-vertical_amount) );
      updateGL();
    }
  }
  else
  {
    if( event->buttons() & Qt::LeftButton )
    {
      m_orthographic_controller.translateView( dx, dy );
      updateGL();
    }
    else if( event->buttons() & Qt::RightButton )
    {
      const GLdouble horizontal_amount = GLdouble( dx );
      const GLdouble vertical_amount =   GLdouble( dy );
      m_orthographic_controller.addToDistFromCenter( 0.1 * ( horizontal_amount - vertical_amount ) );
      updateGL();
    }
  }

  assert( checkGLErrors() );
}

void GLWidget::wheelEvent( QWheelEvent* event )
{
  if( m_use_perspective_camera )
  {
    m_camera_controller.addToDistFromCenter( -0.002 * m_camera_controller.getDistFromCenter() * event->delta() );
  }
  else
  {
    m_orthographic_controller.addToDistFromCenter( -0.05 * event->delta() );
  }

  updateGL();

  assert( checkGLErrors() );
}

void GLWidget::displayAxes() const
{
  glPushAttrib(GL_COLOR);
  glPushAttrib(GL_LIGHTING);
  glPushAttrib(GL_LINE_WIDTH);

  glDisable( GL_LIGHTING );

  glLineWidth(2.0);

  // Draw the positive x, y, and z axis
  glColor3d(1.0,0.0,0.0);
  glBegin(GL_LINES);
  glVertex4d(0.0,0.0,0.0,1.0);
  glVertex4d(1.0,0.0,0.0,0.0);
  glEnd();
  glColor3d(0.0,1.0,0.0);
  glBegin(GL_LINES);
  glVertex4d(0.0,0.0,0.0,1.0);
  glVertex4d(0.0,1.0,0.0,0.0);
  glEnd();
  glColor3d(0.0,0.0,1.0);
  glBegin(GL_LINES);
  glVertex4d(0.0,0.0,0.0,1.0);
  glVertex4d(0.0,0.0,1.0,0.0);
  glEnd();

  // Draw the negative x, y, and z axis
  glLineStipple(1,0x00FF);
  glEnable(GL_LINE_STIPPLE);
  glColor3d(1.0,0.0,0.0);
  glBegin(GL_LINES);
  glVertex4d(-1.0,0.0,0.0,0.0);
  glVertex4d(0.0,0.0,0.0,1.0);
  glEnd();
  glColor3d(0.0,1.0,0.0);
  glBegin(GL_LINES);
  glVertex4d(0.0,-1.0,0.0,0.0);
  glVertex4d(0.0,0.0,0.0,1.0);
  glEnd();
  glColor3d(0.0,0.0,1.0);
  glBegin(GL_LINES);
  glVertex4d(0.0,0.0,-1.0,0.0);
  glVertex4d(0.0,0.0,0.0,1.0);
  glEnd();
  glDisable(GL_LINE_STIPPLE);

  glPopAttrib();
  glPopAttrib();
  glPopAttrib();
}

void GLWidget::renderInputMesh() const
{
  renderTriangleMesh( m_input_triangles, m_input_vertices );
}

// TODO: Push mesh rendering code into display lists
void GLWidget::renderTriangleMesh( const Matrix3Xuc& triangles, const Matrix3Xsc& vertices ) const
{
  glBegin( GL_TRIANGLES );
  for( int tri_num = 0; tri_num < triangles.cols(); ++tri_num )
  {
    assert( triangles( 0, tri_num ) < unsigned( vertices.cols() ) );
    const Eigen::Matrix<GLdouble,3,1> v0 = vertices.col( triangles( 0, tri_num ) ).cast<GLdouble>();
    assert( triangles( 1, tri_num ) < unsigned( vertices.cols() ) );
    const Eigen::Matrix<GLdouble,3,1> v1 = vertices.col( triangles( 1, tri_num ) ).cast<GLdouble>();
    assert( triangles( 2, tri_num ) < unsigned( vertices.cols() ) );
    const Eigen::Matrix<GLdouble,3,1> v2 = vertices.col( triangles( 2, tri_num ) ).cast<GLdouble>();
    assert( fabs( ( v1 - v0 ).cross( v2 - v0 ).norm() ) != 0.0 );
    const Eigen::Matrix<GLdouble,3,1> n = ( v1 - v0 ).cross( v2 - v0 ).normalized();
    assert( fabs( n.norm() - 1.0 ) <= 1.0e-6 );
    glNormal3d( n.x(), n.y(), n.z() );
    glVertex3d( v0.x(), v0.y(), v0.z() );
    glVertex3d( v1.x(), v1.y(), v1.z() );
    glVertex3d( v2.x(), v2.y(), v2.z() );
  }
  glEnd();
}

void GLWidget::renderSignedDistanceField() const
{
  if( m_sdf.numVals() == 0 ) return;

  glPushAttrib( GL_POINT_SIZE );
  glPushAttrib( GL_LIGHTING );

  glPointSize(2.0);
  glDisable( GL_LIGHTING );

  glBegin(GL_POINTS);
  for( unsigned k = 0; k < m_sdf.nz(); ++k ) for( unsigned j = 0; j < m_sdf.ny(); ++j ) for( unsigned i = 0; i < m_sdf.nx(); ++i )
  {
    if( m_sdf.getValue( i, j, k ) < 0 )
    {
      glColor3d( 1.0, 0.0, 0.0 );
    }
    else
    {
      continue;
      //glColor3d(0.0,0.0,1.0);
    }

    const Vector3s grid_point = m_sdf.gridStart() + Vector3s( i * m_sdf.dx(), j * m_sdf.dy(), k * m_sdf.dz() );
    glVertex3d( grid_point.x(), grid_point.y(), grid_point.z() );
  }
  glEnd();

  glPopAttrib();
  glPopAttrib();
}

void GLWidget::renderGrid() const
{
  glPushAttrib( GL_LIGHTING );
  glPushAttrib( GL_COLOR );
  glPushAttrib( GL_LINE_WIDTH );

  glDisable( GL_LIGHTING );

  glColor3d( 0.0, 0.0, 0.0 );
  glLineWidth( 1.0 );
  
  glBegin(GL_LINES);

  for( unsigned j = 0; j < m_sdf.ny(); ++j )
  {
    const Vector3s line_start = m_sdf.gridStart() + j * m_sdf.dy() * Vector3s::UnitY();
    const Vector3s line_end = m_sdf.gridStart() + j * m_sdf.dy() * Vector3s::UnitY() + ( m_sdf.nz() - 1 ) * m_sdf.dz() * Vector3s::UnitZ();
    glVertex3d(line_start.x(),line_start.y(),line_start.z());
    glVertex3d(line_end.x(),line_end.y(),line_end.z());
  }

  for( unsigned k = 0; k < m_sdf.nz(); ++k )
  {
    const Vector3s line_start = m_sdf.gridStart() + k * m_sdf.dz() * Vector3s::UnitZ();
    const Vector3s line_end = m_sdf.gridStart() + k * m_sdf.dz() * Vector3s::UnitZ() + ( m_sdf.ny() - 1 ) * m_sdf.dy() * Vector3s::UnitY();
    glVertex3d( line_start.x(), line_start.y(), line_start.z() );
    glVertex3d( line_end.x(), line_end.y(), line_end.z() );
  }

  for( unsigned i = 0; i < m_sdf.nx(); ++i )
  {
    const Vector3s line_start = m_sdf.gridStart() + i * m_sdf.dx() * Vector3s::UnitX();
    const Vector3s line_end = m_sdf.gridStart() + i * m_sdf.dx() * Vector3s::UnitX() + ( m_sdf.nz() - 1 ) * m_sdf.dz() * Vector3s::UnitZ();
    glVertex3d( line_start.x(), line_start.y(), line_start.z() );
    glVertex3d( line_end.x(), line_end.y(), line_end.z() );
  }

  for( unsigned k = 0; k < m_sdf.nz(); ++k )
  {
    const Vector3s line_start = m_sdf.gridStart() + k * m_sdf.dz()*Vector3s::UnitZ();
    const Vector3s line_end = m_sdf.gridStart() + k * m_sdf.dz()*Vector3s::UnitZ() + ( m_sdf.nx() - 1 ) * m_sdf.dx() * Vector3s::UnitX();
    glVertex3d( line_start.x(), line_start.y(), line_start.z() );
    glVertex3d( line_end.x(), line_end.y(), line_end.z() );
  }

  for( unsigned i = 0; i < m_sdf.nx(); ++i )
  {
    const Vector3s line_start = m_sdf.gridStart() + i * m_sdf.dx()*Vector3s::UnitX();
    const Vector3s line_end = m_sdf.gridStart() + i * m_sdf.dx() * Vector3s::UnitX() + ( m_sdf.ny() - 1 ) * m_sdf.dy() * Vector3s::UnitY();
    glVertex3d( line_start.x(), line_start.y(), line_start.z() );
    glVertex3d( line_end.x(), line_end.y(), line_end.z() );
  }

  for( unsigned j = 0; j < m_sdf.ny(); ++j )
  {
    const Vector3s line_start = m_sdf.gridStart() + j * m_sdf.dy() * Vector3s::UnitY();
    const Vector3s line_end = m_sdf.gridStart() + j * m_sdf.dy() * Vector3s::UnitY() + ( m_sdf.nx() - 1 ) * m_sdf.dx() * Vector3s::UnitX();
    glVertex3d( line_start.x(), line_start.y(), line_start.z() );
    glVertex3d( line_end.x(), line_end.y(), line_end.z() );
  }

  glEnd();

  glPopAttrib();
  glPopAttrib();
  glPopAttrib();
}

void GLWidget::renderSurfaceSampling() const
{
  if( m_samples_vs_convex && m_surface_sampling.empty() ) { return; }
  if( !m_samples_vs_convex && m_convex_hull.cols() == 0 ) { return; }

  const scalar grid_width = ( m_sdf.gridEnd() - m_sdf.gridStart() ).mean();

  // Draw a sphere for each sample
  const Matrix3Xsc& samples = m_samples_vs_convex ? m_surface_sampling.samples() : m_convex_hull;
  const Eigen::Matrix<GLfloat,3,1> primary_color( 0.0, 0.0, 0.5 );
  for( unsigned sample_num = 0; sample_num < samples.cols(); ++sample_num )
  {
    glPushMatrix();
    glTranslated( GLdouble( samples( 0, sample_num ) ), GLdouble( samples( 1, sample_num ) ), GLdouble( samples( 2, sample_num ) ) );
    // TODO: Might work better to scale by min edge length here
    glScaled( 0.01 * m_mesh_max_bbox_width, 0.01 * m_mesh_max_bbox_width, 0.01 * m_mesh_max_bbox_width );
    m_sphere_renderer.drawVertexArray( primary_color, primary_color );
    glPopMatrix();
  }

  // Draw the sphere tree
  if( m_draw_tree )
  {
    std::vector<Vector3s> centers;
    std::vector<scalar> radii;
    assert( m_level_to_draw <= m_sphere_tree.depth() );
    m_sphere_tree.spheres( m_level_to_draw, centers, radii );
    glPushAttrib( GL_LIGHTING );
    glPushAttrib( GL_COLOR );
    glPushAttrib( GL_POINT_SIZE );
    glDisable( GL_LIGHTING );
    glColor3d( 0.0, 0.0, 0.0 );
    glPointSize( 5.0 );
    for( std::vector<Vector3s>::size_type node_idx = 0; node_idx < centers.size(); ++node_idx )
    {
      glPushMatrix();
      glTranslated( centers[node_idx].x(), centers[node_idx].y(), centers[node_idx].z() );
      glScaled( radii[node_idx], radii[node_idx], radii[node_idx] );
      m_sphere_renderer.drawVertices();
      glPopMatrix();
    }
    glPopAttrib();
    glPopAttrib();
    glPopAttrib();
  }

  if( !m_samples_vs_convex && m_draw_normals ) { return; }
  if( m_sdf.empty() || !m_draw_normals ) { return; }

  assert( m_sample_normals.rows() == samples.rows() );
  assert( m_sample_normals.cols() == samples.cols() );

  // Draw a line for each normal
  glPushAttrib( GL_LIGHTING );
  glPushAttrib( GL_COLOR );
  glPushAttrib( GL_LINE_WIDTH );

  glDisable( GL_LIGHTING );

  glColor3d( 0.0, 0.0, 0.0 );
  glLineWidth( 1.0 );

  glBegin( GL_LINES );
  for( unsigned sample_num = 0; sample_num < samples.cols(); ++sample_num )
  {
    assert( fabs( m_sample_normals.col( sample_num ).norm() - 1.0 ) <= 1.0e-6 );
    glVertex3d( GLdouble( samples( 0, sample_num ) ), GLdouble( samples( 1, sample_num ) ), GLdouble( samples( 2, sample_num ) ) );
    glVertex3d( GLdouble( samples( 0, sample_num ) + 0.1 * m_sample_normals( 0, sample_num ) * grid_width ), GLdouble( samples( 1, sample_num ) + 0.1 * m_sample_normals( 1, sample_num ) * grid_width ), GLdouble( samples( 2, sample_num ) + 0.1 * m_sample_normals( 2, sample_num ) * grid_width ) );
  }
  glEnd();

  glPopAttrib();
  glPopAttrib();
  glPopAttrib();
}
