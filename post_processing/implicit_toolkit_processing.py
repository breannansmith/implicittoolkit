'''Tools for loading 3D rigid body simulations.'''

import sys
import numpy


def matrixIsRotation(R):
    '''Check if a 3x3 matrix is orthonormal and orientation preserving.'''
    assert R.shape == (3, 3)
    idresid = numpy.amax(numpy.absolute(numpy.transpose(numpy.matrix(R)) * numpy.matrix(R) - numpy.identity(3)))
    if idresid > 1.0e-9:
        return False
    if abs(numpy.linalg.det(R) - 1.0) > 1.0e-9:
        return False
    return True


class ImplicitGeometry(object):
    '''A container for an implicit representation of geometry.'''
    def __init__(self, h5_file):
        try:
            # Input settings to the distance field generator
            self.cell_width = h5_file['/settings/cell_width'][:]
            assert self.cell_width.shape == (1, 1)
            self.cell_width = self.cell_width[0, 0]
            self.edge_subsamples = h5_file['/settings/edge_subsamples'][:]
            assert self.edge_subsamples.shape == (1, 1)
            self.edge_subsamples = self.edge_subsamples[0, 0]
            self.face_subsamples = h5_file['/settings/face_subsamples'][:]
            assert self.face_subsamples.shape == (1, 1)
            self.face_subsamples = self.face_subsamples[0, 0]
            self.grid_padding = h5_file['/settings/grid_padding'][:]
            assert self.grid_padding.shape == (1, 1)
            self.grid_padding = self.grid_padding[0, 0]
            # Signed distance field settings
            self.cell_delta = h5_file['/sdf/cell_delta'][:, 0]
            assert self.cell_delta.shape == (3,)
            self.grid_dimensions = h5_file['/sdf/grid_dimensions'][:, 0]
            assert self.grid_dimensions.shape == (3,)
            self.grid_origin = h5_file['/sdf/grid_origin'][:, 0]
            assert self.grid_origin.shape == (3,)
            self.signed_distances = h5_file['/sdf/signed_distance'][:, 0]
            assert self.signed_distances.shape == (numpy.product(self.grid_dimensions),)
            # Mesh data
            self.mesh_verts = h5_file['/mesh/vertices'][:]
            assert self.mesh_verts.shape[0] == 3
            self.mesh_faces = h5_file['/mesh/faces'][:]
            assert self.mesh_faces.shape[0] == 3
            self.mesh_edges = h5_file['/mesh/edges'][:]
            assert self.mesh_edges.shape[0] == 2
            # The mesh's moments
            self.I_on_rho = h5_file['/moments/I_on_rho'][:, 0]
            assert self.I_on_rho.shape[0] == 3
            assert numpy.all(self.I_on_rho > 0.0)
            self.R = h5_file['/moments/R'][:]
            assert self.R.shape == (3, 3)
            assert matrixIsRotation(self.R)
            self.volume = h5_file['/moments/volume'][0, 0]
            assert self.volume > 0.0
            self.x = h5_file['/moments/x'][:, 0]
            assert self.x.shape[0] == 3
            # Vertices in the convex hull
            self.convex_hull_verts = h5_file['/convex_hull/vertices'][:]
            assert self.convex_hull_verts.shape[0] == 3
            # Surface samples
            self.surface_samples = h5_file['/surface_samples/samples'][:]
            assert self.surface_samples.shape[0] == 3
            # Sphere tree
            self.sphere_tree_depth = h5_file['/sphere_tree/depth'][0, 0]
            assert self.sphere_tree_depth >= 0
            self.sphere_tree_leaf_contents = h5_file['/sphere_tree/leaf_contents'][:]
            assert self.sphere_tree_leaf_contents.shape[0] == 2
            self.sphere_tree_leaf_node_ids = h5_file['/sphere_tree/leaf_node_ids'][:, 0]
            self.sphere_tree_node_centers = h5_file['/sphere_tree/node_centers'][:]
            assert self.sphere_tree_node_centers.shape[0] == 3
            self.sphere_tree_node_radii = h5_file['/sphere_tree/node_radii'][:, 0]
        except KeyError as key_exception:
            sys.exit('HDF5 Key Error: ' + key_exception.message)
