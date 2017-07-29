'''Simple example of ImplicitGeometry usage that prints the contents of a precomputed implicit representation.'''

import sys
import os
import h5py
import numpy
import implicit_toolkit_processing

if len(sys.argv) != 2:
    sys.exit('Usage: python {} hdf5_file'.format(sys.argv[0]))

hdf5_file_name = sys.argv[1]

if not os.path.isfile(hdf5_file_name):
    sys.exit('Error, file {} does not exist.'.format(hdf5_file_name))


def printPaddedArray(msg, A):
    '''Prints an array with spaces padded to the front.'''
    padding_length = len(msg)
    if A.shape == 2:
        print msg, ' '.join(numpy.array2string(A[0, :]).replace('\n', '').split())
        for row in range(1, A.shape[0]):
            print ' ' * padding_length, ' '.join(numpy.array2string(A[row, :]).replace('\n', '').split())
    else:
        print msg, ' '.join(numpy.array2string(A[:]).replace('\n', '').split())


try:
    with h5py.File(hdf5_file_name, 'r') as h5_file:
        implicit_geometry = implicit_toolkit_processing.ImplicitGeometry(h5_file)
    print 'Input settings:'
    print '   Cell width:', implicit_geometry.cell_width
    print '   Edge subsamples:', implicit_geometry.edge_subsamples
    print '   Face subsamples:', implicit_geometry.face_subsamples
    print '   Grid padding:', implicit_geometry.grid_padding
    print 'Implicit state:'
    print '   Cell delta:', implicit_geometry.cell_delta
    print '   Grid dimensions:', implicit_geometry.grid_dimensions
    print '   Grid origin:', implicit_geometry.grid_origin
    print '   Signed distances:', implicit_geometry.signed_distances
    print 'Input mesh:'
    printPaddedArray('   Vertices:', implicit_geometry.mesh_verts)
    printPaddedArray('   Faces:', implicit_geometry.mesh_faces)
    printPaddedArray('   Edges:', implicit_geometry.mesh_edges)
    print 'Moments:'
    print '   I / rho:', implicit_geometry.I_on_rho
    printPaddedArray('   R:', implicit_geometry.R)
    print '   Volume:', implicit_geometry.volume
    print '   x:', implicit_geometry.x
    printPaddedArray('Convex Hull Verts:', implicit_geometry.convex_hull_verts)
    printPaddedArray('Surface Samples:', implicit_geometry.surface_samples)
    print 'Sphere Tree:'
    print '   Depth:', implicit_geometry.sphere_tree_depth
    printPaddedArray('   Leaf Contents:', implicit_geometry.sphere_tree_leaf_contents)
    print '   Lead Node IDs:', implicit_geometry.sphere_tree_leaf_node_ids
    printPaddedArray('   Node Centers:', implicit_geometry.sphere_tree_node_centers)
    printPaddedArray('   Node Radii:', implicit_geometry.sphere_tree_node_radii)
except IOError as io_exception:
    sys.exit(str(io_exception))
