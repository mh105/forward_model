import numpy as np
import os.path as op
import nibabel

"""
This script is written to translate the mri2mesh surfaces in order to align with the coordinate system in T1.mgz.

The source of the offset comes from the mri2mesh pipeline that uses fsl to orient with q/sform of 127 for the last 
dimension and for reasons unclear to me swapping y and z dimensions with sign flips. This, combined with flirt co-
registration between images and mri_convert -conform processing, accumulates some small offsets between the surfaces 
created in the mri2mesh pipeline and a typical T1.mgz RAS coordinate system. 

The solution to this problem is to use the T1fs_resample.nii.gz file created in the mri2mesh pipeline, flip back the y 
and z dimensions, and use flirt to co-register to T1.mgz. Then we ask for the transformed RAS coordinate of (0,0,0) in 
the T1fs_swap.nii.gz in the destination volume called T1fs_swap_flirt.nii.gz. We compare this transformed RAS coordinate 
with the --cras offset of T1.mgz, which corresponds to the RAS value of voxel [128,128,128] in the T1.mgz volume. We 
compute the difference between the two offset values and use this script to translate the mri2mesh surfaces by the 
offset difference. This effectively re-aligns the surface center to the T1.mgz [128,128,128] center and makes sure 
these surfaces are in the same coordinate system as the rest of the surfaces from headreco pipeline and MNE BEM.  

There are several other terminal commands that go before and after this script. See full documentation of meshing 
pipeline for details. 
"""


def translate_RAS(coords, x, y, z):
    xform = np.array([[1, 0, 0, x],
                      [0, 1, 0, y],
                      [0, 0, 1, z],
                      [0, 0, 0, 1 ]])

    newcoords = np.concatenate((coords, np.ones((coords.shape[0], 1))), axis=1)
    newcoords = np.matmul(xform, newcoords.transpose()).transpose()
    newcoords = np.delete(newcoords, 3, 1)

    return newcoords


def fix_mri2mesh_surface(filename, trans_offset):
    surfpath = op.join(filename)
    coords, faces = nibabel.freesurfer.io.read_geometry(surfpath, read_metadata=False, read_stamp=False)

    newcoords = translate_RAS(coords, trans_offset[0], trans_offset[1], trans_offset[2])

    filename = filename.replace('rh.', '')
    filename = filename.replace('.surf', '')
    newfilename = filename + '_fixed.surf'
    nibabel.freesurfer.io.write_geometry(newfilename, newcoords, faces, create_stamp='Edited to T1.mgz coordinate system',
                                         volume_info=None)
    return


# Find the offset value to translate

file = open("tmp_offset.txt", "r")
T1offset = file.readline()
T1offset = T1offset.strip('\n')
T1offset_list = np.asarray([float(x) for x in T1offset.split(' ')])

m2moffset = file.readline()
m2moffset = m2moffset.strip('\n')
m2moffset_list = np.asarray([float(x) for x in m2moffset.split('  ')])

trans_offset = np.subtract(m2moffset_list, T1offset_list)


# Fix all mri2mesh surfaces

fix_mri2mesh_surface('mri2mesh_wm.surf', trans_offset)
fix_mri2mesh_surface('mri2mesh_gm.surf', trans_offset)
fix_mri2mesh_surface('mri2mesh_csf.surf', trans_offset)
fix_mri2mesh_surface('mri2mesh_cerebellum.surf', trans_offset)
fix_mri2mesh_surface('mri2mesh_skull.surf', trans_offset)
fix_mri2mesh_surface('mri2mesh_fs_lh_wm.surf', trans_offset)
fix_mri2mesh_surface('mri2mesh_fs_rh_wm.surf', trans_offset)


print('All mri2mesh surfaces are fixed - Done!')

