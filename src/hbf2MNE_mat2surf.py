import os.path as op
import scipy.io as sio
import nibabel

# Load the edited .mat file containing surface information
mat_contents = sio.loadmat('morphed_skull.mat')
new_coords = mat_contents['vertices']
new_faces = mat_contents['faces']
new_faces = new_faces - 1  # update the index for python indexing scheme that starts from 0

# Save the surface as a Freesurfer .surf file
filename = 'morphed_skull.surf'
nibabel.freesurfer.io.write_geometry(filename, new_coords, new_faces,
                                     create_stamp='Edited by morphing mri2mesh and headreco surfaces', volume_info=None)


# Load the edited .mat file containing surface information
mat_contents = sio.loadmat('headreco_csf_compatible.mat')
new_coords = mat_contents['vertices']
new_faces = mat_contents['faces']
new_faces = new_faces - 1  # update the index for python indexing scheme that starts from 0

# Save the surface as a Freesurfer .surf file
filename = 'headreco_csf_compatible.surf'
nibabel.freesurfer.io.write_geometry(filename, new_coords, new_faces,
                                     create_stamp='Edited headreco csf surface', volume_info=None)


# Load the edited .mat file containing surface information
if op.exists('mri2mesh_gm_fixed_updated.mat'):
    mat_contents = sio.loadmat('mri2mesh_gm_fixed_updated.mat')
    new_coords = mat_contents['vertices']
    new_faces = mat_contents['faces']
    new_faces = new_faces - 1  # update the index for python indexing scheme that starts from 0

    # Save the surface as a Freesurfer .surf file
    filename = 'mri2mesh_gm_fixed_updated.surf'
    nibabel.freesurfer.io.write_geometry(filename, new_coords, new_faces,
                                         create_stamp='Edited mri2mesh pial surface', volume_info=None)
