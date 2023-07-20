import scipy.io as sio
import nibabel
import glob

valid_files = glob.glob('4layer_BEM_model_all_meshes-ico*.mat')
assert len(valid_files) > 0, 'No valid 4layer_BEM_model_all_meshes-ico*.mat to create .surf files'
valid_fn = valid_files[0]  # all ico resolutions use the same bem surfaces

# load in the surfaces in .mat file format and save the 3shell BEM as .surf
mat_contents = sio.loadmat(valid_fn)

# save outer_skin.surf
new_coords = mat_contents['v1']
new_faces = mat_contents['f1']
new_faces = new_faces - 1  # update the index for python indexing scheme that starts from 0

# Save the surface as a Freesurfer .surf file
filename = 'outer_skin.surf'
nibabel.freesurfer.io.write_geometry(filename, new_coords, new_faces,
                                     create_stamp='Created using hbf bmeshes surfaces', volume_info=None)

# save outer_skull.surf
new_coords = mat_contents['v2']
new_faces = mat_contents['f2']
new_faces = new_faces - 1  # update the index for python indexing scheme that starts from 0

# Save the surface as a Freesurfer .surf file
filename = 'outer_skull.surf'
nibabel.freesurfer.io.write_geometry(filename, new_coords, new_faces,
                                     create_stamp='Created using hbf bmeshes surfaces', volume_info=None)

# save inner_skull.surf
new_coords = mat_contents['v3']
new_faces = mat_contents['f3']
new_faces = new_faces - 1  # update the index for python indexing scheme that starts from 0

# Save the surface as a Freesurfer .surf file
filename = 'inner_skull.surf'
nibabel.freesurfer.io.write_geometry(filename, new_coords, new_faces,
                                     create_stamp='Created using hbf bmeshes surfaces', volume_info=None)
