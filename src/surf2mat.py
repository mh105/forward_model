import scipy.io as sio
import nibabel
import glob

"""
This script is used to read in the .surf files and save the triangulation surfaces using vertices and faces into .mat 
files. Subsequent processing will be done in MATLAB.
"""

# Loop through all files in the current directory that ends with .surf

listing = glob.glob('*.surf')
for filename in listing:
    coords, faces = nibabel.freesurfer.io.read_geometry(filename, read_metadata=False, read_stamp=False)
    faces = faces + 1  # update the index for MATLAB indexing scheme that starts from 1

    # Save a .mat file for visualization and editing
    newfilename = filename.replace('.surf', '.mat')
    sio.savemat(newfilename, {'vertices': coords, 'faces': faces})

print('All triangulation surfaces are saved to .mat files - Done!')

