"""
Script for loading Freesurfer .surf file of BEM surfaces
and source object. Saving them as .mat file to be edited MATLAB

- Uses the nibabel package to load and read Freesurfer .surf files

------------------------
Alex He, July 9th, 2019
"""

# First of all, load in the source space on subject's White Matter surface
# and visualize it with Mayavi
import mne
import matplotlib.pyplot as plt
plt.get_backend()
subject = 'SLEEP_TDEL_recon'
subjects_dir = '/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/MNE/subject_data'
raw = mne.io.read_raw_fif('/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/EEG/data/TDelp/set/TDelp_resting_raw.fif')
trans = mne.read_trans('/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/EEG/data/TDelp/set/TDelp_resting-trans.fif')
src = mne.read_source_spaces('/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/EEG/data/TDelp/set/TDelp_resting-ico5-src.fif')
mne.viz.plot_alignment(raw.info, trans=trans, subject=subject,
                       src=src, subjects_dir=subjects_dir, meg=[], dig=True,
                       surfaces=['head', 'white'], coord_frame='head', interaction='terrain')


# break the following lines to a function...


# Now let's load in the BEM inner skull surface .surf file
import nibabel
filepath = '/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/MNE/subject_data/SLEEP_TDEL_recon/bem/flash/inner_skull.surf'
coords, faces = nibabel.freesurfer.io.read_geometry(filepath, read_metadata=False, read_stamp=False)
faces = faces + 1  # update the index for MATLAB indexing scheme that starts from 1


# Save a .mat file for visualization and editing
import scipy.io as sio
# Multiple by 1000 to convert to mm unit
sio.savemat('BEM_coord.mat', {'vertices': coords, 'faces': faces, 'lsrc': src[0]['rr']*1000, 'lidx': src[0]['inuse'],
                              'rsrc': src[1]['rr']*1000, 'ridx': src[1]['inuse']})

# Editing of triangulation happens in Matlab...

# Load the edited .mat file containing surface information
mat_contents = sio.loadmat('new_BEM_surf.mat')
new_coords = mat_contents['vertices']
new_faces = mat_contents['faces']
new_faces = new_faces - 1  # update the index for python indexing scheme that starts from 0

# Save the surface as a Freesurfer .surf file
filepath = '/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/MNE/subject_data/SLEEP_TDEL_recon/bem/flash/inner_skull_test.surf'
nibabel.freesurfer.io.write_geometry(filepath, new_coords, new_faces, create_stamp='Edited to contain all source points', volume_info=None)

