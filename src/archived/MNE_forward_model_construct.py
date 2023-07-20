"""
Script for constructing MNE forward-models from MRI scans
- T1 structural with watershed algorithm
- 5/30 degree FLASH sequences
------------------------
Alex He, May 29th, 2019
"""

###########################
# Configuring MNE python
# https://martinos.org/mne/stable/auto_tutorials/misc/plot_configuration.html
###########################
import os.path as op
import mne
from mne.datasets.sample import data_path

fname = op.join(data_path(), 'MEG', 'sample', 'sample_audvis_raw.fif')
raw = mne.io.read_raw_fif(fname).crop(0, 10)
original_level = mne.get_config('MNE_LOGGING_LEVEL', 'INFO')

# print the path to the configuration file
print(mne.get_config_path())

# Files inside this folder should never be modified manually. Let’s see what the configurations contain.
print(mne.get_config())

# Set the default logging level for the functions
mne.set_config('MNE_LOGGING_LEVEL', 'INFO')
print(mne.get_config(key='MNE_LOGGING_LEVEL'))

# We can set the global logging level for only this session by calling mne.set_log_level() function.
mne.set_log_level('WARNING')
print(mne.get_config(key='MNE_LOGGING_LEVEL'))

# Call a function and see what messages are printed
cov = mne.compute_raw_covariance(raw, verbose=True)

# Revert back to original config value
mne.set_config('MNE_LOGGING_LEVEL', original_level)
print('Config value restored to: %s' % mne.get_config(key='MNE_LOGGING_LEVEL'))






###########################
# Source-level analysis:
# Forward Operator Construction

# - to compute the forward solution, we need 3 things:
# 1) BEM surfaces and BEM solution
# 2) head-MRI coregistration -trans.fif
# 3) a source space
###########################


###########################
# 1) BEM surfaces and BEM solution
###########################

import os.path as op
import mne
from mne.datasets import sample
data_path = sample.data_path()


# the raw file containing the channel location + types
raw_fname = data_path + '/MEG/sample/sample_audvis_raw.fif'
# The paths to Freesurfer reconstructions
subjects_dir = data_path + '/subjects'
subject = 'sample'

mne.bem.make_flash_bem(subject=subject, subjects_dir=subjects_dir, verbose=True)

# Let's tackle 3) BEM surfaces and BEM solution first!
# Primary goal: --- Compute and visualize BEM meshes

# first let's see the computed one for MNE sample data
mne.viz.plot_bem(subject=subject, subjects_dir=subjects_dir,
                 brain_surfaces='white', orientation='coronal')

# new sequence 8 separate images
subject = 'SLEEP_TDEL_recon'
subjects_dir = '/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/MNE/subject_data'
flash_path = '/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/MNE/subject_data/SLEEP_TDEL_M_95_20190507/raw/flash5_8scan/014'

mne.bem.make_flash_bem(subject=subject, subjects_dir=subjects_dir, flash_path=flash_path, verbose=True)


# new sequence RMS image
subject = 'SLEEP_TDEL_newRMS'
subjects_dir = '/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/MNE/subject_data'
flash_path = '/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/MNE/subject_data/SLEEP_TDEL_M_95_20190507/raw/MEFLASH_8e_1mm_iso_5deg/015'

mne.bem.make_flash_bem(subject=subject, subjects_dir=subjects_dir, flash_path=flash_path, verbose=True)

# old sequence RMS image
subject = 'SLEEP_TDEL_oldflash'
subjects_dir = '/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/MNE/subject_data'
flash_path = '/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/MNE/subject_data/SLEEP_TDEL_M_95_20190507/raw/FLASH05deg/013'

mne.bem.make_flash_bem(subject=subject, subjects_dir=subjects_dir, flash_path=flash_path, verbose=True)


###########################
# 2) head-MRI coregistration -trans.fif
###########################
import mne
mne.gui.coregistration(subject='SLEEP_TDEL_recon',
                       subjects_dir='/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/MNE/subject_data',
                       inst='/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/EEG/data/TDelp/set/TDelp_resting_raw.fif',
                       mark_inside=True, interaction='terrain')


###########################
# 3) a source space
###########################
import mne

# we will set-up a surface-based source space to test out the pipeline,
# alternatively, we could set up a volumetric (discrete) source space.
subject = 'SLEEP_TDEL_recon'
subjects_dir = '/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/MNE/subject_data'

# maybe we want to use higher sampling with ico5 and then patching
src = mne.setup_source_space(subject=subject, subjects_dir=subjects_dir,
                             spacing='ico4', add_dist=False, n_jobs=2)

# visualize the source space on BEM surfaces
mne.viz.plot_bem(subject=subject, subjects_dir=subjects_dir,
                 brain_surfaces='white', src=src, orientation='coronal')

# visualize the source space on expanded white matter surface
import numpy as np  # noqa
from mayavi import mlab  # noqa
from surfer import Brain  # noqa
# both hemispheres
brain = Brain(subject_id=subject, hemi='both', surf='inflated', subjects_dir=subjects_dir, interaction='terrain')
surf = brain.geo['lh']
vertidx = np.where(src[0]['inuse'])[0]
mlab.points3d(surf.x[vertidx], surf.y[vertidx],
              surf.z[vertidx], color=(1, 1, 0), scale_factor=1.5)
surf = brain.geo['rh']
vertidx = np.where(src[1]['inuse'])[0]
mlab.points3d(surf.x[vertidx], surf.y[vertidx],
              surf.z[vertidx], color=(1, 1, 0), scale_factor=1.5)

# save the source space to fif file
mne.write_source_spaces(fname='/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/EEG/data/TDelp/set/TDelp_resting-ico5-src.fif',
                        src=src)


##########################################
# Check alignment of electrodes and source
##########################################
import mne
import matplotlib.pyplot as plt
plt.get_backend()
subject = 'SLEEP_TDEL_recon'
subjects_dir = '/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/MNE/subject_data'
raw = mne.io.read_raw_fif('/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/EEG/data/TDelp/set/TDelp_resting_raw.fif')
trans = mne.read_trans('/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/EEG/data/TDelp/set/TDelp_resting-trans.fif')
src = mne.read_source_spaces('/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/EEG/data/TDelp/set/TDelp_resting-ico5-src.fif')
mne.viz.plot_alignment(raw.info, trans=trans, eeg='projected', subject=subject,
                       subjects_dir=subjects_dir, meg=[], dig=True,
                       surfaces=['head'], coord_frame='head', interaction='terrain')
mne.viz.plot_alignment(raw.info, trans=trans, subject=subject,
                       src=src, subjects_dir=subjects_dir, meg=[], dig=True,
                       surfaces=['head', 'white'], coord_frame='head', interaction='terrain')


#########################################
# Finally ready to compute forward model!
#########################################
import mne
import matplotlib.pyplot as plt
subject = 'SLEEP_TDEL_recon'
subjects_dir = '/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/MNE/subject_data'

# make bem solution
conductivity = (0.3, 0.006, 0.3)  # for three layers
model = mne.make_bem_model(subject=subject, ico=4,
                           conductivity=conductivity,
                           subjects_dir=subjects_dir)
bem = mne.make_bem_solution(model)

# read in the trans file
trans = mne.read_trans('/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/EEG/data/TDelp/set/TDelp_resting-trans.fif')

# read in the source space
src = mne.read_source_spaces('/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/EEG/data/TDelp/set/TDelp_resting-ico5-src.fif')

# compute the forward operator (aka lead field matrix)
raw_fname = '/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/EEG/data/TDelp/set/TDelp_resting_raw.fif'
fwd = mne.make_forward_solution(raw_fname, trans=trans, src=src, bem=bem,
                                meg=False, eeg=True, mindist=5.0, n_jobs=2)

leadfield = fwd['sol']['data']
print("Leadfield size : %d sensors x %d dipoles" % leadfield.shape)

# save the forward solution to a fif file
# use mne.read_forward_solution() to read it back in
mne.write_forward_solution(fname='/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/EEG/data/TDelp/set/TDelp_resting-fwd.fif',
                           fwd=fwd, overwrite=True)





