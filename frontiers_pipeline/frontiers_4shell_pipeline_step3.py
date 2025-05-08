"""
Python codes to construct a 3-shell forward model

using and only using the MNE-Python package (no SimNIBS)
"""

import subprocess
import os
import os.path as op
import mne
import sys
import numpy as np
from mayavi import mlab
from surfer import Brain

################################################################################################################
# [1] configure path
print("[1] ...Configure path...")

subject = 'FS6'

if __name__ == '__main__':
    if len(sys.argv) == 2:
        subjects_dir = sys.argv[1]
    else:
        subjects_dir = input("Enter the full path to the folder that contains FS6: ")
else:
    subjects_dir = input("Enter the full path to the folder that contains FS6: ")

ant_interface_dir = input("Enter the full path to the ANT_interface_code folder (hit ENTER for default): ")
if ant_interface_dir == '':
    ant_interface_dir = '../../ANT_interface_code'
sys.path.insert(1, ant_interface_dir)
from ANT_MNE_python_util import *

set_filepath = input("Enter full path to the set file folder (hit ENTER for default): ")
if set_filepath == '':
    set_filepath = subjects_dir.replace('/mri/recons', '/set')

dig_filepath = input("Enter full path to the fastscan file folder (hit ENTER for default): ")
if dig_filepath == '':
    dig_filepath = subjects_dir.replace('/mri/recons', '/fastscan')

print('')
print(os.listdir(set_filepath))
print('')
setfn = input("Enter a resting state set filename (e.g. XXXX_M_00_night1_Resting_ds500_Z3.set): ")

print('')
print(os.listdir(dig_filepath))
print('')
csvfn = input("Enter a fastscan_dig.csv filename (e.g. XXXX_M_00_night1_fastscan_dig.csv): ")

assert op.exists(op.join(set_filepath, setfn)), op.join(set_filepath, setfn) + " does not exist!"
assert op.exists(op.join(dig_filepath, csvfn)), op.join(dig_filepath, csvfn) + " does not exist!"

verbose_plot = input("Do you want to skip extra sanity check plotting? (y/n): ")
if verbose_plot == 'y':
    verbose = False
else:
    verbose = True

################################################################################################################
# [2] extract BEM surfaces using two different methods
print("[2] ...Construct BEM surfaces...")

# create watershed BEM surfaces
mne.bem.make_watershed_bem(subject=subject, subjects_dir=subjects_dir,
                           overwrite=True, verbose=True, show=False)

# create FLASH BEM surfaces
flash_path = op.join(subjects_dir, subject, 'mri', 'flash', 'flash5_reg.mgz')
mne.bem.make_flash_bem(subject=subject, subjects_dir=subjects_dir, flash5_img=flash_path, register=False,
                       overwrite=True, verbose=True, show=False)

# create a new symbolic link to point the outer_skin.surf file to the watershed surface
# NOTE: there is a bug with MNE 0.20.0 on how it handles the watershed algorithm such that
# the produced brain is upside-down. Follow the tutorial slides to fix this bug if using MNE 0.20.0.
cwd = os.getcwd()
os.chdir(op.join(subjects_dir, subject))
subprocess.run(['rm', '-f', 'bem/outer_skin.surf'])
outer_skin_link = op.join('watershed', subject+'_outer_skin_surface')
subprocess.run(['ln', '-s', outer_skin_link, 'bem/outer_skin.surf'])
os.chdir(cwd)

# visualize the updated BEM surfaces
# NOTE: there's a bug in mne/utils/misc.py, so you may get an AttributeError: 'NoneType' object has no attribute
# 'readline' error when you close the freeview window. It does not interfere with the display so not a fatal bug.
if verbose:
    subprocess.run(['mne', 'freeview_bem_surfaces', '-s', subject, '-d', subjects_dir])

    checkpoint = input("Do the BEM surfaces look ok? (y/n): ")
    if checkpoint == 'n':
        raise(Exception('Creation of BEM surfaces failed.'))

################################################################################################################
# [3] construct a source space
print("[3] ...Construct source spaces...")

# loop through all of ico3, ico4, and ico5 resolutions
all_src = []
for ico in ('ico3', 'ico4', 'ico5'):
    src = mne.setup_source_space(subject=subject, subjects_dir=subjects_dir,
                                 spacing=ico, add_dist=True, n_jobs=4)
    # save the source space to fif file
    fname = op.join(subjects_dir, subject, 'bem', ico + '-src.fif')
    mne.write_source_spaces(fname=fname, src=src, overwrite=True)
    all_src.append(src)

# src = mne.read_source_spaces(op.join(subjects_dir, subject, 'bem', 'ico5-src.fif'))

################################################################################################################
# [4] Optional: correct BEM inner_skull surface
print("[4] ...Correct inner skull surface...")

import scipy.io as sio
import shutil
import nibabel

src = all_src[-1]  # based on ico-5 sources

# Now let's load in the BEM inner skull surface .surf file
surfpath = op.join(subjects_dir, subject, 'bem', 'flash', 'inner_skull.surf')
shutil.copy2(surfpath, op.join(subjects_dir, subject, 'bem', 'flash', 'inner_skull_original.surf'))
coords, faces = nibabel.freesurfer.io.read_geometry(surfpath, read_metadata=False, read_stamp=False)
faces = faces + 1  # update the index for MATLAB indexing scheme that starts from 1

# Save a .mat file for visualization and editing
# Multiply by 1000 to convert to mm unit
savepath = op.join(subjects_dir, subject, 'bem', 'flash', 'BEM_coord.mat')
sio.savemat(savepath, {'vertices': coords, 'faces': faces, 'lsrc': src[0]['rr']*1000, 'lidx': src[0]['inuse'],
                       'rsrc': src[1]['rr']*1000, 'ridx': src[1]['inuse']})

# Editing of triangulation happens in MATLAB...
subprocess.run(['matlab', '-nodesktop', '-r', 'cd ../src; MNE_editbem("BEM_coord.mat",' +
                '"' + op.join(subjects_dir, subject, 'bem', 'flash') + '"); exit;'])

# Load the edited .mat file containing surface information
loadpath = op.join(subjects_dir, subject, 'bem', 'flash', 'BEM_coord_edited.mat')
mat_contents = sio.loadmat(loadpath)
new_coords = mat_contents['vertices']
new_faces = mat_contents['faces']
new_faces = new_faces - 1  # update the index for python indexing scheme that starts from 0

# Save the surface as a Freesurfer .surf file
filepath = op.join(subjects_dir, subject, 'bem', 'flash', 'inner_skull.surf')
nibabel.freesurfer.io.write_geometry(filepath, new_coords, new_faces,
                                     create_stamp='Edited to contain all source points', volume_info=None)

# visualize the source space on expanded white matter surface
if verbose:
    brain = Brain(subject_id=subject, hemi='both', surf='inflated', subjects_dir=subjects_dir, interaction='terrain')
    surf = brain.geo['lh']
    vertidx = np.where(src[0]['inuse'])[0]
    mlab.points3d(surf.x[vertidx], surf.y[vertidx],
                  surf.z[vertidx], color=(1, 1, 0), scale_factor=1.5)
    surf = brain.geo['rh']
    vertidx = np.where(src[1]['inuse'])[0]
    mlab.points3d(surf.x[vertidx], surf.y[vertidx],
                  surf.z[vertidx], color=(1, 1, 0), scale_factor=1.5)

    checkpoint = input("Does the white matter surface look ok? (y/n): ")
    if checkpoint == 'n':
        raise(Exception('Source space decimation of white matter surface failed.'))

################################################################################################################
# [5] construct the digitization of electrodes from FASTSCAN
print("[5] ...Construct digitization from FastScan...")
print("dig file already created for all sleep subjects.")

# Due to the requirement of user input, you MUST run this from a separate MATLAB terminal
# Look at Session5_MATLAB_fastscan.html -> should already be completed for all sleep subjects

################################################################################################################
# [6] create a _raw.fif file.
print("[6] ...Create raw fif file...")

# Note: this could be any set file. If you want to source localize overnight sleep or anaesthesia data,
# then you might want to create the raw.fif file from the set file for desired recording. We will use
# the short resting state data here as an example to show how to use the interface codes.
raw = ant_mne_create_raw(setfn, set_filepath, csvfn, dig_filepath, overwrite=None)  # don't save during the function

# Explicitly save the raw fif file outside the function
fif_filepath = subjects_dir.replace('/mri/recons', '/fif')
if op.exists(fif_filepath) is False:
    os.makedirs(fif_filepath)
fif_fname = op.join(fif_filepath, setfn.strip('.set') + '-raw.fif')
raw.save(fif_fname, overwrite=True)

# raw = mne.io.read_raw_fif(op.join(fif_filepath, setfn.strip('.set') + '-raw.fif'))

################################################################################################################
# [7] register electrodes to the MRI head, i.e. generating a trans.fif file
print("[7] ...Electrode registration to MRI head...")

fif_fname = op.join(fif_filepath, setfn.strip('.set') + '-raw.fif')

# depending on which MNE version this is, use different commands for coregistration
# NOTE: If you are using Catalina MacOS, the mne.gui.coregistration() call may not work, use bash call mne coreg.
mne_version = mne.__version__
mne_version_n = [int(x) for x in mne_version.split('.')]
if mne_version_n[0] > 0:  # starting from MNE-Python 1.0.0, coreg uses pyvistaqt 3d backend rendering
    mne.gui.coregistration(subject=subject, subjects_dir=subjects_dir, inst=fif_fname,
                           scale_by_distance=False)
elif mne_version_n[1] < 20:
    mne.gui.coregistration(subject=subject, subjects_dir=subjects_dir, inst=fif_fname,
                           mark_inside=True, interaction='terrain')
else:
    subprocess.run(['mne', 'coreg', '-d', subjects_dir, '-s', subject,
                    '-f', fif_fname, '--mark-inside', '--interaction=terrain'])

# load the -trans.fif file just saved in the GUI - note that this name have to be entered manually!
trans = mne.read_trans(op.join(fif_filepath, csvfn.replace('_fastscan_dig.csv', '-trans.fif')))

################################################################################################################
# [8] check the alignment of electrodes and source space [This is the last validity check before computing G]
print("[8] ...Final manual quality checks...")

# registered and projected electrodes and skin surface
mne.viz.plot_alignment(raw.info, trans=trans, eeg='projected', subject=subject,
                       subjects_dir=subjects_dir, meg=[], dig=True,
                       surfaces=['head'], coord_frame='head', interaction='terrain')

checkpoint = input("Does the electrode registration look ok? (y/n): ")
if checkpoint == 'n':
    raise(Exception('Electrode alignment to head MRI failed.'))

# registered but not projected electrodes and source space
mne.viz.plot_alignment(raw.info, trans=trans, subject=subject,
                       src=src, subjects_dir=subjects_dir, meg=[], dig=True,
                       surfaces=['head', 'white'], coord_frame='head', interaction='terrain')

checkpoint = input("Do the electrodes and sources look ok? (y/n): ")
if checkpoint == 'n':
    raise(Exception('Proper placement of electrodes and source space failed.'))

################################################################################################################
# [9] compute forward model G matrix
print("[9] ...Compute forward model solutions...")

# set up the BEM transfer matrices
conductivity = (0.3, 0.006, 0.3)  # for three layers
model = mne.make_bem_model(subject=subject, ico=4,
                           conductivity=conductivity,
                           subjects_dir=subjects_dir)
bem = mne.make_bem_solution(model)

# loop through all of ico3, ico4, and ico5 resolutions
for src, ico, nsrc in zip(all_src, ('ico3', 'ico4', 'ico5'), (1284, 5124, 20484)):
    # compute the forward operator (aka lead field matrix, G matrix)
    fwd = mne.make_forward_solution(raw.info, trans=trans, src=src,
                                    bem=bem, meg=False, eeg=True,
                                    mindist=0, n_jobs=4)  # set mindist=0 to avoid excluding sources

    assert fwd['sol']['data'].shape[1] == nsrc * 3, 'Some sources are excluded in the forward model under ' + ico + '.'

    # save the free-orientation forward solution to a fif file
    # use mne.read_forward_solution() to read it back in
    fwdfn = op.join(fif_filepath, setfn.replace('_Resting_ds500_Z3.set', '-' + ico + '-fwd.fif'))
    mne.write_forward_solution(fname=fwdfn, fwd=fwd, overwrite=True)

################################################################################################################
print("[Congrats!] ...frontiers_4shell_pipeline_step3.py finished...")
