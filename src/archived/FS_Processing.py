"""
Bash functions to run to process a subject's MRI data

Also constructs watershed and FLASH BEM surfaces
"""
import os.path as op
import mne
import nibabel
import shutil
import scipy.io as sio
import matplotlib.pyplot as plt
import numpy as np  # noqa
from mayavi import mlab  # noqa
from surfer import Brain  # noqa
plt.get_backend()

# Step 1: grab the dicom files from the scan

findsession <subject>
rsync /cluster/archive/305/siemens/TrioTim-35006-20131029-135028-570000/* /path/to/archive/subj/dicom/

# Step 2: extract scan.info

unpacksdcmdir
   -src <srcdir>
   -targ <targdir>
   -scanonly <targdir>/scan.info

# Step 3: use the run numbers from scan.info and
# unpack the desired runs: MPRAGE and 8-scan FLASH05 degree

unpacksdcmdir -src <srcdir> -targ <targdir>
   -run <runnumb> <subdir> <outputformat> <outputname>

e.g.: unpacksdcmdir -src /cluster/archive/328/siemens/TrioTim-35162-20191008-141618-001527
-targ . -run 7 MEMPRAGE mgz MEMPRAGE.mgz -run 12 FLASH05deg mgz flash5.mgz

# Step 4: run recon-all on the subject

setenv SUBJECTS_DIR /full/path/to/project/directory/
recon-all -subjid <subject> -i /full/path/to/your/unpacked/mprage.mgz -all

# Step 5: check the quality of the recon

fvb() { freeview -v "$SUBJECTS_DIR/$1/mri/T1.mgz" \
                    "$SUBJECTS_DIR/$1/mri/wm.mgz:colormap=heat:opacity=0.25:visible=0" \
                    "$SUBJECTS_DIR/$1/mri/aparc+aseg.mgz:colormap=lut:opacity=0.2" \
                    "$SUBJECTS_DIR/$1/mri/brainmask.mgz" \
-f "$SUBJECTS_DIR/$1/surf/lh.white:edgecolor=blue" "$SUBJECTS_DIR/$1/surf/rh.white:edgecolor=blue" \
"$SUBJECTS_DIR/$1/surf/lh.pial:edgecolor=red" "$SUBJECTS_DIR/$1/surf/rh.pial:edgecolor=red"; }

# Step 6: move the flash5.mgz files to appropriate folder within the recon folder

mkdir ./recon/mri/flash
cp ./raw/FLASH05deg/0*/* ./recon/mri/flash

# Step 7: run watershed and FLASH5 to construct BEM surfaces

subject = 'JFIT_F_41'
subjects_dir = '/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/MNE/subject_data'
flash_path = op.join(subjects_dir, subject, 'mri', 'flash')

# create FLASH BEM surfaces
mne.bem.make_flash_bem(subject=subject, subjects_dir=subjects_dir, flash_path=flash_path,
                       overwrite=True, verbose=True)
# visualize the FLASH BEM surfaces
mne freeview_bem_surfaces -s 'recon' -m 'flash'

# create watershed BEM surfaces
mne.bem.make_watershed_bem(subject=subject, subjects_dir=subjects_dir,
                           overwrite=False, verbose=True)
# visualize the watershed BEM surfaces
mne freeview_bem_surfaces -s 'recon' -m 'watershed'

# Step 8: create a new symbolic link to point the outer_skin.surf file to the watershed surface

rm -f ./recon/bem/outer_skin.surf
ln -s /Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/MNE/subject_data/NSS05/recon/bem/watershed/recon_outer_skin_surface ./recon/bem/outer_skin.surf
# visualize the updated BEM surfaces
mne freeview_bem_surfaces -s 'recon'

# Step 9: construct a source space

src = mne.setup_source_space(subject=subject, subjects_dir=subjects_dir,
                             spacing='ico5', add_dist=True, n_jobs=2)
# save the source space to fif file
fname = op.join(subjects_dir, subject, 'bem', '-ico5-src.fif')
mne.write_source_spaces(fname=fname, src=src)

# Step 10: run the correct BEM surface to correct faulty inner_skull surface

# Now let's load in the BEM inner skull surface .surf file
surfpath = op.join(subjects_dir, subject, 'bem', 'flash', 'inner_skull.surf')
shutil.copy2(surfpath, op.join(subjects_dir, subject, 'bem', 'flash', 'inner_skull_original.surf'))
coords, faces = nibabel.freesurfer.io.read_geometry(surfpath, read_metadata=False, read_stamp=False)
faces = faces + 1  # update the index for MATLAB indexing scheme that starts from 1

# Save a .mat file for visualization and editing
# Multiple by 1000 to convert to mm unit
savepath = op.join(subjects_dir, subject, 'bem', 'flash', 'BEM_coord.mat')
sio.savemat(savepath, {'vertices': coords, 'faces': faces, 'lsrc': src[0]['rr']*1000, 'lidx': src[0]['inuse'],
                              'rsrc': src[1]['rr']*1000, 'ridx': src[1]['inuse']})

# Editing of triangulation happens in Matlab...
MNE_editbem('BEM_coord.mat', op.join(subjects_dir, subject, 'bem', 'flash'))

# Load the edited .mat file containing surface information
loadpath = op.join(subjects_dir, subject, 'bem', 'flash', 'BEM_coord_editted.mat')
mat_contents = sio.loadmat(loadpath)
new_coords = mat_contents['vertices']
new_faces = mat_contents['faces']
new_faces = new_faces - 1  # update the index for python indexing scheme that starts from 0

# Save the surface as a Freesurfer .surf file
filepath = op.join(subjects_dir, subject, 'bem', 'flash', 'inner_skull.surf')
nibabel.freesurfer.io.write_geometry(filepath, new_coords, new_faces, create_stamp='Edited to contain all source points', volume_info=None)

# Step 11: visualize the source space on expanded white matter surface

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

# Step 12: find the digitization of electrode locations
# if a raw.fif file was generated by ANT_interface_code toolbox, it will have
# the correct digitization of electrodes so this step can be skipped as well

dig_fif_fn = 'NSS11-2.fif'
dig_fif_filepath = op.join(subjects_dir, subject, 'bem')

# Step 13: create a fake raw object if one isn't generated using ANT_interface_code toolbox
from NSS_MNE_python_util import *

setfn = 'fake_65_setfile.set'
# this set file can be copied and used for all subjects from the NSS data set
set_filepath = '/Users/alexhe/Dropbox (Personal)/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/reference_modeling/util_struct'

# save this raw object in a FIF format to save_filepath directory
save_filepath = op.join(subjects_dir, subject, 'bem')
fif_fname = op.join(subjects_dir, subject, 'bem', '_raw.fif')
# need to manually rename this to _raw.fif file. Some file identifier is kept to avoid overwriting
# from different subjects' data
raw = NSS_mne_create_raw(setfn, set_filepath, dig_fif_fn, dig_fif_filepath, save_filepath)

# Step 14: generate a -trans.fif file
fif_fname = '/Users/alexhe/Dropbox/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/EEG/data/JFIT_F_41/set/JFIT_F_41_night2_Resting_ds500_raw.fif'
mne.gui.coregistration(subject=subject, subjects_dir=subjects_dir,
                       inst=fif_fname,
                       mark_inside=True, interaction='terrain')
# use the command line version, python interaction doesn't work in Catalina OS
mne coreg -d subjects_dir -s subject -f fif_fname --mark-inside --interaction="terrain"

# Step 15: check the alignment of electrodes and source space

raw = mne.io.read_raw_fif(op.join(subjects_dir, subject, 'bem', '_raw.fif'))
trans = mne.read_trans(op.join(subjects_dir, subject, 'bem', 'recon-trans.fif'))
src = mne.read_source_spaces(op.join(subjects_dir, subject, 'bem', '-ico5-src.fif'))
mne.viz.plot_alignment(raw.info, trans=trans, eeg='projected', subject=subject,
                       subjects_dir=subjects_dir, meg=[], dig=True,
                       surfaces=['head'], coord_frame='head', interaction='terrain')
mne.viz.plot_alignment(raw.info, trans=trans, subject=subject,
                       src=src, subjects_dir=subjects_dir, meg=[], dig=True,
                       surfaces=['head', 'white'], coord_frame='head', interaction='terrain')

# Step 16: construct a forward model solution (G matrix)


# Step 17: save the outer_skin surface, trans matrix, and electrode locations to .mat file for
# constructing a realistic scalp surface for displaying sensor information (like spindle
# distributions) -> move to the forward construction step to save with rest of data structure

# # Now let's load in the BEM outer_skin surface .surf file
# surfpath = op.join(subjects_dir, subject, 'bem', 'outer_skin.surf')
# coords, faces = nibabel.freesurfer.io.read_geometry(surfpath, read_metadata=False, read_stamp=False)
# faces = faces + 1  # update the index for MATLAB indexing scheme that starts from 1
#
# # Save a .mat file for visualization and editing
# # Multiple by 1000 to convert to mm unit
# savepath = op.join(subjects_dir, subject, 'bem', 'NSS11_surf_eloc.mat')
# sio.savemat(savepath, {'vertices': coords, 'faces': faces,
#                        'dig': [x['loc'][0:3] for x in raw.info['chs']], 'trans': trans['trans']})

