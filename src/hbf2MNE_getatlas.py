import mne
import os
import sys
import pickle
import pandas as pd
import scipy.io as sio
import numpy as np
from utilities import get_vert_atlas_lobe_info, atlas_matlab_update_idx

ico = 'ico5'

if __name__ == '__main__':
    if len(sys.argv) == 2:
        ico = sys.argv[1]

if ico == '345':
    ico_list = ('ico3', 'ico4', 'ico5')
else:
    ico_list = (ico, )

subject = 'm2m_recon'
subjects_dir = os.path.join(os.getcwd(), '../final_surface')

# read in labels
labels = mne.read_labels_from_annot(subject, parc='aparc', subjects_dir=subjects_dir)
atlas_templates = pd.read_csv('atlas_templates.csv')
atlas_templates = atlas_templates.set_index('atlas')

for ico in ico_list:
    # load in the source space
    src = mne.read_source_spaces(ico + '-src.fif')

    # special note: FREESURFER labels such as lh.aparc and rh.aparc are encoded using the vertex numbers for
    # the white matter surface. Technically, the lh.aparc and rh.aparc we just pulled from fs recon folder
    # in the mri2mesh pipeline are labelling the native lh.white and rh.white. However, our sources are
    # built using the mri2mesh_fs_lh_wm_fixed.surf and mri2mesh_fs_rh_wm_fixed.surf. This is ok because
    # the only thing we changed in running fix_mri2mesh_surface.py is just translating the coordinates
    # rather than changing the vertex indexing. Hence, the same vertex-index based labelling on the original
    # lh.white and rh.white work just as well for our source space. What a blessing!

    # find indices of sources for different ROIs
    atlas_info = get_vert_atlas_lobe_info(labels, src, atlas_templates)

    # save the atlas_info into a pickle object for later usage
    with open('atlas_info-' + ico + '.pickle', 'wb') as openfile:
        pickle.dump(atlas_info, openfile)

    # save the atlas_info containing indices of sources for different ROIs and the wm/gm boundary surface
    atlas_info = atlas_matlab_update_idx(atlas_info)  # update index to MATLAB format starting from 1
    sio.savemat('atlas_info-' + ico + '.mat',
                {'atlas_info': atlas_info,
                 'idx': np.concatenate((src[0]['inuse'], src[1]['inuse']), axis=0),
                 'vtx': np.concatenate((src[0]['rr'], src[1]['rr']), axis=0),
                 'fce': np.concatenate((src[0]['tris']+1,
                                        src[1]['tris']+1+len(src[0]['rr'])), axis=0),
                 'src_face': np.concatenate((src[0]['use_tris']+1,
                                             src[1]['use_tris']+1+len(src[0]['rr'])), axis=0)})
