import mne
import os
import sys
import scipy.io as sio
import numpy as np

ico = 'ico5'

if __name__ == '__main__':
    if len(sys.argv) == 2:
        ico = sys.argv[1]

if ico == '345':
    ico_list = ('ico3', 'ico4', 'ico5')
else:
    ico_list = (ico, )

# We will build the MNE forward model object using the exact same surfaces as in hbf
subject = 'm2m_recon'
subjects_dir = os.path.join(os.getcwd(), '../final_surface')

# mne.viz.plot_alignment(info, trans=trans, eeg='projected', subject=subject,
#                        subjects_dir=subjects_dir, meg=[], dig=True,
#                        surfaces=['head'], coord_frame='head', interaction='terrain')

# using hbf default conductivity values
conductivity = (0.33, 0.0066, 0.33)  # for three layers: Brain/CSF - Skull - Head

# make bem solution
model = mne.make_bem_model(subject=subject, ico=None,
                           conductivity=conductivity,
                           subjects_dir=subjects_dir)
bem = mne.make_bem_solution(model)

# read in trans
trans = mne.read_trans('trans.fif')

# read in raw.info
info = mne.io.read_info('raw.fif')

for ico in ico_list:
    # load in the source space
    src = mne.read_source_spaces(ico + '-src.fif')

    # update the source coordinates
    mat_contents = sio.loadmat('4layer_BEM_model_all_meshes-' + ico + '.mat')
    src[0]['rr'] = mat_contents['lsrc']/1000
    src[1]['rr'] = mat_contents['rsrc']/1000

    # compute the forward operator (aka lead field matrix)
    fwd = mne.make_forward_solution(info, trans=trans, src=src,
                                    bem=bem, meg=False, eeg=True,
                                    mindist=0, n_jobs=4)  # set mindist=0 to avoid excluding sources

    # verify that no source point is excluded
    if ico == 'ico3':
        assert fwd['sol']['data'].shape[1] == 1284 * 3,\
            'Some sources are excluded in the forward model under ' + ico + '.'
    elif ico == 'ico4':
        assert fwd['sol']['data'].shape[1] == 5124 * 3,\
            'Some sources are excluded in the forward model under ' + ico + '.'
    elif ico == 'ico5':
        assert fwd['sol']['data'].shape[1] == 20484 * 3,\
            'Some sources are excluded in the forward model under ' + ico + '.'

    # convert to fixed orientation
    fwd_fixed = mne.convert_forward_solution(fwd, surf_ori=True, force_fixed=True)

    # Let's grab the source coordinates and lead_field_matrix
    # src = fwd_fixed['src'] # let's not grab the src in fwd because it has been shifted to head coordinate system
    sio.savemat('hbf2MNE_hbfBEM-' + ico + '-fwd.mat',
                {'coordinate': np.concatenate((src[0]['rr'], src[1]['rr']), axis=0),
                 'idx': np.concatenate((src[0]['inuse'], src[1]['inuse']), axis=0),
                 'G': fwd_fixed['sol']['data'],
                 'normal': np.concatenate((src[0]['nn'], src[1]['nn']), axis=0),
                 'dig': [x['loc'][0:3] for x in info['chs']],
                 'trans': trans['trans']})

    # save the free-orientation forward solution to a fif file
    # use mne.read_forward_solution() to read it back in
    # to get back the fixed orientation solution, apply convert_forward_solution after loading
    mne.write_forward_solution(fname='hbf2MNE_hbfBEM-' + ico + '-fwd.fif', fwd=fwd, overwrite=True)
