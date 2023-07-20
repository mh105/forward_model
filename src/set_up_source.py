import mne
import os
import sys
import scipy.io as sio

ico = 'ico5'

if __name__ == '__main__':
    if len(sys.argv) == 2:
        ico = sys.argv[1]

if ico == '345':
    ico_list = ('ico3', 'ico4', 'ico5')
else:
    ico_list = (ico, )

subject = 'm2m_recon'
subjects_dir = os.getcwd()

for ico in ico_list:
    src = mne.setup_source_space(subject=subject, subjects_dir=subjects_dir,
                                 spacing=ico, add_dist=True, n_jobs=2)

    fname = os.path.join(subjects_dir, ico + '-src.fif')
    mne.write_source_spaces(fname=fname, src=src, overwrite=True)

    # Save a .mat file for checking bounding by pial surface
    # Multiply by 1000 to convert to mm unit
    savepath = os.path.join(subjects_dir, ico + '-src.mat')
    sio.savemat(savepath, {'lsrc': src[0]['rr']*1000, 'lidx': src[0]['inuse'], 'lnrm': src[0]['nn'],
                           'rsrc': src[1]['rr']*1000, 'ridx': src[1]['inuse'], 'rnrm': src[1]['nn']})
