import mne
import scipy.io as sio

trans = mne.read_trans('trans.fif')
info = mne.io.read_info('raw.fif')

sio.savemat('elc-trans.mat', {'trans': trans['trans'], 'dig': [x['loc'][0:3] for x in info['chs']]})
