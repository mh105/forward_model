import mne
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

for ico in ico_list:
    # load in the G matrix
    mat_contents = sio.loadmat('4shell_hbfBEM_finalLFM-' + ico + '.mat')

    # read a forward model fif file
    fwd = mne.read_forward_solution(fname='hbf2MNE_hbfBEM-' + ico + '-fwd.fif')

    # confirm that dimensions of LFM are correct
    assert fwd['sol']['data'].shape == mat_contents['xyzG'].shape, 'Dimensions of LFM mismatch!'
    assert fwd['_orig_sol'].shape == mat_contents['xyzG'].shape, 'Dimensions of LFM mismatch!'

    # update the leadfield matrix
    fwd['sol']['data'] = mat_contents['xyzG']
    fwd['_orig_sol'] = mat_contents['xyzG']

    # write the updated forward model fif file
    mne.write_forward_solution(fname='four_layer-' + ico + '-fwd.fif', fwd=fwd, overwrite=True)
