%% Example script of how to use prep_forward_model_struct() function

close all
clear all

% Change the current folder to the folder of this m-file.
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
clearvars tmp

datapath = '/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/forward_modeling/4layer_BEM/sample_mesh/pipelines/bem_surfaces/final_structure';

filename = fullfile(datapath, '4shell_hbfBEM-fwd.mat');
fsfn = fullfile(datapath, 'TDelp_resting_fastscan_dig.mat');
atlasfn = fullfile(datapath, 'atlas_info.mat');
savefn = fullfile(datapath, 'example_duke128_4shell_forward_model');

forward_model = prep_forward_model_struct(filename, fsfn, atlasfn, savefn);
