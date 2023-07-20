function [] = frontiers_4shell_pipeline_step5(datapath, codedir, ico)
% Wrapper function to build 4-shell model in the frontiers pipeline

if nargin < 2
    tmp = matlab.desktop.editor.getActive;
    codedir = [fileparts(tmp.Filename), '/..'];  % full directory to the forward_model code folder
    clearvars tmp
    ico = 'ico5';
elseif nargin < 3
    ico = 'ico5';
end

% addpath to the forward modeling code folders
addpath(genpath(fullfile(codedir, 'src')))
addpath(genpath(fullfile(codedir, 'external_toolbox')))

%% Step 1: Create CSF surface from headreco - 6min
mesh_make_csf_surf(datapath)

%% Step 2: Create Skull from MNE_outer_skull and headreco - 10min
mesh_make_skull_surf(datapath)

%% Step 3: Create 4shell BEM models - 125min
mesh_create_4shell(datapath, ico)

%% Step 4: HBF BEM solution - 20min on nodes with RAM > 100GB
hbf_build_fwd(datapath, ico)

%% Step 5: LFM Post-processing - 1min
hbf_clean_LFM(datapath, ico)

end
