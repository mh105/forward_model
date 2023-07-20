function [ ] = hbf_build_fwd(datapath, ico)
%% Helsinki BEM Framework LCISA Solver to build 4-shell Forward Model

% close all
% clear all

% % Change the current folder to the folder of this m-file.
% tmp = matlab.desktop.editor.getActive;
% cd(fileparts(tmp.Filename));
% clearvars tmp

if nargin < 2
    ico = 'ico5';
end

%% Check if enough RAM on the running OS
if ismac  % Mac OS
    [status, cmdout] = system("sysctl hw.memsize | awk '{print $2}'");
    mem_GB = str2double(cmdout) / 10^9;
elseif isunix  % Linux OS
    [status, cmdout] = system("grep MemTotal /proc/meminfo | awk '{print $2}'");
    mem_GB = str2double(cmdout) / 10^6;
else  % Windows OS
    user = memory;
    status = 0;
    mem_GB = user.MemAvailableAllArrays / 10^9;
end

if status ~= 0
    error('Memory check failed during solving HBF BEM.')
end

assert(mem_GB > 100, ['Insufficient amount of memory for the HBF BEM solver: ', num2str(mem_GB), ' GB, need at least 100 GB.'])

%% Path configuration

% add path to this EEG/code/forward_model folder
% addpath(genpath('./'));

% create datapath for our own BEM meshes
% datapath = '/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/MNE/subject_data/SP011_ah/simnibs_pipe/bem_surfaces/final_structure';

warning('off', 'all')

%% readme.txt from hbf package

% The example scripts present all necessary file formats. Regarding boundary meshes, it is essential to know that
% 1) meshes are ordered from inwards to outwards
% 2) meshes must be closed and
% 3) triangle orientation is CCW --- please run hbf_CheckTriangleOrientation to check this
% 4) the resulting models have been verified with meshes of approx 2500 vertices per surface.
% 5) the recommended distance for dipole sources from the inner skull is half of triangle side length.
% 6) all physical measures are in SI units --- spatial variables in meters, dipole moments in ampere meters.
% 7) all sources must be inside the ISA surface --- in the case of three-shell MEG models, inside surface one.

%% Accomodate looping through ico3-5 source spaces

if strcmp(ico, '345')
    ico_list = {'ico3', 'ico4', 'ico5'};
else
    ico_list = {ico};
end

for ii = 1:length(ico_list)
    % iterate through different icosahedron resolutions
    ico = ico_list{ii};
    
    %% Build the four-shell forward model with our own mesh
    
    % load the 4-shell BEM models
    load(fullfile(datapath, ['4shell_hbfcmpt_BEM_model-', ico, '.mat']),...
        'bmeshes', 'dig', 'sources');%bmeshes, eegsens, sources
    
    faces_list = {bmeshes{5}.e, bmeshes{4}.e, bmeshes{3}.e, bmeshes{2}.e, bmeshes{1}.e};
    vertices_list = {bmeshes{5}.p, bmeshes{4}.p, bmeshes{3}.p, bmeshes{2}.p, bmeshes{1}.p};
    mesh_plot_bem_surfaces(faces_list, vertices_list);
    title('BEM surfaces', 'FontSize', 20)
    drawnow
    
    disp('Vertex numbers:')
    totalvertex = 0;
    for j = 1:length(bmeshes)
        disp(['Surface ', num2str(j) ': ', num2str(bmeshes{j}.nop)])
        totalvertex = totalvertex + bmeshes{j}.nop;
    end
    disp(['Total vertex number: ', num2str(totalvertex)])
    
    % concatenate the two hemispheres
    sourcepos=[sources{1}.points; sources{2}.points]; %source locations
    sourcedir=[sources{1}.normals; sources{2}.normals]; %source directions --- for LFM, each vector must be unit length!
    
    elecpos=dig;
    % project electrodes to scalp and make "elecs" struct
    elecs=hbf_ProjectElectrodesToScalp(elecpos,bmeshes);
    
    % set conductivities
    cratio=50;
    c_brain=.33;c_csf=1.79;c_skull=c_brain/cratio;c_skin=.33;
    ci=[c_brain,c_brain,c_csf,c_skull,c_skin];%conductivity inside each surface --- remember the order!
    co=[c_csf,c_csf,c_skull,c_skin,0];%conductivity outside each surface --- remember the order!
    
    %--------------------------------------------------------------------------
    % these check functions are same as before, we should be able to pass them
    % with ease. If not, I'm in huge trouble.
    Nmeshes=length(bmeshes);
    success=zeros(Nmeshes,1);
    for M=1:Nmeshes
        [bmeshes{M},success(M)]=hbf_CorrectTriangleOrientation(bmeshes{M});
    end
    %   check, is the mesh valid for BEM
    status=cell(Nmeshes,1);
    for M=1:Nmeshes
        [status{M}]=hbf_CheckMesh(bmeshes{M});
    end
    
    % ------------------------------------------------------------------------
    % Hey, ho, let's go! - this takes two hours or so for a total of 45000
    % vertices.
    
    starttime=clock;
    % Build double-layer matrices for the BEM
    D=hbf_BEMOperatorsPhi_LC(bmeshes);
    % Build BEM transfer matrix for potential using the isolated source
    % approach, isolation set to surface 3 = inner skull
    Tphi_full=hbf_TM_Phi_LC_ISA2(D,ci,co,3);
    %Extract/interpolate EEG transfer matrix for electrodes only
    Tphi_elecs=hbf_InterpolateTfullToElectrodes(Tphi_full,bmeshes,elecs);
    % compute lead field matrices for directed dipoles
    %   forward solutions for directed dipoles
    LFMphi_dir=hbf_LFM_LC(bmeshes,Tphi_elecs,sourcepos,sourcedir);
    %   forward solutions for xyz-oriented triplets
    LFMphi_xyz=hbf_LFM_LC(bmeshes,Tphi_elecs,sourcepos);
    
    fprintf('Time for LFM computation was %ds.\n',round(etime(clock,starttime)));
    
    save(fullfile(datapath, ['4shell_hbfBEM_rawLFM-', ico, '.mat']), ...
        'LFMphi_dir', 'LFMphi_xyz', 'sourcepos', 'sourcedir', 'elecs', 'bmeshes')

end

%% Completed!

warning('on', 'all')

pause(10)
close all

end
