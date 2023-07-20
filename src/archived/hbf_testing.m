%% Helsinki BEM Framework LCISA Solver (Matti Stenroos) Testing
% this script is used to test the hbf solver by Matti Stenroos. 

% We will first examine the data input structures from the sample meshes
% used in the hbf_example codes. We will then run the solver with these
% example meshes to make sure it the solver works. 

% Then we will configure the BEM meshes I built based off our own data and
% check if we obtain a correct forward model G matrix. 

close all
clear all

% Change the current folder to the folder of this m-file.
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
clearvars tmp

%% Path configuration

% add path to this EEG/code/forward_model folder
addpath(genpath('./'));

% add path to reference modeling folder
addpath('/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/reference_modeling/code/')

% create datapath to matti's example mesh folder
meshpath= '/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/forward_modeling/4layer_BEM/sample_mesh/matti_example_mesh';

% create datapath for our own BEM meshes
datapath = '/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/forward_modeling/4layer_BEM/sample_mesh/pipelines/bem_surfaces/final_structure';

%% readme.txt from hbf package

% The example scripts present all necessary file formats. Regarding boundary meshes, it is essential to know that
% 1) meshes are ordered from inwards to outwards
% 2) meshes must be closed and
% 3) triangle orientation is CCW --- please run hbf_CheckTriangleOrientation to check this
% 4) the resulting models have been verified with meshes of approx 2500 vertices per surface.
% 5) the recommended distance for dipole sources from the inner skull is half of triangle side length.
% 6) all physical measures are in SI units --- spatial variables in meters, dipole moments in ampere meters.
% 7) all sources must be inside the ISA surface --- in the case of three-shell MEG models, inside surface one.

%% hbf_example_EEG_3shell.m
% let's first look at the more basic three shell model to understand the
% formats of input meshes and electrode locations.

% Example on the use of Helsinki BEM Framework kernel:
% 3-shell EEG model.
% 
% The example uses a three-shell head model made by me at MRC Cognition and
% Brain Sciences Unit, Cambridge, UK. The model was made using FLASH MR
% images, FreeSurfer, and MNE toolbox. If you wish to use this model for
% something that you wish to share or publish, please contact me.
%
% (c) Matti Stenroos
% 2 Mar 2016
% clear

% addpath ./hbf_calc
% addpath ./hbf_mesh
%Load / set geometries
% load('~/bioem/exampledata/cbusample/samplehead_cbu.mat',...
%     'bmeshes', 'eegsens', 'sources');%bmeshes, eegsens, sources
load(fullfile(meshpath, 'samplehead_cbu.mat'),...
    'bmeshes', 'eegsens', 'sources');%bmeshes, eegsens, sources
%   the meshes are nested; now make sure that the meshes are in correct order
bmeshes=hbf_SortNestedMeshes(bmeshes);

% Let me plot the bmeshes and see what they look like and what the
% SortNesstedMeshes are doing.
mesh1 = bmeshes{1};
mesh2 = bmeshes{2};
mesh3 = bmeshes{3};

faces_list = {mesh3.e, mesh2.e, mesh1.e};
vertices_list = {mesh3.p, mesh2.p, mesh1.p};
mesh_plot_bem_surfaces(faces_list, vertices_list);

% couple of things to note about bmeshes: 
% 1) it is a column-wise cell, each cell is a structure with 4 fields: p is
% basically vertices, e is basically faces. nop is the number of vertices,
% and noe is the number of faces. 
% 2) these meshes are at extremely low resolution with only 2562 vertices.
% My skin surface alone has 10242 vertices, and my pial surface is at
% 120490 vertices. This is almost 50-fold increase. I don't know how this
% would impact computing time yet. We may need to downsample all meshes.
% 3) everything is in SI unit, so for my BEM meshes we need to modify the
% units from mm to m. 

% hbf_SortNestedMeshes is simply calculating the approximate radius of each
% mesh and sorting them in increasing order. This is fine and should work
% for our meshes as well.




sourcepos=[sources{1}.p; sources{2}.p];%source locations
sourcedir=[sources{1}.nn; sources{2}.nn]; %source directions --- for LFM, each vector must be unit length!

% couple of things to note about sources.
% 1) there are extra information in the source structure that we probably
% won't need to worry about.
% 2) there are only 8196 source points. This is a oct6 resolution. Our
% source space is under ico5 with 20482 source points. Another huge
% increase in dimensionality. 
% 3) source points are divided into left versus right hemispheres. 
sourcep1 = sources{1}.p; % 1 is left
sourcep2 = sources{2}.p; % 2 is right

figure 
hold on
scatter3(sourcep1(:,1), sourcep1(:,2), sourcep1(:,3), 10, 'r', 'filled')
axis equal
scatter3(sourcep2(:,1), sourcep2(:,2), sourcep2(:,3), 10, 'b', 'filled')

faces_list = {sources{1}.e};
vertices_list = {sources{1}.p};
mesh_plot_bem_surfaces(faces_list, vertices_list);
% 4) we need to get source normal directions for our meshes as well. Fine.





elecpos=eegsens.p;
%   project electrodes to scalp and make "elecs" struct
elecs=hbf_ProjectElectrodesToScalp(elecpos,bmeshes);

% ha, just 70 electrodes. We won with 128/9. 
% fine, we just need the p, we don't need e or p2d fields.
faces_list = {eegsens.e};
vertices_list = {eegsens.p};
mesh_plot_bem_surfaces(faces_list, vertices_list);




%Set conductivities
ci=[1 1/50 1]*.33; %conductivity inside each surface --- remember the order!
co=[ci(2:3) 0]; %conductivity outside each surface

% fine, this is typical.




%--------------------------------------------------------------------------
% This part is not needed, if meshing pipeline is established and correctly
% configured
%If new meshes, perform some checks...
Nmeshes=length(bmeshes);
%   check/correct orientation
success=zeros(Nmeshes,1);
for M=1:Nmeshes
    [bmeshes{M},success(M)]=hbf_CorrectTriangleOrientation(bmeshes{M});
end

% we might need to worry about the triangle faces orientation in CCW. If
% not, let's figure out a way to correct that. 


%   check, is the mesh valid for BEM
status=cell(Nmeshes,1);
for M=1:Nmeshes
    status{M}=hbf_CheckMesh(bmeshes{M});
end

% this should be fine, especially since I've carried out more extensive
% tests and corrections of these meshes. The only thing I'm concerned about
% is the vertice-size issue in our high-resolution meshes.

%--------------------------------------------------------------------------
%Now let's get to business!
starttime=clock;
%   BEM geometry matrices
D=hbf_BEMOperatorsPhi_LC(bmeshes);
%   full transfer matrix
Tphi_full=hbf_TM_Phi_LC_ISA2(D,ci,co,1);
%   extract transfer matrix for electrodes
Tphi_elecs=hbf_InterpolateTfullToElectrodes(Tphi_full,bmeshes,elecs);

%   forward solutions for directed dipoles
LFMphi_dir=hbf_LFM_LC(bmeshes,Tphi_elecs,sourcepos,sourcedir);
%   forward solutions for xyz-oriented triplets
LFMphi_xyz=hbf_LFM_LC(bmeshes,Tphi_elecs,sourcepos);

fprintf('Total time was %ds.\n',round(etime(clock,starttime)));

% Ok cool it finished without error. Looks like the 3shell solver is
% working properly. 

%% Examine LFM in closer details
% Well. Let's try to examine the lead field matrix in closer details and
% visualize the loading plots. This should give us confidence that the
% solver is actually working.

% Let's examine the leadfield matrix
G = LFMphi_dir;

% check the histogram of G matrix entries
figure;
subplot(1,2,1)
histogram(G(:), 100)
title('Before detect_badsource()', 'FontSize', 16, 'interpreter', 'none')

[ G, badsrc, badety] = detect_badsource(G);
disp(['Number of bad sources found: ' num2str(sum(badsrc))])
disp(['Number of bad entries found: ' num2str(sum(badety))])

subplot(1,2,2)
histogram(G(:), 100)
title('After detect_badsource()', 'FontSize', 16, 'interpreter', 'none')

% Forward model variables
forward_model.G = G;
forward_model.source = sourcepos;
forward_model.dig = elecs.pproj;
forward_model.normal = sourcedir;
forward_model.skin_surf_face = bmeshes{3}.e;
forward_model.skin_surf_vertice = bmeshes{3}.p;

% Loading plot of a channel 
plot_loading(26, forward_model, false, true, false);

% visualize the bad sources
badsrc_idx = find(badsrc==1);
scatter3(sourcepos(badsrc_idx,1), sourcepos(badsrc_idx,2), sourcepos(badsrc_idx,3), 300, 'r', 'filled')

% looks reasonable. However, we would like to do some more quantitative
% tests of this BEM solver to be more confident in using it.

%% Build the three-shell forward model with our own mesh
% and we also build the exact same forward model using MNE and we will
% compare the G matrices to see if hbf solver does a good job, at least in
% the 3shell case. 

% addpath ./hbf_calc
% addpath ./hbf_mesh

load(fullfile(datapath, '3shell_hbfcmpt_BEM_model.mat'),...
    'bmeshes', 'dig', 'sources');%bmeshes, eegsens, sources
% the meshes are nested; now make sure that the meshes are in correct order
bmeshes=hbf_SortNestedMeshes(bmeshes);

% concatenate the two hemispheres
sourcepos=[sources{1}.points; sources{2}.points]; %source locations
sourcedir=[sources{1}.normals; sources{2}.normals]; %source directions --- for LFM, each vector must be unit length!

elecpos=dig;
% project electrodes to scalp and make "elecs" struct
elecs=hbf_ProjectElectrodesToScalp(elecpos,bmeshes);
% Set conductivities
ci=[1 1/50 1]*.33; %conductivity inside each surface --- remember the order!
co=[ci(2:3) 0]; %conductivity outside each surface

%--------------------------------------------------------------------------
% This part is not needed, if meshing pipeline is established and correctly
% configured
%If new meshes, perform some checks...
Nmeshes=length(bmeshes);
%   check/correct orientation
success=zeros(Nmeshes,1);
for M=1:Nmeshes
    [bmeshes{M},success(M)]=hbf_CorrectTriangleOrientation(bmeshes{M});
end
%   check, is the mesh valid for BEM
status=cell(Nmeshes,1);
for M=1:Nmeshes
    status{M}=hbf_CheckMesh(bmeshes{M});
end

%--------------------------------------------------------------------------
%Now let's get to business!
starttime=clock;
%   BEM geometry matrices
D=hbf_BEMOperatorsPhi_LC(bmeshes);
%   full transfer matrix
Tphi_full=hbf_TM_Phi_LC_ISA2(D,ci,co,1);
%   extract transfer matrix for electrodes
Tphi_elecs=hbf_InterpolateTfullToElectrodes(Tphi_full,bmeshes,elecs);

%   forward solutions for directed dipoles
LFMphi_dir=hbf_LFM_LC(bmeshes,Tphi_elecs,sourcepos,sourcedir);
%   forward solutions for xyz-oriented triplets
LFMphi_xyz=hbf_LFM_LC(bmeshes,Tphi_elecs,sourcepos);

fprintf('Total time was %ds.\n',round(etime(clock,starttime)));

%% Examine LFM in closer details
% Let's examine the leadfield matrix
G = LFMphi_dir;
xyzG = LFMphi_xyz;

% check the histogram of G matrix entries
figure;
subplot(1,2,1)
histogram(G(:), 100)
title('Before detect_badsource()', 'FontSize', 16, 'interpreter', 'none')

[ G, badsrc, badety, xyzG ] = detect_badsource(G, sourcepos, true, xyzG);
disp(['Number of bad sources found: ' num2str(sum(badsrc))])
disp(['Number of bad entries found: ' num2str(sum(badety))])

subplot(1,2,2)
histogram(G(:), 100)
title('After detect_badsource()', 'FontSize', 16, 'interpreter', 'none')

% Forward model variables
forward_model.G = G;
forward_model.xyzG = xyzG;
forward_model.source = sourcepos;
forward_model.dig = elecs.pproj;
forward_model.normal = sourcedir;
forward_model.skin_surf_face = bmeshes{3}.e;
forward_model.skin_surf_vertice = bmeshes{3}.p;

% Loading plot of a channel 
plot_loading(22, forward_model, false, true, false);

% visualize the bad sources
badsrc_idx = find(badsrc==1);
scatter3(sourcepos(badsrc_idx,1), sourcepos(badsrc_idx,2), sourcepos(badsrc_idx,3), 300, 'r', 'filled')

% Let's save the forward solution at this point.
save(fullfile(datapath, '3shell_hbfBEM-fwd.mat'),...
    'bmeshes', 'dig', 'sources', 'G', 'xyzG', 'forward_model')

%% Build the same forward model using MNE 
% we would like to compare the LFMs computed from the two solvers. MNE
% forwad model leadfield matrices were calculated by the bash script
% hbf2MNE.txt 

% Case 1 - we use the exact same BEM surfaces as used by hbf (outer_skin,
% outer_skull, inner_skull) without changing the vertex resolution of them
% by setting ico=None in MNE mne.make_bem_model() function. This should
% allow for a direct comparison of the BEM solver component between hbf and
% MNE, since everything including sources and electrodes should be the
% same. So we expect the leadfield matrix to be very similar, including the
% problematic ones at the back of the occipital lobe. 

% load the data
load(fullfile(datapath, 'hbf2MNE_hbfBEM-fwd.mat'),...
    'G', 'coordinate', 'idx', 'normal')

% first confirm that the parameters for BEM solvers are identical
% 1) source coordinates 
MNE_sourcepos = coordinate(idx==1, :);
assert(all(MNE_sourcepos == sourcepos, 'all'), 'Not all source points are the same!')
% 2) source normals 
MNE_sourcedir = normal(idx==1, :);
assert(all(MNE_sourcedir == sourcedir, 'all'), 'Not all source normal directions are the same!')
% 3) electrode positions
MNE_elecpos = mesh_getelectrode(fullfile(datapath, 'hbf2MNE_hbfBEM_leadfieldMat.mat'), bmeshes{3}.e, bmeshes{3}.p*1000);
MNE_elecpos = MNE_elecpos/1000;
assert(sum(vecnorm(MNE_elecpos - elecpos, 2, 2))/size(elecpos,1) < 0.005, 'Significant difference between electrode positions. Please check!')


% then we will compare the leadfield matrices 
MNE_G = G;
hbf_G = LFMphi_dir;

% quick check of distribution before correcting for bad entries
figure
histogram(MNE_G, 100)
hold on
histogram(hbf_G, 100)

% detect bad sources and exclude them 
[ MNE_G, MNE_badsrc, MNE_badety] = detect_badsource(MNE_G, MNE_sourcepos, true);
[ hbf_G, hbf_badsrc, hbf_badety] = detect_badsource(hbf_G, MNE_sourcepos, true);
% are they very different?
disp(['Number of bad sources for MNE: ', num2str(sum(MNE_badsrc))])
disp(['Number of bad sources for hbf: ', num2str(sum(hbf_badsrc))])
disp(['Number of bad entries for MNE: ', num2str(sum(MNE_badety))])
disp(['Number of bad entries for hbf: ', num2str(sum(hbf_badety))])

% visualize the bad sources
% set up a forwad_model
forward_model.G = hbf_G;
forward_model.source = sourcepos;
forward_model.dig = elecs.pproj;
forward_model.normal = sourcedir;
forward_model.skin_surf_face = bmeshes{3}.e;
forward_model.skin_surf_vertice = bmeshes{3}.p;
% Loading plot of a channel 
plot_loading(22, forward_model, false, true, false);
% visualize the bad sources
badsrc_idx = find(MNE_badsrc==1);
scatter3(sourcepos(badsrc_idx,1), sourcepos(badsrc_idx,2), sourcepos(badsrc_idx,3), 250, 'r', 'filled')
badsrc_idx = find(hbf_badsrc==1);
scatter3(sourcepos(badsrc_idx,1), sourcepos(badsrc_idx,2), sourcepos(badsrc_idx,3), 250, 'gd', 'filled')

% compare the distributions of entries in leadfield matrices
figure
histogram(MNE_G, 100)
hold on
histogram(hbf_G, 100)
legend('MNE', 'hbf')
title('LFM entries distribution')
set(gca, 'FontSize', 20)
% pretty similar 

% Let's now plot a source topoplot of differences between the G matrices
channel = 86;
[ H,hlink ] = hbf2MNE_compareG(MNE_G, hbf_G, 'MNE', 'hbf', MNE_sourcepos, channel, elecpos);
% great! Looking very similar!


% calculate the difference matrix
diffG = MNE_G - hbf_G;

% plot a histogram
figure
histogram(diffG(:), 100)
title('Distribution of differences in G matrices')
hold on
L = xline(std(diffG(:))*2, 'r--', 'LineWidth', 1);
xline(-std(diffG(:))*2, 'r--', 'LineWidth', 1);
legend([L], {'Two-standard-deviation cut-off'})
set(gca, 'FontSize', 14)

disp(['Total number of entries with large difference: ',...
    num2str(sum(abs(diffG(:)) > std(diffG(:))*2))])
disp(['Different entry percentage out of all entries: ',...
    num2str(sum(abs(diffG(:)) > std(diffG(:))*2) / length(diffG(:))*100), '%'])

% plot the difference in G on sources in 3D space
figure
weight = diffG(channel, 1:size(MNE_sourcepos,1)).^2;
size_weight = ceil((weight-min(weight))/(max(weight)-min(weight))*300)+1;
scatter3(MNE_sourcepos(:,1), MNE_sourcepos(:,2), MNE_sourcepos(:,3), size_weight, weight, 'filled')
axis equal
hold on
scatter3(dig(channel,1), dig(channel,2), dig(channel,3), 300, [1,0,0], 's', 'filled')
colorbar
title('MNE - hbf Weighting', 'FontSize', 16)

%%
% Case 2 - we use the MNE BEM surfaces (outer_skin, outer_skull,
% inner_skull) rather than the MNE ones. This means different BEM
% resolution, MNE by default sets it to ico4 with 2562 vertices BEM
% boundary. This is at reduced resolution from our meshes. There might be
% larger differences between the BEM solvers.

% load the data
load(fullfile(datapath, 'hbf2MNE_MNEBEM-fwd.mat'),...
    'G', 'coordinate', 'idx', 'normal')

% first confirm that the parameters for BEM solvers are identical
% 1) source coordinates 
MNE_sourcepos = coordinate(idx==1, :);
assert(all(MNE_sourcepos == sourcepos, 'all'), 'Not all source points are the same!')
% 2) source normals 
MNE_sourcedir = normal(idx==1, :);
assert(all(MNE_sourcedir == sourcedir, 'all'), 'Not all source normal directions are the same!')
% 3) electrode positions
MNE_elecpos = mesh_getelectrode(fullfile(datapath, 'hbf2MNE_hbfBEM_leadfieldMat.mat'), bmeshes{3}.e, bmeshes{3}.p*1000);
MNE_elecpos = MNE_elecpos/1000;
assert(sum(vecnorm(MNE_elecpos - elecpos, 2, 2))/size(elecpos,1) < 0.005, 'Significant difference between electrode positions. Please check!')


% then we will compare the leadfield matrices 
MNE_G = G;
hbf_G = LFMphi_dir;

% due to the unfortunate downsampling to ico4, some sources got cut-off in
% the MNE BEM surfaces. So this no longer allows for a direct comparison
% since we don't know which sources got cut off. We can look at
% distributional differences though. 

% quick check of distribution before correcting for bad entries
figure
histogram(MNE_G, 100)
hold on
histogram(hbf_G, 100)

% detect bad sources and exclude them 
[ MNE_G, MNE_badsrc, MNE_badety] = detect_badsource(MNE_G, [], false);
[ hbf_G, hbf_badsrc, hbf_badety] = detect_badsource(hbf_G, [], false);
% are they very different?
disp(['Number of bad sources for MNE: ', num2str(sum(MNE_badsrc))])
disp(['Number of bad sources for hbf: ', num2str(sum(hbf_badsrc))])
disp(['Number of bad entries for MNE: ', num2str(sum(MNE_badety))])
disp(['Number of bad entries for hbf: ', num2str(sum(hbf_badety))])

% visualize the bad sources
% set up a forwad_model
forward_model.G = hbf_G;
forward_model.source = sourcepos;
forward_model.dig = elecs.pproj;
forward_model.normal = sourcedir;
forward_model.skin_surf_face = bmeshes{3}.e;
forward_model.skin_surf_vertice = bmeshes{3}.p;
% Loading plot of a channel 
plot_loading(22, forward_model, false, true, false);
% visualize the bad sources
badsrc_idx = find(MNE_badsrc==1);
scatter3(sourcepos(badsrc_idx,1), sourcepos(badsrc_idx,2), sourcepos(badsrc_idx,3), 250, 'r', 'filled')
badsrc_idx = find(hbf_badsrc==1);
scatter3(sourcepos(badsrc_idx,1), sourcepos(badsrc_idx,2), sourcepos(badsrc_idx,3), 250, 'gd', 'filled')

% compare the distributions of entries in leadfield matrices
figure
histogram(MNE_G, 100)
hold on
histogram(hbf_G, 100)
legend('MNE', 'hbf')
title('LFM entries distribution')
set(gca, 'FontSize', 20)
% pretty similar 

% Let's now plot a source topoplot of differences between the G matrices
channel = 22;
[ H,hlink ] = hbf2MNE_compareG(MNE_G, hbf_G, 'MNE native', 'hbf', MNE_sourcepos, channel, elecpos); %#ok<*ASGLU>

% The native MNE surface looks pretty problematic...

% Conclusion: meshing is much much more important than the BEM solver. hbf
% solver seems correct and accurate. We will trust its computations. Now
% we just need to set up the meshes correctly for a 4layer model. 

%% hbf_example_EMEG_Sample4.m

% Example on the use of Helsinki BEM Framework kernel:
% 4-shell MEG+EEG model.
% 
% The example uses meshes from Head 1 of (Stenroos and Nummenmaa, 2016),
% remeshed by me based on sample data of SimNIBS toolbox. If you wish to use the
% meshes or model for something, please contact me first.
%
% (c) Matti Stenroos 2017
% clear
%add BEM framework paths
% addpath ./hbf_calc
% addpath ./hbf_mesh


% load example data for meshes, coils, and
% sources.
loaddir=meshpath;


load(fullfile(loaddir,'bmeshes-m40-50-50-60-70.mat'));%bmeshes
bmeshes=bmeshes.meshes;
% this mesh data set is the Test4 data set in Matti's 2016 paper. It has
% 18479 vertices in total for the 5 meshes taken together. Not very fine
% resolution!

disp('Vertex numbers:')
totalvertex = 0;
for ii = 1:length(bmeshes)
    disp(['Surface ', num2str(ii) ': ', num2str(bmeshes{ii}.nop)])
    totalvertex = totalvertex + bmeshes{ii}.nop;
end
disp(['Total vertex number: ', num2str(totalvertex)])

% let's look at these meshes
faces_list = {bmeshes{5}.e, bmeshes{4}.e, bmeshes{3}.e, bmeshes{2}.e, bmeshes{1}.e};
vertices_list = {bmeshes{5}.p, bmeshes{4}.p, bmeshes{3}.p, bmeshes{2}.p, bmeshes{1}.p};
mesh_plot_bem_surfaces(faces_list, vertices_list);
% fine - looks reasonal, just at very low resolution.


% we probably won't need these coils, since we are not using this solver
% for MEG but only for EEG.
load(fullfile(loaddir,'coils-cmeg306.mat'));%coils

% load(strcat(loaddir,'geometry/sources-sfs35'));%sources
% sourcepos=sources.smeshes{1}.p;sourcedir=sources.smeshes{1}.nn;
load(fullfile(loaddir,'sources-s30d15-p40d15.mat'));%sources
% so... I only have a different source file from the one called by the
% original sample script. But I think I got what's needed. In this new
% source file Matti combined the two hemispheres together.
sourcepos = sources.spos{1}; 
sourcedir = sources.sdir{1};

% Nope... this is actually only the left half of the head. Ok.... fine.
figure
scatter3(sourcepos(:,1), sourcepos(:,2), sourcepos(:,3), 20, 'k', 'filled')
axis equal

% load(fullfile(loaddir,'geometry/elec256_s70.mat'));%electrodes
% so... I don't have this file. But I can guess what the electrodes look
% like, maybe just like in the 3shell example script. So we will instead
% load those electrodes and use them for this one.
load(fullfile(meshpath, 'samplehead_cbu.mat'), 'eegsens');

%Set EEG electrodes
% elecpos=elecs.p;
elecpos=eegsens.p;
%Project electrodes to scalp and build interpolation operator from model
%nodes to projected electrodes.
elecs=hbf_ProjectElectrodesToScalp(elecpos,bmeshes);



% these check functions are same as before, we should be able to pass them
% with ease. If not, I'm in huge trouble.
Nmeshes=length(bmeshes);
success=zeros(Nmeshes,1);
for M=1:Nmeshes
    [bmeshes{M},success(M)]=hbf_CorrectTriangleOrientation(bmeshes{M});
end

status=cell(Nmeshes,1);
for M=1:Nmeshes
    [status{M}]=hbf_CheckMesh(bmeshes{M});
end


% set conductivities 
cratio=50;
c_brain=.33;c_csf=1.79;c_skull=c_brain/cratio;c_skin=.33;
ci=[c_brain,c_brain,c_csf,c_skull,c_skin];%conductivity inside each surface --- remember the order!
co=[c_csf,c_csf,c_skull,c_skin,0];%conductivity outside each surface --- remember the order!


% ------------------------------------------------------------------------
% Hey, ho, let's go!

starttime=clock;
% Build double-layer matrices for the BEM
D=hbf_BEMOperatorsPhi_LC(bmeshes);
% this is the same function as in the 3shell script.


% % turning this off since we won't need MEG
% % Build matrices for integrating the magnetic field due to volume currents 
% DB=hbf_BEMOperatorsB_Linear(bmeshes,coils);


% Build BEM transfer matrix for potential using the isolated source
% approach, isolation set to surface 3 = inner skull
Tphi_full=hbf_TM_Phi_LC_ISA2(D,ci,co,3);
% this is the same function as in the 3shell script.
% which surface do we want to set the ISA to needs some thoughts. I don't
% know enough about ISA approach, maybe I should learn it? 


%Extract/interpolate EEG transfer matrix for electrodes only
Tphi_elecs=hbf_InterpolateTfullToElectrodes(Tphi_full,bmeshes,elecs);
% this is the same function as in the 3shell script.


% % turning this off since we won't need MEG
% % ...and make a transfer matrix for the magnetic field generated by the
% % volume currents.
% TBvol=hbf_TM_Bvol_Linear(DB,Tphi_full,ci,co);


% compute lead field matrices
% LFMstarttime=clock;
% LFM_Edir=hbf_LFM_Phi_LC_dir(bmeshes,Tphi_elecs,sourcepos,sourcedir);
% this is a slightly different function, but it's essentially called by the
% wrapper function used in the 3shell script hbf_LFM_LC()
% Update: somehow calling this function didn't work... Let's try using hbf_LFM_LC()
LFMphi_dir=hbf_LFM_LC(bmeshes,Tphi_elecs,sourcepos,sourcedir);


% % turning this off since we won't need MEG
% LFM_Mdir=hbf_LFM_B_LC_dir(bmeshes,coils,TBvol,sourcepos,sourcedir);
fprintf('Time for LFM computation was %ds.\n',round(etime(clock,starttime)));

%% Examine LFM in closer details
% Well. Let's try to examine the lead field matrix in closer details and
% visualize the loading plots. This should give us confidence that the
% solver is actually working.

% Let's examine the leadfield matrix
G = LFMphi_dir;

% check the histogram of G matrix entries
figure;
subplot(1,2,1)
histogram(G(:), 100)
title('Before detect_badsource()', 'FontSize', 16, 'interpreter', 'none')

[ G, badsrc, badety] = detect_badsource(G);
disp(['Number of bad sources found: ' num2str(sum(badsrc))])
disp(['Number of bad entries found: ' num2str(sum(badety))])

subplot(1,2,2)
histogram(G(:), 100)
title('After detect_badsource()', 'FontSize', 16, 'interpreter', 'none')

% Forward model variables
forward_model.G = G;
forward_model.source = sourcepos;
forward_model.dig = elecs.pproj;
forward_model.normal = sourcedir;
forward_model.skin_surf_face = bmeshes{end}.e;
forward_model.skin_surf_vertice = bmeshes{end}.p;

% Loading plot of a channel 
plot_loading(10, forward_model, false, true, false);

% visualize the bad sources
badsrc_idx = find(badsrc==1);
scatter3(sourcepos(badsrc_idx,1), sourcepos(badsrc_idx,2), sourcepos(badsrc_idx,3), 300, 'r', 'filled')

% looks very reasonable! Maybe it's because the meshes were good? I'm not
% sure. Now let's build it with our own mesh! 

%% Build the four-shell forward model with our own mesh

%add BEM framework paths
% addpath ./hbf_calc
% addpath ./hbf_mesh

load(fullfile(datapath, '4shell_hbfcmpt_BEM_model.mat'),...
    'bmeshes', 'dig', 'sources');%bmeshes, eegsens, sources

disp('Vertex numbers:')
totalvertex = 0;
for ii = 1:length(bmeshes)
    disp(['Surface ', num2str(ii) ': ', num2str(bmeshes{ii}.nop)])
    totalvertex = totalvertex + bmeshes{ii}.nop;
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

%% Examine LFM in closer details
% Well. Let's try to examine the lead field matrix in closer details and
% visualize the loading plots. This should give us confidence that the
% solver is actually working.

% Let's examine the leadfield matrix
G = LFMphi_dir;
xyzG = LFMphi_xyz;

% check the histogram of G matrix entries
figure;
subplot(1,2,1)
histogram(G(:), 100)
title('Before detect_badsource()', 'FontSize', 16, 'interpreter', 'none')

[ G, badsrc, badety, xyzG ] = detect_badsource(G, sourcepos, true, xyzG);
disp(['Number of bad sources found: ' num2str(sum(badsrc))])
disp(['Number of bad entries found: ' num2str(sum(badety))])

subplot(1,2,2)
histogram(G(:), 100)
title('After detect_badsource()', 'FontSize', 16, 'interpreter', 'none')

% Forward model variables
forward_model.G = G;
forward_model.xyzG = xyzG;
forward_model.source = sourcepos;
forward_model.dig = elecs.pproj;
forward_model.normal = sourcedir;
forward_model.skin_surf_face = bmeshes{end}.e;
forward_model.skin_surf_vertice = bmeshes{end}.p;

% Loading plot of a channel 
plot_loading(78, forward_model, false, true, false);

% visualize the bad sources
badsrc_idx = find(badsrc==1);
scatter3(sourcepos(badsrc_idx,1), sourcepos(badsrc_idx,2), sourcepos(badsrc_idx,3), 300, 'r', 'filled')

%%
% Now. What it looks like is that there are some artifact source points
% near the base of the scalp on the medial surface of the temporal lobe
% because that's when the pial surface accuracy is most messed up. It is
% not worth writing an automatic algorithm to detect this becasue it has to
% be judged on a case-by-base basis. What I will propose here instead is a
% manual checking to identify erroneous source poitns and interpolate them
% using mean of near-by-source points. 

% Due to the nature of these source points, they are most easily detectable
% when you do the loading plot of extreme electrodes and at top central
% electrodes. Hence we will iteratively search through these and use
% histograms to guide our marking of source points to interpolate.

% we can clearly see these points using electrode 47 or 1 (mastoid
% electrodes)
plot_loading(47, forward_model, false, true, false);

figure; 
histogram(G(47,:), 50)
figure; 
histogram(G(47,:).^2, 50)

% from the histogram, it looks like some points are abnormal. Let's find
% these points.
[minv, idx] = min(G(47,:));
scatter3(sourcepos(idx,1), sourcepos(idx,2), sourcepos(idx,3), 300, 'r', 'filled')

% use the squared entries to find outliers 
high_src_p = find(G(47,:).^2 > 5*10^3);

scatter3(sourcepos(high_src_p,1), sourcepos(high_src_p,2), sourcepos(high_src_p,3), 300, 'r', 'filled')
% ok great! So high_src_p are the source points that we should interpolate!
% 

[ Gprime, xyzGprime ] = interpolate_sources(G, sourcepos, high_src_p, xyzG);

% let's see how it looks.
forward_model.G = Gprime;
plot_loading(47, forward_model, false, true, false);

% there might still be some issues but looks much better
plot_loading(1, forward_model, false, true, false);

figure; 
histogram(Gprime(1,:).^2, 50)

% use the squared entries to find outliers 
high_src_p = find(Gprime(1,:).^2 > 4700);

scatter3(sourcepos(high_src_p,1), sourcepos(high_src_p,2), sourcepos(high_src_p,3), 300, 'r', 'filled')

[ Gprime, xyzGprime ] = interpolate_sources(Gprime, sourcepos, high_src_p, xyzGprime);

% let's see how it looks.
forward_model.G = Gprime;
plot_loading(1, forward_model, false, true, false);

% let's look through other electrodes
plot_loading(54, forward_model, false, true, false);

% ok - looking good. Update G;
G = Gprime;
xyzG = xyzGprime;
forward_model.G = G;
forward_model.xyzG = xyzG;

% in theory we could have done the same with patch decomposition, but i
% think patch decomposition sacrifices the number of sources so it's better
% to manually identify the very few odd source points and preserve the
% source resolution. If we need to patch decompose later, that's fine.

% Let's save all the data at this point.
save(fullfile(datapath, '4shell_hbfBEM-fwd.mat'),...
    'bmeshes', 'dig', 'sources', 'G', 'xyzG', 'forward_model')

%% Comparison of 3shell versus 4shell hbf BEM models 

% load the forward models 
load(fullfile(datapath, '3shell_hbfBEM-fwd.mat'),...
    'bmeshes', 'dig', 'sources', 'G', 'forward_model')
hbf3shell_G = G;
load(fullfile(datapath, '4shell_hbfBEM-fwd.mat'),...
    'bmeshes', 'dig', 'sources', 'G', 'forward_model')
hbf4shell_G = G;

sourcepos=[sources{1}.points; sources{2}.points];%source locations
elecpos=dig;%electrode positions

% compare the distributions of entries in leadfield matrices
figure
histogram(hbf3shell_G, 100)
hold on
histogram(hbf4shell_G, 100)
legend('hbf3shell', 'hbf4shell')
title('LFM entries distribution')
set(gca, 'FontSize', 20)
% pretty similar 

% Let's now plot a source topoplot of differences between the G matrices
channel = 78;
[ H,hlink ] = hbf2MNE_compareG(hbf3shell_G, hbf4shell_G, 'hbf3shell', 'hbf4shell', sourcepos, channel, elecpos); 
% ok! an overall reduction in magnitude as expected 
% let's visualize the relative distributions as well
[ H,hlink ] = hbf2MNE_compareG(hbf3shell_G, hbf4shell_G, 'hbf3shell', 'hbf4shell', sourcepos, channel, elecpos, false); 


% calculate the difference matrix
diffG = hbf3shell_G - hbf4shell_G;

% plot a histogram
figure
histogram(diffG(:), 100)
title('Distribution of differences in G matrices')
hold on
L = xline(std(diffG(:))*2, 'r--', 'LineWidth', 1);
xline(-std(diffG(:))*2, 'r--', 'LineWidth', 1);
legend([L], {'Two-standard-deviation cut-off'})
set(gca, 'FontSize', 14)

disp(['Total number of entries with large difference: ',...
    num2str(sum(abs(diffG(:)) > std(diffG(:))*2))])
disp(['Different entry percentage out of all entries: ',...
    num2str(sum(abs(diffG(:)) > std(diffG(:))*2) / length(diffG(:))*100), '%'])

% plot the difference in G on sources in 3D space
figure
weight = diffG(channel, 1:size(sourcepos,1)).^2;
size_weight = ceil((weight-min(weight))/(max(weight)-min(weight))*300)+1;
scatter3(sourcepos(:,1), sourcepos(:,2), sourcepos(:,3), size_weight, weight, 'filled')
axis equal
hold on
scatter3(dig(channel,1), dig(channel,2), dig(channel,3), 300, [1,0,0], 's', 'filled')
colorbar
title('hbf3shell - hbf4shell Weighting', 'FontSize', 16)

%% Completed!




