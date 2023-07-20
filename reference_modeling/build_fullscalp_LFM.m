%% Helsinki BEM Framework LCISA Solver (Matti Stenroos) - FULL SCALP LFM
% this script is used to build a full scalp forward model using the tested
% hbf package. The use of the full scalp forwad model will be for
% investigating the effect of re-referencing on scalp potential
% distribution. This is NOT for the use of inverse solution since we never
% have that many channels.

close all
clear all

% Change the current folder to the folder of this m-file.
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
clearvars tmp

%% Path configuration

eeglab

% add path to this EEG/code/forward_model folder
addpath(genpath('./'));

% add path to reference modeling folder
addpath('/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/reference_modeling/code/')
% 
% create datapath for our own BEM meshes
datapath = '/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/forward_modeling/4layer_BEM/sample_mesh/pipelines/bem_surfaces/final_structure';

%% Build a four-shell forward model with our own mesh with 1000 electrodes

load('4shell_hbfcmpt_BEM_model.mat',...
    'bmeshes', 'dig', 'sources');%bmeshes, sources. NO eegsens since using skin!

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

% construct a full scalp electrode distribution
% elecpos=dig;
figure
hold on
P = patch('Faces',bmeshes{end}.e,'Vertices',bmeshes{end}.p,'facecolor',[.5 .5 .5],'edgecolor','none');
camlight('headlight','infinite')
axis equal
rotate3d on
scatter3(dig(:,1), dig(:,2), dig(:,3), 100, 'rs', 'filled')

faces_list = {bmeshes{1}.e};
vertices_list = {bmeshes{1}.p};
mesh_plot_bem_surfaces(faces_list, vertices_list, sources);

% Ok, we need to define a tilted plane based on the electrodes. 
% step 1: find midpoint of the two most anteriror and inferior electrodes:
% 2 and 44.

% defining a tilted plane isn't so easy. Let's first set point 92 to be the
% new origin. Update: rather than setting 92, let's set midpoint of 1 and
% 47. This has the advantage that we can still define the mastoid points in
% the fullscalp electrode space! 
og_dig = dig;
old_92 = mean(dig([1,47],:));
dig = dig - mean(dig([1,47],:));
ant_p = mean(dig([2,44],:));

% figure
% hold on
% scatter3(dig(:,1), dig(:,2), dig(:,3), 100, 'rs', 'filled')
% scatter3(ant_p(:,1), ant_p(:,2), ant_p(:,3), 300, 'g', 'filled')
% scatter3(pos_p(:,1), pos_p(:,2), pos_p(:,3), 300, 'g', 'filled')
% axis equal
% rotate3d on

% then we define two vectors to span a plane 
v1 = ant_p ./ norm(ant_p);
v2 = dig(1,:) - dig(47,:); % between the two mastoids;
v2 = v2 ./ norm(v2);
% 
% mArrow3([0,0,0],v1./10,'color','blue', 'tipWidth', 0.002, 'facealpha',0.5);
% mArrow3([0,0,0],v2./10,'color','blue', 'tipWidth', 0.002, 'facealpha',0.5);

% now define the four ends of a plane 
miny = min(dig(:,2));
lowz = dot([0, miny, 0], v1)*1.5*v1;
% scatter3(lowz(:,1), lowz(:,2), lowz(:,3), 300, 'm', 'filled')
spanv = range(dig(:,1))/2*1.2;
low_left = lowz + v2*spanv; 
low_right = lowz - v2*spanv;
% scatter3(low_left(:,1), low_left(:,2), low_left(:,3), 300, 'm', 'filled')
% scatter3(low_right(:,1), low_right(:,2), low_right(:,3), 300, 'm', 'filled')
maxy = max(dig(:,2));
highz = dot([0,maxy,0],v1)*2*v1;
% scatter3(highz(:,1), highz(:,2), highz(:,3), 300, 'm', 'filled')
high_left = highz + v2*spanv;
high_right = highz - v2*spanv;
% scatter3(high_left(:,1), high_left(:,2), high_left(:,3), 300, 'm', 'filled')
% scatter3(high_right(:,1), high_right(:,2), high_right(:,3), 300, 'm', 'filled')

% translate these points to original coordiante 
low_left = low_left + old_92;
low_right = low_right + old_92;
high_left = high_left + old_92;
high_right = high_right + old_92;
dig = og_dig;
plane_ends = [low_left; low_right; high_right; high_left];

% visualize the plane for cutting 
patch(plane_ends(:,1), plane_ends(:,2), plane_ends(:,3), [0,0,0,0],'facecolor',[0 0 1], 'facealpha', 1)

% looks good! Now we need to decide which vertices to include!

% to do that, we need to define a normal vector that is point up
plane_normal = cross(v1, v2);

mArrow3(low_right,low_right+plane_normal./10,'color','blue', 'tipWidth', 0.002, 'facealpha',0.5);

up_index = zeros(length(bmeshes{end}.p), 1);
for ii = 1:length(up_index)
    if dot(bmeshes{end}.p(ii,:) - mean(dig([1,47],:)), plane_normal) >= -0.005 % we over grab just a little bit
        up_index(ii) = 1;
    end
end

% % get all vertex of the skin surface above the plane 
% elecpos = bmeshes{end}.p(logical(up_index),:);

% % visualize these vertices
% scatter3(elecpos(:,1), elecpos(:,2), elecpos(:,3), 10, 'y', 'filled')

% looks good! Now we want to remesh this to a denser resolution.
% first let's grab all the faces 
faces_index = [];
vertices_index = find(up_index == 1);
for ii = 1:length(bmeshes{end}.e)
    if any(ismember(bmeshes{end}.e(ii,:), vertices_index))
        faces_index = [faces_index, ii];
    end
end
elec_faces = bmeshes{end}.e(faces_index,:);
vertices_index = unique(elec_faces);

% figure
% hold on
% P = patch('Faces',elec_faces,'Vertices',bmeshes{end}.p,'facecolor',[.5 .5 .5],'edgecolor','none');
% camlight('headlight','infinite')
% axis equal
% rotate3d on
% scatter3(dig(:,1), dig(:,2), dig(:,3), 100, 'rs', 'filled')
% scatter3(elecpos(:,1), elecpos(:,2), elecpos(:,3), 10, 'y', 'filled')

% update face_indices
for ii = 1:length(elec_faces)
    [~, elec_faces(ii,:)] = ismember(elec_faces(ii,:), vertices_index);
end
elec_vertices = bmeshes{end}.p(vertices_index,:);

% figure
% hold on
% P = patch('Faces',elec_faces,'Vertices',elec_vertices,'facecolor',[.5 .5 .5],'edgecolor','none');
% camlight('headlight','infinite')
% axis equal
% rotate3d on
% scatter3(dig(:,1), dig(:,2), dig(:,3), 100, 'rs', 'filled')
% scatter3(elecpos(:,1), elecpos(:,2), elecpos(:,3), 10, 'y', 'filled')

% now we want to remesh, first we need to upsample
%subdivide to upsample
nsub = 1; % number of subdivision steps
options.sub_type = 'loop';
options.verb = 0;
[vertices,faces] = perform_mesh_subdivision(elec_vertices',elec_faces',nsub,options);
vertices = vertices'; faces = faces';

% figure
% hold on
% P = patch('Faces',faces,'Vertices',vertices,'facecolor',[.5 .5 .5],'edgecolor','none');
% camlight('headlight','infinite')
% axis equal
% rotate3d on
% scatter3(dig(:,1), dig(:,2), dig(:,3), 100, 'rs', 'filled')
% scatter3(elecpos(:,1), elecpos(:,2), elecpos(:,3), 10, 'y', 'filled')

% remesh to 1000 equidistant points
indStart=1; %Index of the start point
numSeeds=1000; %number of vertices in the new mesh
optionStruct.toleranceLevel=0; %Tolerance for convergence
optionStruct.waitBarOn=1; %Turn on/off waitbar
[ faces,vertices ]=remeshTriSurfDistMap(faces,vertices,numSeeds,indStart,optionStruct);

figure
hold on
P = patch('Faces',faces,'Vertices',vertices,'facecolor',[.5 .5 .5],'edgecolor','none');
camlight('headlight','infinite')
axis equal
rotate3d on
scatter3(dig(:,1), dig(:,2), dig(:,3), 100, 'rs', 'filled')
scatter3(vertices(:,1), vertices(:,2), vertices(:,3), 10, 'y', 'filled')

% ok! vertices are our new electrode positions! Great!

% project electrodes to scalp and make "elecs" struct
elecs=hbf_ProjectElectrodesToScalp(vertices,bmeshes);

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

save('4shell_fullscalp_computed_LFM', 'LFMphi_dir', 'LFMphi_xyz')


% load the computed LFM structures
load('4shell_fullscalp_computed_LFM')
LFMphi_dir_fullscalp = LFMphi_dir;
LFMphi_xyz_fullscalp = LFMphi_xyz;
load('4shell_computed_LFM')

%% Examine LFM in closer details
% Well. Let's try to examine the lead field matrix in closer details and
% visualize the loading plots. This should give us confidence that the
% solver is actually working.

% Let's examine the leadfield matrix
G = LFMphi_dir;
xyzG = LFMphi_xyz;
G_fullscalp = LFMphi_dir_fullscalp;
xyzG_fullscap = LFMphi_xyz_fullscalp;

% check the histogram of G matrix entries
% let's check that the same bad sources will be picked up in the two
% forward models. Detecting badentries is extremely time consuming on
% fullscalp because we have 2032x20484 = 41623488 data points to run
% iterative algorithms on. Instead, we will detect and interpolate bad
% entries and then manually feed in the badentry index to interpolate in
% the fullscalp matrix
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

% Now, let's fix G_fullscalp
figure;
subplot(1,2,1)
histogram(G_fullscalp(:), 100)
title('Before detect_badsource()', 'FontSize', 16, 'interpreter', 'none')

[ G_fullscalp, badsrc_fullscalp, badety_fullscalp, xyzG_fullscap ] = detect_badsource(G_fullscalp, sourcepos, true, xyzG_fullscap, false);
disp(['Number of bad sources found: ' num2str(sum(badsrc))])
assert(all(badsrc == badsrc_fullscalp), 'different bad sources were detected for the two scalp resolutions!')
assert(isempty(badety_fullscalp), 'somehow badety_fullscalp is not empty!')

% let's build cutoff values based off G entries on G_fullscalp entries
cutoff_range = [min(G(:)), max(G(:))];
badety_fullscalp = logical(G_fullscalp(:)<cutoff_range(1)) | logical(G_fullscalp(:)>cutoff_range(2));
disp(['Number of bad entries found: ' num2str(sum(badety_fullscalp))])

% now let's manually interpolate these bad entries
[ G_fullscalp, xyzG_fullscap ] = interpolate_entries(G_fullscalp, sourcepos, badety_fullscalp, xyzG_fullscap);

subplot(1,2,2)
histogram(G_fullscalp(:), 100)
title('After detect_badsource()', 'FontSize', 16, 'interpreter', 'none')

% compare the distributions
figure;
hold on
histogram(G(:), 100, 'Normalization', 'pdf')
histogram(G_fullscalp(:), 100, 'Normalization', 'pdf')
title('Compare distributions of fullscalp and duke128', 'FontSize', 16)
% looking good!!



% Forward model variables
forward_model = struct;
forward_model.G = G_fullscalp;
forward_model.xyzG = xyzG_fullscap;
forward_model.source = sourcepos;
forward_model.dig = elecs.pproj;
forward_model.normal = sourcedir;
forward_model.skin_surf_face = bmeshes{end}.e;
forward_model.skin_surf_vertice = bmeshes{end}.p;

% Loading plot of a channel 
plot_loading(580, forward_model, false, true, false);

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

% Again, since G_fullscalp has too many electrodes to handle, we will work
% on G and then interpolate to G_fullscalp. 
forward_model.G = G;
dig128=hbf_ProjectElectrodesToScalp(dig,bmeshes);
forward_model.dig = dig128.pproj;

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
[ G, xyzG ] = interpolate_sources(G, sourcepos, high_src_p, xyzG);
% we will then use the same source indices to fix G_fullscalp and xyzG_fullscap
[ G_fullscalp, xyzG_fullscap ] = interpolate_sources(G_fullscalp, sourcepos, high_src_p, xyzG_fullscap);

% let's see how it looks.
forward_model.G = G;
plot_loading(47, forward_model, false, true, false);




% there might still be some issues but looks much better
% we can use the other mastoid to check for strange sources.
plot_loading(1, forward_model, false, true, false);

figure; 
histogram(G(1,:).^2, 50)

% use the squared entries to find outliers 
high_src_p = find(G(1,:).^2 > 4700);

scatter3(sourcepos(high_src_p,1), sourcepos(high_src_p,2), sourcepos(high_src_p,3), 300, 'r', 'filled')

[ G, xyzG ] = interpolate_sources(G, sourcepos, high_src_p, xyzG);
% we will then use the same source indices to fix G_fullscalp and xyzG_fullscap
[ G_fullscalp, xyzG_fullscap ] = interpolate_sources(G_fullscalp, sourcepos, high_src_p, xyzG_fullscap);

% let's see how it looks.
forward_model.G = G;
plot_loading(1, forward_model, false, true, false);




% when mastoid is selected on one side, there should be minimal
% contribution from sources on the contralateral side. We use this
% knowledge to further identify and interpolate abnormal sources 
% Note - this should probably be written as a function and be called
% consistently for all forward model lead field matrices
%   LM: 1
channel = 1;
plot_loading(channel, forward_model, false, true, false);
right_src_idx = find(sourcepos(:,1) > 0);
scatter3(sourcepos(right_src_idx, 1), sourcepos(right_src_idx,2), sourcepos(right_src_idx,3), 10,'g', 'filled')
[ ~, max_idx ] = sort(G(channel,right_src_idx).^2, 'descend');
first20 = right_src_idx(max_idx(1:20));
scatter3(sourcepos(first20, 1), sourcepos(first20,2), sourcepos(first20,3), 50, 'r', 'filled')
% interpolate these points 
[ G, xyzG ] = interpolate_sources(G, sourcepos, first20, xyzG);
% we will then use the same source indices to fix G_fullscalp and xyzG_fullscap
[ G_fullscalp, xyzG_fullscap ] = interpolate_sources(G_fullscalp, sourcepos, high_src_p, xyzG_fullscap);
% look at it after interpolation
forward_model.G = G;
plot_loading(channel, forward_model, false, true, false);

%   RM: 47
channel = 47;
plot_loading(channel, forward_model, false, true, false);
left_src_idx = find(sourcepos(:,1) < 0);
scatter3(sourcepos(left_src_idx, 1), sourcepos(left_src_idx,2), sourcepos(left_src_idx,3), 10,'g', 'filled')
[ sorted_list, max_idx ] = sort(G(channel,left_src_idx).^2, 'descend');
first20 = left_src_idx(max_idx(1:20));
scatter3(sourcepos(first20, 1), sourcepos(first20,2), sourcepos(first20,3), 50, 'r', 'filled')
% interpolate these points 
[ G, xyzG ] = interpolate_sources(G, sourcepos, first20, xyzG);
% we will then use the same source indices to fix G_fullscalp and xyzG_fullscap
[ G_fullscalp, xyzG_fullscap ] = interpolate_sources(G_fullscalp, sourcepos, high_src_p, xyzG_fullscap);
% look at it after interpolation
forward_model.G = G;
plot_loading(channel, forward_model, false, true, false);

%% ok - everything looking good. 
% Ok, now that we have fixed sources according to mastoid electrodes on G,
% we will try to look around through other scalp vertices to see if we see
% weird spots that are worth fixing.
forward_model.G = G_fullscalp;
forward_model.xyzG = xyzG_fullscap;
forward_model.dig = elecs.pproj;

% let's look through other electrodes
plot_loading(315, forward_model, false, true, false);

% fine, looks ok
dig = elecs.pproj; % update dig with 1000 electrodes

% Let's save all the data at this point.
save(fullfile(datapath, '4shell_hbfBEM_fullscalp-fwd.mat'),...
    'bmeshes', 'dig', 'sources', 'G_fullscalp', 'xyzG_fullscap', 'forward_model')

% Ok great. We have the lead field matrix corrected. 

%% Now we just need to save everything into a struct
clearvars forward_model

% load the 128-duke configuration as a forward model struct
filename = fullfile(datapath, '4shell_hbfBEM-fwd.mat');
fsfn = fullfile(datapath, 'TDelp_resting_fastscan_dig.mat');
atlasfn = fullfile(datapath, 'atlas_info.mat');

forward_model = prep_forward_model_struct(filename, fsfn, atlasfn);

% now, the only remaining obstacle is that we don't have the 1000
% electrodes on the scalp in the coordinate system of EEGLAB topoplot().
% Hence we don't really have a fsfn for '4shell_hbfBEM_fullscalp-fwd.mat'
% that we just saved. Now we need to somehow figure out the transformation
% between the two coordinate systems and apply that to the 1000 scalp
% electrodes in the _fullscalp forward model. 

% let's first visualize the problem
elc_fn = 'TDelp_resting-elc-trans.mat';
% not projecting since we are trying to find the transformation between coordinate systems
MRI_dig = mesh_getelectrode(fullfile(datapath,elc_fn), [], [], false); 
MRI_dig = MRI_dig/1000;

EEGLAB_dig = [];
for ii = 1:length(MRI_dig)
    EEGLAB_dig(ii,1) = forward_model.eloc(ii).X/1000;
    EEGLAB_dig(ii,2) = forward_model.eloc(ii).Y/1000;
    EEGLAB_dig(ii,3) = forward_model.eloc(ii).Z/1000;
end

figure
hold on
scatter3(MRI_dig(:,1), MRI_dig(:,2), MRI_dig(:,3), 50, 'r', 'filled')
scatter3(EEGLAB_dig(:,1), EEGLAB_dig(:,2), EEGLAB_dig(:,3), 50, 'b', 'filled')
P = patch('Faces',forward_model.skin_surf_face,'Vertices',forward_model.skin_surf_vertex,...
    'facecolor',[.5 .5 .5],'edgecolor','none');
camlight('headlight','infinite')
axis equal
rotate3d on

% NOPE! Ok I give up on trying to back-engineering what MNE did to the
% coordinate system. Let me use ICP to find the exact transformation
% matrix. Ok - ICP algorithms suck. Let me write my own. 

% this function does a OK job in registering the two sets of electrodes not
% perfect though.
[ TR, trans_v, final_cost ] = register_electrode_aprx(MRI_dig, EEGLAB_dig, pi/2, 5);

% let's apply the same transformation to the skin surface points 
skin_v_rot = forward_model.skin_surf_vertex * TR + trans_v;
figure
hold on
scatter3(MRI_dig(:,1), MRI_dig(:,2), MRI_dig(:,3), 50, 'r', 'filled')
scatter3(EEGLAB_dig(:,1), EEGLAB_dig(:,2), EEGLAB_dig(:,3), 50, 'b', 'filled')
P = patch('Faces',forward_model.skin_surf_face,'Vertices',skin_v_rot,...
    'facecolor',[.5 .5 .5],'edgecolor','none');
camlight('headlight','infinite')
axis equal
rotate3d on


% looking good!
% let's apply the transformation to the 1000 electrodes 
dig_rot = dig * TR + trans_v;

figure
hold on
scatter3(dig_rot(:,1), dig_rot(:,2), dig_rot(:,3), 20, 'k', 'filled')
scatter3(EEGLAB_dig(:,1), EEGLAB_dig(:,2), EEGLAB_dig(:,3), 50, 'b', 'filled')
axis equal
rotate3d on


% Now let's create an eloc struct for these 1000 electrodes
eloc = struct;
for i = 1:size(dig_rot, 1)
    eloc(i).X = dig_rot(i, 1)*1000; % convert to mm unit
    eloc(i).Y = dig_rot(i, 2)*1000;
    eloc(i).Z = dig_rot(i, 3)*1000;
    eloc(i).labels = sprintf('EEG%04.f', i);
end
eloc = convertlocs(eloc, 'cart2all');

figure
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1:2) pos(3:4)*1.5])
topoplot([],eloc,'style','both','electrodes','on','emarker', {'.', 'k', 10, 1});
% looks good

%% Find corresponding points on 1000-electrode config for referencing 
% we need to define mastoid points and Z3 for creating the reference
% matrices. Also we need to define a neighbor list for laplacian
% referencing. We have to do this individually, as remeshing of the skin
% surface may result in different electrode order

duke128_dig = forward_model.dig;
e1000_dig = dig;

figure
hold on
scatter3(e1000_dig(:,1), e1000_dig(:,2), e1000_dig(:,3), 20, 'k', 'filled')
scatter3(duke128_dig(:,1), duke128_dig(:,2), duke128_dig(:,3), 50, 'b', 'filled')
axis equal
rotate3d on

% define e1000_85 
[ ~, e1000_85 ] = min(vecnorm(e1000_dig - duke128_dig(85,:), 2, 2));

scatter3(e1000_dig(e1000_85,1), e1000_dig(e1000_85,2), e1000_dig(e1000_85,3), 100, 'r', 'filled')
scatter3(duke128_dig(85,1), duke128_dig(85,2), duke128_dig(85,3), 100, 'g', 'filled')

% define e1000_1 and e1000_47
[ ~, e1000_1 ] = min(vecnorm(e1000_dig - duke128_dig(1,:), 2, 2));
[ ~, e1000_47 ] = min(vecnorm(e1000_dig - duke128_dig(47,:), 2, 2));

scatter3(e1000_dig(e1000_1,1), e1000_dig(e1000_1,2), e1000_dig(e1000_1,3), 100, 'r', 'filled')
scatter3(duke128_dig(1,1), duke128_dig(1,2), duke128_dig(1,3), 100, 'g', 'filled')
scatter3(e1000_dig(e1000_47,1), e1000_dig(e1000_47,2), e1000_dig(e1000_47,3), 100, 'r', 'filled')
scatter3(duke128_dig(47,1), duke128_dig(47,2), duke128_dig(47,3), 100, 'g', 'filled')

% define channelneighbors
channelneighbors = zeros(length(e1000_dig));
for ii = 1:size(channelneighbors, 1)
    neighbor_face = find(any([faces(:,1)==ii, faces(:,2)==ii, faces(:,3)==ii], 2));
    neighbors = unique(faces(neighbor_face,:));
    neighbors(neighbors==ii) = []; % get rid of oneself;
    channelneighbors(ii, neighbors) = 1;
end
assert(issymmetric(channelneighbors), 'channelneighbors is not symmetric! Something went wrong!')

% define all duke 128 electrodes
e1000_2_duke = zeros(1, size(duke128_dig,1));
for ii = 1:length(e1000_2_duke)
    [ ~, e1000_2_duke(ii) ] = min(vecnorm(e1000_dig - duke128_dig(ii,:), 2, 2));
end

% summarize everything into a struct
e1000_info = struct;
e1000_info.e1000_to_duke128 = e1000_2_duke;
e1000_info.channelneighbors = channelneighbors;
e1000_info.e1000_85 = e1000_85;
e1000_info.e1000_1 = e1000_1;
e1000_info.e1000_47 = e1000_47;
e1000_info.faces = faces;


%% We are now ready to prepare the forward_model struct
% we can't use the prep_forward_model_struct() function due to different
% inputs are required. We will implement everything manually here.

%% Compute reference matrices
% we will compute the different R matrices for different referencing
% schemes
% Unipolar references:
% R1 - common average
% R2 - single electrode (recording reference at Z3)
% R3 - single mastoid (left mastoid)
% R4 - linked mastoid
% R5 - REST

% Non-unipolar references:
% R6 - contralateral mastoid
% R7~R11 - Laplacian

ref_label = {'Common Average',... 
    'Recording Reference at Z3',...
    'Left Mastoid',...
    'Linked Mastoid',...
    'REST',...
    'Contralateral mastoid',...
    '4 neightbor LP',...
    'sphere radius 50mm LP',...
    'Manual LP',...
    'Extended Manual LP',...
    'Extended sphere raiud 100mm LP'};


%%%%%%%%%%%%%%%%%%%%%%
% Unipolar references:
%%%%%%%%%%%%%%%%%%%%%%
% R1 - Common average
R1 = eye(size(G_fullscalp,1)) - ones(size(G_fullscalp,1))*1/size(G_fullscalp,1);

% R2 - Single electrode (recording reference at 85)
ref_chan = e1000_85;
R2 = eye(size(G_fullscalp,1));
R2(:, ref_chan) = R2(:, ref_chan) - 1;

% R3 - Single mastoid (left mastoid 1)
R3 = eye(size(G_fullscalp,1));
R3(:, e1000_1) = R3(:, e1000_1) - 1;

% R4 - Linked mastoid (left mastoid 1 and right mastoid 47)
R4 = eye(size(G_fullscalp,1));
R4(:, [e1000_1, e1000_47]) = R4(:, [e1000_1, e1000_47]) - 1/2;

% R5 - REST
fREST = pinv(G_fullscalp)'*pinv(G_fullscalp)*ones(size(G_fullscalp,1), 1) ./ (ones(size(G_fullscalp,1), 1)'*pinv(G_fullscalp)'*pinv(G_fullscalp)*ones(size(G_fullscalp,1), 1));
R5 = eye(size(G_fullscalp,1)) - ones(size(G_fullscalp,1), 1)*fREST';

% Assert that these are indeed unipolar reference operators
assert(testUR(R1) & testUR(R2) & testUR(R3) & testUR(R4) & testUR(R5), 'Failed unipolar reference test')


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-unipolar references:
%%%%%%%%%%%%%%%%%%%%%%%%%%
% R6 - Contralateral mastoid reference
R6 = zeros(size(G_fullscalp,1));
for i = 1:size(R6,1)
    R6(i, i) = 1;
    if e1000_dig(i, 1) > 0 % right-side electrodes
        R6(i, e1000_1) = -1;
    else
        R6(i, e1000_47) = -1;
    end
end
R6(e1000_1,e1000_1) = 0; R6(e1000_1,e1000_47) = 0; 
R6(e1000_47,e1000_47) = 0; R6(e1000_47,e1000_1) = 0; 

% R7 - Laplacian - 4 nearest neighbors
R7 = zeros(size(G_fullscalp,1));
for i = 1:size(R7,1)
    R7(i, i) = 1;
    distance_list = [];
    % search all electrodes
    for j = 1:size(e1000_dig,1)
        distance_list(j) = norm(e1000_dig(i,:) - e1000_dig(j,:));
    end
    
    % get rid of the sensor itself
    distance_list(i) = inf;
    
    % find the smallest 4 neighbor electrodes
    [~, index] = sort(distance_list);
    R7(i, index(1:4)) = R7(i, index(1:4)) - 1/4;
end


% R8 - Laplacian - sphere radius threshold
R8 = zeros(size(G_fullscalp,1));
for i = 1:size(R8,1)
    R8(i, i) = 1;
    neighbor_list = [];
    % search within radius
    for j = 1:size(e1000_dig,1)
        if norm(e1000_dig(i,:) - e1000_dig(j,:)) <= 0.050 % 50mm / 5cm
            neighbor_list = [neighbor_list, j];
        end
    end
    % get rid of the sensor itself
    neighbor_list(neighbor_list == i) = [];
    R8(i, neighbor_list) = -1 / length(neighbor_list);
    
%     figure
%     hold on
%     scatter3(e1000_dig(:,1),e1000_dig(:,2),e1000_dig(:,3),10,'k','filled')
%     scatter3(e1000_dig(i,1),e1000_dig(i,2),e1000_dig(i,3),30,'r','filled')
%     scatter3(e1000_dig(neighbor_list,1),e1000_dig(neighbor_list,2),e1000_dig(neighbor_list,3),20,'g','filled')
%     axis equal
%     rotate3d on
%     pause()
%     close all
end


% R9 - Laplacian - manual definition based on channelneighbors that we just
% created using face triangulation relationship
R9 = zeros(size(G_fullscalp,1));
for i = 1:size(R9,1)
    R9(i, i) = 1;
    neighbor_list = find(channelneighbors(i,:) == 1);
    R9(i, neighbor_list) = -1 / length(neighbor_list);
end


% R10 - Laplacian - manual neighbors extended by second order
R10 = zeros(size(G_fullscalp,1));
for i = 1:size(R10,1)
    R10(i, i) = 1;
    fo_neighbor_list = find(channelneighbors(i,:) == 1);
    all_neighbor_list = fo_neighbor_list;
    for j = 1:length(fo_neighbor_list)
        all_neighbor_list = [all_neighbor_list, find(channelneighbors(fo_neighbor_list(j),:) == 1)];
    end
    all_neighbor_list = unique(all_neighbor_list);
    % get rid of the sensor itself
    all_neighbor_list(all_neighbor_list==i) = [];
    R10(i, all_neighbor_list) = -1 / length(all_neighbor_list);
end


% R11 - Laplacian - extended sphere radius threshold
R11 = zeros(size(G_fullscalp,1));
for i = 1:size(R11,1)
    R11(i, i) = 1;
    neighbor_list = [];
    % search within radius
    for j = 1:size(e1000_dig,1)
        if norm(e1000_dig(i,:) - e1000_dig(j,:)) <= 0.100 % 100mm / 10cm
            neighbor_list = [neighbor_list, j];
        end
    end
    % get rid of the sensor itself
    neighbor_list(neighbor_list == i) = [];
    R11(i, neighbor_list) = -1 / length(neighbor_list);
    
%     figure
%     hold on
%     scatter3(dig(:,1),dig(:,2),dig(:,3),10,'k','filled')
%     scatter3(dig(i,1),dig(i,2),dig(i,3),30,'r','filled')
%     scatter3(dig(neighbor_list,1),dig(neighbor_list,2),dig(neighbor_list,3),20,'g','filled')
%     axis equal
%     rotate3d on
%     pause()
%     close all
end

% Assert that these are indeed non-unipolar reference operators
assert(testUR(R6, true) & testUR(R7, true) & testUR(R8, true) & testUR(R9, true) & testUR(R10, true) & testUR(R11, true), 'Failed non-unipolar reference test')

% Construct augmented forward model matrices
modified_G = cell(1,11);
Reference_mat = [];
for i = 1:length(modified_G)
    % Save the augmented forward model matrices
    eval(sprintf('modified_G{i} = R%d*G_fullscalp;', i))
    % Create a large matrix with reference transformation matrices
    eval(sprintf('Reference_mat(i,:,:) = R%d;', i))
end


%% Load the anatomical ROI based on Desikan atlas
% load atlas_info structure
% this also loads the cortical surface
atlasfn = fullfile(datapath, 'atlas_info.mat'); % this is the same as the one we used for duke128 layout
load(atlasfn, 'atlas_info', 'fce', 'vtx', 'idx', 'src_face');

% check the atlas_info for accuracy
totalvertex = 0;
for ii = 1:length(atlas_info.atlas_vertidx)
    totalvertex = totalvertex + length(atlas_info.atlas_vertidx{ii});
end
assert(totalvertex < size(G_fullscalp,2), 'Atlas labelled more sources than available! Something is wrong!')

% convert the source face to indexing on the source indices rather than
% vertex indices of the whole cortical surface 
assert(length(idx) == length(vtx), 'Size of select index mismatches the vertex number of cortical surface!')
vertno = find(idx == 1);
[ logic_a, source_face ] = ismember(src_face, vertno);
assert(all(logic_a, 'all'), 'Not all vertices contained in src_face are being used as a source!')


%% Create output structure

fwd = struct;

% lead field matrix variables
fwd.G = G_fullscalp;
fwd.xyzG = xyzG_fullscap;
fwd.source_face = source_face;
fwd.source = sourcepos;
fwd.normal = sourcedir;
fwd.dig = dig;

% cortical surface variables
fwd.cortex_surf_face = double(fce);
fwd.cortex_surf_vertex = vtx;
fwd.select_idx = logical(idx);
% Anatomical ROI source indices 
fwd.anato_ROI = atlas_info;

% skin surface
fwd.skin_surf_face = bmeshes{end}.e;
fwd.skin_surf_vertex = bmeshes{end}.p;

% also save the BEM surfaces
fwd.bmeshes = bmeshes;

% Sensor space variables for EEGLAB topoplot()
fwd.eloc = eloc;

% Augmented LFM with Referencing
fwd.modified_G = modified_G;
fwd.ref_matrix = Reference_mat;
fwd.ref_label = ref_label';

% special electrode configuration info
fwd.dig_type = 'e1000';
fwd.e1000_info = e1000_info;

%% Package this with the forward_model from duke128 configuration
% load the 128-duke configuration as a forward model struct
filename = fullfile(datapath, '4shell_hbfBEM-fwd.mat');
fsfn = fullfile(datapath, 'TDelp_resting_fastscan_dig.mat');
atlasfn = fullfile(datapath, 'atlas_info.mat');

forward_model1 = prep_forward_model_struct(filename, fsfn, atlasfn);
forward_model1.dig_type = 'duke128';

forward_model = struct;
forward_model.duke128 = forward_model1;
forward_model.e1000 = fwd;

% sanity checks for equivalence
assert(all(forward_model.duke128.cortex_surf_face == forward_model.e1000.cortex_surf_face, 'all'), 'equivalence check failed')
assert(all(forward_model.duke128.cortex_surf_vertex == forward_model.e1000.cortex_surf_vertex, 'all'), 'equivalence check failed')
assert(all(forward_model.duke128.select_idx == forward_model.e1000.select_idx, 'all'), 'equivalence check failed')
assert(all(forward_model.duke128.skin_surf_face == forward_model.e1000.skin_surf_face, 'all'), 'equivalence check failed')
assert(all(forward_model.duke128.skin_surf_vertex == forward_model.e1000.skin_surf_vertex, 'all'), 'equivalence check failed')

% save the forward_model struct
save('example_4shell_forward_model.mat', 'forward_model', '-v7.3')

%% Done!


