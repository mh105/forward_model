%% mesh_bem_meshing script
% This script is used to check the quality of different BEM surfaces and
% correct intersections of surfaces if necessary.

% We also check the inclusion of source space within the pial surface and
% check for minimal cortex thickness to be >= 0.5mm. We will point-wise
% shift the source away from pial surface to obtain a minimum thickness of
% 0.5mm.

close all
clear all

% Change the current folder to the folder of this m-file.
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
clearvars tmp

procstart = tic;

%% Path configuration
addpath(genpath('./'));

datapath = '/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/forward_modeling/4layer_BEM/sample_mesh/pipelines/bem_surfaces/final_structure';

%% Mesh loading
% also smoothing with iso2mesh smoothsurf()
% we need to downsample to a coarser resolution after smoothing as well
% with GIBBON codes. We specify the number of vertices we want in the
% output triangulation mesh.

% if the smoothed and downsampled surfaces are already saved, we load them
% directly rather than using the mesh_loadmesh() function, which takes a
% long time!

if isfile(fullfile(datapath, 'archive_4layer_BEM_model.mat'))
    % if the surfaces are already saved to the disk
    load(fullfile(datapath, 'archive_4layer_BEM_model.mat'),...
        'f1','f2','f3','f4','f5','f6',...
        'v1','v2','v3','v4','v5','v6')
    
else
    % for new meshes, we need to use the mesh_loadmesh() function to load.
    % Caution: due to the downsampling step, this section could take up to 4h!
    
    % load the surfaces with smoothing
    skin_fn = 'MNE_outer_skin.mat';
    skull_fn = 'morphed_skull.mat';
    csf_fn = 'mri2mesh_csf_fixed.mat';
    pial_fn = 'mri2mesh_gm_fixed.mat';
    cerebll_fn = 'mri2mesh_cerebellum_fixed.mat';
    wm_fn = 'mri2mesh_wm_fixed.mat';
    
    tic
    [ f1,v1,hlink1 ] = mesh_loadmesh(skin_fn, datapath, 0, [], true);
    [ f2,v2,hlink2 ] = mesh_loadmesh(skull_fn, datapath, 0, 16384, true);
    [ f3,v3,hlink3 ] = mesh_loadmesh(csf_fn, datapath, 100, 16384, true);
    [ f4,v4,hlink4 ] = mesh_loadmesh(pial_fn, datapath, 100, 32768, true);
    [ f5,v5,hlink5 ] = mesh_loadmesh(cerebll_fn, datapath, 50, [], true);
    [ f6,v6,hlink6 ] = mesh_loadmesh(wm_fn, datapath, 100, [], true);
    disp('Total time taken in loading, smoothing and downsampling the meshes:')
    toc
    
    close all
end

%% Remesh BEM surfaces to different vertex resolutions
% sometiems we may want to remesh the surfaces into a smaller resolution to
% test out BEM solvers. we can use the mesh_remeshsurface() function.
% wm/gm boundary (f6, v6) is kept as is since it's not used for any
% computation in MATLAB.

% for skin, skull, csf, pial, cerebellum in this order
vertex_resolution_list = [4096, 8192, 8192, 16384, 8192];

meshcell = {};
for ii = 1:5
    eval(sprintf('surf1.p=v%s; surf1.e=f%s; surf1.nop=size(v%s,1);', num2str(ii), num2str(ii), num2str(ii)))
    meshcell{ii} = surf1;
end

% remesh the surfaces
meshcell = mesh_remeshsurface(meshcell, vertex_resolution_list);

for ii = 1:5
    surf1 = meshcell{ii};
    eval(sprintf('v%s=surf1.p; f%s=surf1.e;', num2str(ii), num2str(ii)))
end

%% Visualize BEM meshes function: mesh_plot_bem_surfaces()
% faces_list = {f4};
% vertices_list = {v4};
% mesh_plot_bem_surfaces(faces_list, vertices_list);

%% Test out the mesh_quality_check function
% let's manually create a hole in the surface and see if
% mesh_quality_check will be able to fix it.
while false
    F = f1; V = v1;
    
    figure
    hold on
    P = patch('Faces',F,'Vertices',V,'facecolor',[.5 .5 .5],'edgecolor','none');
    set(P, 'facealpha', 0.5)
    camlight('headlight','infinite')
    axis equal
    rotate3d on
    
    node_to_plot = [1];
    scatter3(V(node_to_plot,1), V(node_to_plot,2), V(node_to_plot,3), 50, 'r', 'filled')
    
    [vertex_index, face_index] = mesh_findneighbor(F, 1, 3);
    
    face_to_delete = F(face_index,:);
    
    P = patch('Faces',face_to_delete,'Vertices',V,'facecolor','y','edgecolor','none');
    
    F_delete = F;
    F_delete(face_index,:) = [];
    
    [V_new, F_new] = removeisolatednode(V,F_delete);
    
    figure
    hold on
    P = patch('Faces',F_new,'Vertices',V_new,'facecolor',[.5 .5 .5],'edgecolor','none');
    % set(P, 'facealpha', 0.5)
    camlight('headlight','infinite')
    axis equal
    rotate3d on
    
    % ok, now F_new and V_new is a triangulation surface with a hole. Let's see
    % how does mesh_quality_check perform.
    
    mesh_quality_check(F_new, V_new, true); % check for defects
    [ F_fix,V_fix ] = mesh_quality_check(F_new, V_new, false);
    
    figure
    ax1 = subplot(1,2,1);
    P = patch('Faces',F_new,'Vertices',V_new,'facecolor',[.5 .5 .5],'edgecolor','none');
    camlight('headlight','infinite')
    axis equal
    rotate3d on
    ax2 = subplot(1,2,2);
    P = patch('Faces',F_fix,'Vertices',V_fix,'facecolor',[.5 .5 .5],'edgecolor','none');
    camlight('headlight','infinite')
    axis equal
    rotate3d on
    
    hlink = linkprop([ax1,ax2],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
    
end
% ok - it does work as expected. Fine!

%% Mesh quality check
% Check the boundary surfaces for erroneous intersections

% I'm not sure whether iso2mesh meshcheckrepair does a good job checking
% for defects and fixing them. It is based on the meshfix package.
% Perhaps Freesurfer has some other tools that we can use?

[ f1,v1 ] = mesh_quality_check(f1, v1, false);
[ f2,v2 ] = mesh_quality_check(f2, v2, false);
[ f3,v3 ] = mesh_quality_check(f3, v3, false);
[ f4,v4 ] = mesh_quality_check(f4, v4, false);
[ f5,v5 ] = mesh_quality_check(f5, v5, false);

%% Mesh inclusion order check
% Check that the meshes bound one another in the correct order

% Skin[1] > Skull[2] > CSF[3] > Pial[4] = Cerebellum[5] > WM/GM source space[6]

assert(mesh_surfaceinclusion(f1, v1, v2), 'Skull[2] is NOT bounded by Skin[1]!')
assert(mesh_surfaceinclusion(f2, v2, v3), 'CSF[3] is NOT bounded by Skull[2]!')
assert(mesh_surfaceinclusion(f3, v3, v4), 'Pial[4] is NOT bounded by CSF[3]!')
assert(mesh_surfaceinclusion(f3, v3, v5), 'Cerebellum[4] is NOT bounded by CSF[3]!')

%% Test out the mesh_surfaceintersect function
% let's manually shift the surface forward and see if the code detects it
while false
    
    % before shifting, we can visually tell skin and skull do not intersect
    faces_list = {f1, f2};
    vertices_list = {v1, v2};
    mesh_plot_bem_surfaces(faces_list, vertices_list);
    
    % test for intersection
    intersected_pre = mesh_surfaceintersect(f1, v1, f2, v2);
    
    % now we shift the skull forward to manually intersect with the skin
    v2_shifted = v2;
    v2_shifted(:,2) = v2(:,2) + 10;
    vertices_list = {v1, v2_shifted};
    mesh_plot_bem_surfaces(faces_list, vertices_list);
    
    % test for intersection
    intersected_post = mesh_surfaceintersect(f1, v1, f2, v2_shifted);
    
    % visualize the intersected surfaces with green coloring
    for ii = 1:length(intersected_post.intx_surf)
        current_surf = intersected_post.intx_surf{ii};
        patch(current_surf,'facecolor','g','edgecolor','none');
    end
    
end
% ok - it does work as expected. Fine!

%% Mesh intersection check
% Check the pial surface / skull surface doesn't intersect with cerebellum.
% As reported by Matti Stenroos' s 2016 paper, there seems to be an
% artifact of csf surface intersecting with the cerebellum at the base of
% the skull. Let's check that.

% check that the csf mesh does not intersect with the cerebellum
tic
disp('Checking csf mesh and cerebellum...')
intersected = mesh_surfaceintersect(f3, v3, f5, v5);
assert(~intersected.intersected, 'csf mesh and cerebellum intersected!')
toc

% check that the pial surface does not intersect with the cerebellum
tic
disp('Checking pial surface and cerebellum...')
intersected = mesh_surfaceintersect(f4, v4, f5, v5);
assert(~intersected.intersected, 'pial surface and cerebellum intersected!')
toc

%% Source space check
% load the source space constructed using MNE
src_fn = 'ico5-src.mat';
load(fullfile(datapath, src_fn))

% check the normal vectors are unit vectors
assert(all(vecnorm(lnrm,2,2)-1 < 10^-6), 'Left hemisphere source-point normals are not unit vectors!')
assert(all(vecnorm(rnrm,2,2)-1 < 10^-6), 'Right hemisphere source-point normals are not unit vectors!')

% check the number of sources matches specified in the source space
% subdivision in the filename
if contains(src_fn, 'ico5')
    assert(sum(lidx) == 10242, 'Number of sources incorrect for ico5 sub-division');
    assert(sum(ridx) == 10242, 'Number of sources incorrect for ico5 sub-division');
end

% make sure source coordinates are of double type
assert(isa(lsrc, 'double') && isa(rsrc, 'double'), 'Source coordinates are not of double type')

% create the sources cell
left_sources = struct;
left_sources.hemisphere = 'left';
left_sources.points = lsrc(lidx==1,:);
left_sources.normals = double(lnrm(lidx==1,:));

right_sources = struct;
right_sources.hemisphere = 'right';
right_sources.points = rsrc(ridx==1,:);
right_sources.normals = double(rnrm(ridx==1,:));

sources = {left_sources; right_sources};

% % quick visual check
% faces_list = {f4};
% vertices_list = {v4};
% mesh_plot_bem_surfaces(faces_list, vertices_list, sources);

%% Test out that the sources are actually subsampling the wm surfaces
% meaning, did MNE use vertices on the wm surface to construct sources?
while false
    lhwm_fn = 'mri2mesh_fs_lh_wm_fixed.mat';
    load(fullfile(datapath, lhwm_fn))
    lh_f = faces; lh_v = vertices;
    rhwm_fn = 'mri2mesh_fs_rh_wm_fixed.mat';
    load(fullfile(datapath, rhwm_fn))
    rh_f = faces; rh_v = vertices;
    
    faces_list = {lh_f, rh_f};
    vertices_list = {lh_v, rh_v};
    mesh_plot_bem_surfaces(faces_list, vertices_list, sources)
end
% Yes and no! Only so on the unsmoothed lh.white and rh.white. This is
% ok, because I believe the wm surface is only used for constructing the
% source space, maybe for defining the normal direction of sources. But
% other than these two functions, it's not used explicitly as a boundary
% surface so smoothing or not / on the surface or not doesn't really
% matter. We load it here just for the sake of visualization.

%% Adjust source space to ensure inclusion and minimal cortex thickness
% Check whether all source points are within the pial surface and we
% correct source point locations to obtain minimal thickness of 0.5mm.

tic
[ sources{1}.points, left_altered_idx ] = mesh_sourceinclusion(sources{1}.points, f4, v4, 0.5);
[ sources{2}.points, right_altered_idx ] = mesh_sourceinclusion(sources{2}.points, f4, v4, 0.5);
disp('Time taken in adjusting source points')
toc

% % Visualize the final source points with pial surface on
% faces_list = {f4};
% vertices_list = {v4};
% mesh_plot_bem_surfaces(faces_list, vertices_list, sources);

% looking good!!

%% Load electrodes and project to skin surface
% we generate this -elc-trans file using generate_elctrans bash script

elc_fn = 'TDelp_resting-elc-trans.mat';
dig = mesh_getelectrode(fullfile(datapath,elc_fn), f1, v1);

% % visualize the electrodes on skin surface
% faces_list = {f1};
% vertices_list = {v1};
% mesh_plot_bem_surfaces(faces_list, vertices_list, [], dig);

%% Mesh preparation completed
disp('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
disp('Total time taken in running mesh_bem_meshing.m script:')
toc(procstart)
disp('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')

%% Visualize BEM head model

faces_list = {f1, f2, f3, f4, f5};
vertices_list = {v1, v2, v3, v4, v5};
mesh_plot_bem_surfaces(faces_list, vertices_list, sources, dig);
title('BEM Head Model')
xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)')
set(gca, 'FontSize', 16)

%% Save all processed variables
% also save modified source points to a mat file to create forward model in MNE
lsrc(lidx==1,:) = sources{1}.points;
rsrc(ridx==1,:) = sources{2}.points;

savefn = fullfile(datapath, '4layer_BEM_model.mat');
save(savefn, 'f1','f2','f3','f4','f5','f6',...
    'v1','v2','v3','v4','v5','v6',...
    'lsrc', 'rsrc',...
    'sources', 'dig')

%% Prepare BEM head model .mat files for hbf_bem solver
% hbf_bem solver by Matti Stenroos expects the meshes to be inputed as
% ordered structure variables in column cell entries.

% convert sources and dig to SI unit in meter
sources{1}.points = sources{1}.points / 1000;
sources{2}.points = sources{2}.points / 1000;
dig = dig / 1000;

% bmeshes are ordered from inwards to outwards

% prepare for the 3shell BEM model:
skin = struct;
skin.p = v1 / 1000; % convert to SI unit meter
skin.e = f1;
skin.nop = size(v1,1);
skin.noe = size(f1,1);
skin.surface = 'outer_skin';

skull = struct;
skull.p = v2 / 1000; % convert to SI unit meter
skull.e = f2;
skull.nop = size(v2,1);
skull.noe = size(f2,1);
skull.surface = 'outer_skull';

csf = struct;
csf.p = v3 / 1000; % convert to SI unit meter
csf.e = f3;
csf.nop = size(v3,1);
csf.noe = size(f3,1);
csf.surface = 'inner_skull';

bmeshes = {csf; skull; skin};

% save the variables for hdf_bem solver
savefn = fullfile(datapath, '3shell_hbfcmpt_BEM_model.mat');
save(savefn, 'bmeshes', 'sources', 'dig')

%%
% prepare for 4shell BEM model:
pial = struct;
pial.p = v4 / 1000; % convert to SI unit meter
pial.e = f4;
pial.nop = size(v4,1);
pial.noe = size(f4,1);
pial.surface = 'pial';

cerebellum = struct;
cerebellum.p = v5 / 1000; % convert to SI unit meter
cerebellum.e = f5;
cerebellum.nop = size(v5,1);
cerebellum.noe = size(f5,1);
cerebellum.surface = 'cerebellum';

bmeshes = {pial; cerebellum; csf; skull; skin};

% save the variables for hdf_bem solver
savefn = fullfile(datapath, '4shell_hbfcmpt_BEM_model.mat');
save(savefn, 'bmeshes', 'sources', 'dig')


