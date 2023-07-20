function [ ] = mesh_create_4shell(datapath, ico)
%% mesh_bem_meshing script -> modified into a function
% This function is based on the mesh_bem_meshing script.
% The purpose of this function is to check the quality of different BEM
% surfaces and correct intersections of surfaces if necessary.

% We also check the inclusion of source space within the pial surface and
% check for minimal cortex thickness to be >= 0.5 mm. We will point-wise
% shift the source away from pial surface to obtain a minimum thickness of
% 0.5 mm.

% close all
% clear all
%
% % Change the current folder to the folder of this m-file.
% tmp = matlab.desktop.editor.getActive;
% cd(fileparts(tmp.Filename));
% clearvars tmp

if nargin < 2
    ico = 'ico5';
end

procstart = tic;

%% Path configuration
% addpath(genpath('./'));

% datapath should point to the final_structure folder
% datapath = '/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/MNE/subject_data/SP011_ah/simnibs_pipe/bem_surfaces/final_structure';

%% Mesh loading and downsampling
% also smoothing with iso2mesh smoothsurf()
% we need to downsample to a coarser resolution after smoothing as well
% with GIBBON codes. We specify the number of vertices we want in the
% output triangulation mesh.

disp('******************************')
disp('Step 1: mesh loading and downsampling')
disp('******************************')

% if the smoothed and downsampled surfaces are already saved, we load them
% directly rather than using the mesh_loadmesh() function, which takes a
% long time!

if isfile(fullfile(datapath, 'archive_4layer_BEM_model.mat'))
    disp('Found archived BEM meshes. Loading from archive_4layer_BEM_model.mat directly. If you would like to reload, delete this file first.')
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
    
    if isfile(fullfile(datapath, 'headreco_csf_compatible.mat'))
        csf_fn = 'headreco_csf_compatible.mat';
    else
        csf_fn = 'mri2mesh_csf_fixed.mat';
    end
    
    if isfile(fullfile(datapath, 'mri2mesh_gm_fixed_updated.mat'))
        pial_fn = 'mri2mesh_gm_fixed_updated.mat';
    else
        pial_fn = 'mri2mesh_gm_fixed.mat';
    end
    
    if isfile(fullfile(datapath, 'mri2mesh_cerebellum_fixed_updated.mat'))
        cerebll_fn = 'mri2mesh_cerebellum_fixed_updated.mat';
    else
        cerebll_fn = 'mri2mesh_cerebellum_fixed.mat';
    end
    
    wm_fn = 'mri2mesh_wm_fixed.mat';
    
    tic
    [ f1,v1,hlink1 ] = mesh_loadmesh(skin_fn, datapath, 0, 0, [], true);
    [ f2,v2,hlink2 ] = mesh_loadmesh(skull_fn, datapath, 0, 10, 16384, true);
    [ f3,v3,hlink3 ] = mesh_loadmesh(csf_fn, datapath, 0, 0, 16384, true);
    [ f4,v4,hlink4 ] = mesh_loadmesh(pial_fn, datapath, 50, 10, 65536, true);
    [ f5,v5,hlink5 ] = mesh_loadmesh(cerebll_fn, datapath, 50, 0, [], true);
    [ f6,v6,hlink6 ] = mesh_loadmesh(wm_fn, datapath, 100, 0, [], true);
    disp('Total time taken in loading, smoothing and downsampling the meshes:')
    toc
    disp(' ')
    
    close all
    
    save(fullfile(datapath, 'archive_4layer_BEM_model.mat'),...
        'f1','f2','f3','f4','f5','f6',...
        'v1','v2','v3','v4','v5','v6')
end

%% Remesh BEM surfaces to smaller vertex resolutions
% sometimes we may want to remesh the surfaces into a smaller resolution to
% test out BEM solvers. we can use the mesh_remeshsurface() function.
% wm/gm boundary (f6, v6) is kept as is since it's not used for any
% computation in MATLAB.

disp('******************************')
disp('Step 2: remesh to target vertex resolutions')
disp('******************************')

% Note that even though we try to aim for a resolution that is compatible
% with icosahedron decimation, it is not always possible to do that while
% maintaining Delauny using a downsampled resolution to begin with for some
% of the convoluted surfaces, such as pial surface. The way Freesurfer does
% it is by significantly upsampling the spatial resolution and derive
% vertex set that allows for icosahedron resolution. However, mri2mesh
% surfaces do not enjoy that property. This is not necessary unless you
% want to use MNE's built-in icosahedron downsampling, which I would
% recommend against because under the ico4 resolution pial surface is
% represented with merely 4096 vertices, losing dramatically the intricate
% properties of the pial surface, sometimes cutting out certain source
% points. So if you want to use these surfaces but still build 3-shell
% model in MNE, you can set ico='none' in MNE make_bem_model() calls. 

% for skin, skull, csf, pial, cerebellum in this order
vertex_resolution_list = [4096, 8192, 8192, 32768, 8192];

meshcell = {};
for ii = 1:5
    eval(sprintf('surf1.p=v%s; surf1.e=f%s; surf1.nop=size(v%s,1);', num2str(ii), num2str(ii), num2str(ii)))
    meshcell{ii} = surf1; %#ok<*AGROW,*NODEF>
end

% remesh the surfaces
meshcell = mesh_remeshsurface(meshcell, vertex_resolution_list);

for ii = 1:5
    surf1 = meshcell{ii}; %#ok<*NASGU>
    eval(sprintf('v%s=surf1.p; f%s=surf1.e;', num2str(ii), num2str(ii)))
end

%% Address pial/cerebellum and csf inclusion issue
% Updated on 06/13/2020: on some young brains, after smoothing and
% downsampling, some pial surface vertices end up being outside of the
% downsampled CSF surface. This is understandable, but we would like to
% address it by dilating the CSF a little.

disp('******************************')
disp('Step 3.1: checking pial inclusion by CSF')
disp('******************************')

% find out pial surface points (if any) that are outside of the csf surface 
in = intriangulation(v3,f3,v4);
out_vpoint = v4(in==0,:);

if ~all(in==1)
    disp(['Number of pial surface points outside of CSF after smoothing and downsampling: ', num2str(size(out_vpoint,1))])
    
%     figure;
%     hold on
%     P = patch('Faces',f3,'Vertices',v3,'facecolor',[.5 .5 .5],'edgecolor','none');
%     P = patch('Faces',f4,'Vertices',v4,'facecolor',[1 0 0],'edgecolor','none');
%     scatter3(out_vpoint(:,1), out_vpoint(:,2), out_vpoint(:,3), 50, 'm', 'filled')
%     camlight('headlight','infinite')
%     axis equal
%     camorbit(0, 180)
%     camlight('headlight')
%     camorbit(0, 180)
%     camorbit(0, 270)
%     rotate3d on
%     title('Pial surface points outside CSF', 'FontSize', 20)

    disp('Fixing the CSF surface to accommodate these outside pial surface vertices.')
    [ distances, surface_points ] = ...
        point2trimesh('Faces', f3, 'Vertices', v3, 'QueryPoints', out_vpoint, 'Algorithm', 'parallel');
    v3 = mesh_dilate_out(v3, f3, out_vpoint, surface_points, 0.3);
    disp('Fixing the CSF surface with 0.3mm margin - Done.')
    disp(' ')

%     figure;
%     hold on
%     P = patch('Faces',f3,'Vertices',v3,'facecolor',[.5 .5 .5],'edgecolor','none');
%     P = patch('Faces',f4,'Vertices',v4,'facecolor',[1 0 0],'edgecolor','none');
%     scatter3(out_vpoint(:,1), out_vpoint(:,2), out_vpoint(:,3), 50, 'm', 'filled')
%     camlight('headlight','infinite')
%     axis equal
%     camorbit(0, 180)
%     camlight('headlight')
%     camorbit(0, 180)
%     camorbit(0, 270)
%     rotate3d on
%     title('CSF surface after correction for Pial surface', 'FontSize', 20)

else
    disp('All pial surface points are inside the CSF surface.')
    disp(' ')
end

% Updated on 03/02/2022: the same issue could occur for cerebellum surface.
disp('******************************')
disp('Step 3.2: checking cerebellum inclusion by CSF')
disp('******************************')

% find out pial surface points (if any) that are outside of the csf surface 
in = intriangulation(v3,f3,v5);
out_vpoint = v5(in==0,:);

if ~all(in==1)
    disp(['Number of cerebellum surface points outside of CSF after smoothing and downsampling: ', num2str(size(out_vpoint,1))])

%     figure;
%     hold on
%     P = patch('Faces',f3,'Vertices',v3,'facecolor',[.5 .5 .5],'edgecolor','none');
%     P = patch('Faces',f5,'Vertices',v5,'facecolor',[1 0 0],'edgecolor','none');
%     scatter3(out_vpoint(:,1), out_vpoint(:,2), out_vpoint(:,3), 50, 'm', 'filled')
%     camlight('headlight','infinite')
%     axis equal
%     camorbit(0, 180)
%     camlight('headlight')
%     camorbit(0, 180)
%     camorbit(0, 270)
%     rotate3d on
%     title('Cerebellum surface points outside CSF', 'FontSize', 20)

    disp('Fixing the CSF surface to accommodate these outside cerebellum surface vertices.')
    [ distances, surface_points ] = ...
        point2trimesh('Faces', f3, 'Vertices', v3, 'QueryPoints', out_vpoint, 'Algorithm', 'parallel');
    v3 = mesh_dilate_out(v3, f3, out_vpoint, surface_points, 0.3);
    disp('Fixing the CSF surface with 0.3mm margin - Done.')
    disp(' ')

%     figure;
%     hold on
%     P = patch('Faces',f3,'Vertices',v3,'facecolor',[.5 .5 .5],'edgecolor','none');
%     P = patch('Faces',f5,'Vertices',v5,'facecolor',[1 0 0],'edgecolor','none');
%     scatter3(out_vpoint(:,1), out_vpoint(:,2), out_vpoint(:,3), 50, 'm', 'filled')
%     camlight('headlight','infinite')
%     axis equal
%     camorbit(0, 180)
%     camlight('headlight')
%     camorbit(0, 180)
%     camorbit(0, 270)
%     rotate3d on
%     title('CSF surface after correction for Cerebellum surface', 'FontSize', 20)

else
    disp('All cerebellum surface points are inside the CSF surface.')
    disp(' ')
end

%% Address csf and skull inclusion issue
% Updated on 07/20/2022: in rare cases, at the skull base around the spinal
% cord where we made an incision, the CSF surface could protrude out from
% the smoothed and downsampled skull surface. We will again address it by
% dilating the skull surface a little. 

disp('******************************')
disp('Step 3.3: checking CSF inclusion by skull')
disp('******************************')

% find out CSF surface points (if any) that are outside of the skull surface 
in = intriangulation(v2,f2,v3);
out_vpoint = v3(in==0,:);

if ~all(in==1)
    disp(['Number of CSF surface points outside of skull after smoothing and downsampling: ', num2str(size(out_vpoint,1))])

    disp('Fixing the skull surface to accommodate these outside CSF surface vertices.')
    [ distances, surface_points ] = ...
        point2trimesh('Faces', f2, 'Vertices', v2, 'QueryPoints', out_vpoint, 'Algorithm', 'parallel');
    v2 = mesh_dilate_out(v2, f2, out_vpoint, surface_points, 0.5, [], 10, 7);
    disp('Fixing the Skull surface with 0.5mm margin - Done.')
    disp(' ')

else
    disp('All CSF surface points are inside the skull surface.')
    disp(' ')
end

%% Visualize BEM meshes function: mesh_plot_bem_surfaces()
% faces_list = {f4};
% vertices_list = {v4};
% mesh_plot_bem_surfaces(faces_list, vertices_list);

%% Mesh quality check
% Check the boundary surfaces for erroneous intersections

disp('******************************')
disp('Step 4: checking mesh quality')
disp('******************************')

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

disp('******************************')
disp('Step 5: checking mesh inclusions')
disp('******************************')

% Skin[1] > Skull[2] > CSF[3] > Pial[4] = Cerebellum[5] > WM/GM source space[6]

assert(mesh_surfaceinclusion(f1, v1, v2), 'Skull[2] is NOT bounded by Skin[1]!')
assert(mesh_surfaceinclusion(f2, v2, v3), 'CSF[3] is NOT bounded by Skull[2]!')
assert(mesh_surfaceinclusion(f3, v3, v4), 'Pial[4] is NOT bounded by CSF[3]!')
assert(mesh_surfaceinclusion(f3, v3, v5), 'Cerebellum[5] is NOT bounded by CSF[3]!')

%% Mesh intersection check
% Check that pial surface / skull surface doesn't intersect with
% cerebellum. As reported by Matti Stenroos' s 2016 paper, there seems to
% be an artifact of csf surface intersecting with the cerebellum at the
% base of the skull. Let's check that.

disp('******************************')
disp('Step 6: checking mesh intersection')
disp('******************************')

% check that the csf mesh does not intersect with the cerebellum
% Updated on 07/29/2022: in rare cases, the cerebellum surface also has
% overgrabbing of the dura mater but we have no way of editing it through
% freesurfer. Hence we will try to adjust one of the two surfaces a little
% if an intersection is detected.
tic
disp('Checking csf mesh and cerebellum...')
% skip remeshing on this first round since both surfaces are of similar
% resolutions. We also turn on edit_V to potentially edit intersections by
% dilating the CSF surface out a bit
[ intersected, v3, v5 ] = mesh_surfaceintersect(f3, v3, f5, v5, [], true, true);
if intersected.intersected  % intersection found and edited
    disp('Intersection between CSF and cerebellum found and corrected.')
    intersected = mesh_surfaceintersect(f3, v3, f5, v5);
    assert(~intersected.intersected, 'csf mesh and cerebellum intersected!')
end
toc

% check that the pial surface does not intersect with the cerebellum
tic
disp('Checking pial surface and cerebellum...')
intersected = mesh_surfaceintersect(f4, v4, f5, v5);
assert(~intersected.intersected, 'pial surface and cerebellum intersected!')
toc

%% Accomodate looping through ico3-5 source spaces
if strcmp(ico, '345')
    ico_list = {'ico3', 'ico4', 'ico5'};
else
    ico_list = {ico};
end

for ii = 1:length(ico_list)
    % iterate through different icosahedron resolutions
    ico = ico_list{ii};
    
    %% Source space check
    % make sure the source space is valid for forward modeling
    
    disp('******************************')
    disp('Step 7: checking source space quality')
    disp('******************************')
    
    % load the source space constructed using MNE
    src_fn = [ico, '-src.mat'];
    load(fullfile(datapath, src_fn)) %#ok<*LOAD>
    
    % check the normal vectors are unit vectors
    assert(all(vecnorm(lnrm,2,2)-1 < 10^-6), 'Left hemisphere source-point normals are not unit vectors!') %#ok<*USENS>
    assert(all(vecnorm(rnrm,2,2)-1 < 10^-6), 'Right hemisphere source-point normals are not unit vectors!')
    
    % check the number of sources matches specified in the source space
    % subdivision in the filename
    if contains(src_fn, 'ico5')
        assert(sum(lidx) == 10242, 'Number of sources incorrect for ico5 sub-division');
        assert(sum(ridx) == 10242, 'Number of sources incorrect for ico5 sub-division');
    elseif contains(src_fn, 'ico4')
        assert(sum(lidx) == 2562, 'Number of sources incorrect for ico4 sub-division');
        assert(sum(ridx) == 2562, 'Number of sources incorrect for ico4 sub-division');
    elseif contains(src_fn, 'ico3')
        assert(sum(lidx) == 642, 'Number of sources incorrect for ico3 sub-division');
        assert(sum(ridx) == 642, 'Number of sources incorrect for ico3 sub-division');
    else
        error(['Invalid source space icosahedron sub-division resolution in ', src_fn, '.'])
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
    
    %% Adjust source space to ensure inclusion and minimal cortex thickness
    % Check whether all source points are within the pial surface and we
    % correct source point locations to obtain minimal thickness of 0.5 mm.
    
    disp('******************************')
    disp('Step 8: adjusting thin cortical thickness sources')
    disp('******************************')
    
    tic
    [ sources{1}.points, left_altered_idx ] = mesh_sourceinclusion(sources{1}.points, f4, v4, 0.5); %#ok<*ASGLU>
    [ sources{2}.points, right_altered_idx ] = mesh_sourceinclusion(sources{2}.points, f4, v4, 0.5);
    disp('Time taken in adjusting source points:')
    toc
    
    % % Visualize the final source points with pial surface on
    % faces_list = {f4};
    % vertices_list = {v4};
    % mesh_plot_bem_surfaces(faces_list, vertices_list, sources);
    
    % looking good!!
    
    %% Load electrodes and project to skin surface
    % we generate this -elc-trans file using generate_elctrans bash script
    
    disp('******************************')
    disp('Step 9: projecting electrodes to scalp')
    disp('******************************')
    
    elc_fn = 'elc-trans.mat';
    dig = mesh_getelectrode(fullfile(datapath,elc_fn), f1, v1);
    
    % % Visualize the electrodes on skin surface
    % faces_list = {f1};
    % vertices_list = {v1};
    % mesh_plot_bem_surfaces(faces_list, vertices_list, [], dig);
    % title('Electrode Locations on Skin', 'FontSize', 20)
    
    %% Visualize BEM head model
    faces_list = {f1, f2, f3, f4, f5};
    vertices_list = {v1, v2, v3, v4, v5};
    mesh_plot_bem_surfaces(faces_list, vertices_list, sources, dig);
    title('BEM Head Model')
    xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)')
    set(gca, 'FontSize', 20)
    
    %% Save all processed variables
    % also save modified source points to a mat file to create forward model in MNE
    lsrc(lidx==1,:) = sources{1}.points;
    rsrc(ridx==1,:) = sources{2}.points;
    
    savefn = fullfile(datapath, ['4layer_BEM_model_all_meshes-', ico, '.mat']);
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
    
    % N.B.: bmeshes are ordered from inwards to outwards
    
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
    savefn = fullfile(datapath, ['3shell_hbfcmpt_BEM_model-', ico, '.mat']);
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
    savefn = fullfile(datapath, ['4shell_hbfcmpt_BEM_model-', ico, '.mat']);
    save(savefn, 'bmeshes', 'sources', 'dig')

end

%% Mesh preparation completed
disp('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
disp('Total time taken in running mesh_create_4shell function:')
toc(procstart)
disp('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')

pause(5)
close all

end



function [] = test_meshfunctions()
% Moving test functions of customized mesh quality checks here 

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
    
    node_to_plot = 1;
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

end
