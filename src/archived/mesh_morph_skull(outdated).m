function [ ] = mesh_morph_skull(datapath, dilation_mode, verbose)

%% morph_skull_surf script -> modified into a function
% This script is used to morph the bone surface extracted from the
% headreco pipeline with the skull surface extracted from the mri2mesh
% pipeline. mri2mesh_skull has the advantage of being a relatively smooth
% surface that terminates at the foramen magnum section of the brainstem.
% However, it has inaccurate deviations along the top of the skull, on the
% two sides of superior sagittal sinus. In contrast, headreco_bone is a
% fairly accurate skull surface that accurately captures the boundaries for
% the whole head. However, because the whole head extraction is done, there
% are sharp edges and hollow sections in the spinal cervical sections as
% well as for bony structures for the face.

% Hence, our solution is to morph the two surfaces together to get an
% optimal skull surface for the head portion, discarding the spinal
% cervical sections and the facial bones.

%% Path configuration
addpath(genpath('./'));

% datapath should point to the final_structure folder
% datapath = '/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/MNE/subject_data/SP011_ah/simnibs_pipe/bem_surfaces/final_structure';

if nargin < 2
    dilation_mode = 'top';
end

if nargin < 3
    verbose = false;
end

%% Load the two skull surfaces and visualize

% Load headreco bone surface
load(fullfile(datapath, 'headreco_bone.mat'))

f1 = faces; %#ok<*NODEF>
v1 = vertices;

figure;
hold on

P = patch('Faces',f1,'Vertices',v1,'facecolor',[.5 .5 .5],'edgecolor','none');
camlight('headlight','infinite')
axis equal
% camorbit(180, 0)
camorbit(0, 180)
camlight('headlight')
camorbit(0, 180)
camorbit(0, 270)
rotate3d on

set(P, 'facealpha', 0.25)

% visualize the SurfaceRAS center
scatter3(0, 0, 0, 100, 'r', 'filled')

% Load mri2mesh bone surface
load(fullfile(datapath, 'mri2mesh_skull_fixed.mat'))

f2 = faces;
v2 = vertices;

P = patch('Faces',f2,'Vertices',v2,'facecolor',[1 0 0],'edgecolor','none');
set(P, 'facealpha', 0.5)

%% Separate the surfaces into 4 quadrants with two planes

% find frontal and occipital poles of the mri2mesh surface
[maxy, fpidx] = max(v2(:,2));
[miny, opidx] = min(v2(:,2));

scatter3(v2(fpidx,1), v2(fpidx,2), v2(fpidx,3), 500, 'b', 'filled')
scatter3(v2(opidx,1), v2(opidx,2), v2(opidx,3), 500, 'g', 'filled')

% pick the more superior point of the two poles. We pick the superior one
% to avoid having to deal with sharp transitions around eye sockets
zthresh = max(v2(fpidx,3), v2(opidx,3));

% we now define an axial plane to cut the surfaces into two parts
maxplane = max([max(v2(:,1)), max(v1(:,1)), max(v2(:,2)), max(v1(:,2))]);
minplane = min([min(v2(:,1)), min(v1(:,1)), min(v2(:,2)), min(v1(:,2))]);
vxtvalue = max(abs(maxplane), abs(minplane))+10;

plane_x = [vxtvalue, vxtvalue, -vxtvalue, -vxtvalue];
plane_y = [vxtvalue, -vxtvalue, -vxtvalue, vxtvalue];
plane_z = ones(1,4) .* zthresh;

% visualize the plane for cutting
patch(plane_x, plane_y, plane_z, [0,0,0,0],'facecolor',[0 0 1], 'facealpha', 1)

% We also find the most inferior point of the mri2mesh surface, define a
% coronal plane that cuts the mri2mesh surface for a anterior-inferior
% portion, which contains most of the convoluted surface edges.

[minz, ipidx] = min(v2(:,3));
ythresh = v2(ipidx,2);

scatter3(v2(ipidx,1), v2(ipidx,2), v2(ipidx,3), 500, 'm', 'filled')

% we now define a coronal plane to further cut the surfaces into four parts
maxplane = max([max(v2(:,1)), max(v1(:,1)), max(v2(:,3)), max(v1(:,3))]);
minplane = min([min(v2(:,1)), min(v1(:,1)), min(v2(:,3)), min(v1(:,3))]);
vxtvalue = max(abs(maxplane), abs(minplane))+10;

plane_x = [vxtvalue, vxtvalue, -vxtvalue, -vxtvalue];
plane_y = ones(1,4) .* ythresh;
plane_z = [vxtvalue, -vxtvalue, -vxtvalue, vxtvalue];

% visualize the plane for cutting
patch(plane_x, plane_y, plane_z, [0,0,0,0],'facecolor',[0 1 0], 'facealpha', 1)

% manual sanity check: do the planes and surfaces look reasonable?
checkok =  input('Do the surfaces and 4 quadrants look ok? (y/n): ', 's');
if strcmpi(checkok, 'n')
    save(fullfile(datapath, 'mesh_morph_skull_failed_checkpoint1_workspace'))
    error('Manual check failed. Workspace is saved. Please debug...')
end

%% Overall game plan
% Our overall goal is the following:

% - For mri2mesh surface above the cutting plane, we want to dilate the
% surface into headreco surface to fix those trange indentations caused by
% inaccurate Freesurfer based segmentation.

% - For mri2mesh surface below the cutting plane, we want to smooth the
% surface a little bit to avoid sharp edges and numerical instability in
% BEM solver. The optimal smoothing is probably a distance-based dilation
% towards a sphere. This should bring out the fossa region much more than
% the occipital skull, which is fairly accurate and we hope to keep it
% unchanged.

%% Goal 1: dilate mri2mesh towards headreco above the plane

% This is a difficult task without sophisticated mesh objects and dedicated
% software. The following heuristic will be used for now:

% 1) We find all headreco vertices above zthresh.
% 2) We ask which of these vertices are outside the mri2mesh compartment.
% 3) For each vertex identified, we construct a small sphere on the
% headreco surface.
% 4) We join all spheres to find out all isolated ROI for dilation on the
% headreco surface. We create separate surfaces for each ROI.
% 5) We then project the ROI surface to mri2mesh to identify a slightly
% larger ROI for dilation on the mri2mesh surface.
% 6) For each vertex on the mri2mesh ROI surface, we propose to dilate it
% by changing the vertex coordinate to be that of projected onto headreco.
% If after projection it is outside the original mri2mesh compartment, we
% accept the dilation and change it. If instead it is now inside, we reject
% the dilation and keep it the vertex in its original position on mri2mesh.
% 7) After looping through all ROI, we would have dilated mri2mesh to
% headreco above the plane.
% 8) We run some meshfix script outside MATLAB to check for mesh defects,
% and we save that as the new mri2mesh skull surface.
% 9) Smoothing below the plane needs to be done separately. I haven't
% thought of a good way to do it yet!

%% ***Update note on 05/12/2020:
% it looks like we might also want to do dilation in the
% posterior-inferior portion of the head.

% The easiest way to do this for now is just to make two cases:
% case 1: dilate only the top
% case 2: dilate both top and posterior-inferior

% dilation_mode = 'top_posinf';

switch dilation_mode
    case 'top'
        m2m_vindex = v2(:,3) > zthresh;
        vindex = v1(:,3) > zthresh;
    case 'top_posinf'
        m2m_vindex = v2(:,3) > zthresh | v2(:,2) < ythresh;
        vindex = v1(:,3) > zthresh | v1(:,2) < ythresh;
end
%%
% % 1) We find all headreco vertices above zthresh.
% vindex = v1(:,3) > zthresh;
% above_vpoint = v1(vindex,:);
%
% % quick visual check
% scatter3(above_vpoint(:,1),above_vpoint(:,2), above_vpoint(:,3), 10, 'g')
%
% %%
% % 2) We ask which of these vertices are outside the mri2mesh compartment.
% in = intriangulation(v2,f2,above_vpoint);
% out_vpoint = above_vpoint(in==0,:);
%
% % quick visual check
% scatter3(out_vpoint(:,1),out_vpoint(:,2), out_vpoint(:,3), 10, 'm')

%%
% 3) It turns out that most of the headreco surface vertices are outside
% the mri2mesh. Hence rather than going through the computationally
% intensive process of constructing ROI patches, we will simply propose to
% dilate ALL mri2mesh points to headreco. Then we ask whether the projected
% points are outside original mri2mesh surface for acceptance criterion.

% candidate points on mri2mesh
%         m2m_vindex = v2(:,3) > zthresh;
m2m_above_vpoint = v2(m2m_vindex,:);

% quick visual check
if verbose
    scatter3(m2m_above_vpoint(:,1),m2m_above_vpoint(:,2), m2m_above_vpoint(:,3), 10, 'g')
end

%%
% 4) create a new triangulation surface for headreco surface above the plane
%         vindex = v1(:,3) > zthresh;
vcount = 1:length(vindex);
vin = vcount(vindex);

face_idx = zeros(size(f1,1), 1);
for ii = 1:size(f1,1)
    if any(ismember(f1(ii,:), vin))
        face_idx(ii) = 1;
    end
end

f1_above = f1(logical(face_idx), :);

% we also need to clean up the indexing of this new triangulation surface
v_include = double(unique(f1_above));
sub_v1 = v1(v_include,:);

% clean the indices contained in f1_above to correspond to sub_v1
for ii = 1:size(f1_above,1)
    [logiout, idxout] = ismember(f1_above(ii,:), v_include);
    assert(all(logiout), 'Some vertices called by f1_above not in v_include!')
    f1_above(ii,:) = idxout;
end

% quick visual check
if verbose
    P = patch('Faces',f1_above,'Vertices',sub_v1,'facecolor',[0 .5 .5],'edgecolor','none');
end

%%
% 5) propose to project the mri2mesh vertices to this new headreco surface
[ distances, surface_points ] = ...
    point2trimesh('Faces', f1_above, 'Vertices', sub_v1, 'QueryPoints', m2m_above_vpoint, 'Algorithm', 'parallel');

% quick visual check
if verbose
    figure;
    hold on
    P = patch('Faces',f1_above,'Vertices',sub_v1,'facecolor',[0 .5 .5],'edgecolor','none');
    camlight('headlight','infinite')
    axis equal
    
    scatter3(m2m_above_vpoint(:,1),m2m_above_vpoint(:,2), m2m_above_vpoint(:,3), 5, 'g')
    scatter3(surface_points(:,1),surface_points(:,2), surface_points(:,3), 5, 'b')
end

%%
% 6) now we deal with acceptance criterion of asking whether the projected
% points are outside the original mri2mesh surface
in = intriangulation(v2,f2,surface_points);

% let's figure out the correct update indices
m2m_vcount = 1:length(m2m_vindex);
m2m_vin = m2m_vcount(m2m_vindex);
update_index = m2m_vin(in==0);

% we update the v2 coordinates
v2_dilate = v2;
v2_dilate(update_index,:) = surface_points(in==0,:);

% quick visual check
if verbose
    figure;
    hold on
    camlight('headlight','infinite')
    axis equal
    P = patch('Faces',f2,'Vertices',v2_dilate,'facecolor',[.5 .5 .5],'edgecolor','none');
    set(P, 'facealpha', 0.5)
    P = patch('Faces',f2,'Vertices',v2,'facecolor',[1 0 0],'edgecolor','none');
    set(P, 'facealpha', 1)
end


% Dilation step done!
disp('Dilation step done!')

%% Goal 2: smooth mri2mesh below the plane
% Our strategy is going to be the following.

% 1) We find the most inferior point of the mri2mesh surface, define a
% coronal plane that cuts the mri2mesh surface for a anterior-inferior
% portion, which contains most of the convoluted surface edges.

% 2) We define a sub-surface using the MNE_outer_skull surface.

% 3) We dilate mri2mesh surface to MNE_outer_skull if projection results in
% points outside the original mri2mesh skull compartment.

% 4) We update the original mri2mesh surface.

%%
% 1) We find the most inferior point of the mri2mesh surface, define a
% coronal plane that cuts the mri2mesh surface for a anterior-inferior
% portion, which contains most of the convoluted surface edges.

% moved to an earlier section - we only smooth the anterior-inferior
% portion using MNE_outer_skull

%%
% 2) We define a sub-surface using the MNE_outer_skull surface.
load(fullfile(datapath, 'MNE_outer_skull.mat'))
f3 = faces;
v3 = vertices;

% Let's construct a sub surface for MNE_outer_skull in the
% anterior-inferior quadrant
vindex = v3(:,3) < zthresh & v3(:,2) > ythresh;
vcount = 1:length(vindex);
vin = vcount(vindex);

face_idx = zeros(size(f3,1), 1);
for ii = 1:size(f3,1)
    if any(ismember(f3(ii,:), vin))
        face_idx(ii) = 1;
    end
end

f3_quad = f3(logical(face_idx), :);

% we also need to clean up the indexing of this new triangulation surface
v_include = double(unique(f3_quad));
sub_v3 = v3(v_include,:);

% clean the indices contained in f1_above to correspond to sub_v1
for ii = 1:size(f3_quad,1)
    [logiout, idxout] = ismember(f3_quad(ii,:), v_include);
    assert(all(logiout), 'Some vertices called by f3_quad not in v_include!')
    f3_quad(ii,:) = idxout;
end

% quick visual check
if verbose
    figure
    P = patch('Faces',f3_quad,'Vertices',sub_v3,'facecolor',[0.5 0.5 0.5],'edgecolor','none');
    hold on
    camlight('headlight','infinite')
    axis equal
    camorbit(0, 180)
    camlight('headlight')
end

%%
% 3) We dilate mri2mesh surface to MNE_outer_skull if projection results in
% points outside the original mri2mesh skull compartment.

m2m_quad_vindex = v2(:,3) < zthresh & v2(:,2) > ythresh;
m2m_quad_vpoint = v2(m2m_quad_vindex,:);

% quick visual check
if verbose
    scatter3(m2m_quad_vpoint(:,1),m2m_quad_vpoint(:,2), m2m_quad_vpoint(:,3), 10, 'g')
end

% We propose to dilate these quadrant points to MNE_outer_skull
[ distances_1, surface_points_1 ] = ...
    point2trimesh('Faces', f3_quad, 'Vertices', sub_v3, 'QueryPoints', m2m_quad_vpoint, 'Algorithm', 'parallel');

% quick visual check
if verbose
    scatter3(surface_points_1(:,1),surface_points_1(:,2), surface_points_1(:,3), 5, 'b')
end

% Based on whether outside the original mri2mesh compartment, we accept the
% dilation
in = intriangulation(v2,f2,surface_points_1);

% visual check of accepted points
accept_p = surface_points_1(in==0,:);
if verbose
    scatter3(accept_p(:,1),accept_p(:,2), accept_p(:,3), 5, 'r')
end

% let's figure out the correct update indices
m2m_quad_vcount = 1:length(m2m_quad_vindex);
m2m_quad_vin = m2m_quad_vcount(m2m_quad_vindex);
update_index_quad = m2m_quad_vin(in==0);

% we update the v2_dilate coordinates
v2_dilate(update_index_quad,:) = surface_points_1(in==0,:);

% quick visual check
figure;
hold on
camlight('headlight','infinite')
axis equal
P = patch('Faces',f2,'Vertices',v2_dilate,'facecolor',[.5 .5 .5],'edgecolor','none');
set(P, 'facealpha', 0.5)
title('Pre-artifact Fixes', 'FontSize', 20)

%%
% 4) There is an artifact at the base of the brainstem. We need to
% somehow correct it. Hmm.. let's see what we can do.

%%%%%%%%%%%%
% Update on 05/12/2020:
% I dont think this step is needed, since we just projected all mri2mesh
% points to the MNE_outer_skull surface. I don't think any more points will
% still be inside MNE_outer_skull other than due to small numerical errors.
%%%%%%%%%%%%

% % Let's visualize the current mri2mesh surface vertices. My intuition is
% % that we can use MNE_outer_skull to fix this artifact.
% m2m_quad_vpoint_new = v2_dilate(m2m_quad_vindex,:);
%
% figure
% P = patch('Faces',f3_quad,'Vertices',sub_v3,'facecolor',[0.5 0.5 0.5],'edgecolor','none');
% hold on
% camlight('headlight','infinite')
% axis equal
% camorbit(0, 180)
% camlight('headlight')
% scatter3(m2m_quad_vpoint_new(:,1),m2m_quad_vpoint_new(:,2), m2m_quad_vpoint_new(:,3), 5, 'b')
%
% % Let's project all points that are inside MNE_outer_skull to MNE_outer_skull
% in = intriangulation(v3,f3,m2m_quad_vpoint_new);
%
% p_inside = m2m_quad_vpoint_new(in, :);
%
% % quick visual check
% scatter3(p_inside(:,1),p_inside(:,2), p_inside(:,3), 5, 'r')
%
% % we project these points to MNE_outer_skull
% [ ~, p_inside_pulled ] = ...
%     point2trimesh('Faces', f3_quad, 'Vertices', sub_v3, 'QueryPoints', p_inside, 'Algorithm', 'parallel');
%
% % figure out the correct update index
% update_index_new = m2m_quad_vin(in==1);
%
% % update the surface vertex coordinates
% v2_dilate(update_index_new,:) = p_inside_pulled;
%
% % quick visual check
% figure;
% hold on
% camlight('headlight','infinite')
% axis equal
% P = patch('Faces',f2,'Vertices',v2_dilate,'facecolor',[.5 .5 .5],'edgecolor','none');
% set(P, 'facealpha', 0.5)


% Smoothing step done!
disp('Smoothing step done!')

%% Goal 3: MeshFix, Smoothing, and Resample
% Now that we have fixed the skull surface and stitched things together, we
% need to clean the surface a bit to get a good triangulation mesh. There
% are three operations we need to do, MeshFix, Smooth, and Resample.

% close all

% visualize initial surface
if verbose
    figure;
    hold on
    camlight('headlight','infinite')
    axis equal
    P = patch('Faces',f2,'Vertices',v2_dilate,'facecolor',[.5 .5 .5],'edgecolor','none');
    camorbit(0, 180)
    camlight('headlight')
    set(P, 'facealpha', 0.5)
    title('Pre-artifact Fixes', 'FontSize', 20)
end

% Step 1 - Gentle smooth
% We do smooth first because there are sharp transitions that we don't
% want to resample on
[conn,connnum,count]=meshconn(f2, length(v2_dilate));
p=smoothsurf(v2_dilate,[],conn,5,0.5,'lowpass',0.5);

% quick visual check
if verbose
    figure
    P = patch('Faces',f2,'Vertices',p,'facecolor',[.5 .5 .5],'edgecolor','none');
    set(P, 'facealpha', 0.5)
    hold on
    camlight('headlight','infinite')
    axis equal
    camorbit(0, 180)
    camlight('headlight')
    title('Slight smoothing', 'FontSize', 20)
end

%%
% Step 2 - Repair
[node,elem]=meshcheckrepair(p,f2,'dup');
[node,elem]=meshcheckrepair(node,elem,'isolated');
[node,elem]=meshcheckrepair(node,elem,'deep');
[node,elem]=meshcheckrepair(node,elem,'meshfix');
[node,elem]=meshcheckrepair(node,elem,'open');
[node,elem]=meshcheckrepair(node,elem,'intersect');
elem=surfaceclean(elem,node);

% Step 3 - More smoothing
[conn,connnum,count]=meshconn(elem, length(node));
snode=smoothsurf(node,[],conn,200,0.8,'lowpass',0.5);

% quick visual check
if verbose
    figure
    P = patch('Faces',elem,'Vertices',snode,'facecolor',[.5 .5 .5],'edgecolor','none');
    set(P, 'facealpha', 0.5)
    hold on
    % scatter3(snode(:,1), snode(:,2), snode(:,3), 5, 'b', 'filled')
    camlight('headlight','infinite')
    axis equal
    camorbit(0, 180)
    camlight('headlight')
    title('Quick fix and more smoothing', 'FontSize', 20)
end

%%
% Step 4 - Fix sharp transition problems
% due to some sharp transitions and smoothing, we have some weird sharp
% transition artifacts at the brainstem position of the skull surface.

% How do we detect and fix these artifacts?
if verbose
    figure
    P = patch('Faces',elem,'Vertices',snode,'facecolor',[.5 .5 .5],'edgecolor','none');
    set(P, 'facealpha', 0.5)
    hold on
    scatter3(snode(:,1), snode(:,2), snode(:,3), 5, 'b', 'filled')
    camlight('headlight','infinite')
    axis equal
    camorbit(0, 180)
    camlight('headlight')
end

% We can try to use curvature around vertices as a method to detect sharp
% transition artifacts.

v4 = snode; f4 = elem;

% Let's first compute a subsurface of the bottom of the skull
% vindex = v4(:,3) < zthresh & v4(:,2) > ythresh;
vindex = v4(:,3) < zthresh;
vcount = 1:length(vindex);
vin = vcount(vindex);

face_idx = zeros(size(f4,1), 1);
for ii = 1:size(f4,1)
    if any(ismember(f4(ii,:), vin))
        face_idx(ii) = 1;
    end
end

f4_quad = f4(logical(face_idx), :);

% we also need to clean up the indexing of this new triangulation surface
v_include = double(unique(f4_quad));
sub_v4 = v4(v_include,:);

% clean the indices contained in f1_above to correspond to sub_v1
for ii = 1:size(f4_quad,1)
    [logiout, idxout] = ismember(f4_quad(ii,:), v_include);
    assert(all(logiout), 'Some vertices called by f4_quad not in v_include!')
    f4_quad(ii,:) = idxout;
end

% quick visual check
if verbose
    figure
    P = patch('Faces',f4_quad,'Vertices',sub_v4,'facecolor',[0.5 0.5 0.5],'edgecolor','none');
    hold on
    camlight('headlight','infinite')
    axis equal
    camorbit(0, 180)
    camlight('headlight')
end

FV = struct;
FV.faces = f4_quad;
FV.vertices = sub_v4;

getderivatives=0;
[PrincipalCurvatures,PrincipalDir1,PrincipalDir2,FaceCMatrix,VertexCMatrix,Cmagnitude]=GetCurvatures(FV ,getderivatives);
GausianCurvature=PrincipalCurvatures(1,:).*PrincipalCurvatures(2,:);

if verbose
    figure('name','Triangle Mesh Curvature','numbertitle','off','color','w');
    hold on
    colormap cool
    caxis([min(GausianCurvature) max(GausianCurvature)]); % color overlay the gaussian curvature
    mesh_h=patch(FV,'FaceVertexCdata',GausianCurvature','facecolor','interp','edgecolor','interp','EdgeAlpha',0.2);
    %set some visualization properties
    set(mesh_h,'ambientstrength',0.35);
    axis off
    axis equal
    camlight();
    lighting phong
    colorbar();
end

% awesome! This correctly identifies the peak of the artifact. We now just
% need to iteratively get the vertices, delete the none-boundary ones and
% create a central node to re-connect them.

%% Well, this is annoying. What if we have more than one topological error?
% % artifact_v_idx = v_include(GausianCurvature>70); % this threshold might need to change!
% % assert(length(artifact_v_idx) == 1, 'More than one artifact center detected! Please check!')

% we need to figure out how many times we need to do this topological
% correction.
possible_artifact_v_idx = v_include(GausianCurvature>50);

% quick visual check
if verbose
    scatter3(v4(possible_artifact_v_idx, 1), v4(possible_artifact_v_idx, 2), v4(possible_artifact_v_idx, 3), 100, 'b', 'filled')
end

% construct a list of topological errors
topolist = {};

for ii = 1:length(possible_artifact_v_idx)
    topo = struct;
    
    current_v = possible_artifact_v_idx(ii);
    
    already_fixed = 0;
    
    for jj = 1:length(topolist)
        if ismember(current_v, topolist{jj}.involved_v)
            already_fixed = jj;
        end
    end
    
    if already_fixed == 0
        [vertex_index, face_index] = mesh_findneighbor(f4, current_v, 6);
        topo.high_artifact_v = current_v;
        topo.face_index = face_index;
        topo.involved_v = unique(f4(face_index,:));
        topolist{length(topolist)+1} = topo;
    else % one of previous topological error already includes it
        [vertex_index, face_index] = mesh_findneighbor(f4, current_v, 6);
        topo = topolist{already_fixed};
        topo.high_artifact_v = [topo.high_artifact_v, current_v];
        topo.face_index = unique([topo.face_index; face_index]);
        topo.involved_v = unique(f4(topo.face_index,:));
        topolist{already_fixed} = topo;
    end
end
% make sure the topolist vertices do not overlap
for ii = 1:length(topolist)
    for jj = 1:length(topolist)
        if ii ~= jj
            assert(~any(ismember(topolist{ii}.involved_v, topolist{jj}.involved_v)), 'some overlaps in vertices in topological error lists')
        end
    end
end

%%
% Now we need to find a way to get the neighboring faces and vertices of
% each artifact center vertex, and we need to handle it for multiple
% topological errors.

% loop through all topological errors and correct them
final_F = f4;
final_V = v4;
delete_face_index = [];
for tt = 1:length(topolist)
    
    current_topo = topolist{tt};
    face_index = current_topo.face_index;
    
    artifact_f = f4(face_index,:);
    
    % quick visual check
    if verbose
        figure
        P = patch('Faces',artifact_f,'Vertices',v4,'facecolor',[0.5 0.5 0.5],'edgecolor','none');
        hold on
        camlight('headlight','infinite')
        axis equal
        camorbit(0, 180)
        camlight('headlight')
    end
    
    % get the edge vertices
    openedge=surfedge(artifact_f);
    edge_vertex = v4(unique(openedge), :);
    if verbose
        scatter3(edge_vertex(:,1), edge_vertex(:,2), edge_vertex(:,3), 50, 'r', 'filled')
    end
    
    %% GIBBON seems to be a nice toolbox that we can use!
    
    % % 1) let's grab the vertices and create as a separate surface
    % sub_v5 = v4(vertex_index, :);
    % f5_artifact = artifact_f;
    %
    % % clean the indices contained in f1_above to correspond to sub_v1
    % for ii = 1:size(f5_artifact,1)
    %     [logiout, idxout] = ismember(f5_artifact(ii,:), vertex_index);
    %     assert(all(logiout), 'Some vertices called by f5_artifact not in vertex_index!')
    %     f5_artifact(ii,:) = idxout;
    % end
    %
    % % quick visual check
    % figure
    % P = patch('Faces',f5_artifact,'Vertices',sub_v5,'facecolor',[0.5 0.5 0.5],'edgecolor','none');
    % hold on
    % camlight('headlight','infinite')
    % axis equal
    % camorbit(0, 180)
    % camlight('headlight')
    
    % 2) create a flat delaunay surface using boundary points
    % add in centroid of boundary vertices to make it flat
    edge_vertex = [edge_vertex; mean(edge_vertex)];
    
    % create delaunay 2D surface on the edge_vertices
    DT = delaunay(edge_vertex(:,1), edge_vertex(:,2));
    % figure
    % triplot(DT,edge_vertex(:,1), edge_vertex(:,2));
    
    if verbose
        figure
        patch('Faces', DT, 'Vertices',edge_vertex,'facecolor',[.5 .5 .5], 'edgecolor', 'k');
        rotate3d on
        axis equal
        camlight('headlight','infinite')
        hold on
        scatter3(edge_vertex(:,1), edge_vertex(:,2), edge_vertex(:,3), 50, 'r', 'filled')
    end
    
    % there are some extra edges created in forming the flat surface. Let's get
    % rid of them to get back the original edges
    original_edge_index = unique(openedge);
    [~, original_edge] = ismember(openedge, original_edge_index);
    
    % we use a while loop to trim the edge faces until we get back the orignal
    % one!
    deleteidx = 0;
    while ~isempty(deleteidx)
        new_edge=surfedge(DT);
        
        % visualize the new edges
        for ii = 1:length(new_edge)
            plotV(edge_vertex(new_edge(ii,:),:),'r.-','MarkerSize',25,'LineWidth',5);
        end
        
        % visualize the old edges
        for ii = 1:length(original_edge)
            plotV(edge_vertex(original_edge(ii,:),:),'b.-','MarkerSize',25,'LineWidth',5);
        end
        
        % find faces that contain these edges
        deleteidx = [];
        for ii = 1:size(DT,1)
            local_edge = surfedge(DT(ii,:));
            isedge_face = false;
            local_edge_idx = [];
            for jj = 1:size(local_edge,1)
                e1 = new_edge(:,1) == local_edge(jj,1) & new_edge(:,2) == local_edge(jj,2);
                e2 = new_edge(:,1) == local_edge(jj,2) & new_edge(:,2) == local_edge(jj,1);
                if any([e1;e2])
                    isedge_face = true;
                    local_edge_idx = [local_edge_idx; jj];
                end
            end
            if isedge_face % is an edge in the new flat surface
                % now we need to check if it is also an edge in the original edge list
                for jj = 1:size(local_edge_idx,1)
                    current_local_edge = local_edge(local_edge_idx(jj),:);
                    e1 = original_edge(:,1) == current_local_edge(1) & original_edge(:,2) == current_local_edge(2);
                    e2 = original_edge(:,1) == current_local_edge(2) & original_edge(:,2) == current_local_edge(1);
                    if ~any([e1;e2]) % doesn't exist in original edge list
                        deleteidx = [deleteidx, ii];
                    end
                end
            end
        end
        DT(deleteidx, :) = [];
        % patch('Faces', DT(deleteidx, :), 'Vertices',edge_vertex,'facecolor',[1 .5 .5], 'edgecolor', 'k');
    end
    
    % visualize the trimmed result - looking good!
    if verbose
        figure
        patch('Faces', DT, 'Vertices',edge_vertex,'facecolor',[.5 .5 .5], 'edgecolor', 'k');
        rotate3d on
        axis equal
        camlight('headlight','infinite')
        hold on
        scatter3(edge_vertex(:,1), edge_vertex(:,2), edge_vertex(:,3), 50, 'r', 'filled')
    end
    
    % refine this flat surface
    [F,V]=subtri(DT,edge_vertex,1);
    for ii = 1:3
        [F,V]=subtri(F,V,1);
    end
    
    % remove the centroid from the edge_vertex list
    edge_vertex(end, :) = [];
    
    % quick visual check
    if verbose
        figure
        P = patch('Faces',artifact_f,'Vertices',v4,'facecolor',[0.5 0.5 0.5],'edgecolor','none');
        set(P, 'facealpha', 0.5)
        hold on
        P = patch('Faces', F, 'Vertices',V,'facecolor',[0 0 1], 'edgecolor', 'none');
        set(P, 'facealpha', 0.5)
        scatter3(edge_vertex(:,1), edge_vertex(:,2), edge_vertex(:,3), 50, 'r', 'filled')
        camlight('headlight','infinite')
        axis equal
        camorbit(0, 180)
        camlight('headlight')
    end
    
    % 3) geodesic remesh with boundary points preserved
    % Control parameters
    pointSpacing=2; %Desired point spacing
    % Estimate number of points required given point spacing
    numPointsInput=size(V,1); %Number of points in the original data
    [A]=patch_area(F,V); %Areas of current faces
    totalArea=sum(A(:)); %Total area
    l=sqrt(totalArea); %Width or length of square with same size
    np=round((l./pointSpacing).^2); %Point spacing for mesh in virtual square
    indListSelect = 1:size(edge_vertex,1);
    
    % geodesic remesh
    clear optionStruct
    optionStruct.toleranceLevel=0; %Tolerance for convergence
    optionStruct.waitBarOn=1; %Turn on/off waitbar
    [Fn,Vn,seedIndex,indSeeds,d]=remeshTriSurfDistMap(F,V,numel(indListSelect)+np,indListSelect,optionStruct); %distance based marching
    
    % quick visual check
    if verbose
        figure
        hold on
        P = patch('Faces', Fn, 'Vertices',Vn,'facecolor',[1 0 0], 'edgecolor', 'k');
        set(P, 'facealpha', 0.5)
        scatter3(edge_vertex(:,1), edge_vertex(:,2), edge_vertex(:,3), 50, 'r', 'filled')
        camlight('headlight','infinite')
        axis equal
        camorbit(0, 180)
        camlight('headlight')
    end
    
    %%
    % Now that we have a more evenly sampled flat surface, we need to
    % incorporate it back into the original skull surface
    
    % [1] Let's make sure the edge vertices are all present. If not, one might
    % need to adjust the desired point spacing parameter during geodesic
    % resampling
    for ii = 1:size(edge_vertex, 1)
        assert(abs(norm(Vn(ii,1) - edge_vertex(ii,1)) - 0) < 10^-4)
    end
    
    % [2] We need to update faces and vertices in the skull surface
    % remove the original faces
    % final_F = f4;
    % final_F(face_index,:) = [];
    assert(~any(ismember(face_index, delete_face_index)), 'a face is deleted before. Something is wrong!')
    delete_face_index = [delete_face_index; face_index];
    
    % append new vertices to the list
    % orig_vn = size(v4,1);
    % final_V = [v4; Vn(size(edge_vertex,1)+1:end, :)];
    orig_vn = size(final_V,1);
    final_V = [final_V; Vn(size(edge_vertex,1)+1:end, :)];
    
    % adjust face indices for these new vertices
    updated_Fn = Fn;
    for ii = size(edge_vertex,1)+1:size(Vn,1)
        updated_Fn(Fn(:,1)==ii, 1) = orig_vn + ii - size(edge_vertex,1);
        updated_Fn(Fn(:,2)==ii, 2) = orig_vn + ii - size(edge_vertex,1);
        updated_Fn(Fn(:,3)==ii, 3) = orig_vn + ii - size(edge_vertex,1);
    end
    
    % adjust face indices for boundary vertices that already exist in the old
    % vertex list
    for ii = 1:size(original_edge_index,1)
        updated_Fn(Fn(:,1)==ii, 1) = original_edge_index(ii);
        updated_Fn(Fn(:,2)==ii, 2) = original_edge_index(ii);
        updated_Fn(Fn(:,3)==ii, 3) = original_edge_index(ii);
    end
    
    % append Fn to final_F
    final_F = [final_F; updated_Fn];
    
    % quick visual check
    if verbose
        figure
        P = patch('Faces',updated_Fn,'Vertices',final_V,'facecolor',[0.5 0.5 0.5],'edgecolor','none');
        hold on
        set(P, 'facealpha', 0.5)
        scatter3(edge_vertex(:,1), edge_vertex(:,2), edge_vertex(:,3), 50, 'r', 'filled')
        camlight('headlight','infinite')
        axis equal
        camorbit(0, 180)
        camlight('headlight')
    end
    
end

% now we need to delete all the useless faces prior to topological error
% correction
final_F(delete_face_index,:) = [];

% And we get rid of un-used vertices
[final_V, final_F]=removeisolatednode(final_V,final_F);

% quick visual check
if verbose
    figure
    P = patch('Faces',final_F,'Vertices',final_V,'facecolor',[0.5 0.5 0.5],'edgecolor','none');
    hold on
    set(P, 'facealpha', 0.5)
    % scatter3(edge_vertex(:,1), edge_vertex(:,2), edge_vertex(:,3), 50, 'r', 'filled')
    camlight('headlight','infinite')
    axis equal
    camorbit(0, 180)
    camlight('headlight')
end

% Awesome! We've fixed all topological errors!
disp('Topological errors fixing done!')

%%
% Step 5 - A bit more smoothing
[conn,connnum,count]=meshconn(final_F, length(final_V));
final_V=smoothsurf(final_V,[],conn,100,0.8,'lowpass',0.5);

% manual visual check before re-meshing
figure
P = patch('Faces',final_F,'Vertices',final_V,'facecolor',[0.5 0.5 0.5],'edgecolor','none');
hold on
set(P, 'facealpha', 0.5)
% scatter3(edge_vertex(:,1), edge_vertex(:,2), edge_vertex(:,3), 50, 'r', 'filled')
camlight('headlight','infinite')
axis equal
camorbit(0, 180)
camlight('headlight')
title('Last round of smoothing', 'FontSize', 20)

% manual sanity check: do the planes and surfaces look reasonable?
checkok =  input('Does the morphed skull surface look ok? (y/n): ', 's');
if strcmpi(checkok, 'n')
    save(fullfile(datapath, 'mesh_morph_skull_failed_checkpoint2_workspace'))
    error('Manual check failed. Workspace is saved. Please debug...')
end

%%
% Step 6 - Resample
% we can use the geodesic remesh function in GIBBON toolbox to resample the
% fixed water-tight surface into desired iso-distances.

F = final_F; V = final_V;

%Plot settings
cMapDist=flipud(igviridis(250));
[cMapIndices,scrambleIndices]=scramble(viridis(250),1); %Colormap

faceAlpha1=1;
faceAlpha2=0.65;
fontSize=25;
markerSize=50;
lineWidth=4;
scatterSize=65;

% Visualize the original mesh in GIBBON
if verbose
    cFigure; hold on;
    title('Input mesh')
    gpatch(F,V,'gw');
    axisGeom(gca,fontSize);
    camlight headlight;
    drawnow;
end

% remesh to round-up to nearest 10000 vertices

% Compute distances on mesh
numSeeds=roundn(size(V,1), 4);

%Option set
[~,indStart]=min(V(:,1)); %Index of the start point
optionStruct.toleranceLevel=0; %Tolerance for convergence
optionStruct.waitBarOn=1; %Turn on/off waitbar

%Compute distances on mesh description - this will take a LONG LONG LONG
%time due to geodesic distance calculation
[Fn,Vn,seedIndex,indSeeds,d]=remeshTriSurfDistMap(F,V,numSeeds,indStart,optionStruct); %distance based marching
[~,~,ind2]=unique(seedIndex);

% quick visual check
figure
P = patch('Faces',Fn,'Vertices',Vn,'facecolor',[0.5 0.5 0.5],'edgecolor','none');
hold on
set(P, 'facealpha', 0.5)
% scatter3(edge_vertex(:,1), edge_vertex(:,2), edge_vertex(:,3), 50, 'r', 'filled')
camlight('headlight','infinite')
axis equal
camorbit(0, 180)
camlight('headlight')

% manual sanity check: do the planes and surfaces look reasonable?
checkok =  input('Does the final skull surface look ok? (y/n): ', 's');
if strcmpi(checkok, 'n')
    save(fullfile(datapath, 'mesh_morph_skull_failed_checkpoint3_workspace'))
    error('Manual check failed. Workspace is saved. Please debug...')
end

%% Save the hardwork result!
% Save data to file
faces = Fn;
vertices = Vn;

save(fullfile(datapath, 'morphed_skull.mat'), 'faces', 'vertices')

% Done!

end
