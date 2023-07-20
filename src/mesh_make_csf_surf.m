function [ ] = mesh_make_csf_surf(datapath, verbose)
%% Make csf surface 
% Sometimes the mri2mesh extracts a perfect looking csf surface. Other
% times, the surface looks quite crappy at the base of the skull, giving
% rise to a huge volume of extra csf that acts like a volume conductor to
% artifically bridge activities across temporal lobes.

% To avoid this kind of problem for all subjects, and to have a consistent
% algorithm to get an usable csf BEM surface, we will use this function
% to construct a good csf surface.

if nargin < 2
    verbose = false;
end

%% Overall processing plan:
% headreco pipeline reliably produces a usable csf surface, however it has
% not been adjusted according to the pial surface (which we are using the
% one from mri2mesh pipeline), nor does it respect the hard cutoff at the
% bottom of the skull since mri2mesh only treats the top half of the head.
% Hence we will do some editing on the headreco csf surface to make it into
% a good one to use for 4-shell forward modeling. 

% Step 1: Identify the bottom of the spinal cord region of the headreco csf
% surface, back trace until we find a series of vertices intersecting with
% the mri2mesh_csf surface, make a cut and close it off with a flat plane. 


% Step 2: compare with pial surface and cerebellum surface. Locally dilate
% the headreco surface in places where the pial surface lies outside or
% where the pial surface gets ridiculously close to the csf layer < 0.05mm.


% Step 3: Sanity check that this new csf surface should be outside all pial
% surface, but be inside the skull surfaces. If not, correct them. 

%% Path configuration
% addpath(genpath('./'));

% datapath should point to the final_structure folder
% datapath = '/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/MNE/subject_data/SP011_ah/simnibs_pipe/bem_surfaces/final_structure';

tic

%% Load everything first 

% Load headreco bone surface
load(fullfile(datapath, 'headreco_csf.mat'), 'faces', 'vertices')
f1 = faces; %#ok<*NODEF>
v1 = vertices;

% Load mri2mesh gm surface
load(fullfile(datapath, 'mri2mesh_gm_fixed.mat'), 'faces', 'vertices')
f2 = faces;
v2 = vertices;

% Load mri2mesh cerebellum surface
load(fullfile(datapath, 'mri2mesh_cerebellum_fixed.mat'), 'faces', 'vertices')
f3 = faces;
v3 = vertices;

% Load mri2mesh csf surface
load(fullfile(datapath, 'mri2mesh_csf_fixed.mat'), 'faces', 'vertices')
f4 = faces;
v4 = vertices;

% Load mri2mesh bone surface
load(fullfile(datapath, 'mri2mesh_skull_fixed.mat'), 'faces', 'vertices')
f5 = faces;
v5 = vertices;

% Load headreco bone surface
load(fullfile(datapath, 'headreco_bone.mat'), 'faces', 'vertices')
f6 = faces;
v6 = vertices;

%% Step 1: Fix skull base protrusion
% This is a tough problem, as we need to figure out a way to find incision
% plane 

% Here's an easy way, we know headreco protrudes out, so I can just grab
% the out points from headreco - there's an assumption here that the
% headreco csf should only be out of mri2mesh_skull at the spinal cord
% region. If it is also out in additional places, this heuristic simply
% falls apart. 

in = intriangulation(v5,f5,v1);
all_out_vertex_index = find(in==0);

% find all csf faces involving these vertices 
all_nob_face_index = any(ismember(f1, all_out_vertex_index), 2);

% Updated on 06/08/2020: occasionally the mri2mesh csf is fairly good
% around the skull base and therefore there are multiple small places where
% the headreco csf protrudes out from the mri2mesh skull. Now I update to
% find the largest connected components for incision. 

% Find the largest component and update indices for incision
f1_split = splitFV(f1(all_nob_face_index,:), v1);
[~,largest_k] = max(cellfun(@(x) size(x,1), {f1_split.vertices}));
out_vertex_index = f1_split(largest_k).vertIDs;
% keep only the original ones that were determined to be out
out_vertex_index = out_vertex_index(ismember(out_vertex_index, all_out_vertex_index));
nob_face_index = any(ismember(f1, out_vertex_index), 2);

% find the edge faces and vertices 
openedge=surfedge(f1(nob_face_index,:));
edge_vertex_index = unique(openedge);

% we erode upwards for (upto) 5 iterations 
start_iter = 5;
erode_complete = false;
while ~erode_complete
    invovled_face = find(nob_face_index==1);
    for ii = 1:length(edge_vertex_index)
        [vertex_index, face_index] = mesh_findneighbor(f1, edge_vertex_index(ii), start_iter);
        invovled_face = [invovled_face; face_index];
    end
    % we want to make sure all of the edge indices are outside of mri2mesh csf,
    % otherwise it might cut the cerebellum base at a weird angle that's hard
    % to fix.
    all_vertex_index = v1(unique(f1(invovled_face,:)),:);
    in = intriangulation(v4,f4,all_vertex_index);
    if sum(in) == 0
        erode_complete = true;
    else
        start_iter = start_iter-1;  % reduce iteration by one
        if start_iter <= 0
            % Updated on 08/10/2022: this means even eroding by 1 iteration
            % of neighbors puts some vertices into the mri2mesh csf
            % surface. This is acceptable, if all the inpoints are on the
            % anterior side, meaning it doesn't have a chance to cut the
            % cerebellum.
            mean_point = mean(all_vertex_index, 1);
            in_points = all_vertex_index(in==1,:);
            
            % make sure all the points inside the csf surface have larger
            % y-axis values than the mean, suggesting they are on the
            % anterior side. If so, allow exit the while loop
            if all(in_points(:, 2) > mean_point(2))
                erode_complete = true;
            else
                error('Could not find an erosion order without entering mri2mesh csf. Please check.')
            end
        end
    end
end
invovled_face = unique(invovled_face);

% visualize incision plane 
if verbose
    figure;
    hold on
    patch('Faces',f1,'Vertices',v1,'facecolor',[.5 .5 .5],'edgecolor','none');
    p = patch('Faces',f5,'Vertices',v5,'facecolor',[1 0 0],'edgecolor','none');
    set(p, 'facealpha', 0.25)
    patch('Faces',f1(all_nob_face_index,:),'Vertices',v1,'facecolor',[0 1 1],'edgecolor','none');
    patch('Faces',f1(invovled_face,:),'Vertices',v1,'facecolor',[0 0 1],'edgecolor','none');
    patch('Faces',f1(nob_face_index,:),'Vertices',v1,'facecolor',[0 1 0],'edgecolor','none');
    title('Incision plane', 'FontSize', 20)
    camlight('headlight','infinite')
    axis equal
    camorbit(0, 180)
    camlight('headlight')
    camorbit(0, 180)
    camorbit(0, 270)
    rotate3d on
end

%% Start cutting at the incision plane and close off the plane
% find the new edge faces
openedge = surfedge(f1(invovled_face,:));
edge_vertex = v1(unique(openedge),:);

% We form a closed plane using these edge vertices
edge_vertex = [edge_vertex; mean(edge_vertex)];

% create delaunay 2D surface on the edge_vertices
DT = delaunay(edge_vertex(:,1), edge_vertex(:,2));

if verbose
    figure
    surfH = patch('Faces', DT, 'Vertices',edge_vertex,'facecolor',[.5 .5 .5], 'edgecolor', 'k');
    rotate3d on
    axis equal
    camlight('headlight','infinite')
    hold on
    scatter3(edge_vertex(:,1), edge_vertex(:,2), edge_vertex(:,3), 50, 'r', 'filled')
    title('Fixing Edges', 'FontSize', 20)
end

% there are some extra edges created in forming the flat surface. Let's get
% rid of them to get back the original edges
original_edge_index = unique(openedge);
[~, original_edge] = ismember(openedge, original_edge_index);

% Make sure all of the open edges are contained in the DT faces. If not, we
% need to add some faces
missing_edges = [];
for ii = 1:size(original_edge)
    edge_in_face = sum((DT == original_edge(ii, 1)) + (DT == original_edge(ii, 2)), 2);
    if sum(edge_in_face == 2) == 0
        missing_edges = [missing_edges; original_edge(ii, :)];
    end
end
for ii = 1:size(missing_edges, 1)
    edge_vertex_A = missing_edges(ii, 1);
    edge_vertex_B = missing_edges(ii, 2);
    
    % find the other neighbor of A and make a face with B, vice versa for B
    neighbor_edge_index = find(sum((original_edge == edge_vertex_A) + (original_edge == edge_vertex_B), 2) == 1);
    assert(length(neighbor_edge_index) == 2, 'Not exactly two neighbors found in the edge list. Please double check!')

    % add two new faces with edge_vertex_A/B and the two side neighbors 
    for jj = 1:length(neighbor_edge_index)
        neighbor_edge_vertex = original_edge(neighbor_edge_index(jj), :);
        
        if any(ismember(neighbor_edge_vertex, edge_vertex_A))
            added_face = [neighbor_edge_vertex, edge_vertex_B];
        else
            added_face = [neighbor_edge_vertex, edge_vertex_A];
        end
        
        DT = [DT; added_face];
    end
end

% we use a while loop to trim the edge faces until we get back the orignal
% edge set
deleteidx = 0;
counter = 0;
while ~isempty(deleteidx)
    % count how many times the while loop was run
    counter = counter + 1;
    if verbose; title(['Trimming in progress: loop ', num2str(counter)]); end
    
    new_edge=surfedge(DT);
    
    % visualize the new edges
    new_edge_lines = [];
    for ii = 1:length(new_edge)
        if verbose 
            H = plotV(edge_vertex(new_edge(ii,:),:),'r.-','MarkerSize',25,'LineWidth',5);
            new_edge_lines = [new_edge_lines, H]; %#ok<*AGROW>
        end
    end
    
    % visualize the old edges
    old_edge_lines = [];
    for ii = 1:length(original_edge)
        if verbose
            H = plotV(edge_vertex(original_edge(ii,:),:),'b.-','MarkerSize',25,'LineWidth',5);
            old_edge_lines = [old_edge_lines, H];
        end
    end
    drawnow
    
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
    
    % delete the extra faces 
    % patch('Faces', DT(deleteidx, :), 'Vertices',edge_vertex,'facecolor',[1 .5 .5], 'edgecolor', 'k');
    DT(deleteidx, :) = [];
    
    % Remove the lines for an updated contour visualization 
    if verbose
        pause(2)
        for ii = 1:length(new_edge_lines)
            delete(new_edge_lines(ii))
        end
        for ii = 1:length(old_edge_lines)
            delete(old_edge_lines(ii))
        end
        
        % update the surface patch display as well
        delete(surfH)
        surfH = patch('Faces', DT, 'Vertices',edge_vertex,'facecolor',[.5 .5 .5], 'edgecolor', 'k');
        drawnow
    end
end

% visualize the trimmed result - looking good!
if verbose
    figure
    patch('Faces', DT, 'Vertices', edge_vertex,'facecolor',[.5 .5 .5], 'edgecolor', 'k');
    rotate3d on
    axis equal
    camlight('headlight','infinite')
    hold on
    scatter3(edge_vertex(:,1), edge_vertex(:,2), edge_vertex(:,3), 50, 'r', 'filled')
    title('After Edge Trimming', 'FontSize', 20)
end

% refine this flat surface
[F,V]=subtri(DT,edge_vertex,1);
for ii = 1:4
    [F,V]=subtri(F,V,1);
end

% remove the centroid from the edge_vertex list
edge_vertex(end, :) = [];

if verbose
    figure
    patch('Faces',f1(invovled_face,:),'Vertices',v1,'facecolor',[0 0 1],'edgecolor','none');
    hold on
    P = patch('Faces', F, 'Vertices',V,'facecolor',[1 0.5 0.5], 'edgecolor', 'k');
    set(P, 'facealpha', 0.5)
    scatter3(edge_vertex(:,1), edge_vertex(:,2), edge_vertex(:,3), 50, 'r', 'filled')
    camlight('headlight','infinite')
    axis equal
    camorbit(0, 180)
    camlight('headlight')
    title('Closed off surface', 'FontSize', 20)
end

% geodesic remesh with boundary points preserved
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
    title('Remeshed closing surface', 'FontSize', 20)
end

% excellent, we've successfully closed off the csf surface after cutting 

%% Now we need to incorporate it back into the original csf surface
% [1] Let's make sure the edge vertices are all present. If not, one might
% need to adjust the desired point spacing parameter during geodesic
% resampling
for ii = 1:size(edge_vertex, 1)
    assert(abs(norm(Vn(ii,1) - edge_vertex(ii,1)) - 0) < 10^-4)
end

% Updated on 06/08/2020: there's no need to do this as we will check for gm
% intersection and correct for it later. Right now it's projecting the
% incision plane onto mri2mesh_csf, which often gets it wrong. So removing
% it from this part improves stability.

% [2] if some of these Vn points are inside mri2mesh_csf layer, then there
% is a chance that it will intersect with the gm surface. So let's pull
% them out to mri2mesh_csf surface if they are inside 
% in = intriangulation(v4,f4,Vn); 
% if verbose; scatter3(Vn(in,1), Vn(in,2), Vn(in,3), 50, 'b', 'filled'); end
% % project to the mri2mesh csf surface 
% [ distances, surface_points ] = ...
%     point2trimesh('Faces', f4, 'Vertices', v4, 'QueryPoints', Vn(in,:), 'Algorithm', 'parallel');
% % update these points 
% Vn(in,:) = surface_points;

% [3] We need to update faces and vertices in the csf surface
% remove the original faces
final_F = f1;
final_F(invovled_face,:) = [];
% append new vertices to the list
orig_vn = size(v1,1);
final_V = [v1; Vn(size(edge_vertex,1)+1:end, :)];

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
    title('After correcting indices added to original mesh', 'FontSize', 20)
end

% And we get rid of un-used vertices
[final_V, final_F]=removeisolatednode(final_V,final_F);

%% Update the csf surface 
v1 = final_V; 
f1 = final_F;

% Step 1 Done. 
disp('Step 1 Done.')

%% Step 2: check out-of-bound pial and cerebellum points
% find out pial surface points that are outside of the csf surface 
in = intriangulation(v1,f1,v2);
while_loop_num = 0;

% Visualize the pial and CSF surfaces
% figure;
% hold on
% patch('Faces',f1,'Vertices',v1,'facecolor',[.5 .5 .5],'edgecolor','none');
% patch('Faces',f2,'Vertices',v2,'facecolor',[0 1 1],'edgecolor','none');
% title('Pial surface outside the CSF surface', 'FontSize', 14)
% camlight('headlight','infinite')
% axis equal
% camorbit(0, 180)
% camlight('headlight')
% camorbit(0, 180)
% camorbit(0, 270)
% rotate3d on

while sum(in==0) > 0
    out_vpoint = v2(in==0,:);
    
    while_loop_num = while_loop_num + 1;
    disp(['[Loop #', num2str(while_loop_num), '] Number of pial surface points outside of csf: ', num2str(size(out_vpoint,1))])
    
    % for every one of these out points, project to the csf surface
    [ distances, surface_points ] = ...
        point2trimesh('Faces', f1, 'Vertices', v1, 'QueryPoints', out_vpoint, 'Algorithm', 'parallel'); %#ok<*ASGLU>
    % dilate around each surface_point
    new_v1 = mesh_dilate_out(v1, f1, out_vpoint, surface_points, 0.5);
    
    % visualize the csf surface with out points and compare after dilation
    if verbose
        figure;
        ax1 = subplot(1,2,1);
        hold on
        P = patch('Faces',f1,'Vertices',v1,'facecolor',[.5 .5 .5],'edgecolor','none');
        set(P, 'facealpha', 0.25)
        camlight('headlight','infinite')
        axis equal
        camorbit(0, 180)
        camlight('headlight')
        camorbit(0, 180)
        camorbit(0, 270)
        rotate3d on
        scatter3(out_vpoint(:,1), out_vpoint(:,2), out_vpoint(:,3), 20, 'r', 'filled')
        title('Pial Surface outside of headreco CSF surface', 'FontSize', 20)
        
        ax2 = subplot(1,2,2);
        hold on
        P = patch('Faces',f1,'Vertices',new_v1,'facecolor',[0 1 0],'edgecolor','none');
        set(P, 'facealpha', 0.25)
        camlight('headlight','infinite')
        axis equal
        camorbit(0, 180)
        camlight('headlight')
        camorbit(0, 180)
        camorbit(0, 270)
        rotate3d on
        scatter3(out_vpoint(:,1), out_vpoint(:,2), out_vpoint(:,3), 20, 'r', 'filled')
        title('headreco CSF surface after adjustment', 'FontSize', 20)
        
        hlink = linkprop([ax1, ax2],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'}); %#ok<*NASGU>
    end
    
    % check again if any points are still out - usually pathological geometry
    in = intriangulation(new_v1,f1,v2);
    
    if sum(in==0) > 0
        % This is usually due to intersections of dilated surfaces that create
        % ambiguity regarding which side is out/in for some points, especially
        % in the medial temporal lobe region. We apply mesh fixes to correct.
        [ f1, new_v1 ] = mesh_quality_check(f1, new_v1, false);
        in = intriangulation(new_v1,f1,v2);
    end
    
    % Update the csf surface
    v1 = new_v1;
end

%% Now repeat the same process for the cerebellum surface
% Updated on 09/01/2022: sometimes the cerebellum could protrude out the
% "sealed" csf surface from Step 1 because closing the incision plane is
% agnostic of the underlying skull surface. Therefore, in rare cases, more
% than one dilation is needed to include all the cerebellum points that are
% outside the new csf surface. Thus, we use a while loop here. 
in = intriangulation(v1,f1,v3);
while_loop_num = 0;

while sum(in==0) > 0
    out_vpoint = v3(in==0,:);
    
    while_loop_num = while_loop_num + 1;
    disp(['[Loop #', num2str(while_loop_num), '] Number of cerebellum surface points outside of csf: ', num2str(size(out_vpoint,1))])
    
    % for every one of these out points, project to the csf surface
    [ distances, surface_points ] = ...
        point2trimesh('Faces', f1, 'Vertices', v1, 'QueryPoints', out_vpoint, 'Algorithm', 'parallel');
    % dilate around each surface_point
    v1 = mesh_dilate_out(v1, f1, out_vpoint, surface_points, 0.5);
    
    % check again for inclusion. Possibly more dilations are needed
    in = intriangulation(v1,f1,v3);
    
    if sum(in==0) > 0
        % Update on 09/30/2022: add a mesh geometry check before more
        % dilations
        [ f1, v1 ] = mesh_quality_check(f1, v1, false);
        in = intriangulation(v1,f1,v3);
    end
end

%% Make sure gm surfaces are entirely enclosed in the updated csf surface
% check if the pial surface is enclosed in the updated csf surface 
in = intriangulation(v1,f1,v2);
assert(sum(in==0) == 0, 'Some pial surface points are still outside of headreco csf.')
% check if the cerebellum surface is enclosed in the updated csf surface 
in = intriangulation(v1,f1,v3);
assert(sum(in==0) == 0, 'Some cerebellum surface points are still outside of headreco csf.')

% Step 2 Done.
disp('Step 2 Done.')

%% Step 3: Smoothing and clean up 
% first of all, smooth the surface significantly
[conn,connnum,count]=meshconn(f1, length(v1));
v1=smoothsurf(v1,[],conn,150,0.8,'lowpass',0.5);

% Note on why we check gm inclusion again:
% often, due to smoothing, some gm points will get out of the csf surface
% again, this is especially a problem in young adults where the pial
% surface gets very close to the csf. Let's repeat Step 2 here. 

%% REPEAT Step 2: check out-of-bound pial and cerebellum points
% find out pial surface points that are outside of the csf surface
in = intriangulation(v1,f1,v2);
out_vpoint = v2(in==0,:);
disp(['After 1st smoothing - Number of pial surface points outside of csf: ', num2str(size(out_vpoint,1))])
% for every one of these out points, project to the csf surface 
[ distances, surface_points ] = ...
    point2trimesh('Faces', f1, 'Vertices', v1, 'QueryPoints', out_vpoint, 'Algorithm', 'parallel');
new_v1 = mesh_dilate_out(v1, f1, out_vpoint, surface_points, 0.5);

% % visualize the csf surface with out points
% figure;
% hold on
% P = patch('Faces',f1,'Vertices',v1,'facecolor',[.5 .5 .5],'edgecolor','none');
% camlight('headlight','infinite')
% axis equal
% camorbit(0, 180)
% camlight('headlight')
% camorbit(0, 180)
% camorbit(0, 270)
% rotate3d on
% scatter3(out_vpoint(:,1), out_vpoint(:,2), out_vpoint(:,3), 20, 'r', 'filled')
% title('Pial Surface outside of headreco CSF surface', 'FontSize', 20)

% Now repeat the same process for the cerebellum surface
in = intriangulation(v1,f1,v3);
out_vpoint = v3(in==0,:);
disp(['After 1st smoothing - Number of cerebellum surface points outside of csf: ', num2str(size(out_vpoint,1))])
[ distances, surface_points ] = ...
    point2trimesh('Faces', f1, 'Vertices', v1, 'QueryPoints', out_vpoint, 'Algorithm', 'parallel');
new_v1 = mesh_dilate_out(v1, f1, out_vpoint, surface_points, 1, new_v1);

% % visualize the csf surface with out points
% figure;
% hold on
% P = patch('Faces',f1,'Vertices',v1,'facecolor',[.5 .5 .5],'edgecolor','none');
% camlight('headlight','infinite')
% axis equal
% camorbit(0, 180)
% camlight('headlight')
% camorbit(0, 180)
% camorbit(0, 270)
% rotate3d on
% scatter3(out_vpoint(:,1), out_vpoint(:,2), out_vpoint(:,3), 20, 'r', 'filled')
% title('Cerebellum Surface outside of headreco CSF surface', 'FontSize', 20)

% Make sure gm surface is entirely enclosed in the updated csf surface
% check if the pial surface is enclosed in the updated csf surface 
in = intriangulation(new_v1,f1,v2);
assert(sum(in==0) == 0, 'Some pial surface points are still outside of headreco csf.')
% check if the cerebellum surface is enclosed in the updated csf surface 
in = intriangulation(new_v1,f1,v3);
assert(sum(in==0) == 0, 'Some cerebellum surface points are still outside of headreco csf.')

% Update the csf surface
v1 = new_v1; 

% Step 2 Done.
disp('Repeat of Step 2 after smoothing Done.')

%% One more round of smoothing 
% Bit more smoothing to smooth out the dilation from repeat of step 2
[conn,connnum,count]=meshconn(f1, length(v1));
v1=smoothsurf(v1,[],conn,50,0.8,'lowpass',0.5);

% visual check
figure
P = patch('Faces',f1,'Vertices',v1,'facecolor',[0.5 0.5 0.5],'edgecolor','none');
hold on
set(P, 'facealpha', 0.5)
camlight('headlight','infinite')
axis equal
camorbit(0, 180)
camlight('headlight')
title('Compatible CSF surface after processing', 'FontSize', 20)
% show the gm surfaces as well 
patch('Faces',f2,'Vertices',v2,'facecolor',[1 0 0],'edgecolor','none');
patch('Faces',f3,'Vertices',v3,'facecolor',[1 0 0],'edgecolor','none');

%% Validity checks 
% Step 3: Sanity check that this new csf surface should be outside all pial
% surface, but be inside the skull surfaces. If not, correct the pial 
% surface, as usually it's due to incorrect freesurfer segmentation of the
% pial surface. After manual quality checking, this should rarely occur. 

overproject = 0.25; % pull 0.25mm over

% check if the pial surface is enclosed in the final csf surface 
in = intriangulation(v1,f1,v2);
if sum(in==0) ~= 0 % sometimes it's just numeric inaccuracy
    disp(['After 2nd smoothing - Number of pial surface points outside of csf: ', num2str(sum(in==0))])
    scatter3(v2(in==0,1), v2(in==0,2), v2(in==0,3), 50, 'g', 'filled');
    possible_out = v2(in==0,:);
    [ distances, surface_points ] = ...
    point2trimesh('Faces', f1, 'Vertices', v1, 'QueryPoints', possible_out, 'Algorithm', 'parallel');

    % Overproject for 0.25mm to make csf thickness realistic
    pull_direction = surface_points - possible_out;
    pull_distance = vecnorm(pull_direction, 2, 2) + overproject; 
    pull_direction = pull_direction ./ vecnorm(pull_direction, 2, 2);
    new_position = possible_out + pull_distance .* pull_direction;

    % update pial surface
    v2(in==0,:) = new_position;
    disp('Pial surface fixed and saved to mri2mesh_gm_fixed_updated.mat. Make sure to use this updated pial surface when building 4-shell forward model.')
    faces = f2;
    vertices = v2;
    save(fullfile(datapath, 'mri2mesh_gm_fixed_updated.mat'), 'faces', 'vertices')
end
in = intriangulation(v1,f1,v2);
assert(sum(in==0) == 0, 'Last validity check failed: some pial surface points are still outside of new csf surface.')

% check if the cerebellum surface is enclosed in the final csf surface 
in = intriangulation(v1,f1,v3);
if sum(in==0) ~= 0 % Updated on 06/06/2022: Cerebellum could be out as well. Need to fix it! 
    disp(['After 2nd smoothing - Number of pial surface points outside of csf: ', num2str(sum(in==0))])
    scatter3(v3(in==0,1), v3(in==0,2), v3(in==0,3), 50, 'g', 'filled');
    possible_out = v3(in==0,:);
    [ distances, surface_points ] = ...
    point2trimesh('Faces', f1, 'Vertices', v1, 'QueryPoints', possible_out, 'Algorithm', 'parallel');
    
    % Overproject for 0.25mm to make csf thickness realistic
    pull_direction = surface_points - possible_out;
    pull_distance = vecnorm(pull_direction, 2, 2) + overproject; % pull 0.25mm over
    pull_direction = pull_direction ./ vecnorm(pull_direction, 2, 2);
    new_position = possible_out + pull_distance .* pull_direction;
    
    % update cerebellum surface
    v3(in==0,:) = new_position;
    disp('Cerebellum surface fixed and saved to mri2mesh_cerebellum_fixed_updated.mat. Make sure to use this updated cerebellum surface when building 4-shell forward model.')
    faces = f3;
    vertices = v3;
    save(fullfile(datapath, 'mri2mesh_cerebellum_fixed_updated.mat'), 'faces', 'vertices')
end
in = intriangulation(v1,f1,v3);
assert(sum(in==0) == 0, 'Last validity check failed: some cerebellum surface points are outside of new csf surface.')

% Step 3 Done.
disp('Step 3 Done.')

%% Save the hardwork result!
% Save data to file
faces = f1;
vertices = v1;

save(fullfile(datapath, 'headreco_csf_compatible.mat'), 'faces', 'vertices')

% Done!

disp('Total time taken:') 
toc

pause(5)
close all

end
