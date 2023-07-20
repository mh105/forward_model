function [ newbem_fn ] = MNE_editbem(oldbem_fn, filepath, offset)
%
% MNE CODES - EDITBEM
%
% - this function is used to correct the BEM surface to not exclude source
% points. 
% - we do so by 1) creating a new vertex perpendicularly beneath the out
% points and restore local Delaunay conditions. We use this new surface to
% define 1st order and 2nd order patches on the original surface, 2) we
% radially dilate these 1st and 2nd order vertices to cover the out points.
%
% Last edit: Alex He 05/20/2022
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Inputs:
%           - oldbem_fn:    filename of the .mat file from calling function
%                           in ANT_MNE_python_util.py to export the
%                           freesurfer inner skull surface in triangulation
%                           format as a .mat file. 
%
%           - filepath:     full path to the folder containing
%                           oldbem_fn.mat, which is a .mat file containing
%                           the surface triangulation of the inner skull
%                           surface exported by ANT_MNE_python_util.py from
%                           a freesurfer format outputted from calling
%                           mne.bem.make_flash_bem() in python or calling
%                           mne flash_bem in shell script.
%
%           - offset:       the offset added to dilation of inner skull to
%                           encompass the source points that are initially
%                           out. In units of mm, usually 0.2-0.5mm are
%                           sufficient to avoid exclusion by MNE-Python
%                           call of mne.make_forward_solution().
%
% Output:
%           - newbem_fn:    newbem_fn = [oldbem_fn '_edited.mat'];
%                           we save the new triangulation file as an
%                           newbem_fn.mat file to be read back into
%                           python functions to create the freesurfer BEM
%                           .surf files within the freesurfer '/bem/flash'
%                           directory. 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if nargin < 3
    offset = 0.5;  % in mm
end

% Addpath to the CurvatureEstimation toolbox using a relative path
addpath('../external_toolbox/CurvatureEstimation')

%% Load the triangulation of inner skull BEM surface
load(fullfile(filepath, oldbem_fn)) %#ok<LOAD>
close all

%% Visualize BEM inner skull surface and out points

figure;
% visualize the BEM inner skull surface
% trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3)); axis equal
patch('Faces',faces,'Vertices',vertices,'facecolor',[.5 .5 .5],'edgecolor','none') %#ok<*NODEF>
camlight('headlight','infinite')
% camorbit(180, 0)
camorbit(0, 180)
camlight('headlight')
camorbit(0, 180)
camorbit(0, 270)
rotate3d on
axis equal

hold on

% visualize the source locations as scatter points
lidx = logical(lidx);
left_x = lsrc(lidx,1);
left_y = lsrc(lidx,2);
left_z = lsrc(lidx,3);
scatter3(left_x,left_y,left_z, 'filled', 'yellow')

ridx = logical(ridx);
right_x = rsrc(ridx,1);
right_y = rsrc(ridx,2);
right_z = rsrc(ridx,3);
scatter3(right_x,right_y,right_z, 'filled', 'yellow')

% use intriangulation.m to check whether source points are outside BEM
% inner skull surface
src_points = [lsrc(lidx,:); rsrc(ridx,:)];
in = intriangulation(vertices,faces,src_points);

out_points = src_points(in==0,:);
left_out = out_points(find(in==0)<length(lsrc(lidx,:)), :);
right_out = out_points(find(in==0)>length(lsrc(lidx,:)), :);

% visualize the out points
scatter3(left_out(:,1),left_out(:,2),left_out(:,3), 'filled', 'red')
scatter3(right_out(:,1),right_out(:,2),right_out(:,3), 'filled', 'red')

title('Raw inner skull surface with out points')

%% Add vertices to the BEM triangulation surface based on these out points
% We will loop through all the out points individually rather than adding
% new surface vertices for all points at the same time.

meta_fo_faces = [];
meta_so_faces = [];
meta_to_faces = [];

% store a copy of the original faces and vertices from .surf file
raw_faces = faces;
raw_vertices = vertices;

for k = 1:size(out_points, 1)
    %%
    % Step 1: create new vertices directly below the out points using the
    % point2trimesh.m function.
    
    % Define the current out point
    testp = out_points(k,:);
    
    % Find surface point and add new vertices and faces to include the
    % surface point as a vertex
    [ ~, ~, faces2, vertices2, corresponding_vertices_ID, ~ ] = point2trimesh('Faces', faces, 'Vertices', vertices, 'QueryPoints', testp);
    
    %%
    % Step 2: identify the neighbor patches 
    
    % For the new surface_point, which is now a new vertex, grab 1st order
    % neighbor vertices and faces
    fo_faces = [];
    for i = 1:length(faces2)
        if ismember(corresponding_vertices_ID, faces2(i,:))
            fo_faces = [fo_faces, i]; %#ok<*AGROW>
        end
    end
    fo_vertices = unique(faces2(fo_faces,:));
    
    % Let's also grab 2nd order neighbor vertices and faces
    so_faces = [];
    for j = 1:length(fo_vertices)
        for i = 1:length(faces2)
            if ismember(fo_vertices(j), faces2(i,:))
                so_faces = [so_faces, i];
            end
        end
    end
    so_vertices = unique(faces2(so_faces,:));
    
    % Let's further grab 3rd order neighbor vertices and faces 
    to_faces = [];
    for j = 1:length(so_vertices)
        for i = 1:length(faces2)
            if ismember(so_vertices(j), faces2(i,:))
                to_faces = [to_faces, i];
            end
        end
    end
    % to_vertices = unique(faces2(to_faces,:));
    
    % store the first and second order faces indices 
    meta_fo_faces = [meta_fo_faces; faces2(fo_faces,:)];
    meta_so_faces = [meta_so_faces; faces2(so_faces,:)];
    meta_to_faces = [meta_to_faces; faces2(to_faces,:)];
    
    %%
    % Step 3: update the vertices and faces 
    vertices = vertices2;
    faces = faces2;
    
end

% %% Visualize the new BEM surface and source points
% 
% figure;
% % visualize the BEM inner skull surface
% % trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3)); axis equal
% patch('Faces',faces,'Vertices',vertices,'facecolor',[.5 .5 .5],'edgecolor','none')
% camlight('headlight','infinite')
% % camorbit(180, 0)
% camorbit(0, 180)
% camlight('headlight')
% camorbit(0, 180)
% camorbit(0, 270)
% rotate3d on
% axis equal
% 
% hold on
% 
% % visualize the meta second-order faces as a different color 
% patch('Faces', meta_so_faces, 'Vertices',vertices,'facecolor',[0.4660 0.6740 0.1880],'edgecolor','none')
% 
% % visualize the source locations as scatter points
% lidx = logical(lidx);
% left_x = lsrc(lidx,1);
% left_y = lsrc(lidx,2);
% left_z = lsrc(lidx,3);
% scatter3(left_x,left_y,left_z, 'filled', 'yellow')
% 
% ridx = logical(ridx);
% right_x = rsrc(ridx,1);
% right_y = rsrc(ridx,2);
% right_z = rsrc(ridx,3);
% scatter3(right_x,right_y,right_z, 'filled', 'yellow')
% 
% % visualize the out points that have been previously identified 
% scatter3(left_out(:,1),left_out(:,2),left_out(:,3), 'filled', 'red')
% scatter3(right_out(:,1),right_out(:,2),right_out(:,3), 'filled', 'red')
% 
% title('Inner skull surface with added vertices')

%% Dilation of original patches defined using 1st and 2nd order neighbors
% Rather than adding more vertices and then smoothing, now that we have the
% patch of second neighbors, let's simply capitalize on this information
% and dilate the original patches locally! This is much much better than
% smoothing, which risk destroying the local geometry?

%% Step 1: Define 1st, 2nd, and 3rd order local patches 
% (this is actually 2nd, 3rd, and 4th order, because we are using patches to
% define vertex indices)

% Here are all the vertices involved in first-order patches
fo_patch_vertices = unique(meta_fo_faces);

% Let's find the ones that were originally there in raw vertices copy
raw_fo_patch_vertices = fo_patch_vertices(fo_patch_vertices < length(raw_vertices)+1);

% Let's find the faces that contain these vertices in the original surface
raw_fo_faces = [];
for j = 1:length(raw_fo_patch_vertices)
    for i = 1:length(raw_faces)
        if ismember(raw_fo_patch_vertices(j), raw_faces(i,:))
            raw_fo_faces = [raw_fo_faces, i];
        end
    end
end

% Here are all the vertices involved in second-order patches
so_patch_vertices = unique(meta_so_faces);

% Let's find the ones that were originally there in raw vertices copy
raw_so_patch_vertices = so_patch_vertices(so_patch_vertices < length(raw_vertices)+1);

% Let's find the faces that contain these vertices in the original surface
raw_so_faces = [];
for j = 1:length(raw_so_patch_vertices)
    for i = 1:length(raw_faces)
        if ismember(raw_so_patch_vertices(j), raw_faces(i,:))
            raw_so_faces = [raw_so_faces, i];
        end
    end
end

% Here are all the vertices involved in third-order patches
to_patch_vertices = unique(meta_to_faces);

% Let's find the ones that were originally there in raw vertices copy
raw_to_patch_vertices = to_patch_vertices(to_patch_vertices < length(raw_vertices)+1);

% Let's find the faces that contain these vertices in the original surface
raw_to_faces = [];
for j = 1:length(raw_to_patch_vertices)
    for i = 1:length(raw_faces)
        if ismember(raw_to_patch_vertices(j), raw_faces(i,:))
            raw_to_faces = [raw_to_faces, i];
        end
    end
end

%% Step 2: Visualize local patches on the BEM inner skull surface

% Visualize the original patches!
figure;
% visualize the BEM inner skull surface
% trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3)); axis equal
patch('Faces',raw_faces,'Vertices',raw_vertices,'facecolor',[.5 .5 .5],'edgecolor','none')
camlight('headlight','infinite')
% camorbit(180, 0)
camorbit(0, 180)
camlight('headlight')
camorbit(0, 180)
camorbit(0, 270)
rotate3d on
axis equal

hold on
% visualize the 3rd, 2nd and 1st order faces in different colors 
patch('Faces', raw_faces(raw_to_faces, :), 'Vertices',raw_vertices,'facecolor',[0.8500, 0.3250, 0.0980],'edgecolor','none')
patch('Faces', raw_faces(raw_so_faces, :), 'Vertices',raw_vertices,'facecolor',[0.4660 0.6740 0.1880],'edgecolor','none')
patch('Faces', raw_faces(raw_fo_faces, :), 'Vertices',raw_vertices,'facecolor',[0.3010 0.7450 0.9330],'edgecolor','none')

% visualize the source locations as scatter points
lidx = logical(lidx);
left_x = lsrc(lidx,1);
left_y = lsrc(lidx,2);
left_z = lsrc(lidx,3);
scatter3(left_x,left_y,left_z, 'filled', 'yellow')

ridx = logical(ridx);
right_x = rsrc(ridx,1);
right_y = rsrc(ridx,2);
right_z = rsrc(ridx,3);
scatter3(right_x,right_y,right_z, 'filled', 'yellow')

% visualize original out points
scatter3(out_points(:,1),out_points(:,2),out_points(:,3), 'filled', 'red')

title('Raw 1st and 2nd order patches before dilation with out points')

%% Step 3: Dilate the local patch 

% Identify the indices of vertices that need to be dilated
dilate_faces = [raw_faces(raw_so_faces, :); raw_faces(raw_fo_faces, :); raw_faces(raw_to_faces, :)];
dilate_vertices = unique(dilate_faces);

% restore the vertices and faces
vertices = raw_vertices;
faces = raw_faces;

% Find out the maximum dilation distance
[ distances, ~] = point2trimesh('Faces', faces, 'Vertices', vertices, 'QueryPoints', out_points);

max_dilate = max(distances) + offset; % add an offset baseline to make dilation smoother

% Find out the normal vectors at each of the vertices so we can dilate in
% the normal direction
FV = struct;
FV.faces = faces;
FV.vertices = vertices;
vertexNormals = CalcVertexNormals(FV, CalcFaceNormals(FV));

% Now we will dilate the involved vertices depending on whether they are
% 1st order, 2nd order, or 3rd order
for i = 1:length(dilate_vertices)
    current_vertex_idx = dilate_vertices(i);
    current_normal = vertexNormals(current_vertex_idx, :);
    
    if ismember(current_vertex_idx, fo_patch_vertices)
        update_vertex_coordinate = vertices(current_vertex_idx,:) + max_dilate .* current_normal;
    elseif ismember(current_vertex_idx, so_patch_vertices)
        update_vertex_coordinate = vertices(current_vertex_idx,:) + max_dilate .* current_normal;
    elseif ismember(current_vertex_idx, to_patch_vertices)
        update_vertex_coordinate = vertices(current_vertex_idx,:) + 1/2*max_dilate .* current_normal;
    else
        update_vertex_coordinate = vertices(current_vertex_idx,:) + 1/4*max_dilate .* current_normal;
    end
    vertices(current_vertex_idx,:) = update_vertex_coordinate;
end

%% Step 4: Visualize the dilated patches on BEM surface

% Visualize the original patches!
figure;
% visualize the BEM inner skull surface
% trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3)); axis equal
patch('Faces',faces,'Vertices',vertices,'facecolor',[.5 .5 .5],'edgecolor','none')
camlight('headlight','infinite')
% camorbit(180, 0)
camorbit(0, 180)
camlight('headlight')
camorbit(0, 180)
camorbit(0, 270)
rotate3d on
axis equal

hold on
% visualize the 3rd, 2nd and 1st order faces in different colors 
patch('Faces', raw_faces(raw_to_faces, :), 'Vertices',vertices,'facecolor',[0.8500, 0.3250, 0.0980],'edgecolor','none')
patch('Faces', raw_faces(raw_so_faces, :), 'Vertices',vertices,'facecolor',[0.4660 0.6740 0.1880],'edgecolor','none')
patch('Faces', raw_faces(raw_fo_faces, :), 'Vertices',vertices,'facecolor',[0.3010 0.7450 0.9330],'edgecolor','none')

% visualize the source locations as scatter points
lidx = logical(lidx);
left_x = lsrc(lidx,1);
left_y = lsrc(lidx,2);
left_z = lsrc(lidx,3);
scatter3(left_x,left_y,left_z, 'filled', 'yellow')

ridx = logical(ridx);
right_x = rsrc(ridx,1);
right_y = rsrc(ridx,2);
right_z = rsrc(ridx,3);
scatter3(right_x,right_y,right_z, 'filled', 'yellow')

% use intriangulation.m again to check whether source points are outside BEM
% inner skull surface. There shouldn't be any but let's visualize and check
src_points = [lsrc(lidx,:); rsrc(ridx,:)];
in = intriangulation(vertices,faces,src_points);

out_points = src_points(in==0,:);
left_out = out_points(find(in==0)<length(lsrc(lidx,:)), :);
right_out = out_points(find(in==0)>length(lsrc(lidx,:)), :);

% visualize any remaining out points, which shouldn't be there
scatter3(left_out(:,1),left_out(:,2),left_out(:,3), 'filled', 'red')
scatter3(right_out(:,1),right_out(:,2),right_out(:,3), 'filled', 'red')

title({'1st, 2nd, 3rd order patches after dilation', '(should have no out point)'})

% Make sure all source points are now inside the inner skull surface
assert(isempty(out_points), 'There are still source points outside the inner skull surface after dilation.')

%%
% Save the edited surface as a .mat file and generate a freesurfer .surf
% file by calling functions in ANT_MNE_python_util.py
a = strsplit(oldbem_fn, '.mat');
oldbem_fn = a{1};
newbem_fn = [oldbem_fn '_edited.mat'];
save(fullfile(filepath, newbem_fn), 'faces', 'vertices')

close all
disp('MATLAB MNE_editbem() completed without error.')

end

