function [ intersected, edited_V1, edited_V2 ] = mesh_surfaceintersect(F1, V1, F2, V2, max_array_size, skip_remesh, edit_V)
%This function is used to check whether two surfaces intersect with each
%other. It utilizes the SurfaceIntersection() function but has to implement
%subdivision of one surface first to preserve memory space. If both
%surfaces are at high resolution, this function could scale up to a really
%long time. For instance, checking the full-resolution pial surface and csf
%mesh has unrealistic runtime.

if nargin < 5
    max_array_size = 10^8;
    skip_remesh = false;
    edit_V = false;
end

if (size(F1,1) < size(F2,1)) || skip_remesh
    % we need to subdivide F2 and modify F1 if needed
    small_f = F1;
    small_v = V1;
    remeshed_surf = '2';
    large_f = F2;
    large_v = V2;
else
    % we need to subdivide F1 and modify F2 if needed
    small_f = F2;
    small_v = V2;
    remeshed_surf = '1';
    large_f = F1;
    large_v = V1;
end

intersected = struct;
intersected.intersected = false;

surface1 = struct;
surface1.faces = small_f; surface1.vertices = small_v;
intersected.intx_small_surf = surface1;

intersected_f_idx = {};
intersected_v_idx = {};
intersected_intx_surf = {};

%% Optional remeshing before checking for intersections between surfaces

if skip_remesh
    Fn = large_f; %#ok<*NASGU>
    Vn = large_v;
    seedIndex = 1:size(large_v, 1);
    indSeeds = 1:size(large_v, 1);
    indSeeds = indSeeds';
    
else
    % % quick visual check before remesh
    % faces_list = {large_f};
    % vertices_list = {large_v};
    % mesh_plot_bem_surfaces(faces_list, vertices_list);
    
    % configure the number of subdivisions
    patch_size = ceil(max_array_size / size(small_v,1)); % number of faces allowed for each patch
    new_num_vertice = ceil(size(large_v,1) / (size(large_f,1) / patch_size));
    
    % remesh the larger surface using this many vertices
    indStart=1; %Index of the start point
    numSeeds=new_num_vertice; % number of vertices in the new mesh
    optionStruct.toleranceLevel=0; %Tolerance for convergence
    optionStruct.waitBarOn=0; %Turn on/off waitbar
    %distance based marching
    [Fn,Vn,seedIndex,indSeeds,~]=remeshTriSurfDistMap(large_f,large_v,numSeeds,indStart,optionStruct); %#ok<ASGLU>
    
    % % quick visual check after remesh
    % faces_list = {Fn};
    % vertices_list = {Vn};
    % mesh_plot_bem_surfaces(faces_list, vertices_list);
    
end

%%
% we then use this remeshed resolution to define patches of the large
% surface for calculating surface intersection with the small surface

for ii = indSeeds' % loop through the patches
    patch_v_idx = find(seedIndex == ii);
    select_f_idx = [];
    for j = 1:size(large_f,1)
        if any(ismember(large_f(j,:), patch_v_idx))
            select_f_idx = [select_f_idx, j];
        end
    end
    new_faces = large_f(select_f_idx,:);
    new_vertices_list = unique(new_faces);
    % we need to update the index
    for j = 1:size(new_faces,1)
        [sel, new_face_idx] = ismember(new_faces(j,:), new_vertices_list);
        assert(all(sel), 'Some of the indices are not in the new_vertices_list!')
        new_faces(j,:) = new_face_idx;
    end
    % construct new vertice list as well
    new_vertices = large_v(new_vertices_list, :);
    
    %     % quick visual check
    %     faces_list = {large_f, new_faces};
    %     vertices_list = {large_v, new_vertices};
    %     mesh_plot_bem_surfaces(faces_list, vertices_list);
    %     scatter3(new_vertices(:,1), new_vertices(:,2), new_vertices(:,3), 20, 'g', 'filled')
    
    % now, we test for surface intersection of this patch with the
    % other surface (surface1) that was kept intact
    surface2 = struct;
    surface2.faces = new_faces; surface2.vertices = new_vertices;
    
    [intersect12, ~] = SurfaceIntersection(surface1, surface2);
    
    if any(intersect12, 'all')
        if edit_V
            warning('Surface intersection was found! edit_V is set to True. Will try to correct the intersection.')
        else
            warning('Surface intersection was found! Intersected surfaces are saved, please investigate!')
        end
        intersected_f_idx{length(intersected_f_idx)+1} = select_f_idx; %#ok<*AGROW>
        intersected_v_idx{length(intersected_v_idx)+1} = new_vertices_list;
        intersected_intx_surf{length(intersected_intx_surf)+1} = surface2;
        intersected.intersected = true;
        
        %         % visualize the intersected surfaces
        %         faces_list = {surface1.faces, surface2.faces};
        %         vertices_list = {surface1.vertices, surface2.vertices};
        %         mesh_plot_bem_surfaces(faces_list, vertices_list);
        %         scatter3(large_v(patch_v_idx,1), large_v(patch_v_idx,2), large_v(patch_v_idx,3), 50, 'g', 'filled')
        
        % try to correct the surface intersection
        if edit_V 
            % try to dilate surface1 a tiny bit in the opposite direction
            center_point = large_v(patch_v_idx, :);  % center point of the patch
            [ ~, surface_points ] = ...
                point2trimesh('Faces', surface1.faces, 'Vertices', surface1.vertices, 'QueryPoints', center_point, 'Algorithm', 'parallel');
            out_vpoint = 2*surface_points - center_point;  % on the opposite side of surface1
            surface1.vertices = mesh_dilate_out(surface1.vertices, surface1.faces, out_vpoint, surface_points, 0.25);  % 0.25mm margin
            
            % verify the patch no longer intersects with the surface
            intersect12 = SurfaceIntersection(surface1, surface2);
            assert(~any(intersect12, 'all'), 'Patch and surface still intersect after dilation.')
        end
    end
end

% update the output structure
intersected.intx_f_dix = intersected_f_idx;
intersected.intx_v_dix = intersected_v_idx;
intersected.intx_surf = intersected_intx_surf;

% return vertices, which are not modified unless edit_V was set to True and
% intersections between surfaces were found
if remeshed_surf == '2'  % then V1 was small_v and maybe modified
    edited_V1 = surface1.vertices;
    edited_V2 = V2;  % V2 is not modified
else  % then V2 was small_v and maybe modified
    edited_V1 = V1; % V1 is not modified
    edited_V2 = surface1.vertices;
end

end
