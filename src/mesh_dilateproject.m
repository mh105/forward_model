function [ new_V ]  = mesh_dilateproject(f3_above, sub_v3, F, V_rot, dilate_index, project_index)
% Used to dilate & project MNE_outer_skull surface onto the headreco skull
% surface in the upper 3 quadrants. 

sigma = 4;

% 1) Project all the dilate points to the headreco surface 
[ distances, surface_points ] = ...
    point2trimesh('Faces', f3_above, 'Vertices', sub_v3, 'QueryPoints', V_rot(dilate_index,:), 'Algorithm', 'parallel');

% 2) Local dilation for each dilate points 
new_V = V_rot;
% dilate local surface for each point, only keep the max dilation 
for ii = 1:length(dilate_index)
    % extract 5-step neighbors 
    [vertex_index, face_index] = mesh_findneighbor(F, dilate_index(ii), 6);
    % propose to dilate
    dilate_direction = surface_points(ii,:) - V_rot(dilate_index(ii),:);
    dilate_distance = norm(dilate_direction);
    dilate_direction = dilate_direction / norm(dilate_direction); % make into unit vector 
    % calculate the distance from each vertex to the seed
    local_map_distance = vecnorm(V_rot(vertex_index,:) - V_rot(dilate_index(ii),:), 2, 2);
    % calculate dilation multiplier using Gaussian filter
    dilate_multi = normpdf(local_map_distance,0,sigma) / normpdf(0,0,sigma);
    % calculate vertex positions after dilation
    V_after_dilation = V_rot(vertex_index,:) + dilate_distance*dilate_multi*dilate_direction;
    % compare with the current distance from original v1 vertex set
    new_distance = vecnorm(V_after_dilation - V_rot(vertex_index,:), 2, 2);
    assert(all(dilate_distance*dilate_multi - new_distance < 10^-10), 'Something is wrong with the dilation.')
    current_distance = vecnorm(new_V(vertex_index,:) - V_rot(vertex_index,:), 2, 2);
    % find out dilations that pushed the points further than already been 
    accepted_dilation = new_distance > current_distance;
    % update the accepted dilations 
    new_V(vertex_index(accepted_dilation), :) = V_after_dilation(accepted_dilation, :);
end

% 3) Project all project points onto the headreco surface 
[ distances, surface_points ] = ...
    point2trimesh('Faces', f3_above, 'Vertices', sub_v3, 'QueryPoints', new_V(project_index,:), 'Algorithm', 'parallel');
% update the positions of these project points 
new_V(project_index,:) = surface_points;

% % 4) sanity visual check
% figure;
% hold on
% P = patch('Faces',f3_above,'Vertices',sub_v3,'facecolor',[0 .5 .5],'edgecolor','none');
% set(P, 'facealpha', 0.5)
% P = patch('Faces',F,'Vertices',new_V,'facecolor',[.5 .5 .5],'edgecolor','none');
% set(P, 'facealpha', 0.5)
% camlight('headlight','infinite')
% axis equal
% camorbit(0, 180)
% camlight('headlight')
% camorbit(0, 180)
% camorbit(0, 270)
% rotate3d on
% view(180, 0)
% % scatter3(V_rot(dilate_index,1), V_rot(dilate_index,2), V_rot(dilate_index,3), 5, 'y', 'filled')

end

