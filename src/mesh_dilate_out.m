function [ new_v1 ] = mesh_dilate_out(v1, f1, out_vpoint, surface_points, dilate_margin, new_v1, sigma, neighbor_n)
% Used to dilate a surface to incorporate points from an inner surface that
% protrudes out. For instance, to dilate the csf surface to accommodate all
% the gm surface points.

if nargin < 5
    dilate_margin = 0.5; % mm
end
if nargin < 6 || isempty(new_v1)
    new_v1 = v1;
end
if nargin < 7
    sigma = 3;
end
if nargin < 8
    neighbor_n = 5;
end

% dilate local surface for each point, only keep the max dilation
for ii = 1:size(out_vpoint, 1)
    % find the closest vertex on csf surface
    [~, min_seed] = min(vecnorm(v1 - surface_points(ii,:), 2, 2));
    % extract 5-step neighbors
    [vertex_index, face_index] = mesh_findneighbor(f1, min_seed, neighbor_n);
    % propose to dilate
    dilate_direction = out_vpoint(ii,:) - surface_points(ii,:);
    dilate_distance = norm(dilate_direction) + dilate_margin; % over dilate by a certain margin in mm
    dilate_direction = dilate_direction / norm(dilate_direction); % make into unit vector
    % calculate the distance from each vertex to the seed
    local_map_distance = vecnorm(v1(vertex_index,:) - v1(min_seed,:), 2, 2);
    % calculate dilation multiplier using Gaussian filter
    dilate_multi = normpdf(local_map_distance,0,sigma) / normpdf(0,0,sigma);
    % calculate vertex positions after dilation
    v1_after_dilation = v1(vertex_index,:) + dilate_distance*dilate_multi*dilate_direction;
    % compare with the current distance from original v1 vertex set
    new_distance = vecnorm(v1_after_dilation - v1(vertex_index,:), 2, 2);
    assert(all(dilate_distance*dilate_multi - new_distance < 10^-10), 'Something is wrong with the dilation.')
    current_distance = vecnorm(new_v1(vertex_index,:) - v1(vertex_index,:), 2, 2);
    % find out dilations that pushed the points further than already been
    accepted_dilation = new_distance > current_distance;
    % update the accepted dilations
    new_v1(vertex_index(accepted_dilation), :) = v1_after_dilation(accepted_dilation, :);
    
    %     % these plotting are for debugging only
    %     x1=scatter3(v1(min_seed,1), v1(min_seed,2), v1(min_seed,3), 20, 'b', 'filled');
    %     x2= patch('Faces',f1(face_index,:),'Vertices',v1,'facecolor',[1 0 0],'edgecolor','none');
    %     set(x2, 'facealpha', 0.25)
    %     x3= scatter3(out_vpoint(ii,1), out_vpoint(ii,2), out_vpoint(ii,3), 20, 'm', 'filled');
    %     x4= scatter3(new_v1(vertex_index,1), new_v1(vertex_index,2), new_v1(vertex_index,3), 20, 'k', 'filled');
    %     x5= scatter3(v1_after_dilation(~accepted_dilation,1), v1_after_dilation(~accepted_dilation,2), v1_after_dilation(~accepted_dilation,3), 20, 'g', 'filled');
    %
    %     X = [x1,x2,x3,x4,x5];
    %     delete(X)
end

end

