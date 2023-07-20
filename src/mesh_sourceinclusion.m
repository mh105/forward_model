function [ src_pts, altered_idx ] = mesh_sourceinclusion(src_pts, pial_f, pial_v, min_distance)
%This function is used to check whether all source points are included
%within the pial surface. If not, we project these points inside the pial
%surface.
%
% In addition, we check for the minimum distances from source points to
% the pial surface to be at least min_distance(mm). If not, we shift the
% source points, that are too close, away from the pial surface to ensure
% the source space respects a minimal cortex thickness condition.

if nargin < 4
    min_distance = 0.5;
end

%% Use a while loop to repeat checking until both criteria are met
% Added on 07/05/2022
end_check = false;
altered_idx = [];
check_loop_num = 0;

while ~end_check
    %% Check if all source points are within the pial surface
    
    % % % quick visual check
    % figure
    % hold on
    % scatter3(src_pts(:,1), src_pts(:,2), src_pts(:,3), 5, 'y', 'filled')
    % P = patch('Faces',pial_f,'Vertices',pial_v,'facecolor',[.5 .5 .5],'edgecolor','none');
    % set(P, 'facealpha', 0.25)
    % camlight('headlight','infinite')
    % axis equal
    % rotate3d on
    
    in = intriangulation(pial_v,pial_f,src_pts);
    
    % highlight the out points
    outp_index = find(in==0);
    altered_idx = [altered_idx; outp_index]; %#ok<*AGROW>
    
    % % visualize these points
    % scatter3(sources(outp_index,1), sources(outp_index,2), sources(outp_index,3), 50, 'r', 'filled')
    
    if ~isempty(outp_index)
        % let's find out the shortest distance by projection of these out points to
        % the pial surface
        [ ~, surface_points ] = ...
            point2trimesh('Faces', pial_f, 'Vertices', pial_v, 'QueryPoints', src_pts(outp_index,:), 'Algorithm', 'parallel');
        
        % We move them in the direction of the projection and overproject by 0.1mm
        % compute the projection vector
        projection_vector = surface_points - src_pts(outp_index,:);
        projection_length = vecnorm(projection_vector, 2, 2);
        projection_scalef = (projection_length+0.1) ./ projection_length;
        scaled_projection = projection_vector .* projection_scalef;
        assert(all(vecnorm(scaled_projection, 2, 2) - projection_length - 0.1*ones(size(projection_length)) < 10^-6), 'Overprojection is not uniformly 0.1mm, please check!')
        % now actually project
        new_points = src_pts(outp_index,:) + scaled_projection;
        % update the source points
        src_pts(outp_index,:) = new_points;
    end
    
    % let's check again for inclusion of all source points
    in = intriangulation(pial_v,pial_f,src_pts);
    assert(sum(in) == size(src_pts,1), 'After projection, still not all source points are included in the pial surface!')
    
    %% Check if cortex thickness is at least [min_distance]mm for all source points
    % compute the distance to the pial surface for all source points
    [ all_distances, all_surface_points ] = ...
        point2trimesh('Faces', pial_f, 'Vertices', pial_v, 'QueryPoints', src_pts, 'Algorithm', 'parallel');
    
    % find the points with thin cortex
    thinp_index = find(abs(all_distances) < min_distance);
    altered_idx = [altered_idx; thinp_index];
    
    % % visualize these points
    % scatter3(src_pts(thinp_index,1), src_pts(thinp_index,2), src_pts(thinp_index,3), 100, 'g', 'filled')
    
    if ~isempty(thinp_index)
        % move them away from the pial surface towards 0.51mm
        projection_vector = src_pts(thinp_index,:) - all_surface_points(thinp_index,:);
        projection_length = vecnorm(projection_vector, 2, 2);
        projection_scalef = (min_distance+0.01) ./ projection_length;
        scaled_projection = projection_vector .* projection_scalef;
        assert(all(vecnorm(scaled_projection, 2, 2) - (min_distance+0.01)*ones(size(projection_length)) < 10^-6),...
            ['Shifting is not uniformly up to ' num2str(min_distance+0.01) 'mm, please check!'])
        % now actually shift
        new_points = all_surface_points(thinp_index,:) + scaled_projection;
        % update the source_points
        src_pts(thinp_index,:) = new_points;
        
        % % visualize the shifted points
        % scatter3(sources(thinp_index,1), sources(thinp_index,2), sources(thinp_index,3), 100, 'b', 'filled')
        
        % now we enter an iterative algorithm that keeps pushing these points until
        % they are at least 0.5mm away from the the pial surface
        fullsrcindex = 1:size(src_pts,1);
        shift_list = fullsrcindex(thinp_index);
        
        % check the shifted distances
        [ shift_distances, shift_surface_points ] = ...
            point2trimesh('Faces', pial_f, 'Vertices', pial_v, 'QueryPoints', src_pts(thinp_index,:), 'Algorithm', 'parallel');
        
        still_thin_index = abs(shift_distances) < min_distance;
        
        while any(still_thin_index==1)
            shift_list = shift_list(still_thin_index);
            
            % move them away from the pial surface towards 0.51mm
            projection_vector = src_pts(shift_list,:) - shift_surface_points(still_thin_index,:);
            projection_length = vecnorm(projection_vector, 2, 2);
            projection_scalef = (min_distance+0.01) ./ projection_length;
            scaled_projection = projection_vector .* projection_scalef;
            assert(all(vecnorm(scaled_projection, 2, 2) - (min_distance+0.01)*ones(size(projection_length)) < 10^-6),...
                ['Shifting is not uniformly up to ' num2str(min_distance+0.01) 'mm, please check!'])
            % now actually shift
            new_points = shift_surface_points(still_thin_index,:) + scaled_projection;
            % update the source_points
            src_pts(shift_list,:) = new_points;
            
            % check the shifted distances
            [ shift_distances, shift_surface_points ] = ...
                point2trimesh('Faces', pial_f, 'Vertices', pial_v, 'QueryPoints', src_pts(shift_list,:), 'Algorithm', 'parallel');
            
            % update the still_thin_index
            still_thin_index = abs(shift_distances) < min_distance;
        end
    end
    
    %% Check again whether criteria are met
    % let's check again for inclusion of all source points.
    in = intriangulation(pial_v,pial_f,src_pts);
    
    % also check again for minimal distance
    all_distances = point2trimesh('Faces', pial_f, 'Vertices', pial_v, 'QueryPoints', src_pts, 'Algorithm', 'parallel');
    thinp_index = find(abs(all_distances) < min_distance, 1);
    
    % ready to exit the while loop
    check_loop_num = check_loop_num + 1;
    if sum(in) == size(src_pts,1) && isempty(thinp_index)
        end_check = true;
    end
    
end

%% Report the total number of source points altered in this processing
altered_idx = unique(altered_idx);
disp(['Total number of source points altered: ' num2str(length(altered_idx)) ' points'])
if check_loop_num == 1
    disp(['Total checking loops used: ' num2str(check_loop_num) ' time'])
else
    warning(['Total checking loops used: ' num2str(check_loop_num) ' times'])
end

end

