function [ H ] = plot_loading(channel, forward_model, skin_surf, max_normal, all_sensor, topoplt)
%% Generates a loading plot

if ~exist('skin_surf', 'var')
    skin_surf = true;
end
if ~exist('max_normal', 'var')
    max_normal = true;
end
if ~exist('all_sensor', 'var')
    all_sensor = true;
end
if nargin < 6
    topoplt = false;
end

if ~topoplt
    % grab weights from the channel electrode
    weight = forward_model.G(channel,1:size(forward_model.source,1)).^2;
    
    % show the normal vector of the max:
    [~, max_idx] = max(weight);
    
    normal_v = forward_model.normal(max_idx, :)/20;
    
    p1 = forward_model.source(max_idx, :);
    p2 = p1 + normal_v;
    
    % plot the sources in 3D spaceh
    H = figure;
    
    size_weight = ceil((weight-min(weight))/(max(weight)-min(weight))*300)+1;
    scatter3(forward_model.source(:,1), forward_model.source(:,2), forward_model.source(:,3),...
        size_weight, weight, 'filled')
    axis equal
    rotate3d on
    
    hold on
    
    if skin_surf
        % Visualize the outer skin surface
        P = patch('Faces',forward_model.skin_surf_face,'Vertices',forward_model.skin_surf_vertex,...
            'facecolor',[.5 .5 .5],'edgecolor','none');
        set(P, 'facealpha', 0.5)
    end
    
    if max_normal
        % Show normal vector
        mArrow3(p1,p2,'color','red', 'tipWidth', 0.002, 'facealpha',0.5);
    end
    
    if all_sensor
        % Plot all channel locations
        scatter3(forward_model.dig(:,1), forward_model.dig(:,2), forward_model.dig(:,3), 300, [1,0,0], 's', 'filled')
    end
    % Plot the loading channel location
    scatter3(forward_model.dig(channel,1), forward_model.dig(channel,2), forward_model.dig(channel,3), 300, [1,0.5,1], 's', 'filled')
    
    % Format the figure
    title(['Loading Plot for Channel ' num2str(channel)], 'FontSize', 16)
    colorbar

else
    figure
    if length(forward_model.eloc) > 200
        topoplot([], forward_model.eloc, 'style','both','electrodes','on','emarker', {'.', 'k', 15, 1});
    else
        topoplot([], forward_model.eloc, 'style','both','electrodes','ptsnumbers','emarker', {'.', 'k', 15, 1});
        
        L = findobj(gcf, 'type', 'Text');
        for ind = 1:length(forward_model.eloc)
            set(L(length(forward_model.eloc)+1-ind), 'FontSize', 14)
        end
    end
    
    set(gca, 'FontSize', 20)
end

