function [ H ] = mesh_plot_bem_surfaces(faces_list, vertices_list, sources, elc)
%Used to plot BEM surfaces for forward modeling

if nargin < 3
    sources = [];
    elc = [];
elseif nargin < 4
    elc = [];
end

%%
H = figure;
hold on

% set up graphical object properties
if length(faces_list) > 1
    transparency_order = linspace(0.25, 0.9, length(faces_list));
else
    transparency_order = 0.5;
end
color_order = distinguishable_colors(length(faces_list) + ~isempty(sources) + ~isempty(elc));

% plot mesh surfaces
for ii = 1:length(faces_list)
    faces = faces_list{ii};
    vertices = vertices_list{ii};
    
    P = patch('Faces',faces,'Vertices',vertices,'facecolor',color_order(ii,:),'edgecolor','none');
    set(P, 'facealpha', transparency_order(ii))
end

% plot sources separately for each hemisphere
if ~isempty(sources)
    lp = sources{1}.points;
    scatter3(lp(:,1), lp(:,2), lp(:,3), 20, color_order(length(faces_list)+1,:), 'filled')
    rp = sources{2}.points;
    scatter3(rp(:,1), rp(:,2), rp(:,3), 20, color_order(length(faces_list)+1,:), 'filled')
end

% plot electrodes
if ~isempty(elc)
    S = scatter3(elc(:,1), elc(:,2), elc(:,3), 150, color_order(length(faces_list) + ~isempty(sources)+1,:), 's', 'filled');
    set(S, 'MarkerFaceAlpha', 0.75)
end

camlight('headlight','infinite')
camorbit(0, 180)
camlight('headlight')
axis equal

end

