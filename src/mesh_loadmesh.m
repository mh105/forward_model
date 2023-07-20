function [ f, v, hlink ] = mesh_loadmesh(filename, datapath, smooth_iter, post_ds_smooth_iter, goal_vertices, ploton)
%Function used to load a triangulation surface saved as .mat file by
%surf2mat.py, with the option to smooth, downsample, and post-downsample
%smooth.

if nargin < 3
    smooth_iter = 0;
    post_ds_smooth_iter = 0;
    goal_vertices = [];
    ploton = false;
elseif nargin < 4
    post_ds_smooth_iter = 0;
    goal_vertices = [];
    ploton = false;
elseif nargin < 5
    goal_vertices = [];
    ploton = false;
elseif nargin < 6
    ploton = false;
end

if isempty(goal_vertices) && post_ds_smooth_iter > 0
    warning('goal_vertices is empty, no downsampling and no post-downsampling smoothing will be done.')
    post_ds_smooth_iter = 0;
end

load(fullfile(datapath, filename), 'faces', 'vertices')
if ~isa(faces, 'double'); faces = double(faces); end
if ~isa(vertices, 'double'); vertices = double(vertices); end

% smoothing
if smooth_iter > 0
    % smooth with smoothsurf()
    [conn,~,~]=meshconn(faces, length(vertices));
    v=smoothsurf(vertices,[],conn,smooth_iter,0.5,'lowpass',0.5);
else
    v=vertices;
end

% % DEBUGGING ONLY
% figure;
% hold on
% P = patch('Faces',faces,'Vertices',v,'facecolor',[.5 .5 .5],'edgecolor','none');
% camlight('headlight','infinite')
% axis equal
% camorbit(0, 180)
% camlight('headlight')
% camorbit(0, 180)
% camorbit(0, 270)
% rotate3d on

% downsampling
if ~isempty(goal_vertices)
    % remesh the high resolution mesh to this many vertices
    indStart=1; %Index of the start point
    numSeeds=goal_vertices; %number of vertices in the new mesh
    optionStruct.toleranceLevel=0; %Tolerance for convergence
    optionStruct.waitBarOn=1; %Turn on/off waitbar
    [ f,v ]=remeshTriSurfDistMap(faces,v,numSeeds,indStart,optionStruct);
else
    f = faces;
end

% post-downsampling smoothing
if post_ds_smooth_iter > 0
    % smooth with smoothsurf()
    [conn,~,~]=meshconn(f, length(v));
    v=smoothsurf(v,[],conn,post_ds_smooth_iter,0.5,'lowpass',0.5);
end

% visualization
if ploton
    figure;

    ax1 = subplot(1,2,1);
    patch('Faces',faces,'Vertices',vertices,'facecolor',[.5 .5 .5],'edgecolor','none');
    camlight('headlight','infinite')
    axis equal
    title(strrep(filename, '.mat', ''), 'FontSize', 16, 'Interpreter','none')

    ax2 = subplot(1,2,2);
    patch('Faces',f,'Vertices',v,'facecolor',[.5 .5 .5],'edgecolor','none');
    camlight('headlight','infinite')
    axis equal
    if ~isempty(goal_vertices)
        titlestr = {['Smoothed for ' num2str(smooth_iter) ' iter'], ['Downsampled to ' num2str(goal_vertices) ' vtx']};
    else
        titlestr = {['Smoothed for ' num2str(smooth_iter) ' iter']};
    end
    title(titlestr, 'FontSize', 16)

    % link camera angles between the two subplots
    hlink = linkprop([ax1,ax2], {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
    rotate3d on

    % adjust the limits so we don't get erroneous holes due to camera
    % cut-off
    xlimit = xlim;
    xlim(xlimit*1.1);
    ylimit = ylim;
    ylim(ylimit*1.1);
    zlimit = zlim;
    zlim(zlimit*1.1);
    view(0, 89)
    drawnow
else
    hlink = [];
end

disp([filename, ' successfully loaded.'])

end
