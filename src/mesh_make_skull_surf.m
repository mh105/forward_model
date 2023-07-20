function [ ] = mesh_make_skull_surf(datapath, verbose)
%% Make skull surface
% It seems that the mri2mesh skull surface in the anterior inferior portion
% is very unreliable from scan to scan. We can't devise a general enough
% algorithm to extract resonable outer skull surface from just mri2mesh.
% Hence, we will adopt the strategy to weave together headreco skull and
% MNE_outer_skull in the anterior inferior portion to put together a good
% skull surface.

if nargin < 2
    verbose = false;
end

%% Path configuration
% addpath(genpath('./'));

% datapath should point to the final_structure folder
% datapath = '/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/MNE/subject_data/SP011_ah/simnibs_pipe/bem_surfaces/final_structure';

tic

%%
%%%%%%%%%%%%%%%%%%%
%%% Preparations 
%%%%%%%%%%%%%%%%%%%

%% Load everything first 
% Load csf surface
if isfile(fullfile(datapath, 'headreco_csf_compatible.mat'))
    load(fullfile(datapath, 'headreco_csf_compatible.mat'))
else
    load(fullfile(datapath, 'mri2mesh_csf_fixed.mat'))
end
f1 = double(faces); %#ok<*NODEF>
v1 = double(vertices);

% Load mri2mesh skull surface
load(fullfile(datapath, 'mri2mesh_skull_fixed.mat'))
f2 = double(faces);
v2 = double(vertices);

% Note on 06/09/2020: now the mri2mesh_skull is not used other than
% finding the frontal and occipital poles of the mri2mesh surface to define
% the quadrant cutting planes. This is because mri2mesh_skull is not very
% reliable overall, so we will upsample MNE skull and project to headreco
% directly in the upper_posterior 3 quadrants, skipping the need for
% mri2mesh_skull all together, producing much much better skull surfaces. 

% Load headreco bone surface
load(fullfile(datapath, 'headreco_bone.mat'))
f3 = double(faces);
v3 = double(vertices);

% Load MNE outer skull surface
load(fullfile(datapath, 'MNE_outer_skull.mat'))
f4 = double(faces);
v4 = double(vertices);

%% Figure out whether a rotation is needed 
% sometimes the head is tilted, making defining planes complicated. We will
% fix explore whether a rotation is needed.
close all

% Do we need a rotation? 
figure(1);
hold on
patch('Faces',f3,'Vertices',v3,'facecolor',[.5 .5 .5],'edgecolor','none');
camlight('headlight','infinite')
axis equal
camorbit(0, 180)
camlight('headlight')
camorbit(0, 180)
camorbit(0, 270)
rotate3d on
view(180, 0)
title('Before Rotation', 'FontSize', 12)

rotation_attempt_count = 0;
stoprotate = false;

figure(2);
hold on
camlight('headlight','infinite')
axis equal
camorbit(0, 180)
camlight('headlight')
camorbit(0, 180)
camorbit(0, 270)
rotate3d on
title('After Rotation', 'FontSize', 12)

while ~stoprotate
    Zdeg =  input('Z rotation (looking left / right) in degree: ', 's');
    Ydeg =  input('Y rotation (tilting left / right) in degree: ', 's');
    Xdeg =  input('X rotation (tilting forward / backward) in degree: ', 's');
    
    Zdeg = str2double(Zdeg);
    t = deg2rad(Zdeg);
    Rz = [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1];
    
    Ydeg = str2double(Ydeg);
    t = deg2rad(Ydeg);
    Ry = [cos(t) 0 -sin(t); 0 1 0; sin(t) 0 cos(t)];
    
    Xdeg = str2double(Xdeg);
    t = deg2rad(Xdeg);
    Rx = [1 0 0; 0 cos(t) sin(t); 0 -sin(t) cos(t)];
    
    v3_rot = v3 * Rz * Ry * Rx;
    
    % visualize the rotated headreco bone
    figure(2);
    P = patch('Faces',f3,'Vertices',v3_rot,'facecolor',[.5 .5 .5],'edgecolor','none');
    view(180, 0)
    
    goodrotate =  input('Is rotation result good? (y/n): ', 's');
    if strcmpi(goodrotate, 'y')
        stoprotate = true;
    else
        delete(P)
    end
    
    rotation_attempt_count = rotation_attempt_count + 1;
end

close all

%% Separate the surfaces into 4 quadrants with two planes
% compute rotated surfaces 
v4_rot = v4 * Rz * Ry * Rx;
v3_rot = v3 * Rz * Ry * Rx;
v2_rot = v2 * Rz * Ry * Rx;

figure(1);
hold on

% visualize headreco bone
P = patch('Faces',f3,'Vertices',v3_rot,'facecolor',[.5 .5 .5],'edgecolor','none');
set(P, 'facealpha', 0.25)
camlight('headlight','infinite')
axis equal
camorbit(0, 180)
camlight('headlight')
camorbit(0, 180)
camorbit(0, 270)
rotate3d on

% visualize MNE_outer_skull
P = patch('Faces',f4,'Vertices',v4_rot,'facecolor',[0 1 0],'edgecolor','none');
set(P, 'facealpha', 0.2)

title({'Zplane should be just above eye sockets', 'Yplane should be just behind spinal cord'}, 'FontSize', 20)

% visualize the SurfaceRAS center origin 
scatter3(0, 0, 0, 100, 'k', 'filled')

% Updated on 03/03/2022: usually the automatic quadrant cutting is fairly
% good, therefore we will make it the default and allow users to adjust if
% the final surfaces and quadrants are not satisfactory rather than
% returning error and having to repeat this function.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find frontal and occipital poles of the mri2mesh surface
[maxy, fpidx] = max(v2_rot(:,2));
[miny, opidx] = min(v2_rot(:,2));

scatter3(v2_rot(fpidx,1), v2_rot(fpidx,2), v2_rot(fpidx,3), 500, 'b', 'filled')
scatter3(v2_rot(opidx,1), v2_rot(opidx,2), v2_rot(opidx,3), 500, 'g', 'filled')

% we now define an axial plane to cut the surfaces into two parts
maxplane = max([max(v2_rot(:,1)), max(v3_rot(:,1)), max(v2_rot(:,2)), max(v3_rot(:,2))]);
minplane = min([min(v2_rot(:,1)), min(v3_rot(:,1)), min(v2_rot(:,2)), min(v3_rot(:,2))]);
vxtvalue_Z = max(abs(maxplane), abs(minplane))+10;

% pick the more superior point of the two poles. We pick the superior one
% to avoid having to deal with sharp transitions around eye sockets
zthresh = max(v2_rot(fpidx,3), v2_rot(opidx,3));
zthresh_keep = zthresh;

plane_x1 = [vxtvalue_Z, vxtvalue_Z, -vxtvalue_Z, -vxtvalue_Z];
plane_y1 = [vxtvalue_Z, -vxtvalue_Z, -vxtvalue_Z, vxtvalue_Z];
plane_z1 = ones(1,4) .* zthresh;

% visualize the plane for cutting
P1 = patch(plane_x1, plane_y1, plane_z1, [0,0,0,0],'facecolor',[0 0 1], 'facealpha', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We also find the most inferior point of the MNE outer skull surface,
% in order to define a coronal plane that cuts the headreco surface for a
% anterior-inferior portion, which contains most of the convoluted surface
% edges. 
[minz, ipidx] = min(v4_rot(:,3));
ythresh = v4_rot(ipidx,2);
ythresh_keep = ythresh;

scatter3(v4_rot(ipidx,1), v4_rot(ipidx,2), v4_rot(ipidx,3), 500, 'y', 'filled')

% we now define a coronal plane to further cut the surfaces into four parts
maxplane = max([max(v4_rot(:,1)), max(v3_rot(:,1)), max(v4_rot(:,3)), max(v3_rot(:,3))]);
minplane = min([min(v4_rot(:,1)), min(v3_rot(:,1)), min(v4_rot(:,3)), min(v3_rot(:,3))]);
vxtvalue_Y = max(abs(maxplane), abs(minplane))+10;

% automatic adjustment of the yplane
stopsearch = false;
while ~stopsearch
    
    [MNE_miny, mneminidx] = min(v4_rot(abs(v4_rot(:,2) - ythresh) < 0.5, 3));
    mne_min = find(abs(v4_rot(:,2) - ythresh) < 0.5);
    mne_min = mne_min(mneminidx);
    S1 = scatter3(v4_rot(mne_min,1), v4_rot(mne_min,2), v4_rot(mne_min,3), 500, 'm', 'filled');
    
    [headreco_miny, minidx] = min(v3_rot(abs(v3_rot(:,2) - ythresh) < 0.5, 3));
    headreco_min = find(abs(v3_rot(:,2) - ythresh) < 0.5);
    headreco_min = headreco_min(minidx);
    S2 = scatter3(v3_rot(headreco_min,1), v3_rot(headreco_min,2), v3_rot(headreco_min,3), 500, 'c', 'filled');
    
    if headreco_miny - MNE_miny < -3
        ythresh = ythresh - 0.01;
        delete(S1)
        delete(S2)
    else
        stopsearch = true;
    end
end

plane_x2 = [vxtvalue_Y, vxtvalue_Y, -vxtvalue_Y, -vxtvalue_Y];
plane_y2 = ones(1,4) .* ythresh;
plane_z2 = [vxtvalue_Y, -vxtvalue_Y, -vxtvalue_Y, vxtvalue_Y];

% visualize the plane for cutting
P2 = patch(plane_x2, plane_y2, plane_z2, [0,0,0,0],'facecolor',[0 0 1], 'facealpha', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% manual sanity check: do the planes and surfaces look reasonable?
view(-90,0)
checkok =  input('Final check: do the surfaces and 4 quadrants look ok? (y/n): ', 's');

if strcmpi(checkok, 'n')

    disp('Enabling manual adjustment of the Z and Y cutting planes...')
    figure(1);
    delete(P1)
    delete(P2)
    delete(S1)
    delete(S2)
    
    % Manual adjustment of the Z-plane 
    Zplane_ready = false;
    zthresh = zthresh_keep; % reset zthresh to initial value
    
    while ~Zplane_ready
        plane_x = [vxtvalue_Z, vxtvalue_Z, -vxtvalue_Z, -vxtvalue_Z];
        plane_y = [vxtvalue_Z, -vxtvalue_Z, -vxtvalue_Z, vxtvalue_Z];
        plane_z = ones(1,4) .* zthresh;
        
        % visualize the plane for cutting
        figure(1);
        P = patch(plane_x, plane_y, plane_z, [0,0,0,0],'facecolor',[0 0 1], 'facealpha', 1);
        
        % manual adjustment
        Zneedshift =  input('Is shifting needed for the Zplane? (y/n): ', 's');
        if strcmpi(Zneedshift, 'y')
            Zshift =  input('How much to shift the Zplane? : ', 's');
            Zshift = str2double(Zshift);
            zthresh = zthresh_keep + Zshift;
            figure(1);
            delete(P)
        elseif strcmpi(Zneedshift, 'n')
            Zplane_ready = true;
        end
    end
    
    % Manual adjustment of the Y-plane
    Yplane_ready = false;
    ythresh = ythresh_keep; % reset ythresh to initial value
    
    while ~Yplane_ready
        plane_x = [vxtvalue_Y, vxtvalue_Y, -vxtvalue_Y, -vxtvalue_Y];
        plane_y = ones(1,4) .* ythresh;
        plane_z = [vxtvalue_Y, -vxtvalue_Y, -vxtvalue_Y, vxtvalue_Y];
        
        % visualize the plane for cutting
        figure(1);
        P = patch(plane_x, plane_y, plane_z, [0,0,0,0],'facecolor',[1 0.5 0.33], 'facealpha', 1);
        
        % manual adjustment
        Yneedshift =  input('Is shifting needed for the Yplane? (y/n): ', 's');
        if strcmpi(Yneedshift, 'y')
            Yshift =  input('How much to shift the Yplane? : ', 's');
            Yshift = str2double(Yshift);
            ythresh = ythresh_keep + Yshift;
            figure(1);
            delete(P)
        elseif strcmpi(Yneedshift, 'n')
            Yplane_ready = true;
        end
    end
    
end

close all

%%
%%%%%%%%%%%%%%%%%%%
%%% Start Processing
%%%%%%%%%%%%%%%%%%%

%% Overall game plan
% Our overall goal is the following:

% We feel pretty good about the headreco surface for the top-back 3
% quadrants. However the first quadrant details of the headreco surface are
% uncessary and problematic. 

% Hence, the overall plan is that we want to make the MNE outer skull
% surface to follow the headreco surface for the top-back 3 quadrants, and
% we will keep the MNE out skull surface as is for the anterior-inferior
% quadrant. Then we need to do some kind of checking to make sure none of
% the csf surface protrudes outside of the skull surface. It is does, then
% we dilate the skull surface to accommodate these positions. 

%% Step 1: Dilate-project MNE to headreco 

% 1) We need to first upsample the MNE outer skull surface 
[F,V]=subtri(f4,v4,1);
for ii = 1:1
    [F,V]=subtri(F,V,1);
end
V_rot = V * Rz * Ry * Rx;

% 2) figure out how many points in the up3quadrants are inside the headreco
% surface 
mne_vindex = V_rot(:,3) >= zthresh | V_rot(:,2) <= ythresh;
testp = V_rot(mne_vindex, :);
in = intriangulation(v3_rot,f3,testp);

% find out the right indices on V_rot
project_index = find(mne_vindex==1);
dilate_index = project_index(in==1);

if verbose
    figure;
    hold on
    P = patch('Faces',F,'Vertices',V_rot,'facecolor',[.5 .5 .5],'edgecolor','none');
    set(P, 'facealpha', 0.25)
    camlight('headlight','infinite')
    axis equal
    camorbit(0, 180)
    camlight('headlight')
    camorbit(0, 180)
    camorbit(0, 270)
    rotate3d on
    view(180, 0)
    scatter3(V_rot(project_index,1), V_rot(project_index,2), V_rot(project_index,3), 10, 'b', 'filled')
    scatter3(V_rot(dilate_index,1), V_rot(dilate_index,2), V_rot(dilate_index,3), 10, 'y', 'filled')
    % overlay the headreco skull surface as well
    P = patch('Faces',f3,'Vertices',v3_rot,'facecolor',[.9 .5 .5],'edgecolor','none');
    set(P, 'facealpha', 0.25)
    
    scatter3(0, 0, 0, 100, 'k', 'filled')
end

% Grab the upper-3-quadrant headreco surface for faster projection.
% instead of using the entire headreco surface which takes forever, we will
% grab a section of the headreco that's relevant and process that section.  

vindex = v3_rot(:,3) > zthresh-5 | v3_rot(:,2) < ythresh+5;
vin = find(vindex==1);

face_idx = zeros(size(f3,1), 1);
for ii = 1:size(f3,1)
    if any(ismember(f3(ii,:), vin))
        face_idx(ii) = 1;
    end
end

f3_above = f3(logical(face_idx), :);

% we also need to clean up the indexing of this new triangulation surface
v_include = double(unique(f3_above));
sub_v3 = v3_rot(v_include,:);

% clean the indices contained in f3_above to correspond to sub_v3
for ii = 1:size(f3_above,1)
    [logiout, idxout] = ismember(f3_above(ii,:), v_include);
    assert(all(logiout), 'Some vertices called by f3_above not in v_include!')
    f3_above(ii,:) = idxout;
end

% Project + Dilate MNE_outer_skull onto headreco skull surface 
V_rot = mesh_dilateproject(f3_above, sub_v3, F, V_rot, dilate_index, project_index);

% smooth the surface 
[conn,connnum,count]=meshconn(F, length(V_rot));
V_rot=smoothsurf(V_rot,[],conn,100,0.8,'lowpass',0.5);

if verbose
    figure;
    hold on
    P = patch('Faces',F,'Vertices',V_rot,'facecolor',[.5 .5 .5],'edgecolor','none');
    set(P, 'facealpha', 0.5)
    camlight('headlight','infinite')
    axis equal
    camorbit(0, 180)
    camlight('headlight')
    camorbit(0, 180)
    camorbit(0, 270)
    rotate3d on
    view(180, 0)
    title('After Step 1', 'FontSize', 20)
end

disp('Step 1 Done!')

%% Step 2: Check for inclusion of csf surface 

% rotate it back 
V = V_rot / Rx / Ry / Rz;

% find out csf surface points that are outside of the skull surface
in = intriangulation(V,F,v1);
out_vpoint = v1(in==0,:);

if verbose
    figure;
    hold on
    P = patch('Faces',F,'Vertices',V,'facecolor',[.5 .5 .5],'edgecolor','none');
    set(P, 'facealpha', 0.5)
    camlight('headlight','infinite')
    axis equal
    camorbit(0, 180)
    camlight('headlight')
    camorbit(0, 180)
    camorbit(0, 270)
    rotate3d on
    view(180, 0)
    % overlay the csf surface as well
    patch('Faces',f1,'Vertices',v1,'facecolor',[.9 .5 .5],'edgecolor','none');
    scatter3(out_vpoint(:,1), out_vpoint(:,2), out_vpoint(:,3), 10, 'b', 'filled')
    title('Pre-fixing for csf surface during Step 2', 'FontSize', 20)
end

% Repeat until the csf surface is entirely enclosed in the updated skull surface
counter = 0;
while sum(in==0) > 0
    % for every one of these out points, project to the skull surface
    [ distances, surface_points ] = ...
        point2trimesh('Faces', F, 'Vertices', V, 'QueryPoints', out_vpoint, 'Algorithm', 'parallel');
    % dilate the skull surface
    V = mesh_dilate_out(V, F, out_vpoint, surface_points, 1.5, [], 3, 6); % use 6-order neighbors to smooth it out with high 1.5mm margin
    
    % smooth the surface a little bit more
    [conn,connnum,count]=meshconn(F, length(V));
    V=smoothsurf(V,[],conn,20,0.8,'lowpass',0.5);
    
    % check again for csf surface points that are outside of skull surface
    in = intriangulation(V,F,v1);
    out_vpoint = v1(in==0,:);
    
    counter = counter + 1;
end

if verbose
    figure;
    hold on
    P = patch('Faces',F,'Vertices',V,'facecolor',[.5 .5 .5],'edgecolor','none');
    set(P, 'facealpha', 0.5)
    camlight('headlight','infinite')
    axis equal
    camorbit(0, 180)
    camlight('headlight')
    camorbit(0, 180)
    camorbit(0, 270)
    rotate3d on
    view(180, 0)
    % overlay the csf surface as well
    patch('Faces',f1,'Vertices',v1,'facecolor',[.9 .5 .5],'edgecolor','none');
    title('Post-fixing for csf surface during Step 2', 'FontSize', 20)
end

disp(['Number of loops of skull surface dilation to include the CSF surface: ', num2str(counter)])
disp('Step 2 Done!')

%% Save the hardwork result!
% Save data to file
faces = F;
vertices = V;

save(fullfile(datapath, 'morphed_skull.mat'), 'faces', 'vertices')

% Done!

disp('Total time taken:') 
toc

pause(5)
close all

end
