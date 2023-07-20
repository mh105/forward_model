function [ forward_model ] = prep_forward_model_struct(filename, fsfn, atlasfn, savefn)
%Function used to set up a forward model structure from a computed BEM
%model using hbf solver pipeline saved as -fwd.mat
if nargin < 4
    savefn = [];
end

%% Load the -fwd.mat file
load(filename, 'bmeshes', 'dig', 'G', 'sources', 'forward_model')

%% Load eloc struct with appropriate coordinate system
% dig loaded in -fwd.mat is in MRI coordinate system and already projected
% onto the scalp. EEGLAB topoplot() function expects the electrode
% positions to be present in a different coordiante system with shifted
% origin and rotation. Without the transformation matrix and additional
% application of rotation, we cannot get it to display correctly with
% topoplot() function. 

% Hence, we will have to load the raw fastscan data created by
% ANT_MNE_fastscan.m that employs an algorithm to match digitization sensor
% points to a template eloc structure provided by ANT that is in the
% correct coordinate system as EEGLAB. In these coordinates, electrodes are
% transformed into the EEGLAB topoplot() coordinate system and therefore
% will be displayed nicely.
load(fsfn)

% we need to reorder the fastscan data to the order of ANT_duke template,
% with electrode 1 being LM. 
eloc = fastscan.fschanlocs;
for ii = 1:length(eloc)
    eloc(fastscan.fschanlocs(ii).ognum) = fastscan.fschanlocs(ii);
end

%% Re-referencing - Matrix Operations
% we will derive everything in a single sequence of matrix
% multiplications to derive a relationship in the form:

% Y = R x G x X
%   Y - channel x time: sensor space time series
%   X - source x time: source space time series
%   G - channel x source: forward model leadfield matrix
%   R - channel x channel: re-ference matrix

% Together, we will denote F = R x G as an augmented forward model that
% captures referencing:
%
%                               Y = F x X
%
%   F - channel x source: augmented forward model matrix

%% Compute reference matrices
% we will compute the different R matrices for different referencing
% schemes
% Unipolar references:
% R1 - common average
% R2 - single electrode (recording reference at Z3)
% R3 - single mastoid (left mastoid)
% R4 - linked mastoid
% R5 - REST

% Non-unipolar references:
% R6 - contralateral mastoid
% R7~R11 - Laplacian

ref_label = {'Common Average',... 
    'Recording Reference at Z3',...
    'Left Mastoid',...
    'Linked Mastoid',...
    'REST',...
    'Contralateral mastoid',...
    '4 neightbor LP',...
    'sphere radius 50mm LP',...
    'Manual LP',...
    'Extended Manual LP',...
    'Extended sphere raiud 100mm LP'};

%%%%%%%%%%%%%%%%%%%%%%
% Unipolar references:
%%%%%%%%%%%%%%%%%%%%%%
% R1 - Common average
R1 = eye(size(G,1)) - ones(size(G,1))*1/size(G,1);

% R2 - Single electrode (recording reference at 85)
ref_chan = 85;
R2 = eye(size(G,1));
R2(:, ref_chan) = R2(:, ref_chan) - 1;

% R3 - Single mastoid (left mastoid 1)
R3 = eye(size(G,1));
R3(:, 1) = R3(:, 1) - 1;

% R4 - Linked mastoid (left mastoid 1 and right mastoid 47)
R4 = eye(size(G,1));
R4(:, [1, 47]) = R4(:, [1, 47]) - 1/2;

% R5 - REST
fREST = pinv(G)'*pinv(G)*ones(size(G,1), 1) ./ (ones(size(G,1), 1)'*pinv(G)'*pinv(G)*ones(size(G,1), 1));
R5 = eye(size(G,1)) - ones(size(G,1), 1)*fREST';

% Assert that these are indeed unipolar reference operators
assert(testUR(R1) & testUR(R2) & testUR(R3) & testUR(R4) & testUR(R5), 'Failed unipolar reference test')


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-unipolar references:
%%%%%%%%%%%%%%%%%%%%%%%%%%
% R6 - Contralateral mastoid reference
R6 = zeros(size(G,1));
for i = 1:size(R6,1)
    R6(i, i) = 1;
    if dig(i, 1) > 0 % right-side electrodes
        R6(i, 1) = -1;
    else
        R6(i, 47) = -1;
    end
end
R6(1,1) = 0; R6(1,47) = 0;
R6(47,47) = 0; R6(47,1) = 0;

% R7 - Laplacian - 4 nearest neighbors
R7 = zeros(size(G,1));
for i = 1:size(R7,1)
    R7(i, i) = 1;
    distance_list = [];
    % search all electrodes
    for j = 1:size(dig,1)
        distance_list(j) = norm(dig(i,:) - dig(j,:));
    end
    
    % get rid of the sensor itself
    distance_list(i) = inf;
    
    % find the smallest 4 neighbor electrodes
    [~, index] = sort(distance_list);
    R7(i, index(1:4)) = R7(i, index(1:4)) - 1/4;
end


% R8 - Laplacian - sphere radius threshold
R8 = zeros(size(G,1));
for i = 1:size(R8,1)
    R8(i, i) = 1;
    neighbor_list = [];
    % search within radius
    for j = 1:size(dig,1)
        if norm(dig(i,:) - dig(j,:)) <= 0.050 % 50mm / 5cm
            neighbor_list = [neighbor_list, j];
        end
    end
    % get rid of the sensor itself
    neighbor_list(neighbor_list == i) = [];
    R8(i, neighbor_list) = -1 / length(neighbor_list);
    
%     figure
%     hold on
%     scatter3(dig(:,1),dig(:,2),dig(:,3),10,'k','filled')
%     scatter3(dig(i,1),dig(i,2),dig(i,3),30,'r','filled')
%     scatter3(dig(neighbor_list,1),dig(neighbor_list,2),dig(neighbor_list,3),20,'g','filled')
%     axis equal
%     rotate3d on
%     pause()
%     close all
end


% R9 - Laplacian - manual definition based on duke waveguard layout
load('duke_128_channelneighbors.mat', 'channelneighbors')
R9 = zeros(size(G,1));
for i = 1:size(R9,1)
    R9(i, i) = 1;
    neighbor_list = find(channelneighbors(i,:) == 1);
    R9(i, neighbor_list) = -1 / length(neighbor_list);
end


% R10 - Laplacian - manual neighbors extended by second order
R10 = zeros(size(G,1));
for i = 1:size(R10,1)
    R10(i, i) = 1;
    fo_neighbor_list = find(channelneighbors(i,:) == 1);
    all_neighbor_list = fo_neighbor_list;
    for j = 1:length(fo_neighbor_list)
        all_neighbor_list = [all_neighbor_list, find(channelneighbors(fo_neighbor_list(j),:) == 1)];
    end
    all_neighbor_list = unique(all_neighbor_list);
    % get rid of the sensor itself
    all_neighbor_list(all_neighbor_list==i) = [];
    R10(i, all_neighbor_list) = -1 / length(all_neighbor_list);
end


% R11 - Laplacian - extended sphere radius threshold
R11 = zeros(size(G,1));
for i = 1:size(R11,1)
    R11(i, i) = 1;
    neighbor_list = [];
    % search within radius
    for j = 1:size(dig,1)
        if norm(dig(i,:) - dig(j,:)) <= 0.100 % 100mm / 10cm
            neighbor_list = [neighbor_list, j];
        end
    end
    % get rid of the sensor itself
    neighbor_list(neighbor_list == i) = [];
    R11(i, neighbor_list) = -1 / length(neighbor_list);
    
%     figure
%     hold on
%     scatter3(dig(:,1),dig(:,2),dig(:,3),10,'k','filled')
%     scatter3(dig(i,1),dig(i,2),dig(i,3),30,'r','filled')
%     scatter3(dig(neighbor_list,1),dig(neighbor_list,2),dig(neighbor_list,3),20,'g','filled')
%     axis equal
%     rotate3d on
%     pause()
%     close all
end

% Assert that these are indeed non-unipolar reference operators
assert(testUR(R6, true) & testUR(R7, true) & testUR(R8, true) & testUR(R9, true) & testUR(R10, true) & testUR(R11, true), 'Failed unipolar reference test')

% Construct augmented forward model matrices
modified_G = cell(1,11);
Reference_mat = [];
for i = 1:length(modified_G)
    % Save the augmented forward model matrices
    eval(sprintf('modified_G{i} = R%d*G;', i))
    % Create a large matrix with reference transformation matrices
    eval(sprintf('Reference_mat(i,:,:) = R%d;', i))
end

%% Load the anatomical ROI based on Desikan atlas
% load atlas_info structure
% this also loads the cortical surface
load(atlasfn, 'atlas_info', 'fce', 'vtx', 'idx', 'src_face');

% check the atlas_info for accuracy
totalvertex = 0;
for ii = 1:length(atlas_info.atlas_vertidx)
    totalvertex = totalvertex + length(atlas_info.atlas_vertidx{ii});
end
assert(totalvertex < size(G,2), 'Atlas labelled more sources than available! Something is wrong!')

% convert the source face to indexing on the source indices rather than
% vertex indices of the whole cortical surface 
assert(length(idx) == length(vtx), 'Size of select index mismatches the vertex number of cortical surface!')
vertno = find(idx == 1);
[ logic_a, source_face ] = ismember(src_face, vertno);
assert(all(logic_a, 'all'), 'Not all vertices contained in src_face are being used as a source!')

%% Create output structure

fwd = struct;

% lead field matrix variables
fwd.G = forward_model.G;
fwd.xyzG = forward_model.xyzG;
fwd.source_face = source_face;
fwd.source = forward_model.source;
fwd.normal = forward_model.normal;
fwd.dig = forward_model.dig;

% cortical surface variables
fwd.cortex_surf_face = double(fce);
fwd.cortex_surf_vertex = vtx;
fwd.select_idx = logical(idx);
% Anatomical ROI source indices 
fwd.anato_ROI = atlas_info;

% skin surface
fwd.skin_surf_face = bmeshes{end}.e;
fwd.skin_surf_vertex = bmeshes{end}.p;

% also save the BEM surfaces
fwd.bmeshes = bmeshes;

% Sensor space variables for EEGLAB topoplot()
fwd.eloc = eloc;

% Augmented LFM with Referencing
fwd.modified_G = modified_G;
fwd.ref_matrix = Reference_mat;
fwd.ref_label = ref_label';



% update struct name 
forward_model = fwd;

%%
% save to disk
if ~isempty(savefn); save(savefn, 'forward_model'); end

end



