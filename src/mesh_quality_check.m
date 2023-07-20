function [ faces,vertices ] = mesh_quality_check(faces, vertices, precheck)
%This function is used to check the quality of a triangulation surface, fix
%errors if present. It relies on the iso2mesh toolbox.

if nargin < 3
    precheck = true;
end

init_f_dim = size(faces);
init_v_dim = size(vertices);

%iso2mesh routines
disp('-------------------------------------')
disp('Running meshcheckrepair() routines...')
disp(' ')

% fix vertex error
[vertices,faces]=meshcheckrepair(vertices,faces,'dup');
[vertices,faces]=meshcheckrepair(vertices,faces,'isolated');

% defect detection before meshfix
if precheck
    [vertices,faces]=meshcheckrepair(vertices,faces,'open');
    [vertices,faces]=meshcheckrepair(vertices,faces,'intersect');
end

% meshfix
[vertices,faces]=meshcheckrepair(vertices,faces,'deep');
[vertices,faces]=meshcheckrepair(vertices,faces,'meshfix');

% defect check again after meshfix
try
    [vertices,faces]=meshcheckrepair(vertices,faces,'open');
    [vertices,faces]=meshcheckrepair(vertices,faces,'intersect');
catch
    % somehow we need to fix again
    [vertices,faces]=meshcheckrepair(vertices,faces,'meshfix');
    [vertices,faces]=meshcheckrepair(vertices,faces,'open');
    [vertices,faces]=meshcheckrepair(vertices,faces,'intersect');
    
    % last check, if wrong gives error
    [vertices,faces]=meshcheckrepair(vertices,faces,'open');
    [vertices,faces]=meshcheckrepair(vertices,faces,'intersect');
end

if ~all(init_f_dim == size(faces)) || ~all(init_v_dim == size(vertices))
    disp(' ')
    disp('Triangulation surface was modified!')
    disp('Initial dimensions: Faces = ')
    disp(init_f_dim)
    disp('Output dimensions: Faces = ')
    disp(size(faces))
    disp('Initial dimensions: Vertices = ')
    disp(init_v_dim)
    disp('Output dimensions: Vertices = ')
    disp(size(vertices))
    
%     if hardnumber
%         % number of vertices is changed. We need to get back to the same number
%         % of vertices we started with.
%         disp('Triangulation surface was modified to a different vertex number.')
%         disp('Vertex number is specified to be hard fixed. We will upsample and remesh again!')
%         attemptn = 1;
%         while ~all(init_v_dim == size(vertices)) && atemptn < 2
%             disp(['Attempt No.' num2str(attemptn), '...'])
%             % subdivide to upsample
%             nsub = 1; % number of subdivision steps
%             options.sub_type = 'loop';
%             options.verb = 0;
%             [vertices,faces] = perform_mesh_subdivision(vertices',faces',nsub,options);
%             vertices = vertices'; faces = faces';
%             
%             % remesh the high resolution mesh to the original number of
%             % vertices
%             indStart=1; %Index of the start point
%             numSeeds=init_v_dim(1); %number of vertices in the new mesh
%             optionStruct.toleranceLevel=0; %Tolerance for convergence
%             optionStruct.waitBarOn=1; %Turn on/off waitbar
%             [ faces,vertices ]=remeshTriSurfDistMap(faces,vertices,numSeeds,indStart,optionStruct);
%             
%             % then we run the check again
%             [ faces,vertices ] = mesh_quality_check(faces, vertices, false, false);
%         end
%     end
end

disp(' ')
disp('mesh_quality_check() completed!')
disp('-------------------------------------')

end

