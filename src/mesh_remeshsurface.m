function [ bmesh ] = mesh_remeshsurface(bmesh, vertice_list)
%Function used to remesh a bmesh cell structure of BEM surfaces into a new
%set of vertice-resolution specified as a list of number of vertices. e.g.
% vertice_list = [14000, 12000, 10000]

edit_index = [];
newresolution = vertice_list;

for ii = 1:length(bmesh)
    if bmesh{ii}.nop ~= vertice_list(ii) && ~isnan(vertice_list(ii))
        disp(['Remeshing the ' num2str(ii) 'th surface'])
        edit_index = [edit_index, ii];
        % remesh the original mesh to a new vertice-resolution
        indStart=1; %Index of the start point
        numSeeds=vertice_list(ii); %number of vertices in the new mesh
        optionStruct.toleranceLevel=0; %Tolerance for convergence
        optionStruct.waitBarOn=1; %Turn on/off waitbar
        [ f,v ]=remeshTriSurfDistMap(bmesh{ii}.e,bmesh{ii}.p,numSeeds,indStart,optionStruct);
        
        %update the mesh structure
        bmesh{ii}.p = v;
        bmesh{ii}.e = f;
        bmesh{ii}.nop = size(v,1);
        bmesh{ii}.noe = size(f,1);
        
        %for printing
        newresolution(ii) = size(v,1);
    end        
end

disp('mesh_remeshsurface() completed.')
if isempty(edit_index)
    disp('No surface was remeshed, already in specified resolution:')
    disp(newresolution)
else
    disp('New surface resolution:')
    disp(newresolution)
end

end

