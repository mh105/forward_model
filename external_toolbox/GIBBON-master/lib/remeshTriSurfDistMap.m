function [varargout]=remeshTriSurfDistMap(varargin)

% function [Fd,Vd,seedIndex,indSeed]=remeshTriSurfDistMap(F,V,numSeeds,startInds,optionStruct)
% ------------------------------------------------------------------------
% Inputs: 
% F = input mesh faces
% V = input mesh vertices
% numSeeds = number of desired vertices (must be much less than length(V))
% startInds = points to start the remeshing from (optional). These points will be included in the output mesh
% W = point weights (optional)
% Outputs: 
% Fd = Output mesh faces
% Vd = Output mesh vertices
% seedIndex = indices of selected output vertices for each index vertex
% indSeed = indices of input vertices selected for the output mesh
%
% 
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% Change log:
% 2018/03/27 updated to use meshDistMarch as alternative to mex functions
% 2018/03/27 Made use of seedIndex2triangulation function which also
% filters unused points. 
%
% ------------------------------------------------------------------------

%% PARSE INPUT

%Default option structure
defaultOptionStruct.toleranceLevel=0;
defaultOptionStruct.waitBarOn=true(1,1);

switch nargin
    case 3
        F=varargin{1};
        V=varargin{2};
        numSeeds=varargin{3};
        startInds=1; %Use first point as only start 
        optionStruct=defaultOptionStruct;
    case 4
        F=varargin{1};
        V=varargin{2};
        numSeeds=varargin{3};
        startInds=varargin{4};         
        optionStruct=defaultOptionStruct;
    case 5
        F=varargin{1};
        V=varargin{2};
        numSeeds=varargin{3};
        startInds=varargin{4};        
        optionStruct=varargin{5};
    otherwise
        error('Wrong number of input arguments!')
end

%% Compute surface seeds based on marching
optionStruct.numSeeds=numSeeds;
[optionStruct]=structComplete(optionStruct,defaultOptionStruct,1);
[d,seedIndex]=meshDistMarch(F,V,startInds,optionStruct);

%% Alex fix:
% for some freaking reason, some seedIndex are NaN. They don't seem to get
% updated in the meshDistMarch function. This is ok since we already
% repeated the process that many number of times to get enough vertices. We
% can manually update seedIndex for these guys.
vertex_idx = unique(seedIndex);
vertex_idx = vertex_idx(~isnan(vertex_idx));
assert(length(vertex_idx) == numSeeds, 'We do not end up with the same number of seeds as desired.')
selected_vertex = V(vertex_idx,:);

still_nan_vertex = find(isnan(seedIndex));
% make sure nothing strange is happening
assert(~any(ismember(still_nan_vertex, vertex_idx)), 'Some selected vertices have seedIndex as NaN!')

% figure out which seed to go to for these still NaN vertices 
for ii = 1:length(still_nan_vertex)
    seed_id = still_nan_vertex(ii);
    [mind, minidx] = min(vecnorm(selected_vertex - V(seed_id,:), 2, 2));
    % update the seedIndex and d 
    seedIndex(seed_id) = vertex_idx(minidx);
    d(seed_id) = mind;
end

% nothing should be NaN now
assert(sum(isnan(seedIndex)) == 0, 'Nope, something is still NaN!')
assert(length(unique(seedIndex)) == numSeeds, 'Unique seed vertices got changed during this fix. Please check.')

%% Derive triangulation
[Fd,Vd,indSeed]=seedIndex2triangulation(F,V,seedIndex);

%% Collect ouput
varargout{1} =  Fd;
varargout{2} =  Vd;
varargout{3} =  seedIndex;
varargout{4} =  indSeed;
varargout{5} =  d;
     
end
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
