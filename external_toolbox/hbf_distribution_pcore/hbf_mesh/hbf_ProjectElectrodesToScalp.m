function elecs=hbf_ProjectElectrodesToScalp(elecpos,meshes,Nmin)
% HBF_PROJECTELECTRODESTOSCALP projects a set of electrodes to the outermost
% mesh, makes some interpolation coefficients, and builds the hbf struct "elecs".
%
% elecs=HBF_PROJECTELECTRODESTOSCALP(elecpos,meshes,Nmin)
% 
% elecpos:  positions of electrodes (somewhat coregistered...)
% bmeshes:  hbf struct bmeshes
% Nmin:     Number of nearest triangles that are used for searching for the
%   nearest position (optional). When projecting nearby electrodes to a
%   typical scalp mesh, Nmin=1 is typically fine. If the point is already
%   on the surface, Nmin==number of triangles that belong to a vertex
%   (typically about 6) is recommended to guarantee finding the absolutely
%   nearest point. The default is Nmin=6.
%
% v161220 Matti Stenroos

if nargin==2 || isempty(Nmin)
    Nmin=6;
end
elecs.porig=elecpos;
[elecs.pproj,elecs.loctype,elecs.locinfo,elecs.projdist]=hbf_ProjectPointsToMesh(meshes{end},elecpos,Nmin);
elecs.NtoE=sparse(hbf_LinearInterpolationMatrix(meshes{end},elecs.pproj,elecs.loctype,elecs.locinfo));


