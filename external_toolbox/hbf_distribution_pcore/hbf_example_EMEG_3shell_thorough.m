% Example on the use of Helsinki BEM Framework kernel:
% 3-shell MEG+EEG model.
% 
% The example uses a three-shell head model made by me at MRC Cognition and
% Brain Sciences Unit, Cambridge, UK. The model was made using FLASH MR
% images, FreeSurfer, and MNE toolbox. If you wish to use this model for
% something that you wish to share or publish, please contact me.
%
% (c) Matti Stenroos
% 4 Nov 2015
clear
%add BEM framework paths
addpath ./hbf_calc
addpath ./hbf_mesh
% load example data for meshes, coils, and
% sources.
load('/proj/bioem/exampledata/cbusample/samplehead_cbu.mat');%bmeshes, eegsens, megsens, sources

% 1. Boundary meshes
%Boundary meshes are described as a cell array Nx1, where each cell is a
%struct that contains fields "p" for points (vertices), and "e" for
%elements (faces, triangles)
% bmeshes: 3 x 1 cell,
% bmeshes{I} =
%   p: [Number of vertices x 3]
%   e: [Number of triangles x 3]
% bmeshes{1}=
%   p: [2562x3 double]
%   e: [5120x3 double]

%The convention of this BEM framework is that the innermost mesh has number
%1 and outermost mesh M; this is mandatory. If the model is nested (like a
%3-shell model), the order can be checked/set with
bmeshes=hbf_SortNestedMeshes(bmeshes);

% 2. MEG sensors
%Set MEG coils. There are some options for describing coils:
% coils = 
%   QP: field computations points [Number of field computation points x 3]
%   QN: sensor orientation        [Number of field computation points x 3];
%       each must be unit length
%   QtoC: integral weights [Number of sensors x Number of field computation points]
%   p:  sensor positions [Number of sensors x 3]
%   n:  sensor normals [Number of sensors x 3]; each must be unit-length
%   name
%   description
% You must give (QP and QN and QtoC) or (p and n). If all these are given,
% the points given in QP & QN & QtoC are used in field computation. In this
% case, we have a very simple coil model: 102 point-like magnetometers.
% These are modelled with 102 points and orientations. I also made a simple
% triangulation for the coilset (element description: 'e'), so you can
% easily plot field topographies on the sensor layout in either 3D or 2D.
coils=megsens;
% coils = 
% 
%       p: [102x3 double]
%       n: [102x3 double]
%       e: [179x3 double]
%     p2d: [102x2 double]

%3. EEG electrodes
%Set EEG electrodes. The electrodes are assumed to lie on scalp; if you
%wish to place electrodes elsewhere, contact me. The electrodes need to be
%roughly coregistered to scalp mesh. In this case, we have 70 electrodes
%and also triangulation like in the case of coils.
% eegsens = 
% 
%       p: [70x3 double], original positions before projection
%       e: [124x3 double], element descriptions
%     p2d: [70x2 double], stereographic projection of electrodes to 2d.
elecpos=eegsens.p;%
%Project electrodes to scalp and build interpolation operator from model
%nodes to projected electrodes.
elecs=hbf_ProjectElectrodesToScalp(elecpos,bmeshes);
%4. Source space
sourcepos=sources{1}.p;%source locations, [Ns x 3] array
%source directions, [Ns x 3] --- if general lead field matrix is desired, 
%each vector must be unit length!
sourcedir=sources{1}.nn;

%5. Conductivities
% Set conductivities for the head model. In this example, we use a (nested)
% 3-layer model, so we need to give conductivities for the brain, skull,
% and scalp compartments. Remember to give the conductivities in the
% correct order!
ci=[1 1/50 1]*.33;%conductivity inside each surface --- remember the order!
co=[ci(2:3) 0];%conductivity outside each surface

%* First time with new meshing: tests...
%Problems with the BEM are typically due to unsuitable or
%ill-specified meshes. Even if meshes are otherwise good, the orientation of the triangles is often
%wrong; this is the most common reason for computations going wrong. This BEM framework assumes To check and flip the orientation,
%you can just type
Nmeshes=length(bmeshes);
success=zeros(Nmeshes,1);
for M=1:Nmeshes;
    [bmeshes{M},success(M)]=hbf_CorrectTriangleOrientation(bmeshes{M});
end
%You can/should also do further quality checks using, e.g., this simple
%tool. If orientation correction failed, this may give you at least some
%tips, where to look for the error.
status=cell(Nmeshes,1);
for M=1:Nmeshes;
    [status{M}]=hbf_CheckMesh(bmeshes{M});
end

% ------------------------------------------------------------------------
% Now we should be ready for computation, unless something went wrong with
% the checks...

starttime=clock;

% 1. BEM double-layer operators for potentials: do this once per meshing...
D=hbf_BEMOperatorsPhi_LC(bmeshes);

% You can check, whether the operator matrices seem correct... this is
% heavily under construction, however!
hbf_CheckBEMOperators_LC(D);

%2. BEM operators for magnetic field due to volume currents: 
% do this, if meshing or MEG sensors change
DB=hbf_BEMOperatorsB_Linear(bmeshes,coils);

%3. Full transfer matrix
% Build BEM transfer matrix for potential using the isolated source
% approach, isolation set to surface 1.
Tphi_full=hbf_TM_Phi_LC_ISA2(D,ci,co,1);

%4. Transfer matrices for potential and volume component of the magnetic field
TBvol=hbf_TM_Bvol_Linear(DB,Tphi_full,ci,co);
%Extract/interpolate EEG transfer matrix for electrodes only
Tphi_elecs=hbf_InterpolateTfullToElectrodes(Tphi_full,bmeshes,elecs);

%5. Forward solutions = leadfield matrices
%Compute directed LFM --- source directions must have unit norm!
%You can use this function as well for computing the magnetic field due to
%any set of dipoles; in that case, just use the dipole moment vectors
%instead of unit-length source direction vectors.
LFM_Mdir=hbf_LFM_B_LC_dir(bmeshes,coils,TBvol,sourcepos,sourcedir);

%LFM for orthogonal unit dipoles (orientations x,y,z in world coordinates).
LFM_Mxyz=hbf_LFM_B_LC_xyz(bmeshes,coils,TBvol,sourcepos);
%EEG...
LFM_Edir=hbf_LFM_Phi_LC_dir(bmeshes,Tphi_elecs,sourcepos,sourcedir);
LFM_Exyz=hbf_LFM_Phi_LC_xyz(bmeshes,Tphi_elecs,sourcepos);

fprintf('Total time was %ds.\n',round(etime(clock,starttime)));

