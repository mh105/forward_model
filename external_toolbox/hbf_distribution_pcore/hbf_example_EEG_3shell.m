% Example on the use of Helsinki BEM Framework kernel:
% 3-shell EEG model.
% 
% The example uses a three-shell head model made by me at MRC Cognition and
% Brain Sciences Unit, Cambridge, UK. The model was made using FLASH MR
% images, FreeSurfer, and MNE toolbox. If you wish to use this model for
% something that you wish to share or publish, please contact me.
%
% (c) Matti Stenroos
% 2 Mar 2016
clear

addpath ./hbf_calc
addpath ./hbf_mesh
%Load / set geometries
load('~/bioem/exampledata/cbusample/samplehead_cbu.mat',...
    'bmeshes', 'eegsens', 'sources');%bmeshes, eegsens, sources
%   the meshes are nested; now make sure that the meshes are in correct order
bmeshes=hbf_SortNestedMeshes(bmeshes);
sourcepos=sources{1}.p; %source locations
sourcedir=sources{1}.nn; %source directions --- for LFM, each vector must be unit length!
elecpos=eegsens.p;
%   project electrodes to scalp and make "elecs" struct
elecs=hbf_ProjectElectrodesToScalp(elecpos,bmeshes);
%Set conductivities
ci=[1 1/50 1]*.33; %conductivity inside each surface --- remember the order!
co=[ci(2:3) 0]; %conductivity outside each surface

%--------------------------------------------------------------------------
% This part is not needed, if meshing pipeline is established and correctly
% configured
%If new meshes, perform some checks...
Nmeshes=length(bmeshes);
%   check/correct orientation
success=zeros(Nmeshes,1);
for M=1:Nmeshes;
    [bmeshes{M},success(M)]=hbf_CorrectTriangleOrientation(bmeshes{M});
end
%   check, is the mesh valid for BEM
status=cell(Nmeshes,1);
for M=1:Nmeshes;
    status{M}=hbf_CheckMesh(bmeshes{M});
end
%--------------------------------------------------------------------------
%Now let's get to business!
starttime=clock;
%   BEM geometry matrices
D=hbf_BEMOperatorsPhi_LC(bmeshes);
%   full transfer matrix
Tphi_full=hbf_TM_Phi_LC_ISA2(D,ci,co,1);
%   extract transfer matrix for electrodes
Tphi_elecs=hbf_InterpolateTfullToElectrodes(Tphi_full,bmeshes,elecs);

%   forward solutions for directed dipoles
LFMphi_dir=hbf_LFM_LC(bmeshes,Tphi_elecs,sourcepos,sourcedir);
%   forward solutions for xyz-oriented triplets
LFMphi_xyz=hbf_LFM_LC(bmeshes,Tphi_elecs,sourcepos);

fprintf('Total time was %ds.\n',round(etime(clock,starttime)));

