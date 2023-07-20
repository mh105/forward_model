% Example on the use of Helsinki BEM Framework kernel:
% 4-shell MEG+EEG model.
% 
% The example uses meshes from Head 1 of (Stenroos and Nummenmaa, 2016),
% remeshed by me based on sample data of SimNIBS toolbox. If you wish to use the
% meshes or model for something, please contact me first.
%
% (c) Matti Stenroos 2017
clear
%add BEM framework paths
addpath ./hbf_calc
addpath ./hbf_mesh
% load example data for meshes, coils, and
% sources.
loaddir='~/wrk/bemlibdata/snsample/';
load(strcat(loaddir,'geometry/bmeshes-m40-50-50-60-70.mat'));%bmeshes
bmeshes=bmeshes.meshes;
load(strcat(loaddir,'geometry/coils-cmeg306.mat'));%coils
load(strcat(loaddir,'geometry/elec256_s70.mat'));%electrodes
load(strcat(loaddir,'geometry/sources-sfs35'));%sources
sourcepos=sources.smeshes{1}.p;sourcedir=sources.smeshes{1}.nn;
Nmeshes=length(bmeshes);
success=zeros(Nmeshes,1);
for M=1:Nmeshes
    [bmeshes{M},success(M)]=hbf_CorrectTriangleOrientation(bmeshes{M});
end

% status=cell(Nmeshes,1);
% for M=1:Nmeshes;
%     [status{M}]=hbf_CheckMesh(bmeshes{M});
% end

%Set EEG electrodes
elecpos=elecs.p;
%Project electrodes to scalp and build interpolation operator from model
%nodes to projected electrodes.
elecs=hbf_ProjectElectrodesToScalp(elecpos,bmeshes);

% ------------------------------------------------------------------------
% Hey, ho, let's go!

starttime=clock;
% Build double-layer matrices for the BEM
D=hbf_BEMOperatorsPhi_LC(bmeshes);
% Build matrices for integrating the magnetic field due to volume currents 
DB=hbf_BEMOperatorsB_Linear(bmeshes,coils);

cratio=50;
c_brain=.33;c_csf=1.79;c_skull=c_brain/cratio;c_skin=.33;
ci=[c_brain,c_brain,c_csf,c_skull,c_skin];%conductivity inside each surface --- remember the order!
co=[c_csf,c_csf,c_skull,c_skin,0];%conductivity outside each surface --- remember the order!

% Build BEM transfer matrix for potential using the isolated source
% approach, isolation set to surface 3 = inner skull
Tphi_full=hbf_TM_Phi_LC_ISA2(D,ci,co,3);
%Extract/interpolate EEG transfer matrix for electrodes only
Tphi_elecs=hbf_InterpolateTfullToElectrodes(Tphi_full,bmeshes,elecs);
% ...and make a transfer matrix for the magnetic field generated by the
% volume currents.
TBvol=hbf_TM_Bvol_Linear(DB,Tphi_full,ci,co);


% compute lead field matrices
LFMstarttime=clock;
LFM_Edir=hbf_LFM_Phi_LC_dir(bmeshes,Tphi_elecs,sourcepos,sourcedir);
LFM_Mdir=hbf_LFM_B_LC_dir(bmeshes,coils,TBvol,sourcepos,sourcedir);
fprintf('Time for LFM computation was %ds.\n',round(etime(clock,LFMstarttime)));
%%
LFMstarttime2=clock;
[LFM_Edir2,LFM2_Mdir]=hbf_LFM_PhiAndB_LC_dir(bmeshes,coils,Tphi_elecs,TBvol,sourcepos,sourcedir);
fprintf('Time for LFM computation was %fs.\n',etime(clock,LFMstarttime2));


fprintf('Total time was %ds.\n',round(etime(clock,starttime)));
