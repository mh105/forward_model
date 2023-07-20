%  HBF_LFM_PHI_LC_XYZ builds electric lead field matrix based on xyz-oriented
%  unit-current dipoles
% 
%  LFM=HBF_LFM_PHI_LC_XYZ(meshes,Tphi,spos,flag_averef)
%    meshes: BEM geometry, cell array of hbf structs
%    Tphi:   Tphi matrix built with the hbf BEM solver
%    spos:   source positions, [M x 3]
%    flag_averef (optional, default value 1): give 0, if you do not want to
%            use average reference
% 
%    LFM:   lead field matrix, [Number of electrodes x 3M]
%        [l_1x l_1y l1_z ... l_Mx l_My l_Mz]
% 
%  This function assumes average reference by default. If 'flag_averef' is
%  selected 0, the potential is computed againts the reference chosen when
%  building the BEM matrix.
% 
%  v160229 Matti Stenroos
%