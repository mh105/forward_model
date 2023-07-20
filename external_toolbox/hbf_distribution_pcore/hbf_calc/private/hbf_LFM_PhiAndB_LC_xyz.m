%  HBF_LFM_PHIANDB_LC_XYZ builds electric and magnetic lead field matrices
%    based on xyz-oriented unit-current dipoles.
% 
%  [LFMphi,LFMm]=HBF_LFM_PHIANDB_LC_XYZ(meshes,coils,Tphi,TB,spos,flag_averef)
%    meshes: BEM geometry, cell array of hbf structs
%    coils:  coil description, hbf struct
%    TB:     TB matrix built with the hbf BEM solver
%    Tphi:   Tphi matrix built with the hbf BEM solver
%    spos:   source positions, [M x 3]
%    flag_averef (optional, default value 1): give 0, if you do not want to
%            use average reference in LFMphi
% 
%    LFMphi: electric lead field matrix, [Number of electrodes x M]
%        [lphi_1x lphi_1y lphi_1z ... lphi_Mx lphi_My lphi_Mz]
%    LFMb: magnetic lead field matrix, [Number of coils x M]
%        [lb_1x lb_1y lb_1z ... lb_Mx lb_My lb_Mz]
% 
% 
%  This function assumes average reference by default. If 'flag_averef' is
%  selected 0, the potential is computed againts the reference chosen when
%  building the BEM matrix.
% 
%  v160303 Matti Stenroos
%