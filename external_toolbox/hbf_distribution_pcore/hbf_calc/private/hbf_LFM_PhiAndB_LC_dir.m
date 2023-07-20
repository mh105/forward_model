%  HBF_LFM_PHIANDB_LC_DIR builds electric and magnetic lead field matrices
%    based on directed current dipoles.
% 
%  [LFMphi,LFMm]=HBF_LFM_PHIANDB_LC_DIR(meshes,coils,Tphi,TB,spos,sdir,flag_averef)
%    meshes: BEM geometry, cell array of hbf structs
%    coils:  coil description, hbf struct
%    TB:     TB matrix built with the hbf BEM solver
%    Tphi:   Tphi matrix built with the hbf BEM solver
%    spos:   source positions, [M x 3]
%    sdir:   source orientations (unit-length), [M x 3]
%    flag_averef (optional, default value 1): give 0, if you do not want to
%            use average reference
% 
%    LFMphi: electric lead field matrix, [Number of electrodes x M]
%        [lphi_1 ... lphi_M]
%    LFMb: magnetic lead field matrix, [Number of coils x M]
%        [lb_1 ... lb_M]
% 
% 
%  You can also compute phi and B due to any set of directed dipoles by
%  giving the dipole moments (with amplitude) in the 'sdir' argument.
% 
%  This function assumes average reference by default. If 'flag_averef' is
%  given 0, the potential is computed againts the reference chosen when
%  building the BEM matrix (default: mean potential over scalp nodes is zero).
% 
%  v160229 Matti Stenroos
%