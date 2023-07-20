function varargout=hbf_LFM_LC(varargin)
%HBF_LFM_LC is a wrapper for building lead field matrices
%varargout=hbf_LFM_LC(varargin)
%
%[LFMphi,LFMb]=HBF_LFM_LC(bmeshes,coils,Tphi,TB,spos,sdir,flag_averef)
%[LFMphi,LFMb]=HBF_LFM_LC(bmeshes,coils,Tphi,TB,spos,flag_averef)
%LFMphi=HBF_LFM_LC(bmeshes,Tphi,spos,sdir,flag_averef)
%LFMphi=HBF_LFM_LC(bmeshes,Tphi,spos,flag_averef)
%LFMb=HBF_LFM_LC(bmeshes,coils,TB,spos,sdir)
%LFMb=HBF_LFM_LC(bmeshes,coils,TB,spos)
%
%   meshes: BEM geometry, cell array of hbf structs
%   coils:  coil description, hbf struct
%   TBvol:     TB matrix built with the hbf BEM solver
%   Tphi:   Tphi matrix built with the hbf BEM solver
%   spos:   source positions, [M x 3]
%   sdir:   source orientations (unit-length), [M x 3]
%   flag_averef (optional, default value 1): give 0, if you do not want to
%           use average reference for LFMphi
%
%   LFMphi: electric lead field matrix, [Number of electrodes x M]
%       [lphi_1 ... lphi_M]
%   LFMb: magnetic lead field matrix, [Number of coils x M]
%       [lb_1 ... lb_M]
%
%
% You can also compute phi and B due to any set of directed dipoles by
% giving the dipole moments (with amplitude) in the 'sdir' argument.
%
% This function assumes average reference for potential by default. If
% 'flag_averef' is given 0, the potential is computed againts the reference
% chosen when building the BEM matrix (default: mean potential over scalp
% nodes is zero). 
% Note: If you are not computing LFMphi, do not give "flag_averef"!
%
% v160303 Matti Stenroos

%parsing what to do ------------------------------------------------------
flag_phi=0;
flag_dir=0;
Nin=length(varargin);

%if the second argument is struct, it must be coil definition -> compute B.
flag_b = isstruct(varargin{2});
%if there is no B, there must be phi then...
if ~flag_b,flag_phi=1;end
%if the last argument does not have 3 dimensions, it must be "flag_averef"
%and there must be phi as well
flag_refgiven = max(size(varargin{Nin}))<3;
if flag_refgiven, flag_phi = 1;end
    
if Nin==7,
    flag_dir=1;
elseif Nin==6, %phi and b but no directions or flag_refgiven;
    flag_phi=1;
    flag_dir=~flag_refgiven;
elseif Nin==5, %phi or b or both
    if flag_b && size(varargin{3},2)==size(varargin{4},2), %both
        flag_phi=1;
    end
    flag_dir = (flag_b && ~flag_phi) || (flag_phi && ~flag_b);
elseif Nin==4, %either phi or b
    flag_phi=~flag_b;
    flag_dir=~(flag_b||flag_refgiven);
elseif Nin==3,
    flag_phi=1;flag_dir=0;
end
% ------
% ...and computing:
if flag_phi && flag_b,
    if flag_dir,
        [varargout{1},varargout{2}]=hbf_LFM_PhiAndB_LC_dir(varargin{:});
    else
        [varargout{1},varargout{2}]=hbf_LFM_PhiAndB_LC_xyz(varargin{:});
    end
elseif flag_phi,
    if flag_dir,
        varargout{1}=hbf_LFM_Phi_LC_dir(varargin{:});
    else
        varargout{1}=hbf_LFM_Phi_LC_xyz(varargin{:});
    end
elseif flag_b,
    if flag_dir,
        varargout{1}=hbf_LFM_B_LC_dir(varargin{:});
    else
        varargout{1}=hbf_LFM_B_LC_xyz(varargin{:});
    end
else
    error('hbf_LFM_LC: Could not parse arguments');
end
