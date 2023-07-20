function [ntoe_full,nodeinds,ntoe_dense]=hbf_LinearInterpolationMatrix(mesh,e_proj,e_loctype,e_locinfo)
% HBF_LINEARINTERPOLATIONMATRIX makes a matrix for interpolating a function
% in electrode positions
% function [ntoe_full,nodeinds,ntoe_dense]=LinearInterpolationMatrix(mesh,e_proj,e_loctype,e_locinfo)
% ntoe_full: [N_elec x mesh.nop]
% nodeinds: [N x 1]; nodes that are actually needed for the electrodes
% ntoe_dense: [N_elec x N]; zero-columns removed
%
% use of ntoe_dense:
% Ttemp=Tmesh(nodeinds,:);
% Telec=ntoe_dense*Ttemp;
%
% v161220 Matti Stenroos

N_elec=size(e_proj,1);
ntoe_full=zeros(N_elec,size(mesh.p,1));
for E=1:N_elec,
    switch e_loctype(E),
        case 1
            tri=e_locinfo(E,1);
            tricoefs=hbf_LinInterpInTriangle(e_proj(E,:),mesh.p(mesh.e(tri,:),:));
            ntoe_full(E,mesh.e(tri,:))=tricoefs;
        case 2
            nodes=e_locinfo(E,:);
            linecoefs=hbf_LinInterpOnLine(e_proj(E,:),mesh.p(nodes,:));
            ntoe_full(E,nodes)=linecoefs;
        case 3
            ntoe_full(E,e_locinfo(E,1))=1;
    end
end
%check, which nodes are actually needed
tris=e_locinfo(e_loctype==1,1);
trinodes=mesh.e(tris,:);
sidenodes=e_locinfo(e_loctype==2,:);
soloists=e_locinfo(e_loctype==3,1);
nodeinds=unique([trinodes(:);sidenodes(:);soloists(:)]);
ntoe_dense=ntoe_full(:,nodeinds);

function coefs=hbf_LinInterpInTriangle(P, varargin)
% HBF_LININTERPINTRIANGLE computes interpolation coefficiens for a point in triangle
% function coefs=LinInterpInTriangle(P,trip);
% A, B, C: corners of the triangle, each [1 x 3];
% trip: [A;B;C]
% P: point inside the triangle (in the plane)
% v161220 Matti Stenroos

% based on LinearInterpolationMatrix.m, a not-too-fast but working old
% function...
if length(varargin)==3,
    A=varargin{1};
    B=varargin{2};
    C=varargin{3};
else
    A=varargin{1}(1,:);
    B=varargin{1}(2,:);
    C=varargin{1}(3,:);
end

AB=B-A;
nAB=sqrt(sum(AB.^2,2));
eAB=AB/nAB;

AC=C-A;
nAC=sqrt(sum(AC.^2,2));
eAC=AC/nAC;

BC=C-B;
nBC=sqrt(sum(BC.^2,2));
eBC=BC/nBC;

%1 2 3 1 2
%A B C A B
AAp=AB+eBC*(dotp(-AB,eBC));
BBp=BC-eAC*(dotp(-BC,-eAC));
CCp=-AC+eAB*(dotp(AC,eAB));

AAn=sqrt(sum(AAp.^2,2));
BBn=sqrt(sum(BBp.^2,2));
CCn=sqrt(sum(CCp.^2,2));

eAAp=AAp/AAn;
eBBp=BBp/BBn;
eCCp=CCp/CCn;

a=dotp(eAAp,P-A)/AAn;
b=dotp(eBBp,P-B)/BBn;
c=dotp(eCCp,P-C)/CCn;

coefs=[1-a 1-b 1-c];

function coefs=hbf_LinInterpOnLine(P, varargin)
% HBF_LININTERPONLINE computes interpolation coefs for a point on a line.
% function coefs=hbf_LinInterpOnLine(P,A,B);
% function coefs=nbf_LinInterpInTriangle(P,linep);
% A, B: endpoints of line, each [1 x 3];
% linep: [A;B]
% P: point on the line (between A and B)
% v161220 Matti Stenroos

% based on LinInterpOnLine.m, an old
% function...
    
if length(varargin)==2,
    A=varargin{1};
    B=varargin{2};
else
    A=varargin{1}(1,:);
    B=varargin{1}(2,:);
end
AB=B-A;
nAB=sqrt(sum(AB.^2));
AP=P-A;
nAP=sqrt(sum(AP.^2));
coefs=[1-nAP/nAB nAP/nAB];

function res=dotp(A,B)
res=A(:,1)*B(:,1)+A(:,2)*B(:,2)+A(:,3)*B(:,3);
