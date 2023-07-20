function [e_proj,e_loctype,e_locinfo,distsurf]=hbf_ProjectPointsToMesh(mesh,epos,Nmin,verboseflag)
%HBF_PROJECTPOINTSTOMESH ... projects a set of points onto a mesh.
%
% [e_proj,e_loctype,e_info,distsurf]=HBF_PROJECTPOINTSTOMESH(mesh,epos,Nmin,verboseflag)
%
% This function projects the pointset 'epos' to nearest positions on the
% mesh 'mesh'. The projection is based on first finding 'Nmin' nearest
% triangle centroids and then finding the nearest points in these triangles
% or on their boundaries.
%
%   mesh: triangle mesh, hbf struct
%   epos: set of points to be projected, [N x 3]
%   Nmin (optional): the number of triangles to use in the search.
%   verboseflag (optional): give a nonzero value to see progress indicator
%
%   e_proj: the projected points, [N x 3]
%   e_loctype: location information of the projection point, [N x 1],
%       1 = in a triangle, 2 = on triangle-side, 3 = in node
%   e_locinfo: further location information depending on loctype, [N x 2], 
%       triangle index, two node indices, or one node index
%   distsurf: distance between the projected point and original point, [N x 1] 
% 
%   The default value for Nmin is 3. If the geometry is complicated or
%   the mesh has large variation in triangle size, higher value is
%   recommended. When projecting nearby electrodes to a typical scalp mesh, 
%   Nmin=1 is typically adequate, but some electrodes might be mis-projected by 1
%   mm or so. If the point is already on surface, Nmin==number of triangles
%   that belong to a vertex is recommended.
%
% v160330 Matti Stenroos

% 4 Nov 2015 (original: ProjectElectrodesToMesh_v4)

if nargin<3 || isempty(Nmin),
    Nmin=3;
end
if nargin<4 || isempty(verboseflag),
    verboseflag=0;
end

p=mesh.p;
e=mesh.e;
mp=TriangleMidpoints(p,e);
un=TriangleNormals(p,e);
noe=size(mesh.e,1);

Ne=size(epos,1);
e_proj=zeros(Ne,3);
e_loctype=zeros(Ne,1);
e_locinfo=zeros(Ne,2);
distsurf=zeros(Ne,1);

dispstep=round(Ne/10);
for E=1:Ne,
    if verboseflag && ~mod(E,dispstep);
        fprintf('%d%% ',round(100*E/Ne));
    end
    fp=epos(E,:);
    distsqr=sum((ones(noe,1)*fp-mp).^2,2);
    
    %Nmin nearest triangle midpoints
    %Matlab sort is slow, so let's do this ourselves...
    if Nmin==1,
        [foo,mintri]=min(distsqr);
    else
        meandistsqr=mean(distsqr);
        minval=zeros(Nmin,1);
        mintri=zeros(Nmin,1);
        for I=1:Nmin,
            [minval(I),mintri(I)]=min(distsqr);
            distsqr(mintri(I))=distsqr(mintri(I))+100*meandistsqr;%make sure that this wont appear again...
        end      
    end
    %triangle descriptions for Nmin nearest triangles
    emin=e(mintri,:);
    %evertices of Nmin nearest triangles
    p1=p(emin(:,1),:);
    p2=p(emin(:,2),:);
    p3=p(emin(:,3),:);
    %project testpoint to surfaces defined by Ntri nearest triangles
    psurf=ProjectPointToSurface(fp,un(mintri,:),p1);
    %are projection points inside triangles?
    intri=ProjectionInTriangle(p1,p2,p3,psurf);
    notintri=~intri;
    intri=find(intri);
    %distances from fieldpoint to projection points that are inside
    %triangles
    dvec=ones(size(intri,1),1)*fp-psurf(intri,:);
    d=sqrt(sum(dvec.*dvec,2));
    %nearest projection point that is inside the triangle
    [mind,minind]=min(d);
    
    %and fill the result fields...
    projpos=psurf(intri(minind),:);
    projinfo=mintri(intri(minind));
    projtype=1;
    projdist=mind;
    
    %now look at triangles, for which the projection point is outside..
    notintri=find(notintri);
    if any(notintri) %no projection inside triangle...
        %find all node pairs = sides of these triangles
        nodeind=emin(notintri,:);
        pairsu=[nodeind(:,[1 2]);nodeind(:,[1 3]);nodeind(:,[2 3])];
        p1=p(pairsu(:,1),:);
        p2=p(pairsu(:,2),:);
        %project fieldpoint to these sides & check is the projection point
        %on a triangle side or outside of it.
        [pedge,onedge]=ProjectPointToLine(fp,p1,p2);
        notonedge=~onedge;
    
        if any(onedge) %projection on sides
            onedge=find(onedge);
            dvec=ones(size(onedge,1),1)*fp-pedge(onedge,:);
            d=sqrt(sum(dvec.*dvec,2));
            [mind,minind]=min(d);
         
            %if the nearest point on the side was closer to fieldpoint
            %than the nearest projection point inside a triangle...
            if isempty(projdist) || mind<projdist ,
                projpos=pedge(onedge(minind),:);
                projinfo=pairsu(onedge(minind),:);
                projtype=2;
                projdist=mind;
            end   
        end
        
        %Now look at those triangle-sides, where the projection point was
        %not inside the corresponding triangles... meaning, the nearest
        %points are in nodes...
        notonedge=find(notonedge);
        if any(notonedge),
            %which nodes were "lonely" = no projection "inside" a line
            solonodes=pairsu(notonedge,:);
            solonodes=unique(solonodes(:));
            dvec=ones(size(solonodes,1),1)*fp-p(solonodes,:);
            d=sqrt(sum(dvec.*dvec,2));
            [mind,minind]=min(d);
            %if it is closer to the surface than any of the points found so
            %far, so be it...
            if isempty(projdist) || mind<projdist
                projinfo=solonodes(minind);
                projpos=p(projinfo,:);
                projtype=3;
                projdist=mind;               
            end
        end
    end
    if isempty(projpos)
        disp(E)
    end
    e_proj(E,:)=projpos;
    e_loctype(E,:)=projtype;
    e_locinfo(E,:)=projinfo;
    distsurf(E,:)=projdist;   
end

%hasn't gone wrong so far, but...
if isempty(projpos)
        fprintf('Something went wrong... Projection not found for point %d\n',E);
end
function unormals=TriangleNormals(nodes,elements)
% function unormals=TriangleNormals(nodes,elements)
%unormals = unit normals
p1=nodes(elements(:,1),:);
p2=nodes(elements(:,2),:);
p3=nodes(elements(:,3),:);
normals=cross(p2-p1,p3-p1);
areas=sqrt(sum(normals.*normals,2))/2;
unormals=normals./(2*areas*[1 1 1]);

function midpoints=TriangleMidpoints(nodes,elements)
% function midpoints=TriangleMidpoints(nodes,elements)
% Calculates midpoints of the mesh triangles.
p1=nodes(elements(:,1),:);
p2=nodes(elements(:,2),:);
p3=nodes(elements(:,3),:);
midpoints=(p1+p2+p3)/3;
    
function res=ProjectionInTriangle(p1,p2,p3,testp)
% function res=ProjectionInTriangle(p1,p2,p3,testp)
% p1, p2, p3 = points 1, 2, 3 of triangles
% p1 = [N x 3] etc
% testp = [1 x 3] or [N x 3]
test1=PointsOnTheSameSideOfTheLine(p1,p2,testp,p3);
test2=PointsOnTheSameSideOfTheLine(p2,p3,testp,p1);
test3=PointsOnTheSameSideOfTheLine(p1,p3,testp,p2);
res=test1 & test2 & test3;

function [pproj,between_p1andp2]=ProjectPointToLine(fp,p1,p2)
% function [pproj,between_p1andp2]=ProjectPointToLine(fp,p1,p2)
% fp = Point to be projected, [1 x 3]
% p1, p2 = points defining a line; [1 x 3] or (for many lines) [N x 3]
% pproj = projected points
% between_p1andp2 = is projection between p1 and p2
% 4 Nov 2010
Np=size(p1,1);
k=p2-p1;
klen=sqrt(sum(k.*k,2));
k0=k./(klen*[1 1 1]);
r=ones(Np,1)*fp-p1;
proj=sum(k0.*r,2);
pproj=p1+proj*[1 1 1].*k0;
between_p1andp2=proj>0 & proj<klen;

function res=ProjectPointToSurface(rf,n,r0)
%function res=ProjectPointToSurface(rf,n,r0)
%res is a projection of rf on the surface defined by unit normal n and point r0
%rf: one point, 1 x 3
%n:  N x 3
%r0: N x 3

rfm=ones(size(n,1),1)*rf;
tau = dots(n,(rfm-r0));
res=rfm-[tau tau tau].*n;

function res=PointsOnTheSameSideOfTheLine(p1,p2,p,q)
% function res=PointsOnTheSameSideOfTheLine(p1,p2,p,q)
% p1, p2 define the lines
% p1, p2 = [N x 3]
% p and q are tested points; either [1 x 3] or [N x 3];
N=size(p1,1);
if size(p,1)~=N,
    p=ones(N,1)*p;
end
if size(q,1)~=N,
    q=ones(N,1)*q;
end
cp1 = cross(p2-p1, p-p1);
cp2 = cross(p2-p1, q-p1);
proj=sum(cp1.*cp2,2);
res=proj>=0;

function dot=dots(R1,R2)
dot=R1(:,1).*R2(:,1)+R1(:,2).*R2(:,2)+R1(:,3).*R2(:,3);
function R=cross(R1,R2)
R(:,1)=R1(:,2).*R2(:,3)-R1(:,3).*R2(:,2);
R(:,2)=R1(:,3).*R2(:,1)-R1(:,1).*R2(:,3);
R(:,3)=R1(:,1).*R2(:,2)-R1(:,2).*R2(:,1);
