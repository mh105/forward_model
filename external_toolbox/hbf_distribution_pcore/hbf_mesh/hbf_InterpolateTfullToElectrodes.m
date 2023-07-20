function Tphi_elecs=hbf_InterpolateTfullToElectrodes(Tfull,bmeshes,elecs)
% function Tphi_elecs=hbf_InterpolateTfullToElectrodes(Tfull,bmeshes,elecs);
[startinds,endinds]=NodeIndices(bmeshes);
Tphi_elecs=elecs.NtoE*Tfull(startinds(end):endinds(end),:);

function [startinds,endinds,Nop]=NodeIndices(meshes)
Nsurf=length(meshes);
startinds=zeros(Nsurf,1);
endinds=zeros(Nsurf,1);
Nop=0;
for I=1:Nsurf,
    startinds(I)=Nop+1;
    endinds(I)=startinds(I)+size(meshes{I}.p,1)-1;
    Nop=endinds(I);
end