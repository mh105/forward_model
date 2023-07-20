function [ H,hlink ] = hbf2MNE_compareG(G1, G2, label1, label2, sources, channel, dig, weightogether)
%This function is to visually compare the G matrix entries between two
%leadfield matrices.

if nargin < 8
    weightogether = true;
end

H = figure;
maxnum_src = min(size(G1,2), size(G2,2));
all_weight = [G1(channel,1:maxnum_src).^2, G2(channel,1:maxnum_src).^2];

% G matrix No. 1
% grab weights from the channel electrode
weight = G1(channel,1:maxnum_src).^2;
if weightogether
    size_weight = ceil((weight-min(all_weight))/(max(all_weight)-min(all_weight))*300)+1;
else
    size_weight = ceil((weight-min(weight))/(max(weight)-min(weight))*300)+1;
end
% plot the sources in 3D space
ax1 = subplot(1,2,1);
hold on 
axis equal
title([label1, ' G matrix'], 'FontSize', 16)
scatter3(sources(1:maxnum_src,1), sources(1:maxnum_src,2), sources(1:maxnum_src,3), size_weight, weight, 'filled')
scatter3(dig(channel,1), dig(channel,2), dig(channel,3), 300, [1,0,0], 's', 'filled')
colorbar
cspect = caxis;
xlimit = xlim; ylimit = ylim; zlimit = zlim;

% G matrix No. 2
% grab weights from the channel electrode
weight = G2(channel,1:maxnum_src).^2;
if weightogether
    size_weight = ceil((weight-min(all_weight))/(max(all_weight)-min(all_weight))*300)+1;
else
    size_weight = ceil((weight-min(weight))/(max(weight)-min(weight))*300)+1;
end
% plot the sources in 3D space
ax2 = subplot(1,2,2);
hold on 
axis equal
title([label2, ' G matrix'], 'FontSize', 16)
scatter3(sources(1:maxnum_src,1), sources(1:maxnum_src,2), sources(1:maxnum_src,3), size_weight, weight, 'filled')
scatter3(dig(channel,1), dig(channel,2), dig(channel,3), 300, [1,0,0], 's', 'filled')
colorbar
if weightogether; caxis(cspect); end
xlim(xlimit); ylim(ylimit); zlim(zlimit);

hlink = linkprop([ax1,ax2],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
rotate3d on

end

