function [ G, xyzG ] = interpolate_sources(G, source, srcidx, xyzG)
%Function used for manual identification of outlier source points and
%generate interpolation of these source points using near-by points 

if nargin < 4
    xyzG = [];
end

% compute a distance matrix between sources
src_distance = zeros(size(source,1));
for ii = 1:size(src_distance,1)
    src_distance(:,ii) = vecnorm(source - source(ii,:), 2, 2);
end

% interpolate using nearby source points < 10mm
update_values = zeros(size(G,1), length(srcidx));
xyz_update_values = zeros(size(G,1), length(srcidx)*3);
for ii = 1:length(srcidx)
    neighbor_sources = find(src_distance(srcidx(ii), :) < 0.01);
    neighbor_sources(neighbor_sources==srcidx(ii)) = []; % exclude oneself
    update_values(:,ii) = mean(G(:, neighbor_sources), 2);
    if ~isempty(xyzG) % need to fix free-orientation LFM as well
        [~, x,y,z] = fix2freeindexing(neighbor_sources);
        first_column_idx = (ii-1)*3+1;
        xyz_update_values(:, first_column_idx) = mean(xyzG(:, x), 2);
        xyz_update_values(:, first_column_idx+1) = mean(xyzG(:, y), 2);
        xyz_update_values(:, first_column_idx+2) = mean(xyzG(:, z), 2);
    end
end
G(:, srcidx) = update_values;
if ~isempty(xyzG)
    xyzG(:, fix2freeindexing(srcidx)) = xyz_update_values;
end

end

