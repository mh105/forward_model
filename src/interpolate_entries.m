function [ G, xyzG ] = interpolate_entries(G, source, bad_entry, xyzG)
%Function used for manual identification of outlier entries in lead field
%matrix and generate interpolation of these bad entries using near-by
%source points

if nargin < 4
    xyzG = [];
end

% compute a distance matrix between sources
src_distance = zeros(size(source,1));
for ii = 1:size(src_distance,1)
    src_distance(:,ii) = vecnorm(source - source(ii,:), 2, 2);
end

% interpolate using nearby source points < 10mm
bad_entry_idx = find(bad_entry==1);
update_values = zeros(size(bad_entry_idx));
sz = size(G);
xyzsz = size(xyzG);
bad_entry_xyz_idx = zeros(length(bad_entry_idx)*3, 1);
xyz_update_values = zeros(length(bad_entry_idx)*3, 1);
total_size = 0;
for ii = 1:length(bad_entry_idx)
    [row,col] = ind2sub(sz,bad_entry_idx(ii));
    neighbor_sources = find(src_distance(col, :) < 0.01);
    neighbor_sources(neighbor_sources==col) = []; % exclude oneself
    update_values(ii) = mean(G(row, neighbor_sources));
    if ~isempty(xyzG) % need to fix free-orientation LFM as well
        [~, x,y,z] = fix2freeindexing(neighbor_sources);
        first_idx = (ii-1)*3+1;
        xyz_update_values(first_idx) = mean(xyzG(row, x), 2);
        xyz_update_values(first_idx+1) = mean(xyzG(row, y), 2);
        xyz_update_values(first_idx+2) = mean(xyzG(row, z), 2);
        xyz_update_idx = sub2ind(xyzsz, repmat(row, 1, 3), fix2freeindexing(col));
        total_size = total_size + length(xyz_update_idx);
        bad_entry_xyz_idx(first_idx:first_idx+2) = xyz_update_idx;
    end
end
G(bad_entry_idx) = update_values;
if ~isempty(xyzG)
    assert(total_size == length(bad_entry_idx)*3, 'total size of entries fixed in xyzG is incorrect!')
    xyzG(bad_entry_xyz_idx) = xyz_update_values;
end

end

