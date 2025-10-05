function [ G, bad_source, bad_entry, xyzG ] = detect_badsource(G, source, interpl, xyzG, detectentry)
%This function is used to detect bad sources and bad G entries. We also add
%the option to interpolate from near by sources.
if nargin < 2
    source = [];
    interpl = false;
    xyzG = [];
    detectentry = true;
elseif nargin < 3
    interpl = false;
    xyzG = [];
    detectentry = true;
elseif nargin < 4
    xyzG = [];
    detectentry = true;
elseif nargin < 5
    detectentry = true;
end


% exclude if the median of a source is abnormal
src_median = median(G,1);
bad_source = isoutlier(src_median);
bad_src_idx = find(bad_source==1);
G(:, bad_src_idx) = 0;
if ~isempty(xyzG) % if we need to fix the free-orientation LFM as well
    bad_src_free_column = fix2freeindexing(bad_src_idx);
    xyzG(:, bad_src_free_column) = 0;
end


% detect if an entry of G is abnormal
if detectentry
    % exclude bad entries in G if they lie outside distribution
    bad_entry = isoutlier(G(:),  'grubbs');
    bad_entry_idx = find(bad_entry==1);
    G(bad_entry_idx) = 0;
    if ~isempty(xyzG) % if we need to fix the free-orientation LFM as well
        [trow,tcol] = ind2sub(size(G),bad_entry_idx);
        unique_rows = unique(trow);
        total_size = 0;
        for ii = 1:length(unique_rows)
            bad_entry_free_column = fix2freeindexing(tcol(trow==unique_rows(ii))');
            testing = xyzG(unique_rows(ii), bad_entry_free_column);
            total_size = total_size + length(testing(:));
            xyzG(unique_rows(ii), bad_entry_free_column) = 0;
        end
        assert(total_size == length(bad_entry_idx)*3, 'total size of entries deleted in xyzG is incorrect!')
    end
else
    bad_entry = nan;
end


% interpolate the deleted entries using near-by sources within 10mm radius
if interpl
    % compute a distance matrix between sources
    src_distance = zeros(size(source,1));
    for ii = 1:size(src_distance,1)
        src_distance(:,ii) = vecnorm(source - source(ii,:), 2, 2);
    end
    
    % fix bad sources
    update_values = zeros(size(G,1), length(bad_src_idx));
    xyz_update_values = zeros(size(G,1), length(bad_src_idx)*3);
    for ii = 1:length(bad_src_idx)

        threshold = 0.01; % start with 10mm
        neighbor_sources = [];
        
        % adaptively increase threshold until at least one neighbor is found
        while isempty(neighbor_sources)
            neighbor_sources = find(src_distance(bad_src_idx(ii), :) < threshold);
            
            % exclude neighbors that are also bad sources
            neighbor_sources(ismember(neighbor_sources, bad_src_idx)) = [];
            
            % if empty or all NaN, increase threshold
            if isempty(neighbor_sources) 
                threshold = threshold + 0.001; % add 1mm
            end
        end
        
        update_values(:,ii) = mean(G(:, neighbor_sources), 2);
        if ~isempty(xyzG) % need to fix free-orientation LFM as well
            [~, x,y,z] = fix2freeindexing(neighbor_sources);
            first_column_idx = (ii-1)*3+1;
            xyz_update_values(:, first_column_idx) = mean(xyzG(:, x), 2);
            xyz_update_values(:, first_column_idx+1) = mean(xyzG(:, y), 2);
            xyz_update_values(:, first_column_idx+2) = mean(xyzG(:, z), 2);
        end
    end
 
    G(:, bad_src_idx) = update_values;
    if ~isempty(xyzG)
        xyzG(:, bad_src_free_column) = xyz_update_values;
    end
    
    if detectentry
        % fix bad entries
        update_values = zeros(size(bad_entry_idx));
        sz = size(G);
        xyzsz = size(xyzG);
        bad_entry_xyz_idx = zeros(length(bad_entry_idx)*3, 1);
        xyz_update_values = zeros(length(bad_entry_idx)*3, 1);
        total_size = 0;

        % find column indices of all bad entries
        [~, col_all] = ind2sub(sz, bad_entry_idx);

        for ii = 1:length(bad_entry_idx)
            [row,col] = ind2sub(sz,bad_entry_idx(ii));
            threshold = 0.01; % start with 10mm
            neighbor_sources = [];
            
            % adaptively increase threshold until at least one neighbor is found
            while isempty(neighbor_sources)
                neighbor_sources = find(src_distance(col, :) < threshold);
                
                % exclude neighbors that are also bad sources
                neighbor_sources(ismember(neighbor_sources, col_all)) = [];
                
                % if empty or all NaN, increase threshold
                if isempty(neighbor_sources)
                    threshold = threshold + 0.001; % add 1mm
                end
            end
            
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
end

end