function [ TR, trans_v, final_cost ] = register_electrode_aprx(e_init, e_dest, init_rot_Z, iterations)
%Function used to register same electrodes under different coordinate
%systems together - trying to register e1 onto e2. Outputs rotation matrix
%and translation vector separately

% Use this function with caution...

e_init_archive = e_init;

figure
hold on
scatter3(e_init(:,1), e_init(:,2), e_init(:,3), 50, 'r', 'filled')
scatter3(e_dest(:,1), e_dest(:,2), e_dest(:,3), 50, 'b', 'filled')
axis equal
rotate3d on

finalZ = zeros(1,iterations);
finalY = zeros(1,iterations);
finalX = zeros(1,iterations);

for ii = 1:iterations
    
    % first let's search the right Z parameter
    step_size = 0.1; %deg
    X_search_range = deg2rad(0:0);
    Y_search_range = deg2rad(0:0);
    if ii == 1
        Z_search_range = deg2rad(rad2deg(init_rot_Z)-10:step_size/10:rad2deg(init_rot_Z)+10);
    else
        Z_search_range = deg2rad(-1:step_size/100:1);
    end
    
    cost_tally = zeros(length(X_search_range), length(Y_search_range), length(Z_search_range));
    
    for ix = 1:length(X_search_range)
        for iy = 1:length(Y_search_range)
            for iz = 1:length(Z_search_range)
                x = X_search_range(ix);
                y = Y_search_range(iy);
                z = Z_search_range(iz);
                Rx = [1 0 0; 0 cos(x) -sin(x); 0 sin(x) cos(x)];
                Ry = [cos(y) 0 sin(y); 0 1 0; -sin(y) 0 cos(y)];
                Rz = [cos(z) -sin(z) 0; sin(z) cos(z) 0; 0 0 1];
                
                e_init_rot = e_init * Rz * Rx * Ry;
                translation_vec = mean(e_dest) - mean(e_init_rot);
                e_init_trans = e_init_rot + translation_vec;
                %
                %             figure
                %             hold on
                %             scatter3(e_init_trans(1,1), e_init_trans(1,2), e_init_trans(1,3), 100, 'k', 'filled')
                %             scatter3(e_dest(1,1), e_dest(1,2), e_dest(1,3), 70, 'md', 'filled')
                %             scatter3(e_init_trans(:,1), e_init_trans(:,2), e_init_trans(:,3), 50, 'r', 'filled')
                %             scatter3(e_dest(:,1), e_dest(:,2), e_dest(:,3), 50, 'b', 'filled')
                %             axis equal
                %             rotate3d on
                cost_tally(ix,iy,iz) = sum(vecnorm(e_init_trans - e_dest, 2, 2));
            end
        end
    end
    
    [val, lidx] = min(cost_tally(:)); %#ok<*ASGLU>
    [ix,iy,iz] = ind2sub(size(cost_tally),lidx);
    
    x = X_search_range(ix);
    y = Y_search_range(iy);
    z = Z_search_range(iz);
    Rx = [1 0 0; 0 cos(x) -sin(x); 0 sin(x) cos(x)];
    Ry = [cos(y) 0 sin(y); 0 1 0; -sin(y) 0 cos(y)];
    Rz = [cos(z) -sin(z) 0; sin(z) cos(z) 0; 0 0 1];
    
    finalZ(ii) = z;
    
    e_init_rot = e_init * Rz * Rx * Ry;
    translation_vec = mean(e_dest) - mean(e_init_rot);
    e_init_trans = e_init_rot + translation_vec;
    
    % figure
    % hold on
    % scatter3(e_init_trans(:,1), e_init_trans(:,2), e_init_trans(:,3), 50, 'r', 'filled')
    % scatter3(e_dest(:,1), e_dest(:,2), e_dest(:,3), 50, 'b', 'filled')
    % axis equal
    % rotate3d on
    %
    
    % then we search the Y rotation parameter
    e_init = e_init_trans;
    
    step_size = 0.1; %deg
    X_search_range = deg2rad(0:0);
    if ii == 1
        Y_search_range = deg2rad(-20:step_size/10:0);
    else
        Y_search_range = deg2rad(-1:step_size/100:1);
    end
    Z_search_range = deg2rad(0:0);
    
    cost_tally = zeros(length(X_search_range), length(Y_search_range), length(Z_search_range));
    
    for ix = 1:length(X_search_range)
        for iy = 1:length(Y_search_range)
            for iz = 1:length(Z_search_range)
                x = X_search_range(ix);
                y = Y_search_range(iy);
                z = Z_search_range(iz);
                Rx = [1 0 0; 0 cos(x) -sin(x); 0 sin(x) cos(x)];
                Ry = [cos(y) 0 sin(y); 0 1 0; -sin(y) 0 cos(y)];
                Rz = [cos(z) -sin(z) 0; sin(z) cos(z) 0; 0 0 1];
                
                e_init_rot = e_init * Rz * Rx * Ry;
                translation_vec = mean(e_dest) - mean(e_init_rot);
                e_init_trans = e_init_rot + translation_vec;
                %
                %             figure
                %             hold on
                %             scatter3(e_init_trans(1,1), e_init_trans(1,2), e_init_trans(1,3), 100, 'k', 'filled')
                %             scatter3(e_dest(1,1), e_dest(1,2), e_dest(1,3), 70, 'md', 'filled')
                %             scatter3(e_init_trans(:,1), e_init_trans(:,2), e_init_trans(:,3), 50, 'r', 'filled')
                %             scatter3(e_dest(:,1), e_dest(:,2), e_dest(:,3), 50, 'b', 'filled')
                %             axis equal
                %             rotate3d on
                cost_tally(ix,iy,iz) = sum(vecnorm(e_init_trans - e_dest, 2, 2));
            end
        end
    end
    
    [val, lidx] = min(cost_tally(:));
    [ix,iy,iz] = ind2sub(size(cost_tally),lidx);
    
    
    x = X_search_range(ix);
    y = Y_search_range(iy);
    z = Z_search_range(iz);
    Rx = [1 0 0; 0 cos(x) -sin(x); 0 sin(x) cos(x)];
    Ry = [cos(y) 0 sin(y); 0 1 0; -sin(y) 0 cos(y)];
    Rz = [cos(z) -sin(z) 0; sin(z) cos(z) 0; 0 0 1];
    
    finalY(ii) = y;
    
    e_init_rot = e_init * Rz * Rx * Ry;
    translation_vec = mean(e_dest) - mean(e_init_rot);
    e_init_trans = e_init_rot + translation_vec;
    
    % figure
    % hold on
    % scatter3(e_init_trans(:,1), e_init_trans(:,2), e_init_trans(:,3), 50, 'r', 'filled')
    % scatter3(e_dest(:,1), e_dest(:,2), e_dest(:,3), 50, 'b', 'filled')
    % axis equal
    % rotate3d on
    
    
    
    
    % finally we search the x parameter
    e_init = e_init_trans;
    
    step_size = 0.1; %deg
    if ii == 1
        X_search_range = deg2rad(-5:step_size/10:5);
    else
        X_search_range = deg2rad(-1:step_size/100:1);
    end
    Y_search_range = deg2rad(0:0);
    Z_search_range = deg2rad(0:0);
    
    cost_tally = zeros(length(X_search_range), length(Y_search_range), length(Z_search_range));
    
    for ix = 1:length(X_search_range)
        for iy = 1:length(Y_search_range)
            for iz = 1:length(Z_search_range)
                x = X_search_range(ix);
                y = Y_search_range(iy);
                z = Z_search_range(iz);
                Rx = [1 0 0; 0 cos(x) -sin(x); 0 sin(x) cos(x)];
                Ry = [cos(y) 0 sin(y); 0 1 0; -sin(y) 0 cos(y)];
                Rz = [cos(z) -sin(z) 0; sin(z) cos(z) 0; 0 0 1];
                
                e_init_rot = e_init * Rz * Rx * Ry;
                translation_vec = mean(e_dest) - mean(e_init_rot);
                e_init_trans = e_init_rot + translation_vec;
                %
                %             figure
                %             hold on
                %             scatter3(e_init_trans(1,1), e_init_trans(1,2), e_init_trans(1,3), 100, 'k', 'filled')
                %             scatter3(e_dest(1,1), e_dest(1,2), e_dest(1,3), 70, 'md', 'filled')
                %             scatter3(e_init_trans(:,1), e_init_trans(:,2), e_init_trans(:,3), 50, 'r', 'filled')
                %             scatter3(e_dest(:,1), e_dest(:,2), e_dest(:,3), 50, 'b', 'filled')
                %             axis equal
                %             rotate3d on
                cost_tally(ix,iy,iz) = sum(vecnorm(e_init_trans - e_dest, 2, 2));
            end
        end
    end
    
    [val, lidx] = min(cost_tally(:));
    [ix,iy,iz] = ind2sub(size(cost_tally),lidx);
    
    x = X_search_range(ix);
    y = Y_search_range(iy);
    z = Z_search_range(iz);
    Rx = [1 0 0; 0 cos(x) -sin(x); 0 sin(x) cos(x)];
    Ry = [cos(y) 0 sin(y); 0 1 0; -sin(y) 0 cos(y)];
    Rz = [cos(z) -sin(z) 0; sin(z) cos(z) 0; 0 0 1];
    
    finalX(ii) = x;
    
    e_init_rot = e_init * Rz * Rx * Ry;
    translation_vec = mean(e_dest) - mean(e_init_rot);
    e_init_trans = e_init_rot + translation_vec;
    %
    % figure
    % hold on
    % scatter3(e_init_trans(:,1), e_init_trans(:,2), e_init_trans(:,3), 50, 'r', 'filled')
    % scatter3(e_dest(:,1), e_dest(:,2), e_dest(:,3), 50, 'b', 'filled')
    % axis equal
    % rotate3d on
    e_init = e_init_trans;
end

%% prepare the final output

TR = eye(3);

for ii = 1:iterations
    x = finalX(ii);
    y = finalY(ii);
    z = finalZ(ii);
    
    Rx = [1 0 0; 0 cos(x) -sin(x); 0 sin(x) cos(x)];
    Ry = [cos(y) 0 sin(y); 0 1 0; -sin(y) 0 cos(y)];
    Rz = [cos(z) -sin(z) 0; sin(z) cos(z) 0; 0 0 1];
    
    TR = TR * Rz * Ry * Rx;
end

e_init_archive = e_init_archive * TR;
trans_v = mean(e_dest) - mean(e_init_archive); % align centroid for translation
e_init_archive = e_init_archive + trans_v;

final_cost = sum(vecnorm(e_init_archive - e_dest, 2, 2));

figure
hold on
scatter3(e_init_archive(:,1), e_init_archive(:,2), e_init_archive(:,3), 50, 'r', 'filled')
scatter3(e_dest(:,1), e_dest(:,2), e_dest(:,3), 50, 'b', 'filled')
axis equal
rotate3d on


end

