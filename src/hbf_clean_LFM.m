function [ ] = hbf_clean_LFM(datapath, ico, verbose)
%% Clean LFM computed from HBF Solver for 4-shell forward models

% close all
% clear all

% % Change the current folder to the folder of this m-file.
% tmp = matlab.desktop.editor.getActive;
% cd(fileparts(tmp.Filename));
% clearvars tmp

if nargin < 2
    ico = 'ico5';
    verbose = false;
elseif nargin < 3
    verbose = false;
end

%% Path configuration

% add path to this EEG/code/forward_model folder
% addpath(genpath('./'));

% create datapath for our own BEM meshes
% datapath = '/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/MNE/subject_data/SP011_ah/simnibs_pipe/bem_surfaces/final_structure';

warning('off', 'all')

tic

%% Accomodate looping through ico3-5 source spaces

if strcmp(ico, '345')
    ico_list = {'ico3', 'ico4', 'ico5'};
else
    ico_list = {ico};
end

for ii = 1:length(ico_list)
    % iterate through different icosahedron resolutions
    ico = ico_list{ii};
    
    %% Examine LFM in closer details and post processing
    % Well. Let's try to examine the lead field matrix in closer details and
    % visualize the loading plots. This should give us confidence that the
    % solver is actually working.
    
    load(fullfile(datapath, ['4shell_hbfBEM_rawLFM-', ico, '.mat']), ...
        'LFMphi_dir', 'LFMphi_xyz', 'sourcepos', 'sourcedir', 'elecs', 'bmeshes')
    
    % Let's examine the leadfield matrix
    G = LFMphi_dir;
    xyzG = LFMphi_xyz;
    
    % check the histogram of G matrix entries
    figure;
    subplot(1,2,1)
    histogram(G(:), 100)
    title('Before detect_badsource()', 'FontSize', 16, 'interpreter', 'none')
    
    [ G, badsrc, badety, xyzG ] = detect_badsource(G, sourcepos, true, xyzG);
    disp(['Number of bad sources found: ' num2str(sum(badsrc))])
    disp(['Number of bad entries found: ' num2str(sum(badety))])
    
    subplot(1,2,2)
    histogram(G(:), 100)
    title('After detect_badsource()', 'FontSize', 16, 'interpreter', 'none')
    
    % Forward model variables
    forward_model.G = G;
    forward_model.xyzG = xyzG;
    forward_model.source = sourcepos;
    forward_model.dig = elecs.pproj;
    forward_model.normal = sourcedir;
    forward_model.skin_surf_face = bmeshes{end}.e;
    forward_model.skin_surf_vertice = bmeshes{end}.p;
    
    if verbose
        % Loading plot of an example channel (reference channel 85)
        plot_loading(85, forward_model, false, true, false);
        title('Does this look reasonable?', 'FontSize', 20)
        
        % visualize the bad sources
        badsrc_idx = find(badsrc==1);
        scatter3(sourcepos(badsrc_idx,1), sourcepos(badsrc_idx,2), sourcepos(badsrc_idx,3), 300, 'r', 'filled')
    end
    
    %% Further sanity checking of source gain
    % Now. What it looks like is that there are some artifact source points
    % near the base of the scalp on the medial surface of the temporal lobe
    % because that's when the pial surface accuracy is most messed up. It is
    % not worth writing an automatic algorithm to detect this becasue it has to
    % be judged on a case-by-base basis. What I will propose here instead is a
    % manual checking to identify erroneous source poitns and interpolate them
    % using mean of near-by-source points.
    
    % Due to the nature of these source points, they are most easily detectable
    % when you do the loading plot of extreme electrodes and at top central
    % electrodes. Hence we will iteratively search through these and use
    % histograms to guide our marking of source points to interpolate.
    
    % we can clearly see these points using electrode 47 or 1 (mastoid
    % electrodes)
    
    % Right mastoid electrode
    channel = 47;
    high_src_p = hbf_check_LFM_outlier(channel, forward_model, 0.999, G, verbose);
    disp(['Number of bad sources detected according to right mastoid: ', num2str(length(high_src_p))])
    % interpolate these bad sources
    [ Gprime, xyzGprime ] = interpolate_sources(G, sourcepos, high_src_p, xyzG);
    % let's see how it looks after fixing.
    forward_model.G = Gprime;
    if verbose
        plot_loading(channel, forward_model, false, true, false);
        title('Right Mastoid Post-fixing', 'FontSize', 20)
    end
    
    % Repeat the same thing for the left mastoid electrode
    channel = 1;
    high_src_p = hbf_check_LFM_outlier(channel, forward_model, 0.999, Gprime, verbose);
    disp(['Number of bad sources detected according to left mastoid: ', num2str(length(high_src_p))])
    % interpolate these bad sources
    [ Gprime, xyzGprime ] = interpolate_sources(Gprime, sourcepos, high_src_p, xyzGprime);
    % let's see how it looks after fixing.
    forward_model.G = Gprime;
    if verbose
        plot_loading(channel, forward_model, false, true, false);
        title('Left Mastoid Post-fixing', 'FontSize', 20)
    end
    
    % Iterate through all other channels and keep track of which are outlier
    % sources
    outlier_source = [];
    for j = 1:size(G,1)
        if ~ismember(j, [1, 47]) % dont repeat it for mastoids
            high_src_p = hbf_check_LFM_outlier(j, forward_model, 0.9999, Gprime);
            outlier_source = [outlier_source, high_src_p]; %#ok<*AGROW>
        end
    end
    outlier_source = unique(outlier_source); % find the unique ones
    disp(['Total number of bad sources detected pooling across channels: ', num2str(length(outlier_source))])
    
    % visualize all the strange sources
    if verbose
        plot_loading(85, forward_model, false, true, false);
        title({'Reference channel Pre-fixing', 'strange sources across channels marked in Red'}, 'FontSize', 16)
        scatter3(forward_model.source(outlier_source,1), forward_model.source(outlier_source,2), forward_model.source(outlier_source,3), 300, 'r', 'filled')
    end
    % interpolate these bad sources
    [ Gprime, xyzGprime ] = interpolate_sources(Gprime, sourcepos, outlier_source, xyzGprime);
    % let's see how it looks after fixing.
    forward_model.G = Gprime;
    if verbose
        plot_loading(85, forward_model, false, true, false);
        title('Reference channel Post-fixing', 'FontSize', 20)
    end
    
    %%
    % ok - Update G;
    G = Gprime;
    xyzG = xyzGprime;
    forward_model.G = G;
    forward_model.xyzG = xyzG;
    
    % in theory we could have done the same with patch decomposition, but i
    % think patch decomposition sacrifices the number of sources so it's better
    % to manually identify the very few odd source points and preserve the
    % source resolution. If we need to patch decompose later, that's fine.
    
    % Let's save all the data at this point.
    save(fullfile(datapath, ['4shell_hbfBEM_finalLFM-', ico, '.mat']), 'G', 'xyzG', 'forward_model')

end

%% Completed!

warning('on', 'all')

disp('Total time taken:') 
toc

pause(10)
close all

end



function [ high_src_p ] = hbf_check_LFM_outlier(channel, forward_model, pcut, Gprime, ploton)
%Used to detect outliers in the LFM on a specified channel 

if nargin < 5
    ploton = false;
end

% check a channel 
if ploton
    plot_loading(channel, forward_model, false, true, false);
    title(['Channel ', num2str(channel), ' Pre-fixing - bad sourecs marked in Red'], 'FontSize', 20)
end

% compute a magnitude vector weighted by distance from the electrode
distance2elc = vecnorm(forward_model.source - forward_model.dig(channel,:), 2, 2);
mag_vec = Gprime(channel,:).^2 .* distance2elc';

% fit exponential distribution to the mag_vec
pd = fitdist(mag_vec', 'Exponential');
cutoff_value = icdf(pd, pcut);

% use this magnitude vector to find outliers
high_src_p = find(mag_vec > cutoff_value);

if ploton
    % plot the problematic sources
    scatter3(forward_model.source(high_src_p,1), forward_model.source(high_src_p,2), forward_model.source(high_src_p,3), 300, 'r', 'filled')
end

end
