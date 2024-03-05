clearvars; clc; close all;

code_folder = 'C:/Users/drdus/emotional_memory_neuronal_analyses/'; % Make sure it ends with /
data_folder = 'C:/Users/drdus/emotional_memory_neuronal_data/'; % Make sure it ends with / and contains the enc and rec folders

addpath([code_folder,'dPCA_DF']);
addpath([code_folder,'utils']);
folder = [data_folder,'enc/'];

% This is setup to create Figure 3
toA = {'all'};  % Use this optin and includeKnow = false for Fig. 3
% toA = {'all','noHPC', 'noAmyg', 'noEC'};  % Use these options and includeKnow for SFig. 4

includeKnow = false; % Set true for SFig. 4

% This saves the optimalLambda values calculated in the top section.
oL = zeros(1,length(toA)); % optimalLamba = [0.0007    0.0007    0.0007    0.0011] for {'all','noHPC', 'noAmyg', 'noEC'}

for ia = 1 :length(toA)

    load([folder,'_dPCA_input_RKF.mat']); % Reload data to plot different leave-one-out results

    firingRates_sparse = firingRates_sparseRKF;
    firingRates_size = firingRates_sizeRKF;
    firingRatesAverage = firingRatesAverageRKF;
    trialNum = trialNumRKF;
    
    timeK = [-0.2 1.5]; % seconds, analysis window
    timeKeep = (time < timeK(2)) & (time > timeK(1));
    time = time(timeKeep); % Reduce the time axis to improve results
    
    % firingRates array is stored in the file in compressed sparse format. This line de-compresses it.
    firingRates = reshape(full(firingRates_sparse), firingRates_size);

    ayType = toA{ia}; 
    if strcmp(ayType,'noHPC')
        subjUse = areaMask == 4 | areaMask == 3; % areaMask = iAHP .* 1 + iPHP .*2 + iAmg .* 3 + iEC .* 4;
    elseif strcmp(ayType,'noAmyg')
        subjUse = areaMask == 1 | areaMask == 2 | areaMask == 4; % areaMask = iAHP .* 1 + iPHP .*2 + iAmg .* 3 + iEC .* 4;
    elseif strcmp(ayType,'noEC')
        subjUse = areaMask == 1 | areaMask == 2 | areaMask == 3; % areaMask = iAHP .* 1 + iPHP .*2 + iAmg .* 3 + iEC .* 4;
    elseif strcmp(ayType,'all')
        subjUse = true(size(firingRates,1),1)';
        ayType = 'all';
    end

    % Select the  desired brain areas or subjects
    firingRates = firingRates(subjUse,:,:,:,:);
    firingRatesAverage = firingRatesAverage(subjUse,:,:,:,:);
    trialNum = trialNum(subjUse,:,:);
    
    % Neuron selection criteria: AT least 3 trials and firing rate less than 50 Hz
    D = size(trialNum,1);
    minN = min(reshape(trialNum(:,:,:), D, []), [], 2);
    meanFiringRate = mean(reshape(firingRatesAverage, D, []), 2);
    n = find(minN >= 3 & meanFiringRate < 50);
    
    firingRates = firingRates(n,:,:,:,:);
    firingRatesAverage = firingRatesAverage(n,:,:,:,:);
    trialNum = trialNum(n,:,:);
    subjMaskFiltered = subjMask(n);
    
    combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
    margNames = {'Memory', 'Emotion', 'Ind.', 'M x E'};
    if includeKnow
        decodingClasses = {[1 1; 2 2; 3 3], [1 2; 1 2; 1 2], [], [1 2; 3 4; 5 6]};
        yLims = [4.5 4.5 6 4.5];
        firingRates = firingRates(:,:,:,timeKeep,:);
        firingRatesAverage = firingRatesAverage(:,:,:,timeKeep,:);
    else % Remove Know trials
        decodingClasses = {[1 1; 2 2], [1 2; 1 2], [], [1 2; 3 4]};
        yLims = [4.5 6.5 8.5 5.5];
        firingRates = firingRates(:,1:2,:,timeKeep,:);
        firingRatesAverage = firingRatesAverage(:,1:2,:,timeKeep,:);
        trialNum = trialNum(:,1:2,:);
    end
    margColours = [rgb('DarkGreen'); rgb('DarkOrange'); rgb('Gray'); [0.4940 0.1840 0.5560]];
    
    %% Cross-validation to find lambda 
    
    % optimalLambda = 4.9879e-04 (If this gives a different number then results/percentages can be different)
    optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
        'combinedParams', combinedParams, ... 
        'numRep', 500, ...
        'filename', 'tmp_optimalLambdas.mat');
    oL(ia) = optimalLambda; % Save optimalLambda values

    
    %% dPCA (with regularization and noise cov)
    
    Cnoise = dpca_getNoiseCovariance(firingRatesAverage, firingRates, trialNum, ...
        'type', 'pooled'); % or averaged?
      
    [W,V,whichMarg, Z] = dpca_df(firingRatesAverage, 50, ...
        'combinedParams', combinedParams, 'lambda', optimalLambda, 'Cnoise', Cnoise);
    
    if includeKnow % This creates some of the Supplementary figures
        explVarK = dpca_explainedVariance(firingRatesAverage, W, V, ...
         'combinedParams', combinedParams, ...
         'Cnoise', Cnoise, 'numOfTrials', trialNum);

        % Plot the first column for SFig. 4
        dpca_plot_1col_df(firingRatesAverage, W, V, @dpca_plot_iapsEnc_2s_2d_2color_KNOW, ...
            'whichMarg', whichMarg,                 ...
            'time', time,                           ...
            'timeEvents', timeEvents,               ...
            'timeMarginalization', 3,               ...
            'ylims', yLims,                  ...
            'legendSubplot', 16,                    ...
            'marginalizationNames', margNames,      ...
            'marginalizationColours', margColours,  ...
            'explainedVar', explVarK, ...
            'numCompToShow',1); 

        saveas(gcf(), [folder, sprintf('_SFig_4_dPCA_1col_%s_%ineurons_RKF.svg',ayType, length(subjMaskFiltered))]);
        exportgraphics(gcf(), [folder,sprintf('_SFig_4_dPCA_1col_%s_%ineurons_RKF.png',ayType, length(subjMaskFiltered))]);

    else % The main figures are created here
        explVar = dpca_explainedVariance_PCAPlots_ENC_RvF(firingRatesAverage, W, V, time, timeEvents, ...
         'combinedParams', combinedParams, ...
         'Cnoise', Cnoise, 'numOfTrials', trialNum);
        saveas(gcf(), [folder, sprintf('_enc_PCA_%s_%ineurons.svg',ayType, length(subjMaskFiltered))]); % Save regular PCs
        exportgraphics(gcf(), [folder, sprintf('_enc_PCA_%s_%ineurons.png',ayType, length(subjMaskFiltered))]);

        Zfull = dpca_plot_df(firingRatesAverage, W, V, @dpca_plot_iapsEnc_2s_2d_2color, ...
            'whichMarg', whichMarg,                 ...
            'time', time,                           ...
            'timeEvents', timeEvents,               ...
            'timeMarginalization', 3,               ...
            'ylims', yLims,                ...
            'legendSubplot', 16,                    ...
            'marginalizationNames', margNames,      ...
            'marginalizationColours', margColours,  ...
            'explainedVar', explVar);     
        
        styles = {'-','--',':'};
        subplot(4,4,9); hold on;
        numOfDecisions = size(Zfull, 2); % first dim of Zfull is ordered by plot order not component order
        for f=1:numOfDecisions 
            plot(-squeeze(Zfull(7, f, 1, :)), squeeze(Zfull(4, f, 1, :)), styles{f}, 'color', rgb('FireBrick'), 'LineWidth', 1.0)
            plot(-squeeze(Zfull(7, f, 2, :)), squeeze(Zfull(4, f, 2, :)), styles{f}, 'color', rgb('Navy'), 'LineWidth', 1.0)
        end
        xlabel('dPC 2')
        ylabel('dPC 4')
        
        subplot(4,4,13); hold on;
        numOfDecisions = size(Zfull, 2); % first dim of Zfull is ordered by plot order not component order
        for f=1:numOfDecisions 
            plot(-squeeze(Zfull(7, f, 1, :)), squeeze(Zfull(10, f, 1, :)), styles{f}, 'color', rgb('FireBrick'), 'LineWidth', 1.0)
            plot(-squeeze(Zfull(7, f, 2, :)), squeeze(Zfull(10, f, 2, :)), styles{f}, 'color', rgb('Navy'), 'LineWidth', 1.0)
        end
        xlabel('dPC 2')
        ylabel('dPC 3')

        fontsize(gcf(),8,'points');
        fontname(gcf(),"Arial");
        saveas(gcf(), [folder, sprintf('dPCA_%s_.svg',ayType)]);
        exportgraphics(gcf(), [folder,sprintf('dPCA_%s_.png',ayType)]);
    end
end

%% decoding part


% Only perform these analyses if excluding Know trials and using all brain regions
if ~includeKnow && strcmp(toA{ia},'all')

    % with 100 iterations this takes around 10*100/60 = 17 min on my laptop
    accuracy = dpca_classificationAccuracy(firingRatesAverage, firingRates, trialNum, ...
        'lambda', optimalLambda, ...
        'numComps', 3, ...
        'combinedParams', combinedParams, ...
        'decodingClasses', decodingClasses, ...
        'noiseCovType', 'pooled', ...
        'numRep', 1000, ...        % increase to 1000
        'filename', 'tmp_classification_accuracy_1000.mat'); % numShuffles = 1000
    
    % with 100 iterations and 100 shuffles this takes 100 times longer than the
    % above function, i.e. 17*100/60 = 28 hours (on my laptop). Be careful.
    % tic
    [accuracyShuffle, brierShuffle] = dpca_classificationShuffled(firingRates, trialNum, ...
        'lambda', optimalLambda, ...
        'combinedParams', combinedParams, ...
        'decodingClasses', decodingClasses, ...
        'noiseCovType', 'pooled', ...
        'numRep', 100, ...        % increase to 100
        'numShuffles', 1000, ...  % increase to 1000 (takes a lot of time)
        'filename', 'tmp_classification_accuracy_shuffle_1000.mat'); % numShuffles = 1000
    % toc
    
    % Get the component significance
    componentsSignif = dpca_signifComponents(accuracy, accuracyShuffle, whichMarg);
    
    % The rest of the code is based on 1000 shuffles and you need to edit the selection if you use a different number of shuffles
    num = 1; % always the first dPC (first column)
    marg = 1;
    B = squeeze(sort(accuracyShuffle(marg,:,1:end),3)); % Sort the shuffled accuracies so we can get p-values
    C = squeeze(accuracy(marg,num,:)) > B;
    CI95_dPC4 = C(:,end-25)'; % Set the signifiance level where p<0.025 for the first column if using 1000 shuffled (25/1000 = 0.025)
    
    marg = 2;
    B = squeeze(sort(accuracyShuffle(marg,:,1:end),3)); % Sort the shuffled accuracies so we can get p-values
    C = squeeze(accuracy(marg,num,:)) > B;
    CI95_dPC2 = C(:,end-25)'; % Set the signifiance level where p<0.025 for the first column if using 1000 shuffled (25/1000 = 0.025)
    
    marg = 4;
    B = squeeze(sort(accuracyShuffle(marg,:,1:end),3)); % Sort the shuffled accuracies so we can get p-values
    C = squeeze(accuracy(marg,num,:)) > B;
    CI95_dPC3 = C(:,end-25)';  % Set the signifiance level where p<0.025 for the first column if using 1000 shuffled (25/1000 = 0.025)
    
    %% Final manuscript plot: Fig. 3
    
    % Update the significance thresholds for plots of main dPCs
    componentsSignif_95CI = componentsSignif; 
    componentsSignif_95CI(2,:) = CI95_dPC2;
    componentsSignif_95CI(3,:) = CI95_dPC3;
    componentsSignif_95CI(4,:) = CI95_dPC4;
    
    events = [0, 0.5];
    dpca_plot_df(firingRatesAverage, W, V, @dpca_plot_iapsEnc_2s_2d_2color, ...
        'whichMarg', whichMarg,                 ...
        'time', time,                           ...
        'timeEvents', timeEvents,               ...
        'timeMarginalization', 3,               ...
        'marginalizationColours', margColours, ...
        'ylims', yLims,           ...
        'legendSubplot', 16,                    ...
        'marginalizationNames', margNames,      ...
        'explainedVar', explVar, ...
        'componentsSignif', componentsSignif_95CI, ...
        'showNonsignificantComponents', true, ...
        'legendSubplot', 4);
    
    styles = {'-','--',':'};
    subplot(4,4,9); hold on;
    numOfDecisions = size(Zfull, 2); % first dim of Zfull is ordered by plot order not component order
    
    for f=1:numOfDecisions 
        plot(-squeeze(Zfull(7, f, 1, :)), squeeze(Zfull(4, f, 1, :)), styles{f}, 'color', rgb('FireBrick'), 'LineWidth', 1.0)
        plot(-squeeze(Zfull(7, f, 2, :)), squeeze(Zfull(4, f, 2, :)), styles{f}, 'color', rgb('Navy'), 'LineWidth', 1.0)
    end
    xlabel('dPC 2')
    ylabel('dPC 4')
    
    subplot(4,4,13); hold on;
    numOfDecisions = size(Zfull, 2); % first dim of Zfull is ordered by plot order not component order
    for f=1:numOfDecisions 
        plot(-squeeze(Zfull(7, f, 1, :)), squeeze(Zfull(10, f, 1, :)), styles{f}, 'color', rgb('FireBrick'), 'LineWidth', 1.0)
        plot(-squeeze(Zfull(7, f, 2, :)), squeeze(Zfull(10, f, 2, :)), styles{f}, 'color', rgb('Navy'), 'LineWidth', 1.0)
    end
    xlabel('dPC 2')
    ylabel('dPC 3')
    
    % Delete plots to make space for the linear decoder results
    delete(subplot(4,4,7))
    delete(subplot(4,4,8))
    delete(subplot(4,4,11))
    delete(subplot(4,4,12))
    delete(subplot(4,4,15))
    
    si = [6,9,12]; % subplot indices
    mi = [1,2,4]; % marginalization (marg) indices
    num = 1; % always the first dPC in the type
    
    for i = 1:length(mi)
        marg = mi(i);
        
        subplot(4, 3, si(i)); hold on;
        B = squeeze(sort(accuracyShuffle(marg,:,1:end),3)); % Sort the shuffled accuracies so we can get p-values
        C = squeeze(accuracy(marg,num,:)) > B;
        
        sh_time = time(1:size(accuracyShuffle, 2));
        
        maxSh = B(:,end-50)'; % When using 1000 shuffles, this looks for areas greater than 50/1000 or p < 0.05 (two-tailed)
        minSh = B(:,50)';
        h = patch([sh_time fliplr(sh_time)], [maxSh fliplr(minSh)], 'k');
        set(h, 'FaceAlpha', 0.15)
        set(h, 'EdgeColor', 'none')
        
        maxSh = B(:,end-25)'; % When using 1000 shuffles, this looks for areas greater than 25/1000 or p < 0.025 (two-tailed)
        minSh = B(:,25)'; 
        h = patch([sh_time fliplr(sh_time)], [maxSh fliplr(minSh)], 'k');
        set(h, 'FaceAlpha', 0.1)
        set(h, 'EdgeColor', 'none')
        
        maxSh = B(:,end)'; % When using 1000 shuffles, this looks for areas greater than 1/1000 or p < 0.001 (two-tailed)
        minSh = B(:,1)';
        h = patch([sh_time fliplr(sh_time)], [maxSh fliplr(minSh)], 'k');
        set(h, 'FaceAlpha', 0.05)
        set(h, 'EdgeColor', 'none')
        
        plot(time,squeeze(accuracy(marg,num,:))','Color',rgb('DarkOrange')); 
        Cd = double(C);
        Cd(Cd==0) = nan;
        if i < 3
            plot(time, Cd(:,end)*0.175, 'Color',[0, 0, 0, 0.9], 'LineWidth', 1.0);
            plot(time, Cd(:,end-25)*0.15, 'Color',[0, 0, 0, 0.8], 'LineWidth', 1.0);
            plot(time, Cd(:,end-50)*0.125, 'Color',[0, 0, 0, 0.7], 'LineWidth', 1.0); 
            ylim([0,1]);
        else
            plot(time, Cd(:,end)*0.125, 'Color',[0, 0, 0, 0.9], 'LineWidth', 1.0) % Plot the significance lines slightly lower for better illustration
            plot(time, Cd(:,end-25)*0.1, 'Color',[0, 0, 0, 0.8], 'LineWidth', 1.0)
            plot(time, Cd(:,end-50)*0.075, 'Color',[0, 0, 0, 0.7], 'LineWidth', 1.0)
            ylim([0,0.5]); % Set smaller ylimit for the interaction because random chance here is 0.25 not 0.5
            xlabel('Time (s)');
        end
    
        if i == 2
            ylabel('Classification accuracy')
        end
    
        xlim([-0.2,1.5]);
        
        xline(timeEvents(1)); xline(timeEvents(2));
        pos = get(gca, 'Position');
    %     xoffset = -0.1;
        pos(1) = pos(1) - 0.1; % Move the figure to the left
        set(gca, 'Position', pos);
        yticks([0, 0.5, 1]);
        set(gca,'XTick',[0, 0.5, 1, 1.5]);
    end
    
    fontsize(gcf(),8,'points');
    fontname(gcf(),"Arial");
    
    exportgraphics(gcf(), [folder,'fig3_enc_dPCA.png'],'Resolution',300);
    saveas(gcf(), [folder,'fig3_enc_dPCA.svg']);
end