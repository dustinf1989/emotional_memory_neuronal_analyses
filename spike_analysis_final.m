% ==============================================================================================================================================================
% Script written by Dustin Fetterhoff to analyze neuronal spike train data recorded during the IAPS emotional memory task
% ==============================================================================================================================================================

clearvars; clc; close all; % clearvars is better than clear all for performance (see warning for details)

code_folder = 'D:\iaps_data_share\code\'; % Make sure it ends with \
data_folder = 'D:\iaps_data_share\'; % Make sure it ends with \ and contains enc and rec folders with dataset

addpath(code_folder);
addpath([code_folder,'fieldtrip-20220531-DF']);
addpath(genpath([code_folder,'boundedline']));
ft_defaults

type = 'combinato_pn'; % Positive and negative waveforms - used for reading spike train data files

% Uncomment the task you want to analyze
task = 'enc'; % IAPS encoding
% task = 'rec'; % IAPS recognition

excludeStimOnset = true; % Remove artifacts occurring +-3ms around stimulus onset and offset
posOnly = true; % Set true to use only positive waveforms - Gets negative from P13 since it was not inverted

regNames = ["Hippocampus" "Hippocampus" "Amygdala" "Entorhinal Cortex" "Pararhinal cortex"];
patahp = "u" + ("HR"|"HL"|"AHR"|"AHL"); % Decode microwire names
patphp = "u" + ("PHR"|"PHL");
patamg = "u" + ("AR"|"AL");
patec = "u" + ("ECR"|"ECL");
patph = "uPCL";

if strcmp(task,'enc')
    thr_trials = 50; % Number of active trials required to include the neuron
    outsmall = {'R','F','K'};
elseif strcmp(task,'rec')
    thr_trials = 100; % Number of active trials required to include the neuron
    outsmall = {'RHit','Miss','CR','KHit'};
end
mfr_thresh = 0.25; % Mean firing rate threshold during task for neuron to be included.
tl = 0.2; % Lowest time to check for average rate comparison
th = 1.5; % Highest time to check for average rate comparison
nShuffAnova = 100; %10,000 for manuscript

psthbin = 0.01; % PSTH bin size MUST CHANGE THE FIRING RATE VECTORS FOR PCA!
t_pre_stim = 2.5;
t_pst_stim = 5.0; % New trial begins at 4s but we want a buffer for the previous response and removing the edge of the Gaussian filter
start_stop = [-t_pre_stim t_pst_stim]; % start and end the PSTH

% Initialize for dPCA
% firingRates: N x S x D x T x maxTrialNum
% firingRatesAverage: N x S x D x T
if strcmp(task,'rec')
    firingRatesRec = zeros(500, 4, 2, 499, 80);
    firingRatesAverageRec = zeros(500, 4, 2, 499);
elseif strcmp(task,'enc')
    firingRatesRKF = zeros(500, 2, 2, 451, 80);
    firingRatesAverageRKF = zeros(500, 2, 2, 451);
end

if strcmp(task,'enc') 
    paths.spikes = [data_folder,'enc\']; % Folder containing combinato output that was reformatted in Python
    paths.pics = [paths.spikes, sprintf('results_enc_%0.1fs_to_%0.1fs_%0.2fmfr_%dtrials_%0.2fbinsize_%ish_CHECK5',tl,th,mfr_thresh,thr_trials,psthbin,nShuffAnova),'\']; % Pictures and data output stored here
elseif strcmp(task,'rec')
    paths.spikes = [data_folder,'rec\'];
    paths.pics = [paths.spikes, sprintf('results_rec_%0.1fs_to_%0.1fs_%0.2fmfr_%dtrials_%0.2fbinsize_%ish_CHECK5',tl,th,mfr_thresh,thr_trials,psthbin,nShuffAnova),'\'];
end

subjects    = { ...
    'Patient1',  {'uAL'; 'uAR'; 'uAHL'; 'uAHR'; 'uECL'; 'uECR'; 'uPCL'; 'uPHR'}; ... 
    'Patient5',  {'uAL'; 'uAR'; 'uHL';  'uAHR'; 'uECL'; 'uECR'; 'uPHR'}; ...
    'Patient6',  {'uAL'; 'uAR'; 'uHL';  'uAHR'; 'uECR'; 'uPHR'}; ...
    'Patient8',  {'uAL'; 'uAR'; 'uAHL'; 'uAHR'; 'uECL'; 'uECR'; 'uPHL'; 'uPHR'}; ...
    'Patient13', {'uAL'; 'uAR'; 'uAHL'; 'uAHR'}; ...  
    };

if ~exist(paths.pics, 'dir')
    mkdir(paths.pics);
end

sp = 0; % Spike Index for saving all of them in one structure

%% loop through subjects
for iSub = 1:size(subjects, 1)
    ni = 0; % Neuron index
    data_folder = strcat(paths.spikes, subjects{iSub, 1},'\');

    if strcmp(task,'enc')
        load(strcat(data_folder,'enc_onsets_v2.mat')); % load trial type outcomes
        onsets = enc_onsets;
    elseif strcmp(task,'rec')
        load(strcat(data_folder,'rec_onsets_v2.mat')); % load trial type outcomes
        onsets = rec_onsets;
    end
    sub(iSub).onsets = onsets; % Save trial types to subjects' structure
 
    %% loop through regions
    for iReg = 1:size(subjects{iSub, 2}, 1)
        
        % Find the spike files for each brain region
        spike_folders   = dir(strcat(data_folder,subjects{iSub, 2}{iReg}, '*'));

        %% loop through files and read them
        for iFile = 1:size(spike_folders, 1)

            % Electrode name for figures 
            [~, f, ~] = fileparts(spike_folders(iFile).name);
            t = tokenize(f, '_'); % extract channel name
            rowsWithString = contains(t, 'u');
            ii = find(rowsWithString);
            s1=[];
            for i=1:sum(rowsWithString)
                s1 = [s1, t{ii(i)}, '_']; % Electrode string names
            end

            file = strcat(data_folder, spike_folders(iFile).name,'\py_spikes_posneg_',spike_folders(iFile).name,'.mat');
            load([data_folder,'ft_trial_',subjects{iSub},'.mat']); % file contains 4 variables: 'event','hdr','trl','trl_sec'

            %% Align spike data with event and ncs data
            if exist(file,'file') % When no spikes are detected, there is no file in the folder
                spike = ft_read_spike(file,'spikeformat',type); % wave_clus: Timesteps are in ms, starting at zero already?
                if length(spike.label) > 0
                    for i=1:length(spike.unit)
                        spike.times_us{i} = spike.timestamp{i}.*1000 - double(hdr.FirstTimeStamp); % in us
                        spike.time{i} = spike.times_us{i} ./ 1000000; % in seconds
                        spike.timestamp{i} = spike.time{i}; % seconds, default in some fieldtrip functions
                        spike.trial{i} = zeros(size(spike.time{i}));
                    end
    
                    % Remove spikes at almost exact the stimulus presentation time because these are artifacts.
                    if excludeStimOnset == true
                        for i1 = 1:length(spike.time)
                            kai = []; 
                            for wi=1:length(event)
                                ki = find(event(wi).timestamp-0.003 < spike.time{i1} & event(wi).timestamp+0.003 > spike.time{i1}); % spike index to remove
                                if ~isempty(ki)
                                    %wai = [wai, wi]; % To check which events had noise spike
                                    kai = [kai, ki]; % kai is a list of spike indicies to exclude
                                end
                            end
                            spike.id_excluded{i1} = kai; % Save excluded indicies to look at waveforms later
                            spike.time{i1}(kai) = []; % Remove the noise from the spike train
                            spike.timestamp{i1}(kai) = []; % Remove the noise from the spike train
                            spike.trial{i1}(kai) = []; % Remove the noise from the spike train
                        end
                    end
                    
                    % Convert timestamps to seconds & start first recorded timestamp at zero
                    spike_notrials = spike;
                    for i=1:length(spike_notrials.label)
                        spike_notrials.trial{i} = ones(1,length(spike_notrials.time{i}));
                    end
                    
                    for si = 1:length(spike.time) % Number of units
                        for i = 1:length(trl)
                            q = (spike.time{si} > trl_sec(i,1) & spike.time{si} < trl_sec(i,2)); % Used to index
                            spike.trial{si}(q) = i;
                        end
                    end
                
                    % Trial structure for spikes: If trials are overlapping, looks like spikes are duplicated so don't get MFR based on this.
                    cfg          = [];
                    cfg.trl = trl_sec; % Use it in seconds and avoid conversion below
                    cfg.timestampspersecond =  1; %spike.hdr.sr; 
                    spikeTrials = ft_spike_maketrials(cfg,spike);
                  
                    for k=1:length(spike.unit)
        
                        lenTrl = trl_sec(end,2) - trl_sec(1,1); % in seconds
                        expSpikes = spike_notrials.time{k} > trl_sec(1,1) & spike_notrials.time{k} < trl_sec(end,2);
                        mean_fr = sum(expSpikes) / lenTrl; % Mean firing rate during the task (some datasets have extra recording before and after the task)
                        n_trials = length(unique(spikeTrials.trial{k})); % number of trials with action potentials
                        
                        if posOnly && strcmp(subjects{iSub}, 'Patient13')
                            wvt = contains(spike.label{k},'neg'); % Signal was NOT inverted in Patient13 and therefore we should use negative spikes 
                            wtName = 'neg';
                        elseif posOnly && (~strcmp(subjects{iSub}, 'Patient13'))
                            wvt = contains(spike.label{k},'pos'); % waveform type filter
                            wtName = 'pos'; 
                        else % if not posOnly analyze all waveforms
                            wvt = true; % waveform type filter
                            if contains(spike.label{k},'pos')
                                wtName = 'pos'; % Temp variable
                            elseif contains(spike.label{k},'neg')
                                wtName = 'neg';
                            end
                        end
                        perc_less3ms = sum(diff(spike_notrials.time{k}(expSpikes)) < 0.003) / sum(expSpikes) * 100; % percentage of ISIs during experiment less than 3ms
        
                        % Only analyze neurons that meet these criteria
                        if (mean_fr > mfr_thresh) && (n_trials > thr_trials) && wvt && (perc_less3ms < 3)
        
                            ni = ni + 1; % Counter within subject. Used for indexing sub() struct.
                            sp = sp + 1; % Counter over all subjects - reaches max number of neurons. Used for dPCA matrices.
    
                            sub(iSub).meanFR(ni) = mean_fr;
                            sub(iSub).n_trials(ni) = n_trials;
                            sub(iSub).name{ni} = [s1,'n',num2str(k)]; % Add neuron name to structure
                            sub(iSub).neurLabel{ni} = spike.label{k};
                            sub(iSub).perc_less3ms(ni) = perc_less3ms; % %ISIs during experiment less than 3ms
        
                            [e,~] = ismember(spikeTrials.trial{k},onsets.emo_all); % Used for raster plots
                            [n,~] = ismember(spikeTrials.trial{k},onsets.neu_all); 
        
                            Eall = onsets.emo_all; % Shorthand variables for trial indices
                            Nall = onsets.neu_all;
        
                            if strcmp(task,'enc')
                                % Indices used for raster plots
                                [eR,~] = ismember(spikeTrials.trial{k},onsets.eR); 
                                [eF,~] = ismember(spikeTrials.trial{k},onsets.eF); 
                                [nR,~] = ismember(spikeTrials.trial{k},onsets.nR); 
                                [nF,~] = ismember(spikeTrials.trial{k},onsets.nF);
            
                                % Shorthand variables for trial indices
                                Er = onsets.eR;
                                Ef = onsets.eF;
                                Nr = onsets.nR;
                                Nf = onsets.nF;
        
                            elseif strcmp(task,'rec')
                                % Indices used for raster plots
                                [eR,~] = ismember(spikeTrials.trial{k},onsets.eRHit); 
                                [eM,~] = ismember(spikeTrials.trial{k},sort(onsets.eMiss));
                                [eCR,~] = ismember(spikeTrials.trial{k},onsets.eCR);
                                [nR,~] = ismember(spikeTrials.trial{k},onsets.nRHit);
                                [nM,~] = ismember(spikeTrials.trial{k},sort(onsets.nMiss)); 
                                [nCR,~] = ismember(spikeTrials.trial{k},onsets.nCR); 
                                 
                                % Shorthand variables for trial indices
                                Er = onsets.eRHit;
                                Em = sort(onsets.eMiss);
                                Ecr = onsets.eCR;
                                Nr = onsets.nRHit;
                                Nm = sort(onsets.nMiss);
                                Ncr = onsets.nCR;
                            end
        
                            % Compute average waveform to get spike width sued in SFig. 1
                            mean_wav = mean(spike.waveform{k},3); % 3rd dimension to get average waveform                        
                            [sub(iSub).width(ni), sub(iSub).ratio(ni)] = get_spike_width_ratio(mean_wav);
                                
                            %% Computing spike densities and peri-stimulus time histograms (PSTHs)
        
                            % This PSTH uses bigger bins sizes to get spike counts in windows to make comparisons between conditions later
                            cfg             = [];
                            cfg.binsize     =  0.1; % if cfgPsth.binsize = 'scott' or 'sqrt', we estimate the optimal bin size from the data itself
                            cfg.outputunit  = 'spikecount'; % give as an output the firing rate
                            cfg.latency     = start_stop; % between -1 and 3 sec.
                            cfg.vartriallen = 'no'; % variable trial lengths are accepted
                            cfg.keeptrials  = 'yes'; % keep the psth per trial in the output
                            psth_cnt = ft_spike_psth(cfg,spikeTrials);
                            baseT = (psth_cnt.time < -tl) & (psth_cnt.time > -th);
                            stimT = (psth_cnt.time > tl) & (psth_cnt.time < th);
        
                            psth_cnt_trial = squeeze(psth_cnt.trial(:,k,:));
                            sub(iSub).baseCnts(:,ni) = sum(psth_cnt_trial(:,baseT),2);
                            sub(iSub).stimCnts(:,ni) = sum(psth_cnt_trial(:,stimT),2);
                            
                            % Used later in Fig 2A: Comparing spike counts between baseline and stimulus periods
                            [~, sub(iSub).pEmoCnt(ni)] = permutest( sub(iSub).baseCnts(Eall,ni)', sub(iSub).stimCnts(Eall,ni)', true, 0.5, 1000, true ); % Emotional
                            [~, sub(iSub).pNeuCnt(ni)] = permutest( sub(iSub).baseCnts(Nall,ni)', sub(iSub).stimCnts(Nall,ni)', true, 0.5, 1000, true ); % Neutral
                            [~, sub(iSub).peRCnt(ni)] = permutest( sub(iSub).baseCnts(Er,ni)', sub(iSub).stimCnts(Er,ni)', true, 0.5, 1000, true ); % emotional remembered (eR)
                            [~, sub(iSub).pnRCnt(ni)] = permutest( sub(iSub).baseCnts(Nr,ni)', sub(iSub).stimCnts(Nr,ni)', true, 0.5, 1000, true ); % neutral remembered (nR)
        
                            if strcmp(task,'enc')
    
                                % Used later in Fig 2A
                                [~,sub(iSub).peFCnt(ni)] = permutest( sub(iSub).baseCnts(Ef,ni)', sub(iSub).stimCnts(Ef,ni)', true, 0.5, 1000, true );
                                [~,sub(iSub).pnFCnt(ni)] = permutest( sub(iSub).baseCnts(Nf,ni)', sub(iSub).stimCnts(Nf,ni)', true, 0.5, 1000, true );
        
                                % Used later in Fig. 2E - ANOVA@ - Remembered v Forgotten
                                emoneu = char(ones(120,1)*100);
                                emoneu(onsets.emo_all) = 'e';
                                emoneu(onsets.neu_all) = 'n';
        
                                RemFor = char(ones(120,1)*100); % Initialize with d to be removed later
                                RemFor(onsets.eR) = 'R'; % Remembered
                                RemFor(onsets.nR) = 'R';
                                RemFor(onsets.eF) = 'F'; % Forgotten
                                RemFor(onsets.nF) = 'F';
                                di = RemFor == 'd'; % Indices to remove from dataset
                                stmCnts = sub(iSub).stimCnts(:,ni);
                                RemFor = RemFor(~di);
                                stmCnts = stmCnts(~di);
                                emoneu = emoneu(~di);
                                [anovaRes{iSub}.p(:,ni), ~, ~] = anovan(stmCnts,{emoneu RemFor},"Model","interaction","Varnames",["emo","mem"],"display","off");
        
                                % ANOVA@ - Remembered v Know v Forgotten
                                emoneuRKF = char(ones(120,1)*100);
                                emoneuRKF(onsets.emo_all) = 'e';
                                emoneuRKF(onsets.neu_all) = 'n';
        
                                RemForRKF = char(ones(120,1)*100); % Initialize with d to be removed later
                                RemForRKF(onsets.eR) = 'R';
                                RemForRKF(onsets.nR) = 'R';
                                RemForRKF(onsets.eF) = 'F';
                                RemForRKF(onsets.nF) = 'F';
                                RemForRKF(onsets.eK) = 'K';
                                RemForRKF(onsets.nK) = 'K';
                                di = RemForRKF == 'd'; % Indices to remove from dataset
                                stmCntsRKF = sub(iSub).stimCnts(:,ni);
                                RemForRKF = RemForRKF(~di);
                                stmCntsRKF = stmCntsRKF(~di);
                                emoneuRKF = emoneuRKF(~di);
                                [anovaResRKF{iSub}.p(:,ni), ~, ~] = anovan(stmCntsRKF,{emoneuRKF RemForRKF},"Model","interaction","Varnames",["emo","mem"],"display","off");
        
                                % ANOVA - Remembered v Know
                                emoneuRK = char(ones(240,1)*100);
                                emoneuRK(onsets.emo_all) = 'e';
                                emoneuRK(onsets.neu_all) = 'n';
        
                                RemForRK = char(ones(240,1)*100);
                                RemForRK(onsets.eR) = 'R'; % Remember
                                RemForRK(onsets.nR) = 'R';
                                RemForRK(onsets.eK) = 'K'; % Know
                                RemForRK(onsets.nK) = 'K';
                                di = RemForRK == 'd'; % Indices to remove from dataset (eg FAs and Know choices)
                                stmCntsRK = sub(iSub).stimCnts(:,ni); % Spike counts during stimulus periods
                                RemForRK = RemForRK(~di); % Remove them
                                stmCntsRK = stmCntsRK(~di);
                                emoneuRK = emoneuRK(~di);
                                [anovaResRK{iSub}.p(:,ni), ~, ~] = anovan(stmCntsRK,{emoneuRK RemForRK},"Model","interaction","Varnames",["emo","mem"],"display","off");
        
                                % ANOVA@ - Know v Forgotten
                                emoneuKF = char(ones(120,1)*100);
                                emoneuKF(onsets.emo_all) = 'e';
                                emoneuKF(onsets.neu_all) = 'n';
                                
                                RemForKF = char(ones(120,1)*100); % Initialize with d to be removed later
                                RemForKF(onsets.eF) = 'F';
                                RemForKF(onsets.nF) = 'F';
                                RemForKF(onsets.eK) = 'K';
                                RemForKF(onsets.nK) = 'K';
                                di = RemForKF == 'd'; % Indices to remove from dataset
                                stmCntsKF = sub(iSub).stimCnts(:,ni);
                                RemForKF = RemForKF(~di);
                                stmCntsKF = stmCntsKF(~di);
                                emoneuKF = emoneuKF(~di);
                                [anovaResKF{iSub}.p(:,ni),~ ,~] = anovan(stmCntsKF,{emoneuKF RemForKF},"Model","interaction","Varnames",["emo","mem"],"display","off");
    
                                combinedRKF = strcat(emoneuRKF, RemForRKF);
                                combinedRK = strcat(emoneuRK, RemForRK);
                                combinedKF = strcat(emoneuKF, RemForKF);
    
                            elseif strcmp(task,'rec')
                                % Used later in Fig 2B
                                [~, sub(iSub).peMCnt(ni)] = permutest( sub(iSub).baseCnts(Em,ni)', sub(iSub).stimCnts(Em,ni)', true, 0.5, 1000, true );
                                [~, sub(iSub).pnMCnt(ni)] = permutest( sub(iSub).baseCnts(Nm,ni)', sub(iSub).stimCnts(Nm,ni)', true, 0.5, 1000, true );
                                [~, sub(iSub).peCRCnt(ni)] = permutest( sub(iSub).baseCnts(Ecr,ni)', sub(iSub).stimCnts(Ecr,ni)', true, 0.5, 1000, true );
                                [~, sub(iSub).pnCRCnt(ni)] = permutest( sub(iSub).baseCnts(Ncr,ni)', sub(iSub).stimCnts(Ncr,ni)', true, 0.5, 1000, true );
        
                                % ANOVA@ for each neuron in main text - RHit v Miss v CR
                                emoneu = char(ones(240,1)*100);
                                emoneu(onsets.emo_all) = 'e';
                                emoneu(onsets.neu_all) = 'n';
        
                                RemFor = char(ones(240,1)*100);
                                RemFor(onsets.eRHit) = 'R'; % Remember
                                RemFor(onsets.nRHit) = 'R';
                                RemFor(onsets.eCR) = 'C'; % Correct Rejections
                                RemFor(onsets.nCR) = 'C';
                                RemFor(onsets.eMiss) = 'M'; % Misses
                                RemFor(onsets.nMiss) = 'M';
                                di = RemFor == 'd'; % Indices to remove from dataset (eg FAs and Know choices)
                                stmCnts = sub(iSub).stimCnts(:,ni); % Spike counts during stimulus periods
                                RemFor = RemFor(~di); % Remove them
                                stmCnts = stmCnts(~di);
                                emoneu = emoneu(~di);
                                [anovaRes{iSub}.p(:,ni),~,~] = anovan(stmCnts,{emoneu RemFor},"Model","interaction","Varnames",["emo","mem"],"display","off");
                                   
                                % ANOVA - RHit v KHit v Miss v CR
                                emoneuRCMK = char(ones(240,1)*100);
                                emoneuRCMK(onsets.emo_all) = 'e';
                                emoneuRCMK(onsets.neu_all) = 'n';
        
                                RemForRCMK = char(ones(240,1)*100);
                                RemForRCMK(onsets.eRHit) = 'R'; % Remember
                                RemForRCMK(onsets.nRHit) = 'R';
                                RemForRCMK(onsets.eCR) = 'C'; % Correct Rejections
                                RemForRCMK(onsets.nCR) = 'C';
                                RemForRCMK(onsets.eMiss) = 'M'; % Misses
                                RemForRCMK(onsets.nMiss) = 'M';
                                RemForRCMK(onsets.eKHit) = 'K'; % Know
                                RemForRCMK(onsets.nKHit) = 'K';
                                di = RemForRCMK == 'd'; % Indices to remove from dataset (eg FAs and Know choices)
                                stmCntsRCMK = sub(iSub).stimCnts(:,ni); % Spike counts during stimulus periods
                                RemForRCMK = RemForRCMK(~di); % Remove them
                                stmCntsRCMK = stmCntsRCMK(~di);
                                emoneuRCMK = emoneuRCMK(~di);
                                [anovaResRCMK{iSub}.p(:,ni),~,~] = anovan(stmCntsRCMK,{emoneuRCMK RemForRCMK},"Model","interaction","Varnames",["emo","mem"],"display","off");
        
                                % ANOVA - RHit v Miss
                                emoneuRM = char(ones(240,1)*100);
                                emoneuRM(onsets.emo_all) = 'e';
                                emoneuRM(onsets.neu_all) = 'n';
        
                                RemForRM = char(ones(240,1)*100);
                                RemForRM(onsets.eRHit) = 'R'; % Remember
                                RemForRM(onsets.nRHit) = 'R';
                                RemForRM(onsets.eMiss) = 'M'; % Misses
                                RemForRM(onsets.nMiss) = 'M';
                                di = RemForRM == 'd'; % Indices to remove from dataset (eg FAs and Know choices)
                                stmCntsRM = sub(iSub).stimCnts(:,ni); % Spike counts during stimulus periods
                                RemForRM = RemForRM(~di); % Remove them
                                stmCntsRM = stmCntsRM(~di);
                                emoneuRM = emoneuRM(~di);
                                [anovaResRM{iSub}.p(:,ni),~,~] = anovan(stmCntsRM,{emoneuRM RemForRM},"Model","interaction","Varnames",["emo","mem"],"display","off");
    
                                % ANOVA - RHit v CR
                                emoneuRC = char(ones(240,1)*100);
                                emoneuRC(onsets.emo_all) = 'e';
                                emoneuRC(onsets.neu_all) = 'n';
        
                                RemForRC = char(ones(240,1)*100);
                                RemForRC(onsets.eRHit) = 'R'; % Remember
                                RemForRC(onsets.nRHit) = 'R';
                                RemForRC(onsets.eCR) = 'C'; % Correct Rejections
                                RemForRC(onsets.nCR) = 'C';
                                di = RemForRC == 'd'; % Indices to remove from dataset (eg FAs and Know choices)
                                stmCntsRC = sub(iSub).stimCnts(:,ni); % Spike counts during stimulus periods
                                RemForRC = RemForRC(~di); % Remove them
                                stmCntsRC = stmCntsRC(~di);
                                emoneuRC = emoneuRC(~di);
                                [anovaResRC{iSub}.p(:,ni),~,~] = anovan(stmCntsRC,{emoneuRC RemForRC},"Model","interaction","Varnames",["emo","mem"],"display","off");
        
                                % ANOVA - RHit v KHit
                                emoneuRK = char(ones(240,1)*100);
                                emoneuRK(onsets.emo_all) = 'e';
                                emoneuRK(onsets.neu_all) = 'n';
        
                                RemForRK = char(ones(240,1)*100);
                                RemForRK(onsets.eRHit) = 'R'; % Remember
                                RemForRK(onsets.nRHit) = 'R';
                                RemForRK(onsets.eKHit) = 'K'; % Know
                                RemForRK(onsets.nKHit) = 'K';
                                di = RemForRK == 'd'; % Indices to remove from dataset (eg FAs and Know choices)
                                stmCntsRK = sub(iSub).stimCnts(:,ni); % Spike counts during stimulus periods
                                RemForRK = RemForRK(~di); % Remove them
                                stmCntsRK = stmCntsRK(~di);
                                emoneuRK = emoneuRK(~di);
                                [anovaResRK{iSub}.p(:,ni),~,~] = anovan(stmCntsRK,{emoneuRK RemForRK},"Model","interaction","Varnames",["emo","mem"],"display","off");
        
                                % ANOVA@ - KHit v Miss
                                emoneuMK = char(ones(120,1)*100);
                                emoneuMK(onsets.emo_all) = 'e';
                                emoneuMK(onsets.neu_all) = 'n';
        
                                RemForMK = char(ones(120,1)*100); % Initialize with d to be removed later
                                RemForMK(onsets.eMiss) = 'M';
                                RemForMK(onsets.nMiss) = 'M';
                                RemForMK(onsets.eKHit) = 'K';
                                RemForMK(onsets.nKHit) = 'K';
                                di = RemForMK == 'd'; % Indices to remove from dataset
                                stmCntsMK = sub(iSub).stimCnts(:,ni);
                                RemForMK = RemForMK(~di);
                                stmCntsMK = stmCntsMK(~di);
                                emoneuMK = emoneuMK(~di);
                                [anovaResMK{iSub}.p(:,ni),~,~] = anovan(stmCntsMK,{emoneuMK RemForMK},"Model","interaction","Varnames",["emo","mem"],"display","off");
        
                                % ANOVA@ - CR v K
                                emoneuCK = char(ones(120,1)*100);
                                emoneuCK(onsets.emo_all) = 'e';
                                emoneuCK(onsets.neu_all) = 'n';
        
                                RemForCK = char(ones(120,1)*100); % Initialize with d to be removed later
                                RemForCK(onsets.eCR) = 'C';
                                RemForCK(onsets.nCR) = 'C';
                                RemForCK(onsets.eKHit) = 'K';
                                RemForCK(onsets.nKHit) = 'K';
                                di = RemForCK == 'd'; % Indices to remove from dataset
                                stmCntsCK = sub(iSub).stimCnts(:,ni);
                                RemForCK = RemForCK(~di);
                                stmCntsCK = stmCntsCK(~di);
                                emoneuCK = emoneuCK(~di);
                                [anovaResCK{iSub}.p(:,ni),~,~] = anovan(stmCntsCK,{emoneuCK RemForCK},"Model","interaction","Varnames",["emo","mem"],"display","off");
        
                                % ANOVA@ - CR v Miss
                                emoneuCM = char(ones(120,1)*100);
                                emoneuCM(onsets.emo_all) = 'e';
                                emoneuCM(onsets.neu_all) = 'n';
        
                                RemForCM = char(ones(120,1)*100); % Initialize with d to be removed later
                                RemForCM(onsets.eCR) = 'C';
                                RemForCM(onsets.nCR) = 'C';
                                RemForCM(onsets.eMiss) = 'M';
                                RemForCM(onsets.nMiss) = 'M';
                                di = RemForCM == 'd'; % Indices to remove from dataset
                                stmCntsCM = sub(iSub).stimCnts(:,ni);
                                RemForCM = RemForCM(~di);
                                stmCntsCM = stmCntsCM(~di);
                                emoneuCM = emoneuCM(~di);
                                [anovaResCM{iSub}.p(:,ni),~,~] = anovan(stmCntsCM,{emoneuCM RemForCM},"Model","interaction","Varnames",["emo","mem"],"display","off");
    
                                combinedRCMK = strcat(emoneuRCMK,RemForRCMK);
                                combinedRK = strcat(emoneuRK,RemForRK);
                                combinedRC = strcat(emoneuRC,RemForRC);
                                combinedRM = strcat(emoneuRM,RemForRM);
                                combinedCK = strcat(emoneuCK,RemForCK);
                                combinedCM = strcat(emoneuCM,RemForCM);
                                combinedMK = strcat(emoneuMK,RemForMK);
                            end
        
                            % ANOVA@: Computed shuffled distributions by maintaining number of trials per condition
                            combined = strcat(emoneu,RemFor); % Combine to maintain numbers of each trial type
        
                            for iu = 1:nShuffAnova % 10,000 in manuscript
                                anovaRes{iSub}.p_sh(ni,iu,:) = randANOVA(combined,stmCnts);
                                if strcmp(task,'enc')
                                    anovaResRKF{iSub}.p_sh(ni,iu,:) = randANOVA(combinedRKF,stmCntsRKF);
                                    anovaResRK{iSub}.p_sh(ni,iu,:) = randANOVA(combinedRK,stmCntsRK);
                                    anovaResKF{iSub}.p_sh(ni,iu,:) = randANOVA(combinedKF,stmCntsKF);
                                elseif strcmp(task,'rec')
                                    anovaResRCMK{iSub}.p_sh(ni,iu,:) = randANOVA(combinedRCMK,stmCntsRCMK);
                                    anovaResRM{iSub}.p_sh(ni,iu,:) = randANOVA(combinedRM,stmCntsRM);
                                    anovaResRC{iSub}.p_sh(ni,iu,:) = randANOVA(combinedRC,stmCntsRC);
                                    anovaResRK{iSub}.p_sh(ni,iu,:) = randANOVA(combinedRK,stmCntsRK);
                                    anovaResCK{iSub}.p_sh(ni,iu,:) = randANOVA(combinedCK,stmCntsCK);
                                    anovaResCM{iSub}.p_sh(ni,iu,:) = randANOVA(combinedCM,stmCntsCM);
                                    anovaResMK{iSub}.p_sh(ni,iu,:) = randANOVA(combinedMK,stmCntsMK);
                                end
                            end
                            
                            % This PSTH uses 0.01 bin size as input to dPCA later.
                            cfg             = [];
                            cfg.binsize     =  psthbin; % if cfgPsth.binsize = 'scott' or 'sqrt', we estimate the optimal bin size from the data itself
                            cfg.outputunit  = 'rate'; % give as an output the firing rate
                            cfg.latency     = start_stop; % between -1 and 3 sec.
                            cfg.vartriallen = 'no'; % variable trial lengths are accepted
                            cfg.keeptrials  = 'yes'; % keep the psth per trial in the output
                            psth = ft_spike_psth(cfg,spikeTrials);
        
                            psth_trial = squeeze(psth.trial(:,k,:));
                            
                            % Gaussian kernel filter for each trial
                            x = -1000:1000;
                            gaussKernel = 1/sqrt(2*pi)/20 * exp(-x.^2/20^2/2);
        
                            psth_trial_norm = zeros(size(psth_trial));
                            for tii = 1:size(psth_trial,1)
                                psth_trial_norm(tii,:) = conv(psth_trial(tii,:), gaussKernel, 'same');
                            end
                            psth_use = psth_trial_norm; % set psth_use = psth_trial to compute all results without the Gaussian kernel
                            psth_use_time = psth.time;
        
                            % Paired t-test to determine significant difference after v before image
                            baseTime = psth_use_time > -th & psth_use_time < -tl; % Baseline time symmetrical to event time
                            evtTime = psth_use_time < th & psth_use_time > tl;
                            allAvg = mean(psth_use);
        
                            % With Gaussian kernel: Used for dPCA
                            sub(iSub).baseRate(:,ni) = mean(psth_use(:,baseTime),2);
                            sub(iSub).stimRate(:,ni) = mean(psth_use(:,evtTime),2);
                            
                            sub(iSub).allZ(:,ni) = (sub(iSub).stimRate(:,ni) - mean(sub(iSub).baseRate(:,ni))) / std(sub(iSub).baseRate(:,ni)); 
                            if mean(sub(iSub).allZ(:,ni)) < 0 % If the overall z-score is positive, flip it to negative
                                sub(iSub).absallZ(:,ni) = -1*sub(iSub).allZ(:,ni);
                            else
                                sub(iSub).absallZ(:,ni) = sub(iSub).allZ(:,ni);
                            end

                            % Without Gaussian kernel: Requested by reviewer for LMM
                            sub(iSub).rawbaseRate(:,ni) = mean(psth_trial(:,baseTime),2); % firing rate without Gaussian kernel
                            sub(iSub).rawstimRate(:,ni) = mean(psth_trial(:,evtTime),2);

                            sub(iSub).rawallZ(:,ni) = (sub(iSub).rawstimRate(:,ni) - mean(sub(iSub).rawbaseRate(:,ni))) / std(sub(iSub).rawbaseRate(:,ni)); 
                            if mean(sub(iSub).rawallZ(:,ni)) < 0
                                sub(iSub).rawabsallZ(:,ni) = -1*sub(iSub).rawallZ(:,ni);
                            else
                                sub(iSub).rawabsallZ(:,ni) = sub(iSub).rawallZ(:,ni);
                            end
        
                            % Used later for Fig 2A
                            eAvg = mean(psth_use(Eall,:));
                            nAvg = mean(psth_use(Nall,:));
    
                            % Save z-scores to determine direction of firing rate change from baseline to stimulus period
                            sub(iSub).ZEmo(ni) = (mean(eAvg(evtTime)) - mean(allAvg(baseTime))) / std(sub(iSub).baseRate(:,ni));
                            sub(iSub).ZNeu(ni) = (mean(nAvg(evtTime)) - mean(allAvg(baseTime))) / std(sub(iSub).baseRate(:,ni));
        
                            % Used for plotting neuronal examples in Fig. 1
                            eRAvg = mean(psth_use(Er,:)); 
                            eRsem = std(psth_use(Er,:)) / sqrt(length(Er));
                            sub(iSub).eRAvg(ni,:) = eRAvg;
                            sub(iSub).ZeR(ni) = (mean(eRAvg(evtTime)) - mean(allAvg(baseTime))) / std(sub(iSub).baseRate(:,ni));
        
                            nRAvg = mean(psth_use(Nr,:));
                            nRsem = std(psth_use(Nr,:)) / sqrt(length(Nr));
                            sub(iSub).nRAvg(ni,:) = nRAvg;
                            sub(iSub).ZnR(ni) = (mean(nRAvg(evtTime)) - mean(allAvg(baseTime))) / std(sub(iSub).baseRate(:,ni));
        
                            if strcmp(task,'enc')
    
                                eFAvg = mean(psth_use(Ef,:));
                                eFsem = std(psth_use(Ef,:)) / sqrt(length(Ef));
                                sub(iSub).ZeF(ni) = (mean(eFAvg(evtTime)) - mean(allAvg(baseTime))) / std(sub(iSub).baseRate(:,ni));
        
                                nFAvg = mean(psth_use(Nf,:));
                                nFsem = std(psth_use(Nf,:)) / sqrt(length(Nf));
                                sub(iSub).ZnF(ni) = (mean(nFAvg(evtTime)) - mean(allAvg(baseTime))) / std(sub(iSub).baseRate(:,ni));
        
                            elseif strcmp(task,'rec')

                                eMAvg = mean(psth_use(Em,:));
                                nMAvg = mean(psth_use(Nm,:));
                                eMsem = std(psth_use(Em,:)) / sqrt(length(Em));
                                nMsem = std(psth_use(Nm,:)) / sqrt(length(Nm));
                                sub(iSub).ZeM(ni) = (mean(eMAvg(evtTime)) - mean(allAvg(baseTime))) / std(sub(iSub).baseRate(:,ni));
                                sub(iSub).ZnM(ni) = (mean(nMAvg(evtTime)) - mean(allAvg(baseTime))) / std(sub(iSub).baseRate(:,ni));
    
                                eCR_Avg = mean(psth_use(Ecr,:));
                                nCR_Avg = mean(psth_use(Ncr,:));
                                sub(iSub).ZeCR(ni) = (mean(eCR_Avg(evtTime)) - mean(allAvg(baseTime))) / std(sub(iSub).baseRate(:,ni));
                                sub(iSub).ZnCR(ni) = (mean(nCR_Avg(evtTime)) - mean(allAvg(baseTime))) / std(sub(iSub).baseRate(:,ni));
    
                            end
        
                            % Save data in the structure for dPCA
                            dPCA_neur_name{sp} = [subjects{iSub}, '_', spike.label{k}];
                            if strcmp(task,'enc')
    
                                % All 3 trials types RKF
                                fr = 1; % Emotional remember
                                dec = 0;
                                repetitions = psth_use(onsets.eR,:);
                                firingRatesAverageRKF(sp, dec+1, fr, 1:size(repetitions,2)) = mean(repetitions, 1); 
                                trialNumRKF(sp, dec+1, fr) = size(repetitions, 1);
                                firingRatesRKF(sp, dec+1, fr, 1:size(repetitions,2), 1:size(repetitions,1)) = repetitions'; 
                                
                                fr = 2; % Neutral Remember
                                dec = 0;
                                repetitions = psth_use(onsets.nR,:);
                                firingRatesAverageRKF(sp, dec+1, fr, 1:size(repetitions,2)) = mean(repetitions, 1);
                                trialNumRKF(sp, dec+1, fr) = size(repetitions, 1);
                                firingRatesRKF(sp, dec+1, fr, 1:size(repetitions,2), 1:size(repetitions,1)) = repetitions';
        
                                fr = 1; % Emotional forgot
                                dec = 1;
                                repetitions = psth_use(onsets.eF,:);
                                firingRatesAverageRKF(sp, dec+1, fr, 1:size(repetitions,2)) = mean(repetitions, 1);
                                trialNumRKF(sp, dec+1, fr) = size(repetitions, 1);
                                firingRatesRKF(sp, dec+1, fr, 1:size(repetitions,2), 1:size(repetitions,1)) = repetitions';
        
                                fr = 2; % Neutral forgot
                                dec = 1;
                                repetitions = psth_use(onsets.nF,:);
                                firingRatesAverageRKF(sp, dec+1, fr, 1:size(repetitions,2)) = mean(repetitions, 1); 
                                trialNumRKF(sp, dec+1, fr) = size(repetitions, 1);
                                firingRatesRKF(sp, dec+1, fr, 1:size(repetitions,2), 1:size(repetitions,1)) = repetitions';
        
                                fr = 1; % Emotional K
                                dec = 2;
                                repetitions = psth_use(onsets.eK,:);
                                firingRatesAverageRKF(sp, dec+1, fr, 1:size(repetitions,2)) = mean(repetitions, 1);
                                trialNumRKF(sp, dec+1, fr) = size(repetitions, 1);
                                firingRatesRKF(sp, dec+1, fr, 1:size(repetitions,2), 1:size(repetitions,1)) = repetitions';
        
                                fr = 2; % Neutral K
                                dec = 2;
                                repetitions = psth_use(onsets.nK,:);
                                firingRatesAverageRKF(sp, dec+1, fr, 1:size(repetitions,2)) = mean(repetitions, 1);
                                trialNumRKF(sp, dec+1, fr) = size(repetitions, 1);
                                firingRatesRKF(sp, dec+1, fr, 1:size(repetitions,2), 1:size(repetitions,1)) = repetitions';
        
                            elseif strcmp(task,'rec')
    
                                for d=1:2 % decision is either emotional or neutral
                                    for f=1:length(outsmall)
                                        if d == 1
                                            repetitions = psth_use(eval(['onsets.e',outsmall{f}]),:);
                                        elseif d == 2
                                            repetitions = psth_use(eval(['onsets.n',outsmall{f}]),:);
                                        end
                                        firingRatesAverageRec(sp, f, d, 1:size(repetitions,2)) = mean(repetitions, 1); 
                                        trialNumRec(sp, f, d) = size(repetitions, 1);
                                        firingRatesRec(sp, f, d, 1:size(repetitions,2), 1:size(repetitions,1)) = repetitions';
                                    end
                                end
                            end
    
                            if strcmp(task,'enc')
                                % Get the trial indicies to plot sorted rasters
                                [~,er] = ismember(spikeTrials.trial{k}(eR), Er); 
                                [~,ef] = ismember(spikeTrials.trial{k}(eF), Ef);
                                [~,nr] = ismember(spikeTrials.trial{k}(nR), Nr);
                                [~,nf] = ismember(spikeTrials.trial{k}(nF), Nf);

                                curr_fig = figure('Visible','Off'); % Great way to keep figures invisible!
                                set(curr_fig, 'Units', 'centimeters', 'Position', [0 0 4 8]);
                                subtightplot(2,1,1, 0.01,0.077,0.22), hold on
                                set(gca,'Ydir','reverse');
                                scatter(spikeTrials.time{k}(eR), er, 5,rgb('FireBrick'),'|','MarkerEdgeAlpha',0.7)
                                scatter(spikeTrials.time{k}(eF), ef+length(Er), 5,rgb('OrangeRed'),'|','MarkerEdgeAlpha',0.7)
                                scatter(spikeTrials.time{k}(nR), nr+length(Ef)+length(Er), 5,rgb('Navy'),'|','MarkerEdgeAlpha',0.7)
                                scatter(spikeTrials.time{k}(nF), nf+length(Nr)+length(Ef)+length(Er), 5,rgb('Turquoise'),'|','MarkerEdgeAlpha',0.7)

                                if contains(spikeTrials.label{k},patahp)
                                    title(regNames(1))
                                    spikeTrials.region(k) = regNames(1);
                                elseif contains(spikeTrials.label{k},patphp)
                                    title(regNames(2))
                                    spikeTrials.region(k) = regNames(2);
                                elseif contains(spikeTrials.label{k},patamg)
                                    title(regNames(3))
                                    spikeTrials.region(k) = regNames(3);
                                elseif contains(spikeTrials.label{k},patec)
                                    title(regNames(4))
                                    spikeTrials.region(k) = regNames(4);
                                elseif contains(spikeTrials.label{k},patph)
                                    title(regNames(5))
                                    spikeTrials.region(k) = regNames(5);
                                end
                                lastT = max(length(Nf)+length(Nr)+length(Ef)+length(Er));
                                ylim([0 lastT+1]);
                                ylabel('Sorted Trial #');
                                xticklabels([]); yticks([1 50 100]); set(gca,'YMinorTick','on');
                                xlim([-2 2]); xticks([-1.5 0 1.5]);

                                subtightplot(2,1,2, 0.01,0.077,0.22), hold on;

                                % Plot twice for legend
                                plot(psth_use_time, eRAvg,'Color',rgb('FireBrick'),'LineWidth',0.75)
                                plot(psth_use_time, eFAvg,'Color',rgb('OrangeRed'),'LineWidth',0.75)
                                plot(psth_use_time, nRAvg,'Color',rgb('Navy'),'LineWidth',0.75)
                                plot(psth_use_time, nFAvg,'Color',rgb('Turquoise'),'LineWidth',0.75)

                                boundedline(psth_use_time, eRAvg, eRsem,'Color',rgb('FireBrick'),'LineWidth',0.75,'alpha')
                                boundedline(psth_use_time, eFAvg, eFsem,'Color',rgb('OrangeRed'),'LineWidth',0.75,'alpha')
                                boundedline(psth_use_time, nRAvg, nRsem, 'Color',rgb('Navy'),'LineWidth',0.75,'alpha')
                                boundedline(psth_use_time, nFAvg, nFsem, 'Color',rgb('Turquoise'),'LineWidth',0.75,'alpha')

                                y3 = ylim();
                                patch([tl tl, th th], [y3(1) y3(2) y3(2) y3(1)], [0.9 0.9 0.9],'LineStyle','none') % stimulus period shading
                                alpha(0.25)
                                patch([-th -th, -tl -tl], [y3(1) y3(2) y3(2) y3(1)], [0.9 0.9 0.9],'LineStyle','none') % baseline period shading
                                alpha(0.25)
                                set(gca,'Ylim',[0 y3(2)]); % Set the bottom y-axis to zero and keep the top automatic

                                legend({'eR','eF','nR', 'nF'},'Location','northwest','AutoUpdate','off');
                                legend boxoff
        
                                ylabel('Firing Rate (Hz)');
                                xlabel('Time (sec)');
                                xlim([-2 2]); xticks([-1.5 0 1.5]);
                                fontsize(curr_fig, 8, 'points');
                                fontname(curr_fig, "Arial");
                                z = num2str((anovaRes{iSub}.p(:,end) < 0.05)', '%d');
                                exportgraphics(curr_fig,strcat(paths.pics,subjects{iSub,1}, '_enc_manuscript_sem_', spike.label{k},'_',num2str(k),'_',z,'.png'));
                                saveas(curr_fig,strcat(paths.pics,subjects{iSub,1}, '_enc_manuscript_sem_', spike.label{k},'_',num2str(k),'_',z,'.svg'));
                                close(curr_fig)

                            elseif strcmp(task,'rec')

                                % Get the trial indicies to plot sorted rasters
                                [~,er] = ismember(spikeTrials.trial{k}(eR), Er); 
                                [~,em] = ismember(spikeTrials.trial{k}(eM), Em);
                                [~,ecr] = ismember(spikeTrials.trial{k}(eCR), Ecr);
                                [~,nr] = ismember(spikeTrials.trial{k}(nR), Nr);
                                [~,nm] = ismember(spikeTrials.trial{k}(nM), Nm);
                                [~,ncr] = ismember(spikeTrials.trial{k}(nCR), Ncr);

                                curr_fig = figure('Visible','Off'); % Make figures invisible to speed up the code.
                                set(curr_fig, 'Units', 'centimeters', 'Position', [0 0 4 8]);
                                subtightplot(2, 1, 1, 0.01, 0.077, 0.22), hold on;
                                set(gca,'Ydir', 'reverse');
                                scatter(spikeTrials.time{k}(eR), er, 5,rgb('FireBrick'),'|','MarkerEdgeAlpha',0.7)
                                scatter(spikeTrials.time{k}(eM), em+length(Er), 5, rgb('OrangeRed'),'|','MarkerEdgeAlpha',0.7)
                                scatter(spikeTrials.time{k}(nR), nr+length(Em)+length(Er), 5,rgb('Navy'),'|','MarkerEdgeAlpha',0.7)
                                scatter(spikeTrials.time{k}(nM), nm+length(Nr)+length(Em)+length(Er), 5,rgb('DeepSkyBlue'),'|','MarkerEdgeAlpha',0.7)

                                lastT = max(length(Nm)+length(Nr)+length(Em)+length(Er)); % Last trial index for axis limit
                                xlim([-2 2]); ylim([0 lastT+1]); xticklabels([]); 
                                yticks([1 50 100 150 200]); set(gca,'YMinorTick','on');
                                ylabel('Sorted Trial #');
                                if contains(spikeTrials.label{k},patahp)
                                    title(regNames(1))
                                    spikeTrials.region(k) = regNames(1);
                                elseif contains(spikeTrials.label{k},patphp)
                                    title(regNames(2))
                                    spikeTrials.region(k) = regNames(2);
                                elseif contains(spikeTrials.label{k},patamg)
                                    title(regNames(3))
                                    spikeTrials.region(k) = regNames(3);
                                elseif contains(spikeTrials.label{k},patec)
                                    title(regNames(4))
                                    spikeTrials.region(k) = regNames(4);
                                elseif contains(spikeTrials.label{k},patph)
                                    title(regNames(5))
                                    spikeTrials.region(k) = regNames(5);
                                end

                                ax4=subtightplot(2,1,2, 0.01,0.077,0.22); hold on;
                                % Plot each twice for legend
                                plot(psth_use_time,eRAvg,'Color',rgb('FireBrick'),'LineWidth',0.75)
                                plot(psth_use_time,eMAvg,'Color',rgb('OrangeRed'),'LineWidth',0.75)
                                plot(psth_use_time,nRAvg,'Color',rgb('Navy'),'LineWidth',0.75)
                                plot(psth_use_time,nMAvg,'Color',rgb('DeepSkyBlue'),'LineWidth',0.75)

                                boundedline(psth_use_time, eRAvg, eRsem,'Color',rgb('FireBrick'),'LineWidth',0.75,'alpha')
                                boundedline(psth_use_time,eMAvg, eMsem,'Color',rgb('OrangeRed'),'LineWidth',0.75,'alpha')
                                boundedline(psth_use_time,nRAvg, nRsem, 'Color',rgb('Navy'),'LineWidth',0.75,'alpha')
                                boundedline(psth_use_time,nMAvg, nMsem, 'Color',rgb('DeepSkyBlue'),'LineWidth',0.75,'alpha')
                                y3 = ylim();

                                legend({'eR','eM', 'nR', 'nM'},'Location','northwest','AutoUpdate','off');
                                legend boxoff
                                ylabel('Firing Rate (Hz)'); 
                                xlabel('Time (sec)'); xlim([-2 2]); xticks([-1.5 0 1.5]); set(gca,'XMinorTick','on'); 
                                ax4.XAxis.MinorTickValues = -2:0.5:2; 
                                patch(ax4,[tl tl, th th], [y3(1) y3(2) y3(2) y3(1)], [0.9 0.9 0.9],'LineStyle','none'); alpha(0.25)
                                patch(ax4,[-th -th, -tl -tl], [y3(1) y3(2) y3(2) y3(1)], [0.9 0.9 0.9],'LineStyle','none'); alpha(0.25)
                                
                                fontsize(curr_fig,8,'points');
                                fontname(curr_fig,"Arial");
                                z = num2str((anovaRes{iSub}.p(:,end) < 0.05)', '%d'); % Is the Anova significant for this neuron? Binary code in order emotional, memory, interaction.
                                exportgraphics(curr_fig,strcat(paths.pics,subjects{iSub,1}, '_rec_manuscript_final_', spike.label{k},'_',num2str(k),'_',z,'.png'));
                                saveas(curr_fig,strcat(paths.pics,subjects{iSub,1}, '_rec_manuscript_final_', spike.label{k},'_',num2str(k),'_',z,'.svg'));
                                close(curr_fig)

                            end % Fig. 1 example neuron plots
                        end % If this spike train meets inclusion criteria
                    end % Iteration over every spike train in the microwire
                end % If there are any neurons in the microwire
            end % If there are any neurons in the microwire
        end % every microwire in a region (should be 8)
    end % every region
end % every subject

%% Create some vectors for regions and subject names
names = string([sub(:).name]);
iAHP = contains(names,'AH') | contains(names,'uH'); % Anterior Hippocampus
iPHP = contains(names,'PH');                        % Posterior Hippocampus (Both are classified as Hippocampus in this manuscript)
iAmg = contains(names,'AR') | contains(names,'AL'); % Amygdala
iEC = contains(names,'EC');                         % Entorhinal Cortex
areaMask = iAHP .* 1 + iPHP .*2 + iAmg .* 3 + iEC .* 4;

subjMask = [];
for iS=1:length(sub)
    subjMask = [subjMask, ones([1,length(sub(iS).pEmoCnt)])*iS ];
end

%% Supp fig 1: Microwire and neuronal statistics

% Find the number of neurons per microwire
for ip = 1:length(dPCA_neur_name)
    dPCA_microwire{ip} = dPCA_neur_name{ip}(1:end-6); % Remove part of the neuron name, keep subject & microwire name
end
[uniqueID, ~, J]=unique(dPCA_microwire);
nWire = histc(J, 1:numel(uniqueID));

fig = figure('Units','centimeters', 'Position', [0 0 20 5]); hold on;
subplot(1,4,1)
histogram(nWire) 
title(sprintf('%d Neurons',length(dPCA_neur_name)))
xlabel('Neurons per wire')
ylabel('# of microwires')

subplot(1,4,2)
histogram([sub(:).perc_less3ms])
xlabel('% ISI < 3ms')
ylabel('# of neurons')

subplot(1,4,3)
histogram([sub(:).meanFR],0:0.5:30)
xlim([0,12])
xlabel('mean FR (Hz)')
ylabel('# of neurons')

subplot(1,4,4)
histogram([sub(:).width])
xlabel('spike width')
ylabel('# of neurons')

fontsize(fig,8,'points');
fontname(fig,"Arial");
saveas(fig, [paths.pics,sprintf('_SFig1_%s_spike_quality.svg',task)]);
exportgraphics(fig,[paths.pics,sprintf('_SFig1_%s_spike_quality.png',task)],'Resolution',300) 
close(fig)

%% Supp Fig 2. 2x2 ANOVA@ results for each neuron
removeDoubleCounts = true; % If a neuron has a significant interaction, set the main effects to not significant. This was performed for the manuscript.
if removeDoubleCounts
    if strcmp(task,'enc')
        anovaRes = removeANOVAdoubleCounts(anovaRes);
        anovaResKF = removeANOVAdoubleCounts(anovaResKF);
        anovaResRK = removeANOVAdoubleCounts(anovaResRK);
        anovaResRKF = removeANOVAdoubleCounts(anovaResRKF);
    
    elseif strcmp(task,'rec')
    
        anovaRes = removeANOVAdoubleCounts(anovaRes);
        anovaResRCMK = removeANOVAdoubleCounts(anovaResRCMK);
        anovaResRM = removeANOVAdoubleCounts(anovaResRM);
        anovaResRC = removeANOVAdoubleCounts(anovaResRC);
        anovaResRK = removeANOVAdoubleCounts(anovaResRK);
        anovaResCK = removeANOVAdoubleCounts(anovaResCK);
        anovaResCM = removeANOVAdoubleCounts(anovaResCM);
        anovaResMK = removeANOVAdoubleCounts(anovaResMK);
    end
end

anova_p = [];
anova_p_sh = [];
for i = 1:length(anovaRes) 
    anova_p = [anova_p, anovaRes{i}.p]; % Save p-values in order: Emotion, Memory, Interaction. 
    anova_p_sh = [anova_p_sh; anovaRes{i}.p_sh]; % Save p-values in order: Emotion, Memory, Interaction. 
end

% Calculate the number of significant neurons in shuffled distributions
anova_cnts_1 = sum(anova_p_sh(:,:,1) < 0.05, 1); % emotion
anova_cnts_2 = sum(anova_p_sh(:,:,2) < 0.05, 1); % memory
anova_cnts_3 = sum(anova_p_sh(:,:,3) < 0.05, 1); % interaction

anova_sig = anova_p < 0.05;
anova_n = sum(anova_sig,2); % NOTE 2x2 for encoding and 3x2 for recognition
% sprintf(['We performed a 2x2 repeated measures ANOVA for each neuron and found %i neurons with a main effect of emotion, %i neurons with a main effect of memory' ...
%     ', and %i neurons with an interaction between memory and emotion.'],anova_n(1), anova_n(2), anova_n(3))

% Compare shuffled distribution to the real ones ANOVA@
p_emo_ANOVA = sum(anova_cnts_1 >= anova_n(1)) / length(anova_cnts_1);
p_mem_ANOVA = sum(anova_cnts_2 >= anova_n(2)) / length(anova_cnts_1);
p_eXm_ANOVA = sum(anova_cnts_3 >= anova_n(3)) / length(anova_cnts_1);

%% ANOVA@ Stacked bar chart
n_anysig = sum(sum(anova_sig,1) > 0); % Does each neuron have any significant value from the ANOVA?
n_e = sum(anova_sig(1,:) & sum(anova_sig,1)==1);
n_m = sum(anova_sig(2,:) & sum(anova_sig,1)==1);
n_x = sum(anova_sig(3,:) & sum(anova_sig,1)==1);
n_em = sum(anova_sig(1,:) & anova_sig(2,:) & sum(anova_sig,1)==2);
n_ex = sum(anova_sig(1,:) & anova_sig(3,:) & sum(anova_sig,1)==2);
n_mx = sum(anova_sig(3,:) & anova_sig(2,:) & sum(anova_sig,1)==2);
n_emx = sum(anova_sig(1,:) & anova_sig(2,:) & anova_sig(3,:) & sum(anova_sig,1)==3);

% Define the total number of neurons
total_neurons = length(anova_sig);

% Calculate percentages
percent_e = n_e / n_anysig * 100;
percent_m = n_m / n_anysig * 100;
percent_x = (n_x+n_ex+n_mx+n_emx) / n_anysig * 100;
percent_em = n_em / n_anysig * 100;

% Create a stacekd bar chart for n_anysig
data_anysig_small = [percent_e, percent_m, percent_em, percent_x];

% Calculate percentages relative to all neurons
percent_e_all = n_e / total_neurons * 100;
percent_m_all = n_m / total_neurons * 100;
percent_x_all = (n_x+n_ex+n_mx+n_emx) / total_neurons * 100;
percent_em_all = n_em / total_neurons * 100;
percent_none = (total_neurons - n_anysig) / total_neurons * 100;

% Create a stacked bar graph for all neurons
data_all_neurons_small = [percent_none, percent_e_all, percent_m_all, percent_em_all, percent_x_all];
n_all_neurons_small = round(data_all_neurons_small * total_neurons / 100);
labels_all_neurons_small = {'n.s.','E', 'M', 'EM', 'X'}; % Emotional, Memory and X=Interaction

for i=1:length(data_all_neurons_small)
    barNamesAll_small{i} = [labels_all_neurons_small{i}, ' (n=', num2str(n_all_neurons_small(i)) ')'];
end

%% Fig 2: Encoding unit summary
if strcmp(task,'enc')

    eR_absZ = [];
    nR_absZ = [];
    eF_absZ = [];
    nF_absZ = [];

    for iS=1:length(sub)
        n_neur(iS) = length(sub(iS).pEmoCnt); % number of neurons

        eR_absZ = [eR_absZ, mean(sub(iS).absallZ(sub(iS).onsets.eR,:),1)];
        nR_absZ = [nR_absZ, mean(sub(iS).absallZ(sub(iS).onsets.nR,:),1)];
        eF_absZ = [eF_absZ, mean(sub(iS).absallZ(sub(iS).onsets.eF,:),1)];
        nF_absZ = [nF_absZ, mean(sub(iS).absallZ(sub(iS).onsets.nF,:),1)];
        
        % Threshold the z-score based on the signifiance test of spike counts between baseline and stimulus periods
        zemo{iS} = sub(iS).ZEmo([sub(iS).pEmoCnt] < 0.05); 
        zneu{iS} = sub(iS).ZNeu([sub(iS).pNeuCnt] < 0.05);
        zeR{iS} = sub(iS).ZeR([sub(iS).peRCnt] < 0.05);
        znR{iS} = sub(iS).ZnR([sub(iS).pnRCnt] < 0.05);
        znF{iS} = sub(iS).ZnF([sub(iS).pnFCnt] < 0.05);
        zeF{iS} = sub(iS).ZeF([sub(iS).peFCnt] < 0.05);
    end

    % Create matrix for bar graphs of Event Responsiveness: If z-score is positive, firing rate increased and vice-versa if negative
    allstacked = [sum([zemo{:}]>0), sum([zneu{:}]>0), sum([zeR{:}]>0), sum([znR{:}]>0), sum([zeF{:}]>0), sum([znF{:}]>0); ...
                  sum([zemo{:}]<0), sum([zneu{:}]<0), sum([zeR{:}]<0), sum([znR{:}]<0), sum([zeF{:}]<0), sum([znF{:}]<0)];

    if p_emo_ANOVA == 0
        p_emo_ANOVA = 1/length(anova_cnts_1);
    end

    fig = figure('Units','centimeters', 'Position', [0 0 9.5 12.5]);
    subplot(321) % Event responsiveness
    bh = bar(allstacked'./sum(n_neur) * 100,'stacked', 'FaceColor','flat');
    bh(1).CData = [0.25 0.25 0.25];  % Change color to first level
    bh(2).CData = [0.75 0.75 0.75];  % Change color to second level
    legend('Increased','Decreased','Location','northeast');
    legend boxoff
    x = {'emo','neu', 'eR', 'nR','eF','nF'};
    xticks([1 2 3 4 5 6]);  xticklabels(x); xtickangle(90)
    xlim([0.5 6.5]); ylim([0 100]); yticks(0:20:100)
    ylabel('% of neurons')
    title('Event responsiveness')

    Y = [data_all_neurons_small;
        0 data_anysig_small];
    subplot(322); % Stacked bar plot
    bar(Y,'stacked');
    xticklabels({'of total','of selective'})
    ylabel('Percentage of neurons')
    lgd = legend(barNamesAll_small,'location','east');
    legend box off;
    ylim([0, 100])
    xlim([0.5 2.5])
    
    % Do we have more neurons with significant ANOVA values than chance levels? Comparing using shuffled trial labels.
    binsV = 1:max(anova_n)+5;
    subplot(334)
    histogram(anova_cnts_1, binsV,'FaceColor',rgb('MidnightBlue')); hold on;
    xline(anova_n(1),'LineWidth',1.5,'Color',rgb('DarkOrange'));
    text(0.1, 0.95, sprintf('p = %0.4g', p_emo_ANOVA), 'Units', 'normalized', 'FontSize', 10);
    ylabel('Number of runs'); title('Emotion'); ylim([0,1400]); xlim([0,30]);
    
    subplot(335)
    histogram(anova_cnts_2, binsV,'FaceColor',rgb('MidnightBlue')); hold on;
    xline(anova_n(2),'LineWidth',1.5,'Color',rgb('DarkOrange'));
    text(0.1, 0.95, sprintf('p = %0.4f', p_mem_ANOVA), 'Units', 'normalized', 'FontSize', 10);
    xlabel('Number of significant neurons'); title('Memory'); ylim([0,1400]); xlim([0,30]);
    
    subplot(336)
    histogram(anova_cnts_3, binsV,'FaceColor',rgb('MidnightBlue')); hold on;
    xline(anova_n(3),'LineWidth',1.5,'Color',rgb('DarkOrange'));
    text(0.1, 0.95, sprintf('p = %0.3f', p_eXm_ANOVA), 'Units', 'normalized', 'FontSize', 10);
    title('Interaction');  ylim([0,1400]); xlim([0,30]);
    
    s3 = subplot(313); hold on;
    x = [1 2 3 4];
    y =  mean([eR_absZ; nR_absZ; eF_absZ; nF_absZ],2);
    err = std([eR_absZ; nR_absZ; eF_absZ; nF_absZ],0,2) / sqrt(length(eR_absZ));
    errorbar(x,y,err,".")
    xlim([0.5 4.5]);
    xl = {'eR', 'nR','eF','nF'};
    xticks([1 2 3 4]);  xticklabels(xl); xtickangle(90)
    ylabel('Absolute Z-Score')
    s3Pos = get(s3,'position');
    s3Pos = s3Pos + [0.015 0 0 0;];
    set(s3,'position',s3Pos);
    
    fontsize(fig,8,'points');
    fontname(fig,"Arial");
    saveas(fig, [paths.pics,sprintf('_Fig2_encoding_unit_summary.svg')]);
    exportgraphics(fig,[paths.pics,sprintf('_Fig2_encoding_unit_summary.png')],'Resolution',300);
    close(fig)

end

%% ANOVA@ per region for SFig 2
if strcmp(task,'enc')
    anova_types = {anovaRes,anovaResRKF};
    anova_names = {'R vs F', 'R vs K vs F'}; % SFig 2A and 2B
elseif strcmp(task,'rec')
    anova_types = {anovaResRCMK};
    anova_names = {'R M C K'}; % SFig 2E
end

binsV = 1:100;
    
for i0 = 1:length(anova_types)
    anova_temp = anova_types{i0};
    anova_p_temp = [];
    anova_p_sh_temp = [];
    for i = 1:length(anova_temp) 
        anova_p_temp = [anova_p_temp, anova_temp{i}.p]; % Save p-values in order: Emotion, Memory, Interaction.
        anova_p_sh_temp = [anova_p_sh_temp; anova_temp{i}.p_sh]; % Save p-values in order: Emotion, Memory, Interaction. 
    end
    
    anova_sig_temp = anova_p_temp < 0.05;
    
    % Are the number of Amygdala neurons significantly above chance? 
    anova_n_amyg_temp = sum(anova_sig_temp(:,areaMask==3),2);
    anova_p_sh_amyg_temp = anova_p_sh_temp(areaMask==3,:,:);

    anova_cnts_1_amyg_temp = sum(anova_p_sh_amyg_temp(:,:,1) < 0.05, 1); % How many were significant?
    anova_cnts_2_amyg_temp = sum(anova_p_sh_amyg_temp(:,:,2) < 0.05, 1);
    anova_cnts_3_amyg_temp = sum(anova_p_sh_amyg_temp(:,:,3) < 0.05, 1);

    p_emo_ANOVA_amyg_temp = sum(anova_cnts_1_amyg_temp >= anova_n_amyg_temp(1)) / length(anova_cnts_1_amyg_temp);
    p_mem_ANOVA_amyg_temp = sum(anova_cnts_2_amyg_temp >= anova_n_amyg_temp(2)) / length(anova_cnts_2_amyg_temp);
    p_eXm_ANOVA_amyg_temp = sum(anova_cnts_3_amyg_temp >= anova_n_amyg_temp(3)) / length(anova_cnts_3_amyg_temp);
    
    % Are the number of EC neurons significantly above chance? 
    anova_n_ec_temp = sum(anova_sig_temp(:,areaMask==4),2);
    anova_p_sh_ec_temp = anova_p_sh_temp(areaMask==4,:,:);

    anova_cnts_1_ec_temp = sum(anova_p_sh_ec_temp(:,:,1) < 0.05, 1);
    anova_cnts_2_ec_temp = sum(anova_p_sh_ec_temp(:,:,2) < 0.05, 1);
    anova_cnts_3_ec_temp = sum(anova_p_sh_ec_temp(:,:,3) < 0.05, 1);

    p_emo_ANOVA_ec_temp = sum(anova_cnts_1_ec_temp >= anova_n_ec_temp(1)) / length(anova_cnts_1_ec_temp);
    p_mem_ANOVA_ec_temp = sum(anova_cnts_2_ec_temp >= anova_n_ec_temp(2)) / length(anova_cnts_2_ec_temp);
    p_eXm_ANOVA_ec_temp = sum(anova_cnts_3_ec_temp >= anova_n_ec_temp(3)) / length(anova_cnts_3_ec_temp);
    
    % Are the number of HPC neurons significantly above chance? 
    anova_n_hpc_temp = sum(anova_sig_temp(:,areaMask==1 | areaMask==2),2);
    anova_p_sh_hpc_temp = anova_p_sh_temp(areaMask==1 | areaMask==2,:,:);

    anova_cnts_1_hpc_temp = sum(anova_p_sh_hpc_temp(:,:,1) < 0.05, 1);
    anova_cnts_2_hpc_temp = sum(anova_p_sh_hpc_temp(:,:,2) < 0.05, 1);
    anova_cnts_3_hpc_temp = sum(anova_p_sh_hpc_temp(:,:,3) < 0.05, 1);

    p_emo_ANOVA_hpc_temp = sum(anova_cnts_1_hpc_temp >= anova_n_hpc_temp(1)) / length(anova_cnts_1_hpc_temp);
    p_mem_ANOVA_hpc_temp = sum(anova_cnts_2_hpc_temp >= anova_n_hpc_temp(2)) / length(anova_cnts_2_hpc_temp);
    p_eXm_ANOVA_hpc_temp = sum(anova_cnts_3_hpc_temp >= anova_n_hpc_temp(3)) / length(anova_cnts_3_hpc_temp);

    fig = figure('Units','centimeters', 'Position', [0 0 9.0 8]); 
    subplot(331)
    histogram(anova_cnts_1_hpc_temp, binsV,'FaceColor',rgb('MidnightBlue')); hold on;
    xline(anova_n_hpc_temp(1),'LineWidth',1.5,'Color',rgb('DarkOrange')); xlim([0,20]);
    ylabel('Emotion')
    title(sprintf('Hippocampus\n(n=%d)\n p = %0.3f',sum(areaMask==1 | areaMask==2),p_emo_ANOVA_hpc_temp))
    
    subplot(334)
    histogram(anova_cnts_2_hpc_temp, binsV,'FaceColor',rgb('MidnightBlue')); hold on;
    xline(anova_n_hpc_temp(2),'LineWidth',1.5,'Color',rgb('DarkOrange')); xlim([0,20]);
    ylabel(sprintf('Number of runs\n Memory'))
    title(sprintf('p = %0.3f',p_mem_ANOVA_hpc_temp))
    
    subplot(337)
    histogram(anova_cnts_3_hpc_temp, binsV,'FaceColor',rgb('MidnightBlue')); hold on;
    xline(anova_n_hpc_temp(3),'LineWidth',1.5,'Color',rgb('DarkOrange')); xlim([0,20]);
    ylabel('Interaction')
    title(sprintf('p = %0.3f',p_eXm_ANOVA_hpc_temp))
    
    subplot(332)
    histogram(anova_cnts_1_amyg_temp, binsV,'FaceColor',rgb('MidnightBlue')); hold on;
    xline(anova_n_amyg_temp(1),'LineWidth',1.5,'Color',rgb('DarkOrange')); xlim([0,20]);
    title(sprintf('Amygdala\n(n=%d)\n p = %0.3f',sum(areaMask==3),p_emo_ANOVA_amyg_temp))
    
    subplot(335)
    histogram(anova_cnts_2_amyg_temp, binsV,'FaceColor',rgb('MidnightBlue')); hold on;
    xline(anova_n_amyg_temp(2),'LineWidth',1.5,'Color',rgb('DarkOrange')); xlim([0,20]);
    title(sprintf('p = %0.3f',p_mem_ANOVA_amyg_temp))
    
    subplot(338)
    histogram(anova_cnts_3_amyg_temp, binsV,'FaceColor',rgb('MidnightBlue')); hold on;
    xline(anova_n_amyg_temp(3),'LineWidth',1.5,'Color',rgb('DarkOrange')); xlim([0,20]);
    title(sprintf('p = %0.3f',p_eXm_ANOVA_amyg_temp))
    xlabel('Number of significant neurons')
    
    subplot(333)
    histogram(anova_cnts_1_ec_temp, binsV,'FaceColor',rgb('MidnightBlue')); hold on;
    xline(anova_n_ec_temp(1),'LineWidth',1.5,'Color',rgb('DarkOrange')); xlim([0,20]);
    title(sprintf('Entorhinal Cortex\n(n=%d)\n p = %0.3f',sum(areaMask==4),p_emo_ANOVA_ec_temp))
    
    subplot(336)
    histogram(anova_cnts_2_ec_temp, binsV,'FaceColor',rgb('MidnightBlue')); hold on;
    xline(anova_n_ec_temp(2),'LineWidth',1.5,'Color',rgb('DarkOrange')); xlim([0,20]);
    title(sprintf('p = %0.3f',p_mem_ANOVA_ec_temp))
    
    subplot(339)
    histogram(anova_cnts_3_ec_temp, binsV,'FaceColor',rgb('MidnightBlue')); hold on;
    xline(anova_n_ec_temp(3),'LineWidth',1.5,'Color',rgb('DarkOrange')); xlim([0,20]);
    title(sprintf('p = %0.3f',p_eXm_ANOVA_ec_temp))
    fontsize(fig,8,'points');
    fontname(fig,"Arial");
    exportgraphics(fig,strcat(paths.pics,sprintf('_SFig2_%s_%s_ANOVA_byRegion.png',task,anova_names{i0})));
    saveas(fig,strcat(paths.pics,sprintf('_SFig2_%s_%s_ANOVA_byRegion.svg',task,anova_names{i0})));
    close(fig);
end

%% ANOVA - Only memory by region - SFig. 2C (encoding) and 2F (recognition)
if strcmp(task,'enc')
    anova_types = {anovaRes, anovaResRK, anovaResKF};
    anova_names = {'R vs. F','R vs. K','K vs. F'};
elseif strcmp(task,'rec')
    anova_types = {anovaResRK, anovaResRM, anovaResRC, anovaResMK, anovaResCK, anovaResCM};
    anova_names = {'R vs. K', 'R vs. M', 'R vs. CR', 'K vs. M', 'K vs. CR', 'M vs. CR'};
end
binsV = 1:100; 

if strcmp(task,'enc')
    fig = figure('Units','centimeters', 'Position', [0 0 9.0 8]); hold on;
elseif strcmp(task,'rec')
    fig = figure('Units','centimeters', 'Position', [0 0 9.0 15]); hold on; % Bigger during recognition for more plots
end
tiledlayout(length(anova_names),3) % 3 brain regions

for i0 = 1:length(anova_types)
    anova_temp = anova_types{i0};
    anova_p_temp = [];
    anova_p_sh_temp = [];
    for i = 1:length(anova_temp) % Would need to eliminate stubjects here
        anova_p_temp = [anova_p_temp, anova_temp{i}.p]; % Save p-values in order: Emotion, Memory, Interaction. Should be optimized.
        anova_p_sh_temp = [anova_p_sh_temp; anova_temp{i}.p_sh]; % Save p-values in order: Emotion, Memory, Interaction. Should be optimized.    
    end
    
    anova_sig_temp = anova_p_temp < 0.05;
    
    % Are the number of Amygdala neurons significantly above chance? 
    anova_n_amyg_temp = sum(anova_sig_temp(:,areaMask==3),2);
    anova_p_sh_amyg_temp = anova_p_sh_temp(areaMask==3,:,:);
    
    anova_cnts_1_amyg_temp = sum(anova_p_sh_amyg_temp(:,:,1) < 0.05, 1);
    anova_cnts_2_amyg_temp = sum(anova_p_sh_amyg_temp(:,:,2) < 0.05, 1);
    anova_cnts_3_amyg_temp = sum(anova_p_sh_amyg_temp(:,:,3) < 0.05, 1);

    p_emo_ANOVA_amyg_temp = sum(anova_cnts_1_amyg_temp >= anova_n_amyg_temp(1)) / length(anova_cnts_1_amyg_temp);
    p_mem_ANOVA_amyg_temp = sum(anova_cnts_2_amyg_temp >= anova_n_amyg_temp(2)) / length(anova_cnts_2_amyg_temp);
    p_eXm_ANOVA_amyg_temp = sum(anova_cnts_3_amyg_temp >= anova_n_amyg_temp(3)) / length(anova_cnts_3_amyg_temp);
    
    % Are the number of EC neurons significantly above chance? 
    anova_n_ec_temp = sum(anova_sig_temp(:,areaMask==4),2);
    anova_p_sh_ec_temp = anova_p_sh_temp(areaMask==4,:,:);

    anova_cnts_1_ec_temp = sum(anova_p_sh_ec_temp(:,:,1) < 0.05, 1);
    anova_cnts_2_ec_temp = sum(anova_p_sh_ec_temp(:,:,2) < 0.05, 1);
    anova_cnts_3_ec_temp = sum(anova_p_sh_ec_temp(:,:,3) < 0.05, 1);

    p_emo_ANOVA_ec_temp = sum(anova_cnts_1_ec_temp >= anova_n_ec_temp(1)) / length(anova_cnts_1_ec_temp);
    p_mem_ANOVA_ec_temp = sum(anova_cnts_2_ec_temp >= anova_n_ec_temp(2)) / length(anova_cnts_2_ec_temp);
    p_eXm_ANOVA_ec_temp = sum(anova_cnts_3_ec_temp >= anova_n_ec_temp(3)) / length(anova_cnts_3_ec_temp);
    
    % Are the number of HPC neurons significantly above chance? 
    anova_n_hpc_temp = sum(anova_sig_temp(:,areaMask==1 | areaMask==2),2);
    anova_p_sh_hpc_temp = anova_p_sh_temp(areaMask==1 | areaMask==2,:,:);

    anova_cnts_1_hpc_temp = sum(anova_p_sh_hpc_temp(:,:,1) < 0.05, 1);
    anova_cnts_2_hpc_temp = sum(anova_p_sh_hpc_temp(:,:,2) < 0.05, 1);
    anova_cnts_3_hpc_temp = sum(anova_p_sh_hpc_temp(:,:,3) < 0.05, 1);

    p_emo_ANOVA_hpc_temp = sum(anova_cnts_1_hpc_temp >= anova_n_hpc_temp(1)) / length(anova_cnts_1_hpc_temp);
    p_mem_ANOVA_hpc_temp = sum(anova_cnts_2_hpc_temp >= anova_n_hpc_temp(2)) / length(anova_cnts_2_hpc_temp);
    p_eXm_ANOVA_hpc_temp = sum(anova_cnts_3_hpc_temp >= anova_n_hpc_temp(3)) / length(anova_cnts_3_hpc_temp);

    nexttile()
    histogram(anova_cnts_2_hpc_temp, binsV,'FaceColor',rgb('MidnightBlue')); hold on;
    xline(anova_n_hpc_temp(2),'LineWidth',1.5,'Color',rgb('DarkOrange')); xlim([0,20]);
    if i0 == 3
        ylabel(sprintf('Number of runs\n %s', anova_names{i0} ))
    else
        ylabel(sprintf('%s',anova_names{i0} ))
    end

    if i0 == 1
        title(sprintf('Hippocampus\np = %0.3f', p_mem_ANOVA_hpc_temp))
    else
        title(sprintf('p = %0.3f', p_mem_ANOVA_hpc_temp))
    end

    nexttile()
    histogram(anova_cnts_2_amyg_temp, binsV,'FaceColor',rgb('MidnightBlue')); hold on;
    xline(anova_n_amyg_temp(2),'LineWidth',1.5,'Color',rgb('DarkOrange')); xlim([0,20]);
    if i0 == 1
        title(sprintf('Amygdala\np = %0.3f', p_mem_ANOVA_amyg_temp))
    else
        title(sprintf('p = %0.3f', p_mem_ANOVA_amyg_temp))
    end

    nexttile()
    histogram(anova_cnts_2_ec_temp, binsV, 'FaceColor',rgb('MidnightBlue')); hold on;
    xline(anova_n_ec_temp(2), 'LineWidth',1.5,'Color',rgb('DarkOrange')); xlim([0,20]);
    if i0==1
        title(sprintf('Entorhinal Cortex\np = %0.3f',p_mem_ANOVA_ec_temp))
    else
        title(sprintf('p = %0.3f',p_mem_ANOVA_ec_temp))
    end

end

fontsize(fig,8,'points');
fontname(fig,"Arial");
if strcmp(task,'enc')
    exportgraphics(fig,[paths.pics,'_SFig2C_enc_ANOVA_OnlyMemory.png']);
    saveas(fig,[paths.pics,'_SFig2C_enc_ANOVA_OnlyMemory.svg']);
elseif strcmp(task,'rec')
    exportgraphics(fig,[paths.pics,'_SFig2F_rec_ANOVA_OnlyMemory.png']);
    saveas(fig,[paths.pics,'_SFig2F_rec_ANOVA_OnlyMemory.svg']);
end
close(fig);

%% SFig 3: ANOVA@ By subject
fig = figure('Units','centimeters', 'Position', [0 0 9.0 13]); hold on;
tiledlayout(length(sub), size(anovaRes{1}.p,1));
ttypes = {'Emotion', 'Memory', 'Interaction'};
for iS=1:length(subjects)
    for ia = 1:length(ttypes) % 3 types: emotion, memory and interaction
        nexttile()
        anova_sub_sh = squeeze(sum(anovaRes{iS}.p_sh < 0.05,1));
        anova_sub_n = sum(anovaRes{iS}.p < 0.05,2);
        anova_sub_p = sum(anova_sub_sh(:,ia) >= anova_sub_n(ia)) / length(anova_sub_sh);
        if anova_sub_p == 0
            anova_sub_p = 1/size(anova_sub_sh,1);
        end
        histogram(anova_sub_sh(:,ia), binsV,'FaceColor',rgb('MidnightBlue')); hold on; % shuffled distributions
        xline(anova_sub_n(ia,:),'LineWidth',1.5,'Color',rgb('DarkOrange')); xlim([0,25]); % Real number of neurons
        if iS ==1
            title(sprintf('%s\n p=%0.4f',ttypes{ia},anova_sub_p),'FontWeight','Normal')
        else
            title(sprintf('p = %0.4f',anova_sub_p),'FontWeight','Normal')
        end

        % Axes labels
        if ia == 1 && iS ==3
            ylabel(sprintf('Number of runs\n %s\n n = %d',subjects{iS,1},size(anovaRes{iS}.p,2)));
        elseif iS == 5 && ia ==2
            xlabel('Number of significant neurons');
        elseif ia == 1
            ylabel(sprintf('%s\n n = %d',subjects{iS,1},size(anovaRes{iS}.p,2)));
        end
            
    end
end
fontsize(fig,8,'points');
fontname(fig,"Arial");
exportgraphics(fig,strcat(paths.pics,'_SFig3_',task,'_ANOVA_bySubject.png'),'Resolution',300) 
saveas(fig, [paths.pics,'_SFig3_',task,'_ANOVA_bySubject.svg']) 
close(fig);

if strcmp(task,'rec')

    %% Fig.2: Recogntion - Analyze percentages of individual neuronal responses

    eR_absZ = [];
    nR_absZ = [];
    eM_absZ = [];
    nM_absZ = [];
    eCR_absZ = [];
    nCR_absZ = [];

    for iS=1:length(sub)
    
        eR_absZ   = [eR_absZ,   mean(sub(iS).absallZ(sub(iS).onsets.eRHit,:),1)];
        nR_absZ   = [nR_absZ,   mean(sub(iS).absallZ(sub(iS).onsets.nRHit,:),1)];
        eCR_absZ  = [eCR_absZ,  mean(sub(iS).absallZ(sub(iS).onsets.eCR,:),1)];
        nCR_absZ  = [nCR_absZ,  mean(sub(iS).absallZ(sub(iS).onsets.nCR,:),1)];
        eM_absZ   = [eM_absZ,   mean(sub(iS).absallZ(sub(iS).onsets.eMiss,:),1)];
        nM_absZ   = [nM_absZ,   mean(sub(iS).absallZ(sub(iS).onsets.nMiss,:),1)];
        
        zeR{iS} = sub(iS).ZeR([sub(iS).peRCnt] < 0.05);
        zeCR{iS} = sub(iS).ZeCR([sub(iS).peCRCnt] < 0.05);
        zeM{iS} = sub(iS).ZeM([sub(iS).peMCnt] < 0.05);
        znR{iS} = sub(iS).ZnR([sub(iS).pnRCnt] < 0.05);
        znCR{iS} = sub(iS).ZnCR([sub(iS).pnCRCnt] < 0.05);
        znM{iS} = sub(iS).ZnM([sub(iS).pnMCnt] < 0.05);
    
        n_neur(iS) = length(sub(iS).pEmoCnt);
    end

    % stacked is the counts of neurons with significantly increasing or decreasing FR per condition - baseline vs stimulus period
    stacked = [sum([zeR{:}]>0), sum([znR{:}]>0), sum([zeM{:}]>0), sum([znM{:}]>0), sum([zeCR{:}]>0), sum([znCR{:}]>0); ...
               sum([zeR{:}]<0), sum([znR{:}]<0), sum([zeM{:}]<0), sum([znM{:}]<0), sum([zeCR{:}]<0), sum([znCR{:}]<0)];

    %% Fig 2 Recognition

    % For plotting pairwise memory effects
    anova_types = {anovaResRM, anovaResRC, anovaResCM};
    anova_names = {'RHit vs. Miss','RHit vs. CR','Miss vs. CR'};

    fig = figure('Units','centimeters', 'Position', [0 0 9.5 12.5]); hold on;
    subplot(321)
    bh = bar(stacked'./sum(n_neur) * 100,'stacked', 'FaceColor','flat');
    bh(1).CData = [0.25 0.25 0.25];  % Change color to first level
    bh(2).CData = [0.75 0.75 0.75];  % Change color to second level
    legend('Increased','Decreased','Location','northeast');
    legend boxoff
    x = {'eRHit', 'nRHit','eMiss','nMiss', 'eCR', 'nCR'};
    xticks(1:6);  xticklabels(x); xtickangle(90)
    xlim([0.5 6.5]); ylim([0 100]); yticks(0:20:100); 
    ylabel('% of neurons');  
    title('Event responsiveness')

    subplot(3,2,2);
    Y = [data_all_neurons_small;
        0 data_anysig_small];
    subplot(322);
    bar(Y,'stacked');
    xticklabels({'of total','of selective'})
    ylabel('% of neurons')
    title('Differential selectivity')
    lgd = legend(barNamesAll_small,'location','east');
    legend box off;
    ylim([0, 100])
    xlim([0.5 2.5])

    % Do we have more neurons with significant ANOVA values than chance levels? Comparing using shuffled trial labels.
    binsV = 1:max(anova_n)+5;
    subplot(334)
    histogram(anova_cnts_1, binsV,'FaceColor',rgb('MidnightBlue')); hold on;
    xline(anova_n(1),'LineWidth',1.5,'Color',rgb('DarkOrange'));
    text(0.1, 0.95, sprintf('p = %0.4g', p_emo_ANOVA), 'Units', 'normalized', 'FontSize', 10);
    ylabel('Number of runs');  title('Emotion'); ylim([0,1300]); xlim([0,45]); %xticks([0:10:55])
    
    subplot(335)
    histogram(anova_cnts_2, binsV,'FaceColor',rgb('MidnightBlue')); hold on;
    xline(anova_n(2),'LineWidth',1.5,'Color',rgb('DarkOrange'));
    text(0.1, 0.95, sprintf('p = %0.4f', p_mem_ANOVA), 'Units', 'normalized', 'FontSize', 10);
    xlabel('Number of significant neurons');  title('Memory'); ylim([0,1300]); xlim([0,45]); %xticks([0:10:55])
    
    subplot(336)
    histogram(anova_cnts_3, binsV,'FaceColor',rgb('MidnightBlue')); hold on;
    xline(anova_n(3),'LineWidth',1.5,'Color',rgb('DarkOrange'));
    text(0.1, 0.95, sprintf('p = %0.4f', p_eXm_ANOVA), 'Units', 'normalized', 'FontSize', 10);
    title('Interaction'); ylim([0,1300]); xlim([0,45]); %xticks([0:10:55]) % ylim([0,150])

    % Plot the pairwise memory effects
    for i0 = 1:length(anova_types)
        anova_temp = anova_types{i0};
        anova_p_temp = [];
        anova_p_sh_temp = [];
        for i = 1:length(anova_temp) % Would need to eliminate stubjects here
            anova_p_temp = [anova_p_temp, anova_temp{i}.p]; % Save p-values in order: Emotion, Memory, Interaction. Should be optimized.
            anova_p_sh_temp = [anova_p_sh_temp; anova_temp{i}.p_sh]; % Save p-values in order: Emotion, Memory, Interaction. Should be optimized.    
        end

        anova_cnts_2_temp = sum(anova_p_sh_temp(:,:,2) < 0.05, 1); % memory effect only
        anova_n_temp = sum(anova_p_temp < 0.05,2); % Number of neurons with significant memory effects. In order: Emotion, memory, interaction
        
        % Compare shuffled distribution to the real ones ANOVA@
        p_mem_ANOVA_temp = sum(anova_cnts_2_temp >= anova_n_temp(2)) / length(anova_cnts_2_temp);

        subplot(3, 3, 6+i0)
        histogram(anova_cnts_2_temp, binsV,'FaceColor',rgb('MidnightBlue')); hold on;
        xline(anova_n_temp(2),'LineWidth',1.5,'Color',rgb('DarkOrange'));
        text(0.1, 0.95, sprintf('p = %0.4f', p_mem_ANOVA_temp), 'Units', 'normalized', 'FontSize', 10);
        title(anova_names{i0}); ylim([0,1300])
        if i0 ==1 || i0 == 4
            ylabel('Number of runs');
        elseif i0 == 5
            xlabel('Number of significant neurons');
        end
    end
    
    fontsize(fig,8,'points');
    fontname(fig,"Arial");
    saveas(fig, [paths.pics,'_Fig2_recognition_final.svg']);
    exportgraphics(fig,[paths.pics,'_Fig2_recognition_final.png'],'Resolution',300);
    close(fig)

    %% SFig 2D: Recognition ANOVA for all brain regions using KHits
    anova_types = {anovaResRK, anovaResCK, anovaResMK}; 
    anova_names = {'RHit vs. KHit', 'CR vs. KHit', 'Miss vs. KHit'}; 

    binsV2 = 1:25; 
    fig = figure('Units','centimeters', 'Position', [0 0 9.5 3]); hold on;

    for i0 = 1:length(anova_types)
        anova_temp = anova_types{i0};
        anova_p_temp = [];
        anova_p_sh_temp = [];
        for i = 1:length(anova_temp) % Would need to eliminate stubjects here
            anova_p_temp = [anova_p_temp, anova_temp{i}.p]; % Save p-values in order: Emotion, Memory, Interaction. Should be optimized.
            anova_p_sh_temp = [anova_p_sh_temp; anova_temp{i}.p_sh]; % Save p-values in order: Emotion, Memory, Interaction. Should be optimized.
        end
        
        anova_cnts_2_temp = sum(anova_p_sh_temp(:,:,2) < 0.05, 1); % Memory effect
        anova_n_temp = sum(anova_p_temp < 0.05,2); % NOTE 2x2 for encoding and 3x2 for recognition
    
        % Compare shuffled distribution to the real ones ANOVA@
        p_mem_ANOVA_temp = sum(anova_cnts_2_temp >= anova_n_temp(2)) / length(anova_cnts_2_temp);
    
        subplot(1,3,i0)
        histogram(anova_cnts_2_temp, binsV2,'FaceColor',rgb('MidnightBlue')); hold on;
        xline(anova_n_temp(2),'LineWidth',1.5,'Color',rgb('DarkOrange'));
        text(0.1, 0.95, sprintf('p = %0.4f', p_mem_ANOVA_temp), 'Units', 'normalized', 'FontSize', 10);
        title(anova_names{i0}); ylim([0,1300])
        if i0 ==1 
            ylabel('Number of runs');
        elseif i0 == 2 
            xlabel('Number of significant neurons');
        end
        
    end
    fontsize(fig,8,'points');
    fontname(fig,"Arial");
    saveas(fig, [paths.pics,'_SFig2D_recognition_KHit_comp.svg']);
    exportgraphics(fig,[paths.pics,'_SFig2D_recognition_KHit_comp.png'],'Resolution',300);
    close(fig)

end % rec

%% Write to csv files to compute linear mixed effects models in R
neur_name = [sub(:).name];

if strcmp(task,'enc')
    csvfile = fullfile(paths.pics, 'encoding_trial_results_LMM.csv');
    fileID = fopen(csvfile,'a+');
    fprintf(fileID, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,\n', 'subjNeur', 'subject', 'neur_name', 'region', 'trial', 'emoneu', 'remfor', 'baseCnts', ...
                    'stimCnts', 'baseRate', 'stimRate','allZ','absAllZ','rawAbsAllZ'); 

    for si2 = 1:length(subjects)
        emoneu = cell(1,120);
        remfor = cell(1,120);
        emoneu(sub(si2).onsets.emo_all) = {'emo'};
        emoneu(sub(si2).onsets.neu_all) = {'neu'};
        remfor(sub(si2).onsets.eR) = {'eR'};
        remfor(sub(si2).onsets.eF) = {'eF'};
        remfor(sub(si2).onsets.eK) = {'eK'};
        remfor(sub(si2).onsets.nR) = {'nR'};
        remfor(sub(si2).onsets.nF) = {'nF'};
        remfor(sub(si2).onsets.nK) = {'nK'};
        for ni2 = 1:length(sub(si2).name) % for each neuron
            
            if contains(sub(si2).name{ni2},'AH') || contains(sub(si2).name{ni2},'uH')
                region = 'aHP';
            elseif contains(sub(si2).name{ni2},'PH')
                region = 'pHP';
            elseif contains(sub(si2).name{ni2},'AR') || contains(sub(si2).name{ni2},'AL')
                region = 'Amyg';
            elseif contains(sub(si2).name{ni2},'EC')
                region = 'EC';
            end

            for ti2=1:size(sub(1).baseCnts,1) % Now through all trials
                fprintf(fileID, '%s,%s,%s,%s,%d,%s,%s,%d,%d,%f,%f,%f,%f,%f,\n',[subjects{si2}, '_', sub(si2).name{ni2}], subjects{si2}, sub(si2).name{ni2}, region, ti2, emoneu{ti2}, remfor{ti2}, sub(si2).baseCnts(ti2,ni2), ...
                    sub(si2).stimCnts(ti2,ni2), sub(si2).baseRate(ti2,ni2), sub(si2).stimRate(ti2,ni2), sub(si2).allZ(ti2,ni2), sub(si2).absallZ(ti2,ni2), sub(si2).rawabsallZ(ti2,ni2)); 
            end
        end
    end
    fclose(fileID);    

elseif strcmp(task,'rec')
    
    csvfile = fullfile(paths.pics, 'recognition_trial_results_LMM.csv');
    fileID = fopen(csvfile,'a+');
    fprintf(fileID, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,\n', 'subjNeur','subject', 'neur_name', 'region', 'trial', 'emoneu', 'remfor', 'baseCnts', ...
                    'stimCnts', 'baseRate', 'stimRate', 'allZ', 'absAllZ', 'rawAbsAllZ');

    for si2 = 1:length(subjects)
        emoneu = cell(1,240);
        remfor = cell(1,240);
        emoneu(sub(si2).onsets.emo_all) = {'emo'};
        emoneu(sub(si2).onsets.neu_all) = {'neu'};
        remfor(sub(si2).onsets.eRHit) = {'eRHit'};
        remfor(sub(si2).onsets.eKHit) = {'eKHit'};
        remfor(sub(si2).onsets.eCR) = {'eCR'};
        remfor(sub(si2).onsets.eMiss) = {'eM'};
        remfor(sub(si2).onsets.nRHit) = {'nRHit'};
        remfor(sub(si2).onsets.nCR) = {'nCR'};
        remfor(sub(si2).onsets.nMiss) = {'nM'};
        remfor(sub(si2).onsets.nKHit) = {'nKHit'};
        for ni2 = 1:length(sub(si2).name) % for each neuron
            
            if contains(sub(si2).name{ni2},'AH') || contains(sub(si2).name{ni2},'uH')
                region = 'aHP';
            elseif contains(sub(si2).name{ni2},'PH')
                region = 'pHP';
            elseif contains(sub(si2).name{ni2},'AR') || contains(sub(si2).name{ni2},'AL')
                region = 'Amyg';
            elseif contains(sub(si2).name{ni2},'EC')
                region = 'EC';
            end

            for ti2=1:size(sub(si2).baseCnts,1) % Now through all trials
                fprintf(fileID, '%s,%s,%s,%s,%d,%s,%s,%d,%d,%f,%f,%f,%f,%f,\n',[subjects{si2}, '_', sub(si2).name{ni2}],subjects{si2}, sub(si2).name{ni2}, region, ti2, emoneu{ti2}, remfor{ti2}, sub(si2).baseCnts(ti2,ni2), ...
                    sub(si2).stimCnts(ti2,ni2), sub(si2).baseRate(ti2,ni2), sub(si2).stimRate(ti2,ni2), sub(si2).allZ(ti2,ni2), sub(si2).absallZ(ti2,ni2), sub(si2).rawabsallZ(ti2,ni2)); 
            end
        end
    end
    fclose(fileID);
end

%% Save dPCA input structure
time = psth_use_time;
timeEvents = [0 0.5]; % Image presentation at 0 sec and disappears at 0.5 sec.
timeEventsNames = {'ON', 'OFF'};

if strcmp(task,'enc')
    firingRatesRKF = firingRatesRKF(1:sp, :, :, :, 1:max(trialNumRKF(:))); % drop nonused zeros - arrays are initialized larger than necessary to speed up script
    firingRatesAverageRKF = firingRatesAverageRKF(1:sp, :, :, :);
    % sparsifying 
    firingRates_sizeRKF = size(firingRatesRKF);
    firingRates_sparseRKF = sparse(firingRatesRKF(:));
    
    save(strcat(paths.pics,'_dPCA_input_RKF.mat'), 'firingRates_sparseRKF', 'firingRates_sizeRKF','firingRatesAverageRKF', 'trialNumRKF', ...
                         'time', 'timeEvents', 'timeEventsNames', 'subjMask', 'areaMask', 'dPCA_neur_name')
elseif strcmp(task,'rec')
    % drop nonused zeros - arrays are initialized larger than necessary to speed up script
    firingRatesRec = firingRatesRec(1:sp, :, :, :, 1:max(trialNumRec(:)));
    firingRatesAverageRec = firingRatesAverageRec(1:sp, :, :, :);
    % sparsifying 
    firingRates_sizeRec = size(firingRatesRec);
    firingRates_sparseRec = sparse(firingRatesRec(:));
    
    save(strcat(paths.pics,'_dPCA_input_4stim_2dec.mat'), 'firingRates_sparseRec', 'firingRates_sizeRec', 'firingRatesAverageRec', 'trialNumRec', ...
                     'time', 'timeEvents', 'timeEventsNames', 'subjMask', 'areaMask', 'dPCA_neur_name')
end