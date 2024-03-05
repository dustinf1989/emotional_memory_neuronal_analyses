%%% Produce Tables S1 and S2 and accompanying statistics
%%% 
%%% IAPS Memory Test
%%% Day 1 Encode, Day 2 Recognition
%%% In recog test, R 97 - Left ( The first column of log.res
%%%                K 100 - Down
%%%                N 98 - Right
%%% For encoding, which_e/n determines which of the pool of stimuli will be
%%% used for this part. which_e/n is used to write "encode_list"
%%% IAPS_enc is a randomisation of these stimuli (cogent uses IAPS_enc to
%%% index encode_list during encoding)
%%%
%%% Originally written by Manuela Costa, adapted by Dustin Fetterhoff

clearvars; close all; clc;

% Change this path to the data directory - Nothing else should need to be changed
data_folder = 'C:/Users/drdus/emotional_memory_neuronal_data/';

all = load(fullfile(data_folder,'images/all_images.txt')); %ALL STIMULI 1st 80 emotional; IAPS codes
 
% Uncomment one of the below to select the desired subject numbers
% subjects = [1 5 6 8 13]; % To get stats for only subjects with neuronal data (Not reported)
subjects = [1 2 5 6 8 10 11 12 13]; % Get stats for all subjects

numberofsubjects = size(subjects,2);

for isub=1:length(subjects)
    sub = subjects(isub);
    enc_dir=[data_folder,'/enc/Patient',num2str(sub),'/Enc/'];
    rec_dir=[data_folder,'/rec/Patient',num2str(sub),'/Rec/'];

    encode_list = load([enc_dir,'encode_list.txt']); %list of 120 IAPS image numbers
    load([enc_dir,'IAPS_enc.mat']); %% IAPS_enc is length 120; E items are 1:40, N from 41:120... which pertains to their position in Encode_List.txt
    
    what_encoded = encode_list(IAPS_enc); %GIVES THE EXACT ORDER IMAGES DISPLAYED
    
    what_encoded_e = intersect(what_encoded,all(1:80));
    what_encoded_n = intersect(what_encoded,all(81:end));
    
    e_pics_enc = find(ismember(what_encoded,what_encoded_e));%gives position of emotional pics in encoding presentation (40x1)
    n_pics_enc = find(ismember(what_encoded,what_encoded_n));%gives position of neutral pics in encoding presentation (80x1)
    
    %%%%RECOGNITION%%%%
        
    log = load([rec_dir,'log.res']); 
    response = log(:,1); % R 97 - Left, K 100 - Down, N 98 - Right
    
    load([rec_dir,'IAPS_rec.mat']);
    what_recognised = all(IAPS_rec);
    
    old_E = find(ismember(what_recognised,what_encoded_e)); % all the old emotional pictures, listed by number of position they were shown
    what_recognised_e = find(ismember(what_recognised,all(1:80))); % all emotional items shown (80), listed by number they were shown at
    new_E = setdiff(what_recognised_e,old_E); % all new emotional items (40) not shown the previous day, listed by number they were shown at
    
    old_N = find(ismember(what_recognised,what_encoded_n));
    what_recognised_n = find(ismember(what_recognised,all(81:end)));
    new_N = setdiff(what_recognised_n,old_N);
    
    rec_onsets = struct('Old_e_items',old_E,'Old_n_items',old_N,'New_e_items',new_E,'New_n_items',new_N,'emo_all',sort([old_E;new_E]),'neu_all',sort([old_N;new_N]),'old_all',sort([old_E;old_N]),'new_all',sort([new_E;new_N]));
    rec_onsets.eRHit    =    old_E(find(response(old_E)==97));% correctly remembered
    rec_onsets.eKHit    =    old_E(find(response(old_E)==100));% correctly familiar with (known)
    rec_onsets.eMiss =    old_E(find(response(old_E)==98));% already seen but not recognised
    rec_onsets.eCR   =    new_E(find(response(new_E)==98)); % never seen and not recognised
    rec_onsets.eRFA =    new_E(find(response(new_E)==97)); % wrongly remembered
    rec_onsets.eKFA =    new_E(find(response(new_E)==100));% wrongly familiar with
    rec_onsets.eNotResp = vertcat(old_E(find(response(old_E)==0)), new_E(find(response(new_E)==0)));  % No response in Recognition Task
    
    rec_onsets.nRHit    =    old_N(find(response(old_N)==97));
    rec_onsets.nKHit    =    old_N(find(response(old_N)==100));
    rec_onsets.nMiss =    old_N(find(response(old_N)==98));
    rec_onsets.nCR   =    new_N(find(response(new_N)==98));
    rec_onsets.nRFA =    new_N(find(response(new_N)==97));
    rec_onsets.nKFA =    new_N(find(response(new_N)==100));
    rec_onsets.nNotResp = vertcat(old_N(find(response(old_N)==0)), new_N(find(response(new_N)==0)));
%     save('rec_onsets.mat','rec_onsets') % Onsets files are already save in the given datasets
    
    old_eCorrRemID     = what_recognised(rec_onsets.eRHit); %gets the stimulus number, IAPS reference number of the picture
    old_eCorrFamID     = what_recognised(rec_onsets.eKHit);
    old_eMissedID  = what_recognised(rec_onsets.eMiss);
    old_eNotRespID = what_recognised(rec_onsets.eNotResp);
    
    old_nCorrRemID     = what_recognised(rec_onsets.nRHit);
    old_nCorrFamID     = what_recognised(rec_onsets.nKHit);
    old_nMissedID  = what_recognised(rec_onsets.nMiss);
    old_nNotRespID = what_recognised(rec_onsets.nNotResp);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%Now go back to Encoding to do subsequent memory analysis%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    enc_onsets = struct('emo_all',e_pics_enc,'neu_all',n_pics_enc); 
    
    enc_onsets.eR = find(ismember(what_encoded,old_eCorrRemID));
    enc_onsets.eK = find(ismember(what_encoded,old_eCorrFamID));
    enc_onsets.eF = find(ismember(what_encoded,old_eMissedID));
    enc_onsets.eNotResp = find(ismember(what_encoded,old_eNotRespID));
    
    enc_onsets.nR = find(ismember(what_encoded,old_nCorrRemID));
    enc_onsets.nK = find(ismember(what_encoded,old_nCorrFamID));
    enc_onsets.nF = find(ismember(what_encoded,old_nMissedID));
    enc_onsets.nNotResp = find(ismember(what_encoded,old_nNotRespID));
%     save('enc_onsets.mat','enc_onsets') % Onsets files are already save in the given datasets

   %% Now do memory performance analysis

   names = {'Hits', 'FA','eRHit', 'eRFA', 'eKHit', 'eKFA', 'eMiss', 'eCR', 'nRHit', 'nRFA', 'nKHit', 'nKFA', 'nMiss', 'nCR', 'RKHits','RKFA','eRKHits','eRKFA','nRKHits','nRKFA'};
    
   recog(isub,:) = [
        (length(find(response(old_E)==97))+length(find(response(old_N)==97))) / (length(old_E)+length(old_N)); % Hits all 
        (length(find(response(new_E)==97))+length(find(response(new_N)==97))) / (length(new_E)+length(new_N)); % False Alarams all
        length(find(response(old_E)==97))/length(old_E); %E R Hits
        length(find(response(new_E)==97))/length(new_E); %E R FAs
        length(find(response(old_E)==100))/length(old_E); %E K Hits
        length(find(response(new_E)==100))/length(new_E); %E K FAs
        length(find(response(old_E)==98))/length(old_E); %E Miss
        length(find(response(new_E)==98))/length(new_E); %E Crj
        length(find(response(old_N)==97))/length(old_N); %N R Hits 
        length(find(response(new_N)==97))/length(new_N); %N R FAs
        length(find(response(old_N)==100))/length(old_N); %N K Hits
        length(find(response(new_N)==100))/length(new_N); %N K FAs
        length(find(response(old_N)==98))/length(old_N); %N Miss
        length(find(response(new_N)==98))/length(new_N); %N Crj
        (length(find(response(old_E)==97))+length(find(response(old_N)==97))+length(find(response(old_E)==100))+length(find(response(old_N)==100))) / (length(old_E)+length(old_N)); % RKHits all 
        (length(find(response(new_E)==97))+length(find(response(new_N)==97))+length(find(response(new_E)==100))+length(find(response(new_N)==100))) / (length(new_E)+length(new_N)); % RK FAs all
        (length(find(response(old_E)==97))+length(find(response(old_E)==100))) / (length(old_E)); % eRK Hits 
        (length(find(response(new_E)==97))+length(find(response(new_E)==100))) / (length(new_E));  % eRK FA 
        (length(find(response(old_N)==97))+length(find(response(old_N)==100))) / (length(old_N)); % nRK Hits 
        (length(find(response(new_N)==97))+length(find(response(new_N)==100))) / (length(new_N))]; % nRK FA

    namesTr = {'eRHit', 'eRFA', 'eKHit', 'eKFA', 'eMiss', 'eCR', 'nRHit', 'nRFA', 'nKHit', 'nKFA', 'nMiss', 'nCR'};

    nTrials(isub,:) = [length(find(response(old_E)==97)); %E R Hits
        length(find(response(new_E)==97)); %E R FAs
        length(find(response(old_E)==100)); %E K Hits
        length(find(response(new_E)==100)); %E K FAs
        length(find(response(old_E)==98)); %E Miss
        length(find(response(new_E)==98)); %E Crj
        length(find(response(old_N)==97)); %N R Hits
        length(find(response(new_N)==97)); %N R FAs
        length(find(response(old_N)==100)); %N K Hits
        length(find(response(new_N)==100)); %N K FAs
        length(find(response(old_N)==98)); %N Miss
        length(find(response(new_N)==98));]; %N Crj
    
end  

%% Get PR and d'prime values

PR_all = recog(:,1)-recog(:,2); %PR: R - FA
R_dprime = dprime(recog(:,1),recog(:,2));

PR_eR = recog(:,3)-recog(:,4); %PR: eR - eRFA
eR_dprime = dprime(recog(:,3), recog(:,4)); % d'prime: eR - eRFA

PR_nR = recog(:,9)-recog(:,10); %PR: nR - nRFA
nR_dprime = dprime(recog(:,9), recog(:,10)); % d'prime: nR - nRFA

[~,p_PR_nR,~,~] = ttest(PR_nR); % Is the percent correct above chance level (Above zero) for emotional remember vs false alarms?

[~,p_PR_eR,~,~] = ttest(PR_eR); % Is the percent correct above chance level (Above zero) for neutral remember vs false alarms?

[~,p_eRFA_nRFA,~,~] = ttest(recog(:,4), recog(:,10)); % eRFA vs nRFA

%% Write to csv file
csvfile = [data_folder,'/Table_S1_number_of_response_types.csv'];
fileID = fopen(csvfile,'a+');
fprintf(fileID, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,\n', 'subject', string(namesTr));
for i=1:length(subjects)
    fprintf(fileID, '%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,\n', sprintf('Patient%d',subjects(i)), nTrials(i,:));
end
fclose(fileID);

%% Write to csv file
csvfile = [data_folder,'/Table_S2_behavioral_performance.csv'];
fileID = fopen(csvfile,'a+');
fprintf(fileID, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,\n', 'subject', names{1}, names{2},'PR_all','dR-all',names{3}, names{4},'ePR','deR',names{9}, names{10},'nPR','dnR');
for i=1:length(subjects)
    fprintf(fileID, '%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,\n', sprintf('Patient%d',subjects(i)), recog(i,1),recog(i,2),PR_all(i),R_dprime(i),recog(i,3),recog(i,4),PR_eR(i),eR_dprime(i), recog(i,9),recog(i,10),PR_nR(i),nR_dprime(i));
end
fclose(fileID);

%% 2x2 ANOVA - Requires Statistics and Machine Learning Toolbox
% Used for stats in supplementary material
data = table([1:size(recog,1)].',recog(:,3), recog(:,4), recog(:,9), recog(:,10), 'VariableNames', {'id', 'eRHit', 'eRFA', 'nRHit', 'nRFA'});
w = table(categorical([1 1 2 2].'), categorical([1 2 1 2].'), 'VariableNames', {'Stimulus', 'Memory'}); % within-desing
disp(data)
rm = fitrm(data, 'nRFA-eRHit ~ 1', 'WithinDesign', w);
ranova(rm, 'withinmodel', 'Stimulus*Memory')
