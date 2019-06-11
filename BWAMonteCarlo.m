%% Specify number of permutations
clear all;close all;%clc;

simTrial = 100; %Iteration number


%% Load Stuff
folder='C:\Users\Ethan\Documents\MATLAB\CPanalysis\';
subjInds=(3);      % vector
pilotType=3;       % 1:singleSpeaker 2:alternateSpeaker 3:cocktailParty

stimulusType=4;    % select predictor (only pilotType>1) 1:female, 2:male, 3:female-male, 4:[female male] 
removeSpeechENV=0; % use 1:boxcar 0:boxcar+ENV
quickMode=1;       % 1:lambdaOptimization on all trials (fast), 0:lambdaOptimization only on training trials (slow)

%%% select data
nirsTypes=[1:4]; % 1:oxyRaw 2:deoxyRaw 3:oxyGlm 4:deoxyGlm {def:[3 4]}, vector
COIs=[666];      % channels of interest, [] uses all direct-neighbor channels, 0 uses all channels, 666=left&right AC
lpCutoff=[0.0099]; % cutoff frequency of lowpass filter {def:<0.5Hz}, for 1-min blocks use 0.05Hz, 0 applies no filtering
inputType=2;     % input data preprocessed in 1:TSI 2:Homer2 {def:2}; NOTE:TSI data may crash!

%%% set timing parameters 
trialDur=1;              % trial duration (in integer minutes), if >1: append consecutive trials {def:2}
lagRange=0:16;           % lags for TRF analysis (in sec), should span the HRF {def:0:20}
lambdaRange=[0 2.^(linspace(0,50,10))]; % range of lambdas to be sampled 

%%% set labels
pilotTypeLabels={'singleSpeaker','alternateSpeaker','cocktailParty'};
inputTypeLabels={'_TSI',''};
stimulusTypeLabels={'female','male','female-male','[female male]'};
nirsTypeLabels={'oxyRaw','deoxyRaw','oxyGlm','deoxyGlm'};

%%% import all or direct-neighbor channels (if requested)
if isempty(COIs) % import all direct-neighbor channels
    load([folder 'data_after_cleaning\channelMask_directNeighbors_pilot.mat']); % imports 'channelMask'
    channelMask=channelMask';channelMask=channelMask(:)';
    if inputType==1 % TSI (contains all possible channels)
        COIs=find(channelMask);
    else % inputType==2 % Homer2 (contains only direct-neighbor channels)
        COIs=1:sum(channelMask);
    end
    clear channelMask
elseif COIs==0 % import all channels
    load([folder 'data_after_cleaning\channelMask_directNeighbors_pilot.mat']); % imports 'channelMask'
    if inputType==1 % TSI (contains all possible channels)
        COIs=1:numel(channelMask);
    else % inputType==2 % Homer2 (contains only direct-neighbor channels)
        COIs=1:sum(channelMask);
    end
    clear channelMask
elseif COIs==666 % both auditory cortices
    load([folder 'data_after_cleaning\channelMask_directNeighbors_pilot.mat']); % imports 'channelMask'
    if inputType==1 % TSI (contains all possible channels)
        COIs=[78 282]; % 78:right AC, 282:left AC
    else % inputType==2 % Homer2 (contains only direct-neighbor channels)
        channelMask=channelMask';channelMask=channelMask(:)';
        channelMask([78 282])=666; % 78:right AC, 282:left AC
        channelMask(channelMask==0)=[];
        COIs=find(channelMask==666);  
        clear channelMask
    end
end


%% Sim Reality
         
for indSubj=subjInds % loop for subjects
    
    disp([' ']);
    disp(['Subject ' num2str(indSubj)]);
    
    %%% import preprocessed FNIRS data and stimuli
    folderSubj=[folder 'data_after_cleaning\S0' num2str(indSubj,'%3.f') '_pilot_' pilotTypeLabels{pilotType} '\'];
    load([folderSubj 'S0' num2str(indSubj,'%3.f') '_pilot_stimuliNIRSproc' inputTypeLabels{inputType} '.mat']); % imports 'stim{ch}(tr,smp)','nirs{nirsType}(tr,smp,ch)','fsNIRS'
    lags=lagRange(1):1/fsNIRS:lagRange(end); % timing parameters for TRF analysis

    %%% design lp filter (if requested)
    if lpCutoff>0
        filterOrder=4; % {def:4}
        [b,a]=butter(filterOrder,lpCutoff*2/fsNIRS,'low');
    end
    
    %%% select desired nirsTypes and channels
    for nirsType=nirsTypes % loop for nirsTypes
        tic
        
        disp(['    ' nirsTypeLabels{nirsType} ':']);
        if pilotType==1 % singleSpeaker = single audio channel --> duplicate
            if inputType==1 % TSI
                stimN{1}=stim;
                stimN{2}=stim;
            else % inputType==2 % Homer2
                stimN{1}=stim{1};
                stimN{2}=stim{1};
            end
        else % pilotType==2 % alternateSpeaker = two audio channels
            if removeSpeechENV==0 % include speech envelope 
                stimN{1}=stim{1}; % left/female boxcar + ENV
                stimN{2}=stim{2}; % right/male boxcar + ENV
            else % removeSpeechENV==1 % exclude speech envelope 
                stimN{1}=stim{3}; % left/female boxcar
                stimN{2}=stim{4}; % right/male boxcar
            end
        end
        nirsN=nirs{nirsType}(:,:,COIs); % selected channels (tr,smp,ch)
        
        %%% concatenate 1-min trials (if requested)
        if trialDur>1
            nbTrs=floor(size(stimN{1},1)/trialDur);
            disp(['        New number of trials: ' num2str(nbTrs) ' (each lasting ' num2str(trialDur) ' min)']);
            for indTr=1:nbTrs % loop for trial sets
                temp=(stimN{1}(1+trialDur*(indTr-1):trialDur+trialDur*(indTr-1),:))'; % select a set of trialDur trials
                tempStim{1}(indTr,:)=temp(:); % concatenate
                temp=(stimN{2}(1+trialDur*(indTr-1):trialDur+trialDur*(indTr-1),:))'; % select a set of trialDur trials
                tempStim{2}(indTr,:)=temp(:); % concatenate
                for indCh=1:size(nirsN,3) % loop for channels
                    temp=(nirsN(1+trialDur*(indTr-1):trialDur+trialDur*(indTr-1),:,indCh))'; % select a set of trialDur trials
                    tempNirs(indTr,:,indCh)=temp(:); % concatenate
                end % loop for channels
            end % loop for trial sets            
            stimN=tempStim;
            nirsN=tempNirs; 
           % clear temp tempStim tempNirs
        end
        nbTrs=size(stimN{1},1);
        
        %%% lowpass filter stimuli and NIRS
        if lpCutoff>0
            for indTr=1:nbTrs
                if removeSpeechENV~=1 && pilotType~=2
                    for indCh=1:size(stimN,2)
                        stimN{indCh}(indTr,:)=filtfilt(b,a,stimN{indCh}(indTr,:));
                    end
                end
                for indCh=1:size(nirsN,3)
                    nirsN(indTr,:,indCh)=filtfilt(b,a,nirsN(indTr,:,indCh));
                end
            end
        end

        %%% z-score all trials
        for indTr=1:nbTrs
            for indCh=1:size(stimN,2)
                stimN{indCh}(indTr,:)=zscore(stimN{indCh}(indTr,:),0,2);
            end
            for indCh=1:size(nirsN,3) % loop for channels
                nirsN(indTr,:,indCh)=zscore(nirsN(indTr,:,indCh),0,2);
            end % loop for channels
        end
        
        %%% compute backward model
        warning off;
        if pilotType==1;stimulusType=1;end % pilot 1 contained only a single speaker
        if ~quickMode % lambdaOptimization on training trials (slow)
            for indTrTest=1:nbTrs % loop for test trials

                %% lambda optimization (sample a large range of lambdas) on training trials
                indTrTrain=0;
                for indTr=setdiff(1:nbTrs,indTrTest) % loop for training trials
                    indTrTrain=indTrTrain+1;
                    if stimulusType==1 % female
                        STIM{indTrTrain}=stimN{1}(indTr,:)'; % female
                    elseif stimulusType==2 % male
                        STIM{indTrTrain}=stimN{2}(indTr,:)'; % male
                    elseif stimulusType==3 % female-male
                        STIM{indTrTrain}=stimN{1}(indTr,:)' - stimN{2}(indTr,:)'; % female-male
                    else % stimulusType==4 % [female male]
                        STIM{indTrTrain}=[stimN{1}(indTr,:)' stimN{2}(indTr,:)']; % [female male]
                    end
                    DATA{indTrTrain}=squeeze(nirsN(indTr,:,:));
                end % loop for training trials
                [~,~,MSerror]=mTRFcrossval(STIM,DATA,fsNIRS,-1,lagRange(1)*1000,lagRange(end)*1000,lambdaRange); % compute backward model to obtain R
                MSerrorAvg=squeeze(nanmean(MSerror,1)); % residual (averaged across training trials)
                %figure;plot([0 linspace(0,20,10)],MSerrorAvg);hold on;xlabel('lambda');ylabel('MSE'); % plot residual as function of lambda 
                if stimulusType<4;MSerrorAvg=MSerrorAvg';end  % single predictor
                for i=1:size(STIM{1},2) % loop for predictors
                    [~,ind]=min(MSerrorAvg(:,i)); % minimum residual (associated with optimal lambda)
                    lambdaOpt(i)=lambdaRange(ind); % optimal lambda
                end % loop for predictors
               % clear indTrTrain indTr STIM DATA MSerror MSerrorAvg ind i

                %% train backward model
                indTrTrain=0;
                for indTr=setdiff(1:nbTrs,indTrTest) % loop for training trials
                    indTrTrain=indTrTrain+1;
                    if stimulusType==1 % female
                        STIM=stimN{1}(indTr,:)'; % female
                    elseif stimulusType==2 % male
                        STIM=stimN{2}(indTr,:)'; % male
                    elseif stimulusType==3 % female-male
                        STIM=stimN{1}(indTr,:)' - stimN{2}(indTr,:)'; % female-male
                    else % stimulusType==4 % [female male]
                        STIM=[stimN{1}(indTr,:)' stimN{2}(indTr,:)']; % [female male]
                    end
                    DATA=squeeze(nirsN(indTr,:,:));
                    for i=1:size(STIM,2) % loop for predictors
                        [modelTr(:,:,indTrTrain,i),~,biasTr(:,indTrTrain,i)]=mTRFtrain(STIM(:,i),DATA,fsNIRS,-1,lagRange(1)*1000,lagRange(end)*1000,lambdaOpt(i)); % ch,smp,tr
                    end % loop for predictors
                end % loop for training trials
                model(indTrTest,:,:,:)=squeeze(mean(modelTr,3)); % averaged across trials (tr,ch,smp)
                bias=squeeze(mean(biasTr,2)); % averaged across trials (ch)
                %clear indTr indTrTrain STIM DATA modelTr biasTr

                %% validate backward model on test trial
                if stimulusType==1 % female
                    STIM=stimN{1}(indTrTest,:)'; % female
                elseif stimulusType==2 % male
                    STIM=stimN{2}(indTrTest,:)'; % male
                elseif stimulusType==3 % female-male
                    STIM=stimN{1}(indTrTest,:)' - stimN{2}(indTrTest,:)'; % female-male
                else % stimulusType==4 % [female male]
                    STIM=[stimN{1}(indTrTest,:)' stimN{2}(indTrTest,:)']; % [female male]
                end
                DATA=squeeze(nirsN(indTrTest,:,:));
                [~,R(indTrTest,:),~,~]=mTRFpredict(STIM,DATA,squeeze(model(indTrTest,:,:,:)),fsNIRS,-1,lagRange(1)*1000,lagRange(end)*1000,bias);
                %clear STIM DATA bias

            end % loop for test trials
            Rs{nirsType}=nanmean(atanh(R),1); % Fisher-transformed correlation coefficient (averaged across test trials)
            models{nirsType}=squeeze(mean(model,1)); % decoder (averaged across test trials) (ch,smp)
            %clear stimN nirsN R model
            
        else % quickMode % lambdaOptimization on all trials (fast)

        
        %% Monte Carlo Simulation Parameters
        for j = 1:simTrial
            female = stim{1}; male = stim{2}; %get stim from cell array
          
            %shuffle the ROWS of the trialsXsamples matrix in order to keep permutations within each trial
            shuffStimFEM = female(randperm(size(female,1)),:); 
            shuffStimMALE = male(randperm(size(male,1)),:);
          
            %put shuffled stimuli into cell array
            stimSIM{1} = shuffStimFEM;
            stimSIM{2} = shuffStimMALE;
            %% lambda optimization (sample a large range of lambdas) on all trials (fast)
 
            for indTr=1:nbTrs % loop for all trials
                if stimulusType==1 % female
                    STIM{indTr}=stimSIM{1}(indTr,:)'; % female
                elseif stimulusType==2 % male
                    STIM{indTr}=stimSIM{2}(indTr,:)'; % male
                elseif stimulusType==3 % female-male
                    STIM{indTr}=stimSIM{1}(indTr,:)' - stimSIM{2}(indTr,:)'; % female-male
                else % stimulusType==4 % [female male]
                    STIM{indTr}=[stimSIM{1}(indTr,:)' stimSIM{2}(indTr,:)']; % [female male]
                end
                DATA{indTr}=squeeze(nirsN(indTr,:,:));
            end % loop for training trials
            [Rtr,~,MSerror,~,model]=mTRFcrossval(STIM,DATA,fsNIRS,-1,lagRange(1)*1000,lagRange(end)*1000,lambdaRange); % compute forward model to obtain R
            Ravg=squeeze(nanmean(atanh(Rtr),1)); % Fisher-transformed correlation coefficient (averaged across all trials)
            MSerrorAvg=squeeze(nanmean(MSerror,1)); % residual (averaged across all trials)
            %figure;plot([0 linspace(0,20,10)],MSerrorAvg);hold on;xlabel('lambda');ylabel('MSE'); % plot residual as function of lambda 
            if stimulusType<4;Ravg=Ravg';MSerrorAvg=MSerrorAvg';end % single predictor
            for i=1:size(STIM{1},2) % loop for predictors
                [~,ind]=min(MSerrorAvg(:,i)); % minimum residual (associated with optimal lambda)
                Rs{nirsType}(i)=Ravg(ind,i);  % R associated with optimal lambda
                models{nirsType}(i,:)=squeeze(mean(model(:,ind,:,i),1))'; % TRF (averaged across all trials) associated with optimal lambda (ch*smp)
            end % loop for predictors
           % clear indTr STIM DATA Rtr Ravg MSerror MSerrorAvg ind
           % clear stimN nirsN R model

        end
        
        
        %%% plot R 
        disp(['        Average stimulus reconstruction accuracy associated with ' stimulusTypeLabels{stimulusType} ': ' num2str(Rs{nirsType}) '    (Pearsons R)']);
        toc
        
    end % loop for nirsTypes

    %%% export R and model
%     save([folderSubj 'S0' num2str(indSubj,'%3.f') '_backwResults_' num2str(trialDur) 'tr.mat'],'Rs','models','COIs'); 

end % loop for subjects

end

%% Notes
%shuffledStimFEM = randperm(numel(female(:,:))); %shuffle the indeces into a vector
%shuffledStimMALE = randperm(numel(male(:,:)));
%permFEMstim = female(shuffledStimFEM); %apply shuffled indeces to whole time series
%permMALEstim = male(shuffledStimMALE);