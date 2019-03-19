%%% forward modeling of FNIRS data
%%% for choice of lambdaRange, see https://stackoverflow.com/questions/12182063/how-to-calculate-the-regularization-parameter-in-linear-regression

clear all;close all;clc;
folder='/Users/e/Desktop/RDM/';
subjInds=[1 2];
quickMode=1; % 1:lambdaOptimization on all trials (fast), 0:lambdaOptimization on training trials (slow)

%%% select data
nirsTypeLabels={'oxyRaw','deoxyRaw','oxyGlm','deoxyGlm'};
nirsTypes=[1 2]; % 1:oxyRaw 2:deoxyRaw 3:oxyGlm 4:deoxyGlm {def:[1 2]}
COIs=[1 6];       % channels of interest, [] uses all direct-neighbor channels, 0 uses all channels
lpCutoff=0.8;    % cutoff frequency of lowpass filter {def:0.5Hz}, 0 applies no filtering

%%% set timing parameters 
trialDur=1;              % trial duration (in integer minutes), if >1: append consecutive trials {def:2}
lagRange=0:60;           % lags for TRF analysis (in sec), should span the HRF {def:0:20}
lambdaRange=[0 2.^(linspace(0,20,10))]; % range of lambdas to be sampled 

%%% import all or direct-neighbor channels (if requested)
if isempty(COIs) % import all direct-neighbor channels
    load([folder '09_Data_after_cleaning/channelMask_directNeighbors_pilot.mat']); % imports 'channelMask'
    channelMask=channelMask';channelMask=channelMask(:)';
    COIs=find(channelMask);
    clear channelMask
elseif COIs==0 % import all channels
    load([folder '09_Data_after_cleaning/channelMask_directNeighbors_pilot.mat']); % imports 'channelMask'
    COIs=1:numel(channelMask);
    clear channelMask
end

for indSubj=subjInds % loop for subjects
    
    disp([' ']);
    disp(['Subject ' num2str(indSubj)]);
    
    %%% import preprocessed FNIRS data and stimuli
    folderSubj=[folder '09_Data_after_cleaning/S' num2str(indSubj,'%02.f') '_pilot/'];
    load([folderSubj 'S' num2str(indSubj,'%02.f') '_pilot_stimuliNIRSproc.mat']); % imports 'stim(tr,smp)','nirs{nirsType}(tr,smp,ch)','fsNIRS'
    lags=lagRange(1):1/fsNIRS:lagRange(end); % timing parameters for TRF analysis

    %%% design lp filter (if requested)
    if lpCutoff>0
    filterOrder=4; % {def:4}
    [b,a]=butter(filterOrder,lpCutoff*2/fsNIRS,'low');
    end
    
    %%% select desired nirsTypes and channels
    for nirsType=nirsTypes % loop for nirsTypes
        stimN=stim;
        nirsN=mean(nirs{nirsType}(:,:,COIs),3); % average across channels (tr,smp)
        
        %%% concatenate 1-min trials (if requested)
        if trialDur>1
            nbTrs=floor(size(stimN,1)/trialDur);
            disp(['   New number of trials: ' num2str(nbTrs) ' (each lasting ' num2str(trialDur) ' min)']);
            for indTr=1:nbTrs % loop for trial sets
                temp=(stimN(1+trialDur*(indTr-1):trialDur+trialDur*(indTr-1),:))'; % select a set of trialDur trials
                tempStim(indTr,:)=temp(:); % concatenate
                temp=(nirsN(1+trialDur*(indTr-1):trialDur+trialDur*(indTr-1),:))'; % select a set of trialDur trials
                tempNirs(indTr,:)=temp(:); % concatenate
            end % loop for trial sets            
            stimN=tempStim;
            nirsN=tempNirs; 
            clear temp tempStim tempNirs
        end
        nbTrs=size(stimN,1);
        
        %%% lowpass filter stimuli and NIRS
        if lpCutoff>0
            for indTr=1:nbTrs
                stimN(indTr,:)=filtfilt(b,a,stimN(indTr,:));
                nirsN(indTr,:)=filtfilt(b,a,nirsN(indTr,:));
            end
        end

        %%% z-score all trials
        for indTr=1:nbTrs
            stimN(indTr,:)=zscore(stimN(indTr,:),0,2);
            nirsN(indTr,:)=zscore(nirsN(indTr,:),0,2);
        end

        %%% compute forward model
        warning off;
        if ~quickMode % lambdaOptimization on training trials (slow)
            for indTrTest=1:nbTrs % loop for test trials

                %% lambda optimization (sample a large range of lambdas) on training trials
                indTrTrain=0;
                for indTr=setdiff(1:nbTrs,indTrTest) % loop for training trials
                    indTrTrain=indTrTrain+1;
                    STIM{indTrTrain}=stimN(indTr,:)';
                    DATA{indTrTrain}=nirsN(indTr,:)';
                end % loop for training trials
                [~,~,MSerror]=mTRFcrossval(STIM,DATA,fsNIRS,+1,lagRange(1)*1000,lagRange(end)*1000,lambdaRange); % compute forward model to obtain R
                MSerrorAvg=nanmean(MSerror,1); % residual (averaged across training trials)
                %figure;plot([0 linspace(0,20,10)],MSerrorAvg);hold on;xlabel('lambda');ylabel('MSE'); % plot residual as function of lambda 
                [~,ind]=min(MSerrorAvg); 
                lambdaOpt=lambdaRange(ind); % optimal lambda
                clear indTrTrain indTr STIM DATA MSerror MSerrorAvg ind

                %% train forward model
                indTrTrain=0;
                for indTr=setdiff(1:nbTrs,indTrTest) % loop for training trials
                    indTrTrain=indTrTrain+1;
                    STIM=stimN(indTr,:)';
                    DATA=nirsN(indTr,:)';
                    [modelTr(:,indTr),~,biasTr(indTr)]=mTRFtrain(STIM,DATA,fsNIRS,+1,lagRange(1)*1000,lagRange(end)*1000,lambdaOpt);
                end % loop for training trials
                model(indTrTest,:)=mean(modelTr,2); % averaged across trials
                bias=mean(biasTr); % averaged across trials
                clear indTrTrain indTr STIM DATA modelTr biasTr

                %% validate forward model on test trial
                STIM=stimN(indTrTest,:)';
                DATA=nirsN(indTrTest,:)';
                [~,R(indTrTest),~,~]=mTRFpredict(STIM,DATA,model(indTrTest,:),fsNIRS,+1,lagRange(1)*1000,lagRange(end)*1000,bias);
                clear STIM DATA bias

            end % loop for test trials
            Rs{nirsType}=nanmean(atanh(R)); % Fisher-transformed correlation coefficient (averaged across test trials)
            models{nirsType}=mean(model,1); % TRF (averaged across test trials)
            clear stimN nirsN R model
        
        else % quickMode % lambdaOptimization on all trials (fast)
            
            %% lambda optimization (sample a large range of lambdas) on all trials
            for indTr=1:nbTrs % loop for all trials
                STIM{indTr}=stimN(indTr,:)';
                DATA{indTr}=nirsN(indTr,:)';
            end % loop for training trials
            [Rtr,~,MSerror,~,model]=mTRFcrossval(STIM,DATA,fsNIRS,+1,lagRange(1)*1000,lagRange(end)*1000,lambdaRange); % compute forward model to obtain R
            Ravg=nanmean(atanh(Rtr),1); % Fisher-transformed correlation coefficient (averaged across all trials)
            MSerrorAvg=nanmean(MSerror,1); % residual (averaged across all trials)
            %figure;plot([0 linspace(0,20,10)],MSerrorAvg);hold on;xlabel('lambda');ylabel('MSE'); % plot residual as function of lambda 
            [~,ind]=min(MSerrorAvg); % minimum residual (associated with optimal lambda)
            Rs{nirsType}=Ravg(ind);  % R associated with optimal lambda
            models{nirsType}=squeeze(mean(model(:,ind,:),1))'; % TRF (averaged across all trials) associated with optimal lambda
            clear indTr STIM DATA Rtr Ravg MSerror MSerrorAvg ind
            clear stimN nirsN R model
        
        end
        
            
        %%% plot R and TRF (if TRF differs from HRF: fine-tune timing parameters)
        disp(['   Average EEG prediction accuracy (Pearsons R): ' num2str(Rs{nirsType})]);
        figure;plot(linspace(lagRange(1),lagRange(end),size(models{nirsType},2)),models{nirsType});hold on;
        xlabel('lag [s]');ylabel('TRF magnitude');title(['S' num2str(indSubj) ' ' nirsTypeLabels{nirsType}]);
        
    end % loop for nirsTypes

    %%% export R and TRF
    save([folderSubj 'S' num2str(indSubj,'%02.f') '_forwResults_' num2str(trialDur) 'tr.mat'],'Rs','models','COIs'); 

end % loop for subjects