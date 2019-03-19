%%% backward modeling of FNIRS data
%%% for choice of lambdaRange, see https://stackoverflow.com/questions/12182063/how-to-calculate-the-regularization-parameter-in-linear-regression

clear all;close all;clc;
folder='/Users/e/Desktop/RDM/';
subjInds=[1 2];
quickMode=1; % 1:lambdaOptimization on all trials (fast), 0:lambdaOptimization on training trials (slow)

%%% select data
nirsTypeLabels={'oxyRaw','deoxyRaw','oxyGlm','deoxyGlm'};
nirsTypes=[1 2]; % 1:oxyRaw 2:deoxyRaw 3:oxyGlm 4:deoxyGlm {def:[1 2]}
COIs=[1:5];         % channels of interest, [] uses all direct-neighbor channels, 0 uses all channels
lpCutoff=0.5;    % cutoff frequency of lowpass filter {def:0.5Hz}, 0 applies no filtering

%%% set timing parameters 
trialDur=6;              % trial duration (in integer minutes), if >1: append consecutive trials {def:2}
lagRange=0:20;           % lags for TRF analysis (in sec), should span the HRF {def:0:20}
lambdaRange=[0 2.^(linspace(0,50,10))]; % range of lambdas to be sampled 

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
        nirsN=nirs{nirsType}(:,:,COIs); % selected channels (tr,smp,ch)
        
        %%% concatenate 1-min trials (if requested)
        if trialDur>1
            nbTrs=floor(size(stimN,1)/trialDur);
            disp(['   New number of trials: ' num2str(nbTrs) ' (each lasting ' num2str(trialDur) ' min)']);
            for indTr=1:nbTrs % loop for trial sets
                temp=(stimN(1+trialDur*(indTr-1):trialDur+trialDur*(indTr-1),:))'; % select a set of trialDur trials
                tempStim(indTr,:)=temp(:); % concatenate
                for indCh=1:size(nirsN,3) % loop for channels
                    temp=(nirsN(1+trialDur*(indTr-1):trialDur+trialDur*(indTr-1),:,indCh))'; % select a set of trialDur trials
                    tempNirs(indTr,:,indCh)=temp(:); % concatenate
                end % loop for channels
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
                for indCh=1:size(nirsN,3)
                    nirsN(indTr,:,indCh)=filtfilt(b,a,nirsN(indTr,:,indCh));
                end
            end
        end

        %%% z-score all trials
        for indTr=1:nbTrs
            stimN(indTr,:)=zscore(stimN(indTr,:),0,2);
            for indCh=1:size(nirsN,3) % loop for channels
                nirsN(indTr,:,indCh)=zscore(nirsN(indTr,:,indCh),0,2);
            end % loop for channels
        end
        
        %%% compute backward model
        warning off;
        if ~quickMode % lambdaOptimization on training trials (slow)
            for indTrTest=1:nbTrs % loop for test trials

                %% lambda optimization (sample a large range of lambdas) on training trials
                indTrTrain=0;
                for indTr=setdiff(1:nbTrs,indTrTest) % loop for training trials
                    indTrTrain=indTrTrain+1;
                    STIM{indTrTrain}=stimN(indTr,:)';
                    DATA{indTrTrain}=squeeze(nirsN(indTr,:,:));
                end % loop for training trials
                Rtr=mTRFcrossval(STIM,DATA,fsNIRS,-1,lagRange(1)*1000,lagRange(end)*1000,lambdaRange); % compute forward model to obtain R
                Ravg=nanmean(atanh(Rtr),1); % Fisher-transformed correlation coefficient (averaged across training trials)
                %figure;plot(Ravg);hold on;xlabel('lambda');ylabel('R'); % plot R as function of lambda 
                [~,ind]=max(Ravg); 
                lambdaOpt=lambdaRange(ind); % optimal lambda
                clear indTrTrain indTr STIM DATA Rtr Ravg ind

                %% train backward model
                indTrTrain=0;
                for indTr=setdiff(1:nbTrs,indTrTest) % loop for training trials
                    indTrTrain=indTrTrain+1;
                    STIM=stimN(indTr,:)';
                    DATA=squeeze(nirsN(indTr,:,:));
                    [modelTr(:,:,indTr),~,biasTr(:,indTr)]=mTRFtrain(STIM,DATA,fsNIRS,-1,lagRange(1)*1000,lagRange(end)*1000,lambdaOpt); % ch,smp,tr
                end % loop for training trials
                model(indTrTest,:,:)=mean(modelTr,3); % averaged across trials (tr,ch,smp)
                bias=mean(biasTr,2); % averaged across trials (ch)
                clear indTrTrain indTr STIM DATA modelTr biasTr

                %% validate backward model on test trial
                STIM=stimN(indTrTest,:)';
                DATA=squeeze(nirsN(indTrTest,:,:));
                [~,R(indTrTest),~,~]=mTRFpredict(STIM,DATA,squeeze(model(indTrTest,:,:)),fsNIRS,-1,lagRange(1)*1000,lagRange(end)*1000,bias);
                clear STIM DATA bias

            end % loop for test trials
            Rs{nirsType}=nanmean(atanh(R)); % Fisher-transformed correlation coefficient (averaged across test trials)
            models{nirsType}=squeeze(mean(model,1)); % decoder (averaged across test trials) (ch,smp)
            clear stimN nirsN R model
            
        else % quickMode % lambdaOptimization on all trials (fast)

            %% lambda optimization (sample a large range of lambdas) on all trials
            for indTr=1:nbTrs % loop for all trials
                STIM{indTr}=stimN(indTr,:)';
                DATA{indTr}=squeeze(nirsN(indTr,:,:));
            end % loop for training trials
%             [Rtr,~,MSerror,~,model]=mTRFcrossval(STIM,DATA,fsNIRS,-1,lagRange(1)*1000,lagRange(end)*1000,lambdaRange); % compute forward model to obtain R
            [Rtr,~,MSerror,~,model]=mTRFcrossval(STIM,DATA,fsNIRS,-1,lagRange(1)*1000,lagRange(end)*1000,lambdaRange); % compute forward model to obtain R
            Ravg=nanmean(atanh(Rtr),1); % Fisher-transformed correlation coefficient (averaged across all trials)
            MSerrorAvg=nanmean(MSerror,1); % residual (averaged across all trials)
            %figure;plot([0 linspace(0,20,10)],MSerrorAvg);hold on;xlabel('lambda');ylabel('MSE'); % plot residual as function of lambda 
            [~,ind]=min(MSerrorAvg); % minimum residual (associated with optimal lambda)
            Rs{nirsType}=Ravg(ind);  % R associated with optimal lambda
%             models{nirsType}=squeeze(mean(model(:,ind,:),1))'; % TRF (averaged across all trials) associated with optimal lambda (ch*smp)
            clear indTr STIM DATA Rtr Ravg MSerror MSerrorAvg ind
            clear stimN nirsN R model

        end
        
        
        %%% plot R 
        disp(['   Average stimulus reconstruction accuracy (Pearsons R): ' num2str(Rs{nirsType})]);
        
    end % loop for nirsTypes

    %%% export R and model
    save([folderSubj 'S' num2str(indSubj,'%02.f') '_backwResults_' num2str(trialDur) 'tr.mat'],'Rs','models','COIs'); 

end % loop for subjects