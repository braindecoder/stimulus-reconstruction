%%% spectral analysis of stimuli and FNIRS data

clear all;close all;clc;
folder='C:\Users\L.Riecke\Dropbox\FPN\Internship Ethan\';
subjInds=[1 2];

%%% select data
sigTypeLabels={'stim','nirs'};
nirsTypeLabels={'oxyRaw','deoxyRaw','oxyGlm','deoxyGlm'};
nirsTypes=[1 2]; % 1:oxyRaw 2:deoxyRaw 3:oxyGlm 4:deoxyGlm {def:[1 2]}
COIs=[155];        % channels of interest, [] uses all direct-neighbor channels, 0 uses all channels

%%% set timing parameters 
trialDur=60;     % trial duration (in integer minutes), if >1: append consecutive trials {def:2}

%%% import all or direct-neighbor channels (if requested)
if isempty(COIs) % import all direct-neighbor channels
    load([folder '09_Data_after_cleaning\channelMask_directNeighbors_pilot.mat']); % imports 'channelMask'
    channelMask=channelMask';channelMask=channelMask(:)';
    COIs=find(channelMask);
    clear channelMask
elseif COIs==0 % import all channels
    load([folder '09_Data_after_cleaning\channelMask_directNeighbors_pilot.mat']); % imports 'channelMask'
    COIs=1:numel(channelMask);
    clear channelMask
end

for indSubj=subjInds % loop for subjects
    
    disp([' ']);
    disp(['Subject ' num2str(indSubj)]);
    
    %%% import preprocessed FNIRS data and stimuli
    folderSubj=[folder '09_Data_after_cleaning\S' num2str(indSubj,'%02.f') '_pilot\'];
    load([folderSubj 'S' num2str(indSubj,'%02.f') '_pilot_stimuliNIRSproc.mat']); % imports 'stim(tr,smp)','nirs{nirsType}(tr,smp,ch)','fsNIRS'

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
        
        %%% z-score all trials
        for indTr=1:nbTrs
            stimN(indTr,:)=zscore(stimN(indTr,:),0,2);
            nirsN(indTr,:)=zscore(nirsN(indTr,:),0,2);
        end
        
        %%% compute single-trial spectra
        for indTr=1:nbTrs % loop for test trials
            s{1}=stimN(indTr,:);
            s{2}=nirsN(indTr,:);
            for i=1:2
                ss=s{i};
                ss=ss.*hann(numel(ss))';
                nfft=numel(ss);
                sDFT=fft(ss);                        % complex two-sided amplitude spectrum
                ampSpec=abs(sDFT/nfft);              % two-sided amplitude spectrum
                ampSpec=ampSpec(1:nfft/2+1);         % one-sided amplitude spectrum
                ampSpec(2:end-1)=2*ampSpec(2:end-1); % DC and Nyquist freq do not occur twice
                rmsSpec=ampSpec/sqrt(2);             % one-sided energy spectrum
                powSpec{i}(indTr,:)=rmsSpec.^2;      % one-sided power spectrum
                f=0:fsNIRS/nfft:fsNIRS/2;            % frequency vector
                clear ss nfft sDFT ampSpec rmsSpec 
            end
            clear s
        end % loop for test trials
            
        %%% plot trial-averaged spectra
        [~,indF]=min(abs(f-ceil(fsNIRS/2)));
        f(indF+1:end)=[];
        for i=1:2
            powSpec{i}=mean(powSpec{i}(:,1:indF),1);
        end
        if nirsType==nirsTypes(1)
            figure;plot(f,(powSpec{1}));hold on;xlabel('frequency [Hz]');ylabel('power');title(['S' num2str(indSubj) ' ' sigTypeLabels{1}]); % stim
        end
        figure;plot(f,(powSpec{2}));hold on;xlabel('frequency [Hz]');ylabel('power');title(['S' num2str(indSubj) ' ' nirsTypeLabels{nirsType}]); % nirs
        clear f powSpec stimN nirsN
        
    end % loop for nirsTypes
end % loop for subjects