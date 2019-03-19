%%% import stimuli and preprocessed FNIRS data from FNIRS experiment
%%% prepares the data for forward/backward modeling

clear all;close all;clc;

folderIN='/Users/e/Desktop/RDM/07_Raw_data/';
folderOUT='/Users/e/Desktop/RDM/Output';
subjInds=[1 2];
onsetDur=30; % duration of stimulation onset interval (to be rejected) {def:30s}

for indSubj=subjInds % loop for subjects
    
    disp(['Subject ' num2str(indSubj)]);
    
    
    %% import stimuli from RDM server
    load([folderIN 'S' num2str(indSubj,'%02.f') '_pilot/S' num2str(indSubj,'%02.f') '_pilot_stimuli.mat']);
    stim=s;
    nbTr=ceil((numel(stim)/fs)/triggerPeriod)+1; % number of triggers (including end trigger) that were supposed to be sent out
    clear trigger gapDuration s
    
    
    %% import FNIRS parameters from original HDR file on RDM server
    hdrFile=dir([folderIN 'S' num2str(indSubj,'%02.f') '_pilot/FNIRS_data/*_orig.hdr']);
    fid=fopen([folderIN 'S' num2str(indSubj,'%02.f') '_pilot/FNIRS_data/' hdrFile(1).name]); % open HDR file
    lnText='             ';
    while ~strcmp(lnText(1:13),'SamplingRate=')
        lnTextN=fgetl(fid); % skip lines
        if length(lnTextN)>=13;lnText=lnTextN;end
    end
    fsNIRS=str2num(lnText(14:end)); % sampling rate starts at position 14
    while ~strcmp(lnText,'Events="#')
        lnText=fgetl(fid); % skip lines
    end
    indTr=0;
    while ~strcmp(lnText,'#"')
        indTr=indTr+1;
        lnText=fgetl(fid); % get trigger lines
        if ~strcmp(lnText,'#"')
            lnTextCell=strsplit(lnText);
            triggerNumber(indTr)=str2num(lnTextCell{2}); % 2nd column contains trigger value
            triggerOnset(indTr)=str2num(lnTextCell{3});  % 3rd column contains trigger onset in samples
        end
    end
    clear lnText lnTextN lnTextCell
    fclose(fid); % close HDR file
    
    
    %% import preprocessed FNIRS data from RDM server and truncate 
    % Format of NIRS data:
    %   nirs{1}=AllDataOxy;     % oxy raw (time*channel)
    %   nirs{2}=AllDataDeOxy;   % deoxy raw (time*channel)
    %   nirs{3}=residualOxy;    % oxy glm (time*channel)
    %   nirs{4}=residualDeOxy;  % deoxy glm (time*channel)
    load([folderIN 'S' num2str(indSubj,'%02.f') '_pilot/S' num2str(indSubj,'%02.f') '_pilot_nirs.mat']); % imports 'nirs','TSIout'
    for nirsType=1:4 
        smpOfLastTrigger=triggerOnset(end);
        lastSmpBeforeFirstTrigger=triggerOnset(1)-1;
        if nirsType>=3 % glm data include no baseline
            smpOfLastTrigger=smpOfLastTrigger-TSIout.OxyDeOxyBaselineEnd;
            lastSmpBeforeFirstTrigger=lastSmpBeforeFirstTrigger-TSIout.OxyDeOxyBaselineEnd;
        end
        nirs{nirsType}(smpOfLastTrigger:end,:)=[];         % remove samples starting with last trigger
        nirs{nirsType}(1:lastSmpBeforeFirstTrigger,:)=[];  % remove samples before first trigger
    end
    clear smpOfLastTrigger lastSmpBeforeFirstTrigger
    lengthNIRSdata=size(nirs{1},1); % get length of FNIRS data (from onset trigger until end trigger-1)
    
    
    %% check recorded triggers for errors or jitter 
    figure;hist(diff(triggerOnset(1:end-1)),min(diff(triggerOnset(1:end-1))):max(diff(triggerOnset(1:end-1)))); % visualize jitter in trial length 
    hold on;xlabel('Trial length');ylabel('Number of trials'); 
    if triggerNumber(1)~=2 || ...        % missed first trigger 
       triggerNumber(end)~=8 || ...      % missed end trigger (or some other trigger)
       numel(triggerNumber)~=nbTr|| ...  % number of recorded triggers ~= number of triggers that were supposed to be sent out
       any(diff(triggerNumber)==0)       % repeated a trigger
            disp('    ERROR: Missing trigger(s) in FNIRS data!');
    elseif triggerOnset(end)<lengthNIRSdata
            disp('    ERROR: Preprocessed FNIRS data are too short!');
    elseif triggerOnset(end)>lengthNIRSdata
            disp('    ERROR: Preprocessed FNIRS data are too long!');
    end
    
    
    %% extract stimulus envelope and resample to match the length of FNIRS data
    stim=abs(hilbert(stim)); % extract envelope (or onset envelope) 
    stim=[0 diff(stim)];     % extract first derivative (slope)
    stim=max(stim,0);        % half-wave rectify
    stim=resample(stim,lengthNIRSdata,numel(stim)); % resampling includes anti-aliasing filter
    
    
    %% stimuli and FNIRS data: remove stimulation onset interval
    stim(1:round(onsetDur*fsNIRS))=[];   % remove stimulation onset interval
    for nirsType=1:size(nirs,2) % loop for FNIRS signal types
        nirs{nirsType}(1:round(onsetDur*fsNIRS),:)=[]; % remove stimulation onset interval
    end % loop for FNIRS signal types
    
    
%     %% apply dummy HRF to 1st channel
%     t=0:1/fsNIRS:25; % HRF interval {def:25s}
%     T0=0;n=4;lambda=2;
%     hrf=((t-T0).^(n-1)).*exp(-(t-T0)/lambda)/((lambda^n)*factorial(n-1)); % define HRF
%     temp=conv(stim,hrf); % convolve stimulus with HRF
%     for nirsType=1:size(nirs,2) % loop for FNIRS signal types
%         nirs{nirsType}(:,1)=temp(1:size(nirs{nirsType},1)); % insert convolved stimulus in 1st channel 
%     end % loop for FNIRS signal types
%     clear t T0 n lambda hrf temp
    
    
    %% stimuli and FNIRS data: cut into 1-min trials
    trLength=round(triggerPeriod*fsNIRS);
    nbTr=floor(length(stim)/trLength);
    for indTr=1:nbTr
        stim2(indTr,:)=stim(1+(indTr-1)*trLength:trLength+(indTr-1)*trLength);
    end
    for nirsType=1:size(nirs,2) % loop for FNIRS signal types
        for indTr=1:nbTr
            nirs2{nirsType}(indTr,:,:)=nirs{nirsType}(1+(indTr-1)*trLength:trLength+(indTr-1)*trLength,:); % tr,smp,ch
        end
    end % loop for FNIRS signal types
    stim=stim2;nirs=nirs2; 
    clear triggerPeriod trLength nbTr stim2 nirs2

    
    %% save preprocessed stimuli and FNIRS data
    folderOUTsubj=[folderOUT '09_Data_after_cleaning\S' num2str(indSubj,'%02.f') '_pilot\'];
    if ~exist(folderOUTsubj);mkdir(folderOUTsubj);end
    save([folderOUTsubj 'S' num2str(indSubj,'%02.f') '_pilot_stimuliNIRSproc.mat'],'stim','nirs','fsNIRS');
    clear folderOUTsubj
    
end % loop for subjects

