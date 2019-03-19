%%% exports preprocessed FNIRS data from Turbo Satori 

clear all;close all;clc;

folderIN='\\ca-um-nas201\fpn_rdm$\DM0934_LR_SpeechDec\07_Raw_data\';
folderOUT='\\ca-um-nas201\fpn_rdm$\DM0934_LR_SpeechDec\07_Raw_data\';
subjInd=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% fix HDR file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% backup and rename original HDR file
hdrFile=dir([folderIN 'S' num2str(subjInd,'%02.f') '_pilot\FNIRS_data\*.hdr']);
copyfile([folderIN 'S' num2str(subjInd,'%02.f') '_pilot\FNIRS_data\' hdrFile(1).name], ...
         [folderIN 'S' num2str(subjInd,'%02.f') '_pilot\FNIRS_data\' hdrFile(1).name(1:end-4) '_orig.hdr']); 

%%% replace original triggers by single dummy task trigger (to enable GLM in TSI)
%%% TSI takes the the 1st trigger ('1') to be the onset of the baseline
%%% TSI takes the the 2nd trigger ('2') to be the onset of the first condition
fid=fopen([folderIN 'S' num2str(subjInd,'%02.f') '_pilot\FNIRS_data\' hdrFile(1).name]); % open HDR file
lnText='             ';
while ~strcmp(lnText(1:13),'SamplingRate=')
    lnTextN=fgetl(fid); % skip lines
    if length(lnTextN)>=13;lnText=lnTextN;end
end
fsNIRS=str2num(lnText(14:end)); % sampling rate starts at position 14
fclose(fid);
str=fileread([folderIN 'S' num2str(subjInd,'%02.f') '_pilot\FNIRS_data\' hdrFile(1).name]); 
pos(1)=strfind(str,sprintf('Events="#'))+11;       % first digit of first trigger line
pos(2)=strfind(str,sprintf('[DataStructure]'))-9;  % last digit of last trigger line
strN=sprintf('%s\t%s\t%s\n%s\t%s\t%s',num2str(51,'%.2f'),'1',num2str(round(51*fsNIRS)),num2str(51+1/fsNIRS,'%.2f'),'2',num2str(round(51*fsNIRS)+1));
str(pos(1):pos(1)+numel(strN)-1)=strN;
str(pos(1)+numel(strN):pos(2))=[];
fid=fopen([folderIN 'S' num2str(subjInd,'%02.f') '_pilot\FNIRS_data\' hdrFile(1).name],'wt');
fprintf(fid,'%s',str);
fclose(fid);
clear hdrFile fsNIRS str pos strN fid


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% import from TurboSatori
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% do these manual steps:
% -launch TSI
% -set parameters in TSI
%     -Data Read Frequency: max
%     -Hb/Hbo Baseline from 5 to 50 (requires that stimulation was started >60s after FNIRS recording onset)
%     -Preprocessing applied to Hbo/Hb data (not WL data): only check on "linear detrend"
%     -GLM parameters:
%           -Add Confound Predictors: only check on "Channel-confound predictors"
%           -Residuals: checked off
%           -GLM: check on "Use fully preprocessed signal for GLM"
% -import adjusted HDR file from folderIN
% -select all channels

%%% connect Matlab with TSI
configs.TSI_IP = 'localhost'; % run Turbo-Satori on the same Computer
configs.TSI_PORT = 55555;     % default port
tsiNetInt = TSINetworkInterface( TSIClient( configs.TSI_IP, configs.TSI_PORT ) );
tsiNetInt.createConnection();

%%% import variables
timePoint = tsiNetInt.tGetCurrentTimePoint();
if timePoint == 0;warning('TimePoint must be greater zero!');return;end
NrOfChannels = tsiNetInt.tGetNrOfChannels();
OxyDataScaleFactor = tsiNetInt.tGetOxyDataScaleFactor();
OxyDeOxyBaselineEnd = tsiNetInt.tGetOxyDeOxyBaselineEnd(); % last sample of baseline (should be before first stimulus trigger)
SelectedChannels = tsiNetInt.tGetSelectedChannels();       % should be equal to NrOfChannels
SamplingRate = tsiNetInt.tGetSamplingRate();
AllDataOxy = tsiNetInt.tGetAllDataOxy(SelectedChannels, timePoint)'*OxyDataScaleFactor;     % oxy data
AllDataDeOxy = tsiNetInt.tGetAllDataDeOxy(SelectedChannels, timePoint)'*OxyDataScaleFactor; % deOxy data
FullNrOfPredictors = tsiNetInt.tGetFullNrOfPredictors();          % required for GLM
NrOfConfoundPredictors = tsiNetInt.tGetNrOfConfoundPredictors();  % required for GLM
Condition = tsiNetInt.tGetProtocolCondition(timePoint);           % required for GLM

%%% verify that all samples and channels were exported
files=dir([folderIN 'S' num2str(subjInd,'%02.f') '_pilot\FNIRS_data\*.wl1']);
fileID = fopen([folderIN 'S' num2str(subjInd,'%02.f') '_pilot\FNIRS_data\' files(1).name],'r');
formatSpec = '%9f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%f%[^\n\r]';
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
temp = [dataArray{1:end-1}]; % get number of samples and channels from raw light data
if any([size(AllDataOxy,1) size(AllDataDeOxy,1)]~=size(temp,1)) % check for correct number of samples
    error(['ERROR: Not all samples ' num2str(size(AllDataOxy,1)) '/' num2str(size(temp,1)) ' were imported!']);
end
if any([size(AllDataOxy,2) size(AllDataDeOxy,2)]~=size(temp,2)) % check for correct number of channels
    error(['ERROR: Not all channels ' num2str(size(AllDataOxy,2)) '/' num2str(size(temp,2)) ' were imported!']);
end
clear files fileID formatSpec dataArray temp

%%% compute GLMs
AllDataOxyBL = AllDataOxy(OxyDeOxyBaselineEnd:end,:);      % remove baseline samples
AllDataDeOxyBL = AllDataDeOxy(OxyDeOxyBaselineEnd:end,:);  % remove baseline samples

%% GLM for oxy data
DMOxy = zeros(timePoint - OxyDeOxyBaselineEnd,1); % compute design matrix for oxyData
for p = 1:FullNrOfPredictors
   disp(['Oxy: predictor ' num2str(p) '/' num2str(FullNrOfPredictors)]);
   for t = OxyDeOxyBaselineEnd:timePoint
       DMOxy(t-OxyDeOxyBaselineEnd+1,p) = tsiNetInt.tGetValueOfDesignMatrix(p,t,1);
   end
end
betaOxy = inv(DMOxy'*DMOxy)*DMOxy'*AllDataOxyBL;   % compute beta for each oxy condition
predictedOxy = (rot90(betaOxy,3)*rot90(DMOxy,1))'; % compute predicted time course for each oxy condition
residualOxy = AllDataOxyBL-predictedOxy;           % compute residual time course for each oxy condition

%% GLM for deoxy data
DMDeOxy = zeros(timePoint - OxyDeOxyBaselineEnd,1); % compute design matrix for deoxyData
for p = 1:FullNrOfPredictors
   disp(['Deoxy: predictor ' num2str(p) '/' num2str(FullNrOfPredictors)]);
   for t = OxyDeOxyBaselineEnd:timePoint
       DMDeOxy(t-OxyDeOxyBaselineEnd+1,p) = tsiNetInt.tGetValueOfDesignMatrix(p,t,0);
   end
end
betaDeOxy = inv(DMDeOxy'*DMDeOxy)*DMDeOxy'*AllDataDeOxyBL; % compute beta for each deoxy condition
predictedDeOxy = (rot90(betaDeOxy,3)*rot90(DMDeOxy,1))';   % compute predicted time course for each deoxy condition
residualDeOxy = AllDataDeOxyBL-predictedDeOxy;             % compute residual time course for each deoxy condition

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% export imported data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% plot
% figure;hold on;
% plot(AllDataOxy);     % filtered Oxy Data
% plot(AllDataDeOxy);   % filtered DeOxy Data
% plot(residualOxy);    % residual Oxy Data
% plot(residualDeOxy);  % residual DeOxy Data

%%% merge variables into single structure/variable
TSIout.timePoint=timePoint;
TSIout.NrOfChannels=NrOfChannels;
TSIout.OxyDataScaleFactor=OxyDataScaleFactor;
TSIout.OxyDeOxyBaselineEnd=OxyDeOxyBaselineEnd;
TSIout.SelectedChannels=SelectedChannels;
TSIout.FullNrOfPredictors=FullNrOfPredictors;
TSIout.NrOfConfoundPredictors=NrOfConfoundPredictors; 
TSIout.SamplingRate=SamplingRate;
TSIout.DMOxy=DMOxy;
TSIout.betaOxy=betaOxy;
TSIout.DMDeOxy=DMDeOxy;
TSIout.betaDeOxy=betaDeOxy;
nirs{1}=AllDataOxy;     % oxy raw (time*channel)
nirs{2}=AllDataDeOxy;   % deoxy raw (time*channel)
nirs{3}=residualOxy;    % oxy glm (time*channel)
nirs{4}=residualDeOxy;  % deoxy glm (time*channel)

%%% save to RDM server
save([folderOUT 'S' num2str(subjInd,'%02.f') '_pilot\S' num2str(subjInd,'%02.f') '_pilot_nirs.mat'],'nirs','TSIout');
