function [] = get_rich_data()

load('/Users/e/Desktop/RDM/09_Data_after_cleaning/S01_pilot_alternateSpeaker/S01_pilot_stimuliNIRSproc.mat');

xc_fig = figure;
ts_fig = figure;

for i = 1:length(stim)
    
    cur_stim_all = stim{i};
    cur_nirs_all = nirs{i};
    
    k = 61;
    chan = 54;
    
    cur_stim = squeeze(cur_stim_all(k,:,:));
    cur_nirs = squeeze(cur_nirs_all(k,:,chan));
    
    cur_stim = rangenorm(cur_stim);
    cur_nirs = rangenorm(cur_nirs);
    
    cur_stim = smooth(cur_stim,61);
    cur_nirs = smooth(cur_nirs,61);
    
    white_nirs = prewhiten(cur_nirs);
    
    figure(ts_fig);
    subplot(2,2,i);
    plot(cur_stim); hold on; plot(cur_nirs,'r'); hold off;
    
    cur_stim_trim = cur_stim(1:120);
    
    xcorr(cur_nirs,cur_stim_trim,100);
    
    figure(xc_fig);
    subplot(2,2,i);
    xcorr(cur_nirs,cur_stim,100);
    
end

function out = rangenorm(data)

out = (data-min(data))/(max(data)-min(data));

function y1=prewhiten(y,A)
% does smoothing with 
% y is the timeseries, should be in rows
y = y(:);
n=length(y);
if nargin<2,A=0.7;end
A1=eye(n);
for (i=2:n)
    A1(i,i-1)=-A;
end
y1=A1*y;