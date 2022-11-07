addpath(genpath('E:\Box Sync\ECoG_Recon\matlab_code\'));
global DUKEDIR
DUKEDIR = 'E:\microECoG surgeries\S14 surgery';
dLabels = dir(DUKEDIR);
dLabels = dLabels(3:end);

tw = [-2 2]; % time window
etw = [-1.5 1.5];
etwG = [-1 1];% epoch time window
prtw = [-1.5 -1]; % preonset time window
pstw = [0.25 0.75]; % postonset time window
gammaF = [70 150]; % frequency in Hz
load('channelMap.mat');
selectedChannels = sort(chanMap(~isnan(chanMap)))';
%% data Loading
Experiment = loadExperiment('S14');
fsD = Experiment.processing.ieeg.sample_rate;
Trials = dbTrials('S14',Experiment.recording.recording_day,'Speech_OvertMimeMove');
trialFiles = strcat('\',Experiment.recording.recording_day,'\mat\trialInfo.mat');
[ieegSplit1,~,trigOnset]=trialIEEGUpdate(Trials,selectedChannels,'SpikeResponseOnset','spike',tw.*1000);