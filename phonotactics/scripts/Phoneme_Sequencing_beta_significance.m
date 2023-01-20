global BOX_DIR
global DUKEDIR
saveFolder = '\TempDecode\Articulator\sensorimotor_delay_sig_v2\';
BOX_DIR = 'C:\Users\sd355\Box'
DUKEDIR = 'C:\Users\sd355\Box\CoganLab\D_Data'


Task=[];

Task.Name='Phoneme_Sequencing';
TASK_DIR=([BOX_DIR '/CoganLab/D_Data/' Task.Name]);
DUKEDIR=TASK_DIR;
saveFolder = '\TempDecode\PooledSubjects\sensorimotor\syllable\';
%% Loading data
Subject = popTaskSubjectData(Task);
subjectIds2remove = [24 27 31:length(Subject) ];
% removing D18 because of large negative response trials
Subject(subjectIds2remove) = [];

%% Loading Normalized High-Gamma for an ROI and corresponding trialInfo
fieldEpoch = 'ResponseStart';
selectRoi = '';
fieldTime = [0 0.5];
respTimeThresh = 0;

%trialInfoStruct = extractTrialInfo(Subject, remFastResponseTimeTrials=respTimeThresh);

ieegFilterHiFreqBase = extractBandPassDataWithROI(Subject,...
    Epoch = 'Start', roi = selectRoi,Time= [-0.5 0],...
    respTimeThresh=respTimeThresh,...
    fband = [250 500], remWMchannels=false);

ieegFilterHiFreqResponse = extractBandPassDataWithROI(Subject,...
    Epoch = fieldEpoch, roi = selectRoi,Time= fieldTime,...
    respTimeThresh=respTimeThresh,...
    fband = [250 500],remWMchannels=false);

%% Looking for significance beta
for iSubject = 1:length(ieegFilterHiFreqResponse)
    iSubject
    ieegPowerBase = ieegFilterHiFreqBase(iSubject).ieegPower;
    ieegPowerResponse = ieegFilterHiFreqResponse(iSubject).ieegPower;
    assert(size(ieegPowerBase,1)==size(ieegPowerResponse,1),' Channel Dimension mismatch');
    assert(size(ieegPowerBase,2)==size(ieegPowerResponse,2),' Trial Dimension mismatch');
    
    pSubj = []; 
    for iChan = 1:size(ieegPowerResponse,1)
        pSubj(iChan) = permtest_sk(ieegPowerResponse(iChan,:),ieegPowerBase(iChan,:),10000);
    end
    ieegFilterHiFreqResponse(iSubject).pHiFreq = pSubj;
    [pval, pmasked] = fdr(pSubj,0.05);
    ieegFilterHiFreqResponse(iSubject).p_masked = pmasked;
end

%% Loading Normalized High-Gamma for an ROI and corresponding trialInfo
fieldEpoch = 'ResponseStart';
selectRoi = '';
fieldTime = [0 0.5];
respTimeThresh = 0;

%trialInfoStruct = extractTrialInfo(Subject, remFastResponseTimeTrials=respTimeThresh);

ieegFilterLowFreqBase = extractBandPassDataWithROI(Subject,...
    Epoch = 'Start', roi = selectRoi,Time= [-0.5 0],...
    respTimeThresh=respTimeThresh,...
    fband = [2 10], remWMchannels=false);

ieegFilterLowFreqResponse = extractBandPassDataWithROI(Subject,...
    Epoch = fieldEpoch, roi = selectRoi,Time= fieldTime,...
    respTimeThresh=respTimeThresh,...
    fband = [2 10],remWMchannels=false);

%% Looking for significance beta
for iSubject = 1:length(ieegFilterLowFreqResponse)
    iSubject
    ieegPowerBase = ieegFilterLowFreqBase(iSubject).ieegPower;
    ieegPowerResponse = ieegFilterLowFreqResponse(iSubject).ieegPower;
    assert(size(ieegPowerBase,1)==size(ieegPowerResponse,1),' Channel Dimension mismatch');
    assert(size(ieegPowerBase,2)==size(ieegPowerResponse,2),' Trial Dimension mismatch');
    
    pSubj = []; 
    for iChan = 1:size(ieegPowerResponse,1)
        pSubj(iChan) = permtest_sk(ieegPowerResponse(iChan,:),ieegPowerBase(iChan,:),10000);
    end
    ieegFilterLowFreqResponse(iSubject).pLowFreq = pSubj;
    [pval, pmasked] = fdr(pSubj,0.05);
    ieegFilterLowFreqResponse(iSubject).p_masked = pmasked;
end

%% Loading Normalized High-Gamma for an ROI and corresponding trialInfo
fieldEpoch = 'ResponseStart';
selectRoi = '';
fieldTime = [0 1];
respTimeThresh = 0;

%trialInfoStruct = extractTrialInfo(Subject, remFastResponseTimeTrials=respTimeThresh);

ieegFilterBetaBase = extractBandPassDataWithROI(Subject,...
    Epoch = 'Start', roi = selectRoi,Time= [-0.5 0],...
    respTimeThresh=respTimeThresh,...
    fband = [15 30], remWMchannels=false);

ieegFilterBetaResponse = extractBandPassDataWithROI(Subject,...
    Epoch = fieldEpoch, roi = selectRoi,Time= fieldTime,...
    respTimeThresh=respTimeThresh,...
    fband = [15 30],remWMchannels=false);

%% Looking for significance beta
for iSubject = 1:length(ieegFilterBetaResponse)
    iSubject
    ieegPowerBase = ieegFilterBetaBase(iSubject).ieegPower;
    ieegPowerResponse = ieegFilterBetaResponse(iSubject).ieegPower;
    assert(size(ieegPowerBase,1)==size(ieegPowerResponse,1),' Channel Dimension mismatch');
    assert(size(ieegPowerBase,2)==size(ieegPowerResponse,2),' Trial Dimension mismatch');
    
    pSubj = []; 
    for iChan = 1:size(ieegPowerResponse,1)
        pSubj(iChan) = permtest_sk(ieegPowerBase(iChan,:),ieegPowerResponse(iChan,:),10000);
    end
    ieegFilterBetaResponse(iSubject).pBeta = pSubj;
    [pval, pmasked] = fdr(pSubj,0.05);
    ieegFilterBetaResponse(iSubject).p_masked = pmasked;
end
%%

elecNameMuscleArtifact = [];
for iSubj = 1:length(ieegFilterHiFreqResponse)
elecNameMuscleArtifact = [elecNameMuscleArtifact ieegFilterHiFreqResponse(iSubj).channelName(ieegFilterHiFreqResponse(iSubj).p_masked & ~(ieegFilterBetaResponse(iSubj).p_masked) )]
end


