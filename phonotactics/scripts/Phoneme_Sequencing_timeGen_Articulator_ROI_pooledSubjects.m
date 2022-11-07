global BOX_DIR
global DUKEDIR
saveFolder = '\TempDecode\Articulator\sensorimotor_delay_sig_v2\';
BOX_DIR = 'C:\Users\sd355\Box'
DUKEDIR = 'C:\Users\sd355\Box\CoganLab\D_Data'
fDown = 200; %Downsampled Sampling Frequency
timeExtract = [-2 2];

Task=[];

Task.Name='Phoneme_Sequencing';
TASK_DIR=([BOX_DIR '/CoganLab/D_Data/' Task.Name]);
DUKEDIR=TASK_DIR;
saveFolder = '\TempDecode\PooledSubjects\sensorimotor\';
%% Loading data
Subject = popTaskSubjectData(Task);
subjectIds2remove = [1 31 37:length(Subject)];
% removing D18 because of large negative response trials
Subject(subjectIds2remove) = [];

%% Loading Normalized High-Gamma for an ROI and corresponding trialInfo
fieldEpoch = 'Auditory';
selectRoi = 'postcentral';
fieldTime = [0 2];
ieegHGStruct = extractHGDataWithROI(Subject,baseName = 'Start',...
    Epoch = fieldEpoch, roi = selectRoi,Time=fieldTime);
trialInfoStruct = extractTrialInfo(Subject, remNegResponseTimeTrials=true);
% Remove empty subjects
emptyIds = [];
for iSubject = 1:length(Subject)
    if(isempty(ieegHGStruct(iSubject).ieegHGNorm))
        emptyIds = [emptyIds iSubject];
    end
end

ieegHGStruct(emptyIds) = [];
trialInfoStruct(emptyIds) = [];
% Pooling across channels based on minimum trial matching
[ieegStructPooled,phonemeTrialPooled,channelNamePooled] = poolChannelWithMinTrial(ieegHGStruct,trialInfoStruct);
delayAnatChan = find(ismember(channelNamePooled,prodDelayElecs));

%% Decoding analysis on pooled data

numFold = 10; % K-Fold cross-validation
varExplained = 80; % PCA variance

numIter = 1; % Number of decoder iterations for stable results
timeWin = 0.2; % Window length for temporal generalization (s)
timeRes = 0.02; % Window hop for temporal generalization (s)
cvcIds = find(phonemeTrialPooled.syllableUnit(:,1)'==2);
vcvIds = find(phonemeTrialPooled.syllableUnit(:,1)'==1);
% Initialize Decoder object
phonDecode = decoderClass(numFold,varExplained,numIter);
disp('1D temporal generalization');
% 1D temporalGeneralization            
decodeTimeStruct1D_c1 = phonDecode.tempGenClassify1D(ieegStructPooled,...
    phonemeTrialPooled.phonemeUnit(:,1)',timeRes,timeWin,delayAnatChan,...
    cvcIds);
% 2D temporalGeneralization            
disp('2D temporal generalization');
decodeTimeStruct2D_c1 = phonDecode.tempGenClassify2D(ieegStructPooled,...
    phonemeTrialPooled.phonemeUnit(:,1)',timeRes,timeWin,delayAnatChan,...
    cvcIds);

disp('1D temporal generalization');
% 1D temporalGeneralization            
decodeTimeStruct1D_c2 = phonDecode.tempGenClassify1D(ieegStructPooled,...
    phonemeTrialPooled.phonemeUnit(:,2)',timeRes,timeWin,delayAnatChan,...
    vcvIds);
% 2D temporalGeneralization            
disp('2D temporal generalization');
decodeTimeStruct2D_c2 = phonDecode.tempGenClassify2D(ieegStructPooled,...
    phonemeTrialPooled.phonemeUnit(:,2)',timeRes,timeWin,delayAnatChan,...
    vcvIds);

disp('1D temporal generalization');
% 1D temporalGeneralization            
decodeTimeStruct1D_c3 = phonDecode.tempGenClassify1D(ieegStructPooled,...
    phonemeTrialPooled.phonemeUnit(:,3)',timeRes,timeWin,delayAnatChan,...
    cvcIds);
% 2D temporalGeneralization            
disp('2D temporal generalization');
decodeTimeStruct2D_c3 = phonDecode.tempGenClassify2D(ieegStructPooled,...
    phonemeTrialPooled.phonemeUnit(:,3)',timeRes,timeWin,delayAnatChan,...
    cvcIds);


decodeTimeStruct1D_v1 = phonDecode.tempGenClassify1D(ieegStructPooled,...
    phonemeTrialPooled.phonemeUnit(:,1)',timeRes,timeWin,delayAnatChan,...
    vcvIds);
% 2D temporalGeneralization            
disp('2D temporal generalization');
decodeTimeStruct2D_v1 = phonDecode.tempGenClassify2D(ieegStructPooled,...
    phonemeTrialPooled.phonemeUnit(:,1)',timeRes,timeWin,delayAnatChan,...
    vcvIds);

disp('1D temporal generalization');
% 1D temporalGeneralization            
decodeTimeStruct1D_v2 = phonDecode.tempGenClassify1D(ieegStructPooled,...
    phonemeTrialPooled.phonemeUnit(:,2)',timeRes,timeWin,delayAnatChan,...
    cvcIds);
% 2D temporalGeneralization            
disp('2D temporal generalization');
decodeTimeStruct2D_v2 = phonDecode.tempGenClassify2D(ieegStructPooled,...
    phonemeTrialPooled.phonemeUnit(:,2)',timeRes,timeWin,delayAnatChan,...
    cvcIds);

disp('1D temporal generalization');
% 1D temporalGeneralization            
decodeTimeStruct1D_v3 = phonDecode.tempGenClassify1D(ieegStructPooled,...
    phonemeTrialPooled.phonemeUnit(:,3)',timeRes,timeWin,delayAnatChan,...
    vcvIds);
% 2D temporalGeneralization            
disp('2D temporal generalization');
decodeTimeStruct2D_v3 = phonDecode.tempGenClassify2D(ieegStructPooled,...
    phonemeTrialPooled.phonemeUnit(:,3)',timeRes,timeWin,delayAnatChan,...
    vcvIds);


%save('pooledSubject_smc_time_gen_articulator_responseOnset.mat','decodeTimeStruct1D','decodeTimeStruct2D');
channelsUsed = channelNamePooled(delayAnatChan);

disp('Saving..');
if ~exist([DUKEDIR saveFolder]) %#ok<EXIST> 
    mkdir([DUKEDIR saveFolder])
end
save(fullfile([DUKEDIR saveFolder 'pooledSubject_' Task.Name '_' ...
    selectRoi '_' fieldEpoch '_Start_cons_vowel_decoded.mat']),...
    'decodeTimeStruct1D_c1','decodeTimeStruct1D_c2','decodeTimeStruct1D_c3',...
    'decodeTimeStruct2D_c1','decodeTimeStruct2D_c2','decodeTimeStruct2D_c3',...
    'decodeTimeStruct1D_v1','decodeTimeStruct1D_v2','decodeTimeStruct1D_v3',...
    'decodeTimeStruct2D_v1','decodeTimeStruct2D_v2','decodeTimeStruct2D_v3',...
    'channelsUsed');
