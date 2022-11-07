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
saveFolder = '\TempDecode\PooledSubjects\prodElecs\prefrontal\phonotactic\';
%% Loading data
Subject = popTaskSubjectData(Task);
subjectIds2remove = [1 27 31 37:length(Subject)];
% removing D18 because of large negative response trials
Subject(subjectIds2remove) = [];

%% Loading Normalized High-Gamma for an ROI and corresponding trialInfo
fieldEpoch = 'Go';
selectRoi = 'frontal';
fieldTime = [-1.5 1];
respTimeThresh = 0;
ieegHGStruct = extractHGDataWithROI(Subject,baseName = 'Start',...
    Epoch = fieldEpoch, roi = selectRoi,Time=fieldTime, respTimeThresh=respTimeThresh,...
    subsetElec=prodDelayElecs,remWMchannels=false);
trialInfoStruct = extractTrialInfo(Subject, remFastResponseTimeTrials=respTimeThresh);
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
[ieegStructPooled,phonemeTrialPooled,channelNamePooled] = poolChannelWithMaxTrial(ieegHGStruct,trialInfoStruct);
% delayAnatSigChan = find(ismember(channelNamePooled,prodDelayElecs));
% channelNamePooled(delayAnatSigChan)
%% Decoding analysis on pooled data - Blick

numFold = 10; % K-Fold cross-validation
varExplained = 80; % PCA variance
decodeTimeStruct1D_cvc = [];
decodeTimeStruct1D_vcv = [];
decodeTimeStruct2D_cvc = [];
decodeTimeStruct2D_vcv = [];
numIter = 1; % Number of decoder iterations for stable results
timeWin = 0.2; % Window length for temporal generalization (s)
timeRes = 0.02; % Window hop for temporal generalization (s)
cvcIds = find(phonemeTrialPooled.syllableUnit(:,1)'==2);
vcvIds = find(phonemeTrialPooled.syllableUnit(:,1)'==1);
phonDecode = decoderClass(numFold,varExplained,numIter);
disp('CVC - 1D temporal generalization');
% cvc - 1D temporalGeneralization            
decodeTimeStruct1D_cvc = phonDecode.tempGenRegress1D(ieegStructPooled,...
    phonemeTrialPooled.phonotactic(:,10)',  selectTrials= cvcIds);
% cvc - 2D temporalGeneralization            
disp('CVC - 2D temporal generalization');
decodeTimeStruct2D_cvc = phonDecode.tempGenRegress2D(ieegStructPooled,...
    phonemeTrialPooled.phonotactic(:,10)', selectTrials = cvcIds);

disp('VCV - 1D temporal generalization');
% 1D temporalGeneralization            
decodeTimeStruct1D_vcv = phonDecode.tempGenRegress1D(ieegStructPooled,...
    phonemeTrialPooled.phonotactic(:,10)',selectTrials = vcvIds);
% 2D temporalGeneralization            
disp('VCV - 2D temporal generalization');
decodeTimeStruct2D_vcv = phonDecode.tempGenRegress2D(ieegStructPooled,...
    phonemeTrialPooled.phonotactic(:,10)',selectTrials = vcvIds);

disp('VCV - 1D temporal generalization');
% 1D temporalGeneralization            
decodeTimeStruct1D = phonDecode.tempGenRegress1D(ieegStructPooled,...
    shuffle(phonemeTrialPooled.phonotactic(:,10)'));
% 2D temporalGeneralization            
disp('VCV - 2D temporal generalization');
decodeTimeStruct2D = phonDecode.tempGenRegress2D(ieegStructPooled,...
    phonemeTrialPooled.phonotactic(:,10)');

%save('pooledSubject_smc_time_gen_articulator_responseOnset.mat','decodeTimeStruct1D','decodeTimeStruct2D');


disp('Saving..');
if ~exist([DUKEDIR saveFolder]) %#ok<EXIST> 
    mkdir([DUKEDIR saveFolder])
end
save(fullfile([DUKEDIR saveFolder 'pooledSubject_' Task.Name '_' ...
    selectRoi '_' fieldEpoch '_Start_blick_decoded_0ms_resptimethresh.mat']),...
    'decodeTimeStruct1D_cvc','decodeTimeStruct2D_cvc',...
    'decodeTimeStruct1D_vcv','decodeTimeStruct2D_vcv',...
    'decodeTimeStruct2D', 'decodeTimeStruct1D', 'channelNamePooled');
%% Decoding analysis on pooled data - phonotactics

numFold = 10; % K-Fold cross-validation
varExplained = 80; % PCA variance
decodeTimeStruct1D_cvc = [];
decodeTimeStruct1D_vcv = [];
decodeTimeStruct2D_cvc = [];
decodeTimeStruct2D_vcv = [];
numIter = 1; % Number of decoder iterations for stable results
timeWin = 0.2; % Window length for temporal generalization (s)
timeRes = 0.02; % Window hop for temporal generalization (s)
cvcIds = find(phonemeTrialPooled.syllableUnit(:,1)'==2);
vcvIds = find(phonemeTrialPooled.syllableUnit(:,1)'==1);
for ipos = 1:9
% Initialize Decoder object
phonDecode = decoderClass(numFold,varExplained,numIter);
disp('CVC - 1D temporal generalization');
% cvc - 1D temporalGeneralization            
decodeTimeStruct1D_cvc{ipos} = phonDecode.tempGenRegress1D(ieegStructPooled,...
    -log(1+phonemeTrialPooled.phonotactic(:,ipos)'),  selectTrials= cvcIds);
% cvc - 2D temporalGeneralization            
disp('CVC - 2D temporal generalization');
decodeTimeStruct2D_cvc{ipos} = phonDecode.tempGenRegress2D(ieegStructPooled,...
    -log(1+phonemeTrialPooled.phonotactic(:,ipos)'), selectTrials = cvcIds);

disp('VCV - 1D temporal generalization');
% 1D temporalGeneralization            
decodeTimeStruct1D_vcv{ipos} = phonDecode.tempGenRegress1D(ieegStructPooled,...
    -log(1+phonemeTrialPooled.phonotactic(:,ipos)'),selectTrials = vcvIds);
% 2D temporalGeneralization            
disp('VCV - 2D temporal generalization');
decodeTimeStruct2D_vcv{ipos} = phonDecode.tempGenRegress2D(ieegStructPooled,...
    -log(1+phonemeTrialPooled.phonotactic(:,ipos)'),selectTrials = vcvIds);
end
%save('pooledSubject_smc_time_gen_articulator_responseOnset.mat','decodeTimeStruct1D','decodeTimeStruct2D');


disp('Saving..');
if ~exist([DUKEDIR saveFolder]) %#ok<EXIST> 
    mkdir([DUKEDIR saveFolder])
end
save(fullfile([DUKEDIR saveFolder 'pooledSubject_' Task.Name '_' ...
    selectRoi '_' fieldEpoch '_Start_phonotactic_decoded_0ms_resptimethresh.mat']),...
    'decodeTimeStruct1D_cvc','decodeTimeStruct2D_cvc',...
    'decodeTimeStruct1D_vcv','decodeTimeStruct2D_vcv',...
    'channelNamePooled');
