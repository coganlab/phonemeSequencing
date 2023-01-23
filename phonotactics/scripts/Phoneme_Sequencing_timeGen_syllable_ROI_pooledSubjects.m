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
subjectIds2remove = [1 27 31 37:length(Subject)];
% removing D18 because of large negative response trials
Subject(subjectIds2remove) = [];

%% Loading Normalized High-Gamma for an ROI and corresponding trialInfo
fieldEpoch = 'Auditory';
 selectRoi = {'temporal','banks'};
%selectRoi = {'supramarginal','inferiorparietal'};
% selectRoi = {'frontal','opercula','triangular'};
%selectRoi = {'central'};
%selectRoi = '';
fieldTime = [-0.5 2];
respTimeThresh = 0;

[trialInfos,trialInfoStruct] = extractTrialInfo(Subject, remFastResponseTimeTrials=respTimeThresh);

ieegHGStruct = extractHGDataWithROI(Subject,baseName = 'Start',...
    Epoch = fieldEpoch, roi = selectRoi,Time=fieldTime,respTimeThresh=respTimeThresh,...
    subsetElec=elecNameFeedBack,remWMchannels=false);

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

save(fullfile(['pooledSubject_' Task.Name '_temporal_feedback' ...
    '_' fieldEpoch '_delayelecs_data.mat']),...
    'ieegStructPooled','phonemeTrialPooled',...
    'channelNamePooled');

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
decodeTimeStruct1D = phonDecode.tempGenClassify1D(ieegStructPooled,...
    shuffle(phonemeTrialPooled.syllableUnit(:,1))',timeRes=timeRes,timeWin=timeWin);
% 2D temporalGeneralization            
disp('2D temporal generalization');
decodeTimeStruct2D = phonDecode.tempGenClassify2D(ieegStructPooled,...
    phonemeTrialPooled.syllableUnit(:,1)',timeRes=timeRes,timeWin=timeWin);



%save('pooledSubject_smc_time_gen_articulator_responseOnset.mat','decodeTimeStruct1D','decodeTimeStruct2D');


disp('Saving..');
if ~exist([DUKEDIR saveFolder]) %#ok<EXIST> 
    mkdir([DUKEDIR saveFolder])
end
save(fullfile([DUKEDIR saveFolder 'pooledSubject_' Task.Name '_' ...
    selectRoi '_' fieldEpoch '_Start_syllable_resptime_0ms_decoded.mat']),...
    'decodeTimeStruct1D','decodeTimeStruct2D',...
    'channelNamePooled');
