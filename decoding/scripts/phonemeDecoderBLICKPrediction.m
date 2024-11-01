addpath(genpath('C:\Users\sd355\Box\ECoG_Recon\matlab_code\'));
global DUKEDIR
DUKEDIR = 'C:\Users\sd355\Box\CoganLab\Data\Micro\Processed Data';
dLabels = dir(DUKEDIR);
dLabels = dLabels(3:end);

tw = [-2 2]; % time window


baseTW = [-0.5 0]; % preonset time window
activeTW = [-1.5 1.5]; % postonset time window
gammaF = [70 150]; % frequency in Hz
fsDown = 200;
% load('channelMap.mat');
% chanMap = (chanMap');
% selectedChannels = sort(chanMap(~isnan(chanMap)))';
%% data Loading
subjectId = 'S14';
iPhon = 1
Experiment = loadExperiment(subjectId);
chanMap = Experiment.recording.channel_map;
chanMap = (chanMap');
selectedChannels = sort(chanMap(~isnan(chanMap)))';
fsD = Experiment.processing.ieeg.sample_rate;
Trials = dbTrials(subjectId,Experiment.recording.recording_day,'Speech_OvertMimeMove');

trialFiles = strcat('\',Experiment.recording.recording_day,'\mat\trialInfo.mat');
load([DUKEDIR '/' subjectId '/' trialFiles ]);
[ieegResponse,~,trigOnset]=trialIEEGUpdate(Trials,selectedChannels,['ResponseOnset'],'ieeg',tw.*1000);
[ieegAuditory,~,trigAuditory]=trialIEEGUpdate(Trials,selectedChannels,'Auditory','ieeg',tw.*1000);

ieegAuditory = permute(ieegAuditory,[2,1,3]);
ieegResponse = permute(ieegResponse,[2,1,3]);
 
respId = find(~isnan(trigOnset));
ieegResponse = ieegResponse(:,respId,:);
ieegAuditory = ieegAuditory(:,respId,:); 
Trials = Trials(respId);
trialInfo = trialInfo(respId);  
trigOns = trigOnset(respId)./fsD;
trigAud = trigAuditory(respId)./fsD;
respTime  = trigOns-trigAud;
[respTimeSort,sortId] = sort(respTime);
%% Ieeg Class definition

ieegAud = ieegStructMicro(ieegAuditory, fsD, tw, [1 fsD/2], 'Auditory', chanMap);
ieegResponse = ieegStructMicro(ieegResponse, fsD, tw, [1 fsD/2], 'Response', chanMap);

%% Phoneme trial parser
phonemeTrial = phonemeSequenceTrialParser(trialInfo);
%% Common average referencing
load([subjectId '_impedance.mat'])
higImpId = find(log10(impedance)>6);
ieegAudCar = ieegAud.extractCar(higImpId);
clear ieegAud
ieegResponseCar = ieegResponse.extractCar(higImpId);
clear ieegResponse
%% Extract High Gamma
normType =2 ;
ieegHGBase = ieegAudCar.extractHiGamma(fsDown,baseTW);
normFactor = extractHGnormFactor(ieegHGBase);
ieegHGAudNorm = normHiGamma(ieegHGBase,normFactor,2);
ieegHGRespNorm = ieegResponseCar.extractHiGamma(fsDown,activeTW,normFactor, normType);

%% Extracting all trials
load([subjectId '_sigChannel.mat'])
phonemeUnits = phonemeTrial.syllableUnit(:,iPhon)';
cvcTrials = phonemeUnits==1;
vcvTrials = ~cvcTrials;

%% Setting up decoder
numFold = 20;
varExplained = 80;
phonRegress = 'Pfwd1';
timeRes = 0.01;
timeWin = 0.2;
phonDecode = phonemeDecoderClass(numFold,varExplained);
disp('CVC 1D temporal generalization');
% 1D temporalGeneralization            
decodeTimeStruct1D = phonDecode.tempGenRegress1D(ieegHGRespNorm,phonemeTrial,phonRegress,timeRes,timeWin,sigChannel,cvcTrials);
% 2D temporalGeneralization            
disp('CVC 2D temporal generalization');
decodeTimeStruct2D = phonDecode.tempGenRegress2D(ieegHGRespNorm,phonemeTrial,phonRegress,timeRes,timeWin,sigChannel,cvcTrials);
            
save([subjectId '_' phonRegress '_predict_cvc_temp_gen.mat'],'decodeTimeStruct1D','decodeTimeStruct2D');

disp('VCV 1D temporal generalization');
% 1D temporalGeneralization            
decodeTimeStruct1D = phonDecode.tempGenRegress1D(ieegHGRespNorm,phonemeTrial,phonRegress,timeRes,timeWin,sigChannel,vcvTrials);
% 2D temporalGeneralization            
disp('VCV 2D temporal generalization');
decodeTimeStruct2D = phonDecode.tempGenRegress2D(ieegHGRespNorm,phonemeTrial,phonRegress,timeRes,timeWin,sigChannel,vcvTrials);
            
save([subjectId '_' phonRegress '_predict_vcv_temp_gen.mat'],'decodeTimeStruct1D','decodeTimeStruct2D');


