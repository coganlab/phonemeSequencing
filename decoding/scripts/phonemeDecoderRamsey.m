addpath(genpath('C:\Users\sd355\Box\ECoG_Recon\matlab_code\'));
global DUKEDIR
DUKEDIR = 'C:\Users\sd355\Box\CoganLab\Data\Micro\Processed Data';
dLabels = dir(DUKEDIR);
dLabels = dLabels(3:end);

tw = [-1.5 1.5]; % time window


baseTW = [-0.5 0]; % preonset time window
activeTW = [-1 1]; % postonset time window
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
[ieegResponse,~,trigOnset]=trialIEEGUpdate(Trials,selectedChannels,['phon' num2str(iPhon) 'Onset'],'ieeg',tw.*1000);
[ieegAuditory,~,trigAuditory]=trialIEEGUpdate(Trials,selectedChannels,'Auditory','ieeg',tw.*1000);

ieegAuditory = permute(ieegAuditory,[2,1,3]);
ieegResponse = permute(ieegResponse,[2,1,3]);
 
respId = find(~isnan(trigOnset));
ieegResponse = ieegResponse(:,respId,:);
ieegAuditory = ieegAuditory(:,respId,:); 
Trials = Trials(respId);
  
trigOns = trigOnset(respId)./fsD;
trigAud = trigAuditory(respId)./fsD;
respTime  = trigOns-trigAud;
[respTimeSort,sortId] = sort(respTime);
%% Ieeg Class definition

ieegBase = ieegStructMicro(ieegAuditory, fsD, tw, [1 fsD/2], 'Auditory', chanMap);
ieegResponse = ieegStructMicro(ieegResponse, fsD, tw, [1 fsD/2], 'Response', chanMap);

%% Phoneme trial parser
phonemeTrial = phonemeSequenceTrialParser(Trials);
%% Common average referencing
load([subjectId '_impedance.mat'])
higImpId = find(log10(impedance)>6);
ieegBaseCar = ieegBase.extractCar(higImpId);
clear ieegBase
ieegResponseCar = ieegResponse.extractCar(higImpId);
clear ieegResponse
%% Extract High Gamma
normType =2 ;
ieegHGBase = ieegBaseCar.extractHiGamma(fsDown,baseTW);
normFactor = extractHGnormFactor(ieegHGBase);
ieegHGBaseNorm = normHiGamma(ieegHGBase,normFactor,2);
ieegHGRespNorm = ieegResponseCar.extractHiGamma(fsDown,activeTW,normFactor, normType);

%% Extracting Ramsey trials p, k, u, a, rest
phonemeUnits = phonemeTrial.phonemeUnit(:,1)';
rTrials = phonemeUnits == 1|phonemeUnits == 4 |phonemeUnits == 6|phonemeUnits == 9;
cTrials = phonemeTrial.syllableUnit(:,1)'==1;
vTrials = ~cTrials;
rPhonemeUnits = phonemeUnits(rTrials);


%% Combining phoneme & rest trials

% ieegHGComb = ieegHGRespNorm;
% ieegHGComb.data = ieegHGComb.data(:,rTrials,:);
% % selectTrials1 = randsample(1:size(ieegHGBaseNorm.data,2), sum(rTrials)/4);
% % selectTrials2 = randsample(1:size(ieegHGBaseNorm.data,2), sum(rTrials)/4);
% % %ieegHGComb.data = cat(2,ieegHGComb.data(:,rTrials,:), ieegHGBaseNorm.data(:,selectTrials1,:));
% %  ieegHGComb.data = cat(2,ieegHGComb.data(:,rTrials,:), cat(3,ieegHGBaseNorm.data(:,selectTrials1,:),ieegHGBaseNorm.data(:,selectTrials2,:)));
% % rPhonemeUnits = [rPhonemeUnits 10.*ones(1,length(selectTrials1))];

%% decoder object
load([subjectId '_sigChannel.mat'])
dObj = decoderClass(20,80,1);
% ifgChan = [16 12 126 122 105 104 114 108 5 62 37 11 21 48 63 8];
%  ifgChanFp =    [124 101 106 111 100 102 20 22 7 24 9 30 39 46 50 6 61 31 55 70 69 77 68 81 82 72];
% smcChan = setdiff(sigChannel,[ifgChan ifgChanFp]);
decodeResultStruct = dObj.baseClassify(ieegHGRespNorm,phonemeUnits, [-0.5 0.5], sigChannel)


% channels2plot = zeros(1,256);
% channels2plot(ifgChan) = 1;
% channels2plot(smcChan) = 2;
% figure; chanView(channels2plot,chanMap2,"cval",[0 2]); formatTicks(gca); axis equal
% set(gca,'YTick',[])
% set(gca,'XTick',[])
% title('');
% decodeTimeStruct1 = dObj.tempGenClassify1D(ieegHGRespNorm,phonemeTrial.phonemeUnit(:,1)',0.01,0.2,sigChannel);
% decodeTimeStruct2 = dObj.tempGenClassify1D(ieegHGRespNorm,phonemeTrial.phonemeUnit(:,2)',0.01,0.2,sigChannel);
% decodeTimeStruct3 = dObj.tempGenClassify1D(ieegHGRespNorm,phonemeTrial.phonemeUnit(:,3)',0.01,0.2,sigChannel);
% 
% 
% figure;
% hold on;
% [accMax,accMaxId1] = max(decodeTimeStruct1.accTime');
% plot(decodeTimeStruct1.timeRange+0.1, decodeTimeStruct1.accTime.*100,'LineWidth',2);
% 
% 
% [accMax,accMaxId2] = max(decodeTimeStruct2.accTime');
% plot(decodeTimeStruct2.timeRange+0.1, decodeTimeStruct2.accTime.*100,'LineWidth',2);
% hold on;
% 
% [accMax,accMaxId3] = max(decodeTimeStruct3.accTime');
% plot(decodeTimeStruct3.timeRange+0.1, decodeTimeStruct3.accTime.*100,'LineWidth',2);
% legend('P1','P2','P3')
% hold on;
% xline(decodeTimeStruct1.timeRange(accMaxId1)+0.1,':','','LineWidth',2);
% xline(decodeTimeStruct2.timeRange(accMaxId2)+0.1,':','','LineWidth',2);
% xline(decodeTimeStruct3.timeRange(accMaxId3)+0.1,':','','LineWidth',2);
% 
% xlabel('Time from utterance onset (s)');
% ylabel('Decoding Accuracy (%)');
% yline(11.11, '--','chance','LineWidth',2);
% set(gca,'FontSize',15);
% title('Subject 1');
% 
% ylim([0 70]);
% axis square

%% Space grid sampling
chanWindow = [2 4; 3 6; 4 8; 5 10; 6 12; 7 14; 8 16];
accSpaceSamp = [];
for iWindow = 1:size(chanWindow,1)
    iWindow
    elecPoints = matrixSubSample(chanMap',chanWindow(iWindow,:),1);
    accSamp = [];
    phonErrorSamp = [];
    for iSamp = 1:size(elecPoints,1)
        elecPointSamp = elecPoints(iSamp,:);
        elecPointSamp = elecPointSamp(~isnan(elecPointSamp));
        dObj = decoderClass(20,[80]);
        decodeResultStruct = dObj.baseClassify(ieegHGRespNorm,phonemeUnits, [-0.5 0.5],elecPointSamp);       
        accSamp(iSamp) =  decodeResultStruct.accPhoneme;
    end
    accSpaceSamp{iWindow} = accSamp;
end
elecSpaceSampStr = [];
accSpaceSampAll = [];
accSpaceMean = [];
for iWindow = 1:size(chanWindow,1)    
    elecSpaceSampStrTemp = [num2str(chanWindow(iWindow,1)) 'x' num2str(chanWindow(iWindow,2))];
    elecSpaceSampStr = [elecSpaceSampStr; repmat({elecSpaceSampStrTemp},length(accSpaceSamp{iWindow}'),1)];
    accSpaceSampAll = [accSpaceSampAll accSpaceSamp{iWindow}];
    accSpaceMean(iWindow) = median(accSpaceSamp{iWindow});    
end 
save([subjectId '_PhonemeDecodeSpaceGrid_' num2str(iPhon)  '_Phoneme_v1.mat'],'accSpaceSampAll', 'elecSpaceSampStr', 'accSpaceMean');   
%save('PhonemeDecodeSpaceGridFirstPhoneme.mat','accSpaceSampAll','phonErrorSpaceSampAll','elecSpaceSampStr');

%% Poisson disc sampling - All phoneme decoder
nonSigChannel = [];

nSamp = 1.5;
firstElec = [0.025 0.05];
totalDist = 120;
accDist = [];
phonErrorDist = [];

sulcusChannels2Remove = [71 67 115 75 121 51 117 29 63 27 58 46 44 35 62 23 112 109];
chanMapSig = chanMap;
chanMapSig(ismember(chanMap(:),sulcusChannels2Remove)) = nan;
numSamp = [];
for iDist = 1:totalDist
    elecSampDensity = firstElec+(iDist-1)*firstElec(1);
    iDist
    parfor iSamp = 1:round(nSamp*(totalDist/iDist))
        elecPt = [];
        while(isempty(elecPt))
            elecSampClose = elecSampDensity;
             nElec = round((elecSampClose(1) + (elecSampClose(2)-elecSampClose(1)) .* rand(1,1))*length(selectedChannels));
             %elecPtIds = round(poissonDisc([8,16],poisSpace,nElec));
             elecPtIds = ceil(poissonDisc2(size(chanMap),nElec));
                
    
            elecPt = [];
            for iElec = 1:size(elecPtIds,1)
                elecPt(iElec) = chanMapSig(elecPtIds(iElec,1),elecPtIds(iElec,2));
            end
            elecPt = elecPt(~isnan(elecPt)); 
        end
        elecPtcm = ismember(selectedChannels,elecPt);
        dObj = decoderClass(10,[80],1);
        decodeResultStruct = dObj.baseClassify(ieegHGRespNorm,phonemeUnits, [-0.5 0.5], find(elecPtcm));       
        accDist = [accDist decodeResultStruct.accPhoneme];
%         dObj = decoderClass(20,[80]);
%         decodeResultStruct = dObj.baseClassify(ieegHGRespNorm,phonemeUnits, [-0.5 0.5], find(elecPtcm),cTrials);        
%         accDistCons = [accDistCons decodeResultStruct.accPhoneme];
%         dObj = decoderClass(20,[80]);
%         decodeResultStruct = dObj.baseClassify(ieegHGRespNorm,phonemeUnits, [-0.5 0.5], find(elecPtcm),vTrials);        
%         accDistVowel = [accDistVowel decodeResultStruct.accPhoneme];
        numSamp = [numSamp sum(elecPtcm)];
    end
    
    save([subjectId '_PhonemeDecodeHighRes_' num2str(iPhon)  '_Phoneme_no_sulcus_v2.mat'],'accDist', 'numSamp');   
end

%% Poisson disc sampling - N-way
nonSigChannel = [];

nSamp = 1.5;
firstElec = [0.025 0.05];
totalDist = 120;



for iWay = 2:9
    accDist = [];
    numSamp = [];
    labels2choose = nchoosek(1:9,iWay);
    for iDist = 1:totalDist
        elecSampDensity = firstElec+(iDist-1)*firstElec(1);
        iDist
        for iSamp = 1:round(nSamp*(totalDist/iDist))
            
           
            elecSampClose = elecSampDensity;
             nElec = round((elecSampClose(1) + (elecSampClose(2)-elecSampClose(1)) .* rand(1,1))*length(selectedChannels));
             %elecPtIds = round(poissonDisc([8,16],poisSpace,nElec));
             elecPtIds = ceil(poissonDisc2(size(chanMap),nElec));
                
    
            elecPt = [];
            for iElec = 1:size(elecPtIds,1)
                elecPt(iElec) = chanMap(elecPtIds(iElec,1),elecPtIds(iElec,2));
            end
            elecPt = elecPt(~isnan(elecPt)); 
            elecPtcm = ismember(selectedChannels,elecPt);
            accTrial = zeros(1,size(labels2choose,1));
            for iTrial = 1:size(labels2choose,1)
                rTrials= ismember(phonemeUnits,labels2choose(iTrial,:));
                dObj = decoderClass(10,80);
                decodeResultStruct = dObj.baseClassify(ieegHGRespNorm,phonemeUnits, [-0.5 0.5], find(elecPtcm),rTrials);
                accTrial(iTrial) = decodeResultStruct.accPhoneme;
            end
            mean(accTrial)
            accDist = [accDist mean(accTrial)];
            
            numSamp = [numSamp sum(elecPtcm)];
        end
        sum(elecPtcm)
        save([subjectId '_PhonemeDecodeHighResFirstPhoneme_' num2str(iWay) '_way_v2.mat'],'accDist', 'numSamp');   
    end
end


