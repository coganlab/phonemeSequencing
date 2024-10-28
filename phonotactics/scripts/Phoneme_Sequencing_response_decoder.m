
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
% subjectIds2remove = [1 2 18 24 31 37:length(Subject)];
% % removing D18 because of large negative response trials
% Subject(subjectIds2remove) = [];

%% Loading Normalized High-Gamma for an ROI and corresponding trialInfo

 fieldEpoch = 'ResponseStart';
roiIdsAll = false(1,length(elecNameProductionIntersect));
%roiList = {'superiortemporal','ipc','ifg','rmfg','cmfg','smc','insula','opercular','triangular','astg','pstg','heschl','sts'};

%roiList = {'ifg','rmfg','cmfg'};
% roiList = {'smc','ipc'};
%  roiList = {'astg','pstg','heschl','sts'};
% hem = 'rh';
% selectRoi = [];
% for iRoi = 1:length(roiList)
%     selectRoi = [selectRoi extractBNlabels(roiList{iRoi},hem)];
% %     roiIds = extractRoiIds(Subject, elecNameProductionClean,selectRoi);
% %     roiIdsAll = roiIdsAll | roiIds;
% end



 selectRoi ={''};
% roiIds = extractRoiIds(Subject, elecNameAllTimePermClean,selectRoi);
% % numTotalRoi = sum(roiIds);
%  roiIds = extractRoiIds(Subject, elecNameProductionClean,selectRoi);
% numProdActiveRoi = sum(roiIds);


respTimeThresh = 0.05;
respDurThresh = 1.5;
timeEpoch = [-1.5  1.5];

% channelInfoAll = extractChannelLocation(Subject,elecNameProductionClean);
% xyzChannel = reshape(extractfield(channelInfoAll,'xyz'),3,length(elecNameProductionClean));
% preAcIds = xyzChannel(2,:)<0;
% roiIds = roiIds & preAcIds;
% sum(roiIds)
%trialInfoStruct = extractTrialInfo(Subject, remFastResponseTimeTrials=respTimeThresh);

% elecs2remove = elecNameFeedBack_000;
% elecs2remove =  elecs2remove(~cellfun('isempty',elecs2remove));    
    
% elecNameProductionInfo = extractChannelLocation(Subject,elecNameProductionClean);
% elecNameProductionWM = elecNameProductionClean(contains({elecNameProductionInfo.Location},["White"]));
%roiIds = extractRoiIds(Subject, elecNameProductionClean,["White","hypointen","known"]);
subSelectElecs = elecNameProductionIntersect;

ieegHGStruct = extractHGDataWithROI(Subject,baseName = 'Start',...
    Epoch = fieldEpoch, roi = selectRoi, Time=timeEpoch,respTimeThresh=respTimeThresh,...
    subsetElec=subSelectElecs,remWMchannels=true,normType=1,fDown = 200,respDurThresh=respDurThresh);

% Remove empty subjects
emptyIds = [];

trialInfoStruct = [];
for iSubject = 1:length(Subject)
    
    if(isempty(ieegHGStruct(iSubject).ieegHGNorm))
        trialInfoStruct(iSubject).phonemeTrial = [];
        emptyIds = [emptyIds iSubject];
    else
      trialInfoStruct(iSubject).phonemeTrial = phonemeSequenceTrialParser(ieegHGStruct(iSubject).trialInfo);   
    end

end

ieegHGStruct(emptyIds) = [];
trialInfoStruct(emptyIds) = [];

maxTrial = [];
for iSubj = 1:length(ieegHGStruct)
    chanSubj(iSubj) = length(ieegHGStruct(iSubj).channelName);
    trialSubj(iSubj) = length(ieegHGStruct(iSubj).trialInfo);
    if(trialSubj(iSubj)>=156)
        maxTrial = [maxTrial iSubj];
    end
end
% 
% ieegHGStruct = ieegHGStruct(maxTrial);
% trialInfoStruct = trialInfoStruct(maxTrial);
% 
% 
 save(fullfile(['pooledSubject_' Task.Name '_Production_car_z_score' ...
    '_prodelecs_data_lex_intersect.mat']),...
    'ieegHGStruct','trialInfoStruct','maxTrial',...
    'timeEpoch','-v7.3');
  
%% Loading Normalized High-Gamma for an ROI and corresponding trialInfo

% fieldEpoch = 'Auditory';

% selectRoi = {'A21c','A21r', 'A37dl', 'aSTS'}; % middletemporal
%  selectRoi = {'A38m','A41/42','TE1.0','TE1.2','A22c','A38l','A22r','rpSTS','cpSTS'}; % superiortemporal
%  selectRoi = {'A20iv','A37elv','A20r','A37elv','A20r','A20il','A37vl','A20cl', 'A20cv','A20rv','A37mv','A37lv'}; % inferiortemporal combined with fusiform
%  selectRoi = {'A39c','A39rd','A40rd','A40c','A39rv','A40rv'}; % inferior-parietal
%  selectRoi = {'A44d','IFS','A45c','A45r', 'A44op', 'A44v'}; %IFG
%  selectRoi = {'A46','A10l','A9/46v'}; % rostral middle frontal
%  selectRoi = {'IFJ','A8vl','A6vl','A9/46d'}; % caudal middle frontal
%  selectRoi = {'A4hf','A6cdl','A4ul','A4t','A4tl','A6cvl','A1/2/3ll','A4ll', 'A1/2/3ulhf','A1/2/3/tonla','A2','A1/2/3tru'}; % sensorimotor, paracentral
%  selectRoi = {'G','vIa','dIa','vId/vIg','dIg','dId'}; %insula
%  selectRoi = {'A6m','A6dl'}; %sma

% selectRoi = {'supramarginal','inferiorparietal','temporal','sts','opercular','triangular','middlefrontal','central','insula','superiorfrontal','fusiform'};
%selectRoi = [extractBNlabels('middletemporal','bh'),extractBNlabels('superiortemporal','bh'),extractBNlabels('inferiortemporal','bh')];
%selectRoi = strcat(selectRoi,'_R');
hemList = {'lh','rh'}
roiList = {'middletemporal','superiortemporal','inferiortemporal','ipc','ifg','rmfg','cmfg','smc','insula','sma','opercular','triangular','astg','pstg','heschl','sts','ifs'};

for iHem = 1:length(hemList)
    for iRoi = 1:length(roiList)
        selectRoi = extractBNlabels(roiList{iRoi},hemList{iHem})
        roiIds = extractRoiIds(Subject, elecNameAllTimePermClean,selectRoi);
        numTotalRoi = sum(roiIds);
        roiIds = extractRoiIds(Subject, elecNameProductionClean,selectRoi);
        numProdActiveRoi = sum(roiIds);
        
        
        respTimeThresh = 0.05;
        respDurThresh = 1.5;
        timeEpoch = [-1.5  1.5];
        
        channelInfoAll = extractChannelLocation(Subject,elecNameProductionClean);
        xyzChannel = reshape(extractfield(channelInfoAll,'xyz'),3,length(elecNameProductionClean));
        preAcIds = xyzChannel(2,:)<0;
        roiIds = roiIds & preAcIds;
        sum(roiIds)
        %trialInfoStruct = extractTrialInfo(Subject, remFastResponseTimeTrials=respTimeThresh);
        
        % elecs2remove = elecNameFeedBack_000;
        % elecs2remove =  elecs2remove(~cellfun('isempty',elecs2remove));    
            
        % elecNameProductionInfo = extractChannelLocation(Subject,elecNameProductionClean);
        % elecNameProductionWM = elecNameProductionClean(contains({elecNameProductionInfo.Location},["White"]));
        %roiIds = extractRoiIds(Subject, elecNameProductionClean,["White","hypointen","known"]);
        subSelectElecs = elecNameProductionClean;
        
        ieegHGStruct = extractHGDataWithROI(Subject,baseName = 'Start',...
            Epoch = 'ResponseStart', roi = selectRoi, Time=timeEpoch,respTimeThresh=respTimeThresh,...
            subsetElec=subSelectElecs,remWMchannels=true,normType=1,fDown = 200,respDurThresh=respDurThresh);
        
        % Remove empty subjects
        emptyIds = [];
        
        trialInfoStruct = [];
        for iSubject = 1:length(Subject)
            
            if(isempty(ieegHGStruct(iSubject).ieegHGNorm))
                trialInfoStruct(iSubject).phonemeTrial = [];
                emptyIds = [emptyIds iSubject];
            else
              trialInfoStruct(iSubject).phonemeTrial = phonemeSequenceTrialParser(ieegHGStruct(iSubject).trialInfo);   
            end
        
        end
        
        ieegHGStruct(emptyIds) = [];
        trialInfoStruct(emptyIds) = [];
        
        maxTrial = [];
        for iSubj = 1:length(ieegHGStruct)
            chanSubj(iSubj) = length(ieegHGStruct(iSubj).channelName);
            trialSubj(iSubj) = length(ieegHGStruct(iSubj).trialInfo);
            if(trialSubj(iSubj)>=156)
                maxTrial = [maxTrial iSubj];
            end
        end
        % 
        % ieegHGStruct = ieegHGStruct(maxTrial);
        % trialInfoStruct = trialInfoStruct(maxTrial);
        % 
        % 
         save(fullfile(['pooledSubject_' Task.Name '_Production_car_' hemList{iHem} '_' roiList{iRoi} '_z_score' ...
            '_prodelecs_data_v3.mat']),...
            'ieegHGStruct','trialInfoStruct','numTotalRoi','numProdActiveRoi','maxTrial',...
            'timeEpoch','-v7.3');
    end
end
%%
% Pooling across channels based on maximum trial matching; zero pad extra
% trials

maxTrial = [];
for iSubj = 1:length(ieegHGStruct)
    chanSubj(iSubj) = length(ieegHGStruct(iSubj).channelName);
    trialSubj(iSubj) = length(ieegHGStruct(iSubj).trialInfo);
    if(trialSubj(iSubj)>167)
        maxTrial = [maxTrial iSubj];
    end
end

 [ieegStructPooled,phonemeTrialPooled,channelNamePooled] = poolChannelWithMinTrial(ieegHGStruct,trialInfoStruct);
% save(fullfile(['pooledSubject_' Task.Name '_production_other_anat_zscore' ...
%     '_prodelecs_data.mat']),...
%     'ieegStructPooled','phonemeTrialPooled',...
%     'channelNamePooled','timeEpoch','-v7.3');


%% High gamma activation movie
violColor = [ 0.231 0.035 0.502];
brightViol1 = brighten(violColor,0.9);
brightViol2 = brighten(violColor,0.5);
roiIds = extractRoiIds(Subject, elecNameProductionClean,["White"]);
plot_elec_density_v2(channelNamePooled, cLim = [0 5],...
    gaussFwhm=15, diskThresh = 100, stype='pial', ptype='inflated',...
    transparentPoint=0, colorscale = 'parula');

[ieegStructPooled,phonemeTrialPooled,channelNamePooled] = poolChannelWithMaxTrial(ieegHGStruct,trialInfoStruct);

vTrials = phonemeTrialPooled.syllableUnit(:,1)'==1;
cTrials = ~vTrials;

ieegMean = squeeze(nanmean(ieegStructPooled.data,2));
absDiffIeeg = abs(diff(ieegMean'));
%chanNoiseId = (max(absDiffIeeg(1:20,:))>0.15)|(max(absDiffIeeg(245:275,:))>0.15)|(max(absDiffIeeg(65:95,:))>0.15)|(max(absDiffIeeg(580:640,:))>0.15)
chanNoiseId = (max(absDiffIeeg(210:240,:))>0.15)|(max(absDiffIeeg(250:270,:))>0.15)
chanNoise = channelNamePooled(chanNoiseId)
meanChan = (mean(ieegMean(:,timeGamma>=-0.5&timeGamma<=0.5)'))';
selectRoi = {'supramarginal','inferiorparietal','temporal','sts','opercula','triangular','middlefrontal','central','insula','superiorfrontal','fusiform'};

roiIds = extractRoiIds(Subject, channelNamePooled,selectRoi);

ieegMeanCVC = squeeze(mean(ieegStructPooled.data(:,cTrials,:),2));
ieegMeanVCV = squeeze(mean(ieegStructPooled.data(:,vTrials,:),2));
%plot_subjs_on_average_activation(channelNamePooled, (max(ieegMean'))', 'fsaverage', cfg)

plot_activation_density_v2(channelNamePooled(roiIds), meanChan(roiIds), cLim = [0 1],...
    gaussFwhm=15, diskThresh = 100, stype='pial', dtype='inflated', colbarTitle =' HG power',...
    transparentPoint=0.25, colorscale = 'parula');

plot_activation_timeSeries_bh_v2(channelNamePooled(roiIds), ieegMean(roiIds,:), tw = [-1 1.5],...
    cLim = [0 1.5], movTitle = 'High-gamma-utterance_20_mm_fwhm_euclidean_gaussian__z_score',...
    frameRate=10, gaussFwhm=20, diskThresh = 100, stype='pial', ptype='inflated',...
    colbarTitle=' HG (z-score)', transparentPoint=0.5, colorscale = 'parula',showMesial=1);

plot_activation_timeSeries_geodesic_rh_v2(channelNamePooled(roiIds), ieegMean(roiIds,timeGamma>=-0.5), tw = [-0.5 1.5],...
    cLim = [0 1.5], movTitle = 'High-gamma-utterance_20_mm_fwhm_geodesic_gaussian_z_score',...
    frameRate=10, gaussFwhm=20, diskThresh = 100, stype='pial', ptype='inflated',...
    colbarTitle=' HG (z-score)', transparentPoint=0.5, colorscale = 'parula',showMesial=1);

%%
% Decoder model

%load('pooledSubject_Phoneme_Sequencing_production_all_zscore_prodelecs_data.mat')


vTrials = phonemeTrialPooled.syllableUnit(:,1)'==1;
syllableUnit = phonemeTrialPooled.syllableUnit(:,1)';
cTrials = ~vTrials;
% phoneme decoder
phonemeUnits = phonemeTrialPooled.phonemeUnit(:,1)';
phonemeClass = phonemeTrialPooled.phonemeClass(:,1)';

dObj = decoderClass(10,[80],1);
decodeResultStruct = dObj.baseClassify(ieegStructPooled,phonemeUnits, d_time_window = [-0.5 0.5])

% consonant decoder
dObj = decoderClass(20,[80],1);
decodeResultStructConsonant = dObj.baseClassify(ieegStructPooled,phonemeUnits, d_time_window = [-0.5 0.5],selectTrial=find(cTrials))

% vowel decoder
dObj = decoderClass(20,[80],1);
decodeResultStructVowel = dObj.baseClassify(ieegStructPooled,phonemeUnits, d_time_window = [-0.5 0.5],selectTrial=find(vTrials))

% articulator decoder
dObj = decoderClass(20,[80],1);
decodeResultStructArticulator = dObj.baseClassify(ieegStructPooled,phonemeClass, d_time_window = [-0.5 0.5])

% syllable decoder
dObj = decoderClass(20,[80],1);
decodeResultStructSyllable = dObj.baseClassify(ieegStructPooled,syllableUnit, d_time_window = [-0.5 0.5])

save(fullfile(['pooledSubject_' Task.Name '_production_all_zscore' ...
    '_prodelecs_data.mat']),...
    'ieegStructPooled','phonemeTrialPooled',...
    'channelNamePooled','timeEpoch','decodeResultStruct','decodeResultStructConsonant','decodeResultStructSyllable'...
    ,'decodeResultStructVowel','decodeResultStructArticulator','-v7.3');
%%
cmatzero = decodeResultStruct.cmat

accTemp = mean(diag(cmatzero))
labels = {'a','ae','i','u','b','p','v','g','k'};        
%  labels = {'low','high','labials','dorsals'}      
% labels = {'CVC','VCV'}
figure; 
imagesc(cmatzero.*100);
colormap parula
set(gca,'XTick',[1:length(labels)]);
set(gca,'YTick',[1:length(labels)]);
set(gca,'xticklabels',labels)
set(gca,'yticklabels',labels)
cb = colorbar;
caxis([0 max(diag(cmatzero.*100))])
ylabel(cb,'Accuracy (%)');
set(gca,'FontSize',15);
title(['Accuracy : ' num2str(accTemp*100) '%']);
formatTicks(gca)
%% Phoneme decoding only


cMatAll = [];
accIter = [];
for iTer = 1:10
    [ieegStructPooled,phonemeTrialPooled,channelNamePooled] = poolChannelWithMaxTrial(ieegHGStruct(maxTrial),trialInfoStruct(maxTrial));
%      roiIds = extractRoiIds(Subject, channelNamePooled,["supramarginal", "inferiorparietal" ]);
    vTrials = phonemeTrialPooled.syllableUnit(:,1)'==1;
    syllableUnit = phonemeTrialPooled.syllableUnit(:,1)';
    cTrials = ~vTrials;
    % phoneme decoder
    phonemeUnit = phonemeTrialPooled.phonemeUnit(:,1)';
    phonemeClass = phonemeTrialPooled.phonemeClass(:,1)';
    decoderUnit = phonemeUnit;
    dObj = decoderClass(10,[80],1);
    decodeResultStruct = dObj.baseClassify(ieegStructPooled,decoderUnit, d_time_window = [-0.5 0.5]);
    
    cMatAll(iTer,:,:) = decodeResultStruct.cmat;
    accIter(iTer) = decodeResultStruct.accPhoneme;

end
%% Model weight extraction
modelweightAll = [];
modelweightMean = [];

modelweightAllShuffle = [];
cMatAll = [];
accIter = [];
hemList = {'lh','rh'}
roiList = {'superiortemporal','ipc','ifg','rmfg','cmfg','smc','insula','opercular','triangular','astg','pstg','heschl','sts','ifs'};



for iTer = 1:4
    [ieegStructPooled,phonemeTrialPooled,channelNamePooled] = poolChannelWithMaxTrial(ieegHGStruct(maxTrial),trialInfoStruct(maxTrial));
%      roiIds = extractRoiIds(Subject, channelNamePooled,["supramarginal", "inferiorparietal" ]);
    roiIdsAll = true(1,length(channelNamePooled));
    for iHem = 1:length(hemList)
        for iRoi = 1:length(roiList)
            selectRoi = extractBNlabels(roiList{iRoi},hemList{iHem})
            roiIds = extractRoiIds(Subject, channelNamePooled,selectRoi);
            roiIdsAll = roiIdsAll | roiIds;
        end
    end
    channelNameSelect = channelNamePooled(roiIdsAll);

    vTrials = phonemeTrialPooled.syllableUnit(:,1)'==1;
    syllableUnit = phonemeTrialPooled.syllableUnit(:,1)';
    cTrials = ~vTrials;
    % phoneme decoder
    phonemeUnit = phonemeTrialPooled.phonemeUnit(:,1)';
    phonemeClass = phonemeTrialPooled.phonemeClass(:,1)';
    decoderUnit = syllableUnit;
    dObj = decoderClass(10,[80],1);
    decodeResultStruct = dObj.baseClassify(ieegStructPooled,decoderUnit, d_time_window = [-0.5 0.5], selectChannel = find(roiIdsAll));
    
    cMatAll(iTer,:,:) = decodeResultStruct.cmat;
    accIter(iTer) = decodeResultStruct.accPhoneme;
    
%     modelWeightCell = extractPcaLdaModelWeights(decodeResultStruct.modelWeights, length(unique(decoderUnit)), length(channelNameSelect), 200);
%     decodeResultStruct =[];    
%     modelweightAll(iTer,:,:,:) = modelWeightCell.modelweightchantime;
%     
%     phonemeShuffle = (decoderUnit);
%     phonemeShuffle(cTrials) = shuffle(phonemeShuffle(cTrials));
%     phonemeShuffle(vTrials) = shuffle(phonemeShuffle(vTrials));
%     
%     dObj = decoderClass(10,[80],1);
%     decodeResultStructShuffle = dObj.baseClassify(ieegStructPooled,phonemeShuffle, d_time_window = [-0.5 0.5], selectChannel = find(roiIdsAll));
%     modelWeightCell = extractPcaLdaModelWeights(decodeResultStructShuffle.modelWeights, length(unique(decoderUnit)), length(channelNameSelect), 200);
%     modelweightAllShuffle(iTer,:,:,:) = modelWeightCell.modelweightchantime;
%     decodeResultStructShuffle = [];
end

modelweightsub = abs(modelweightAll-modelweightAllShuffle);
% modelweight = squeeze(sum(sum((modelweightAll),2),4));
% modelweightshuffle = squeeze(sum(sum((modelweightAllShuffle),2),4));
modelweightnorm = squeeze(sum(sum((modelweightsub),2),4));

% plot_activation_density_v2(channelNamePooled, modelweightnorm(1,:)', cLim = [1 3],...
%     gaussFwhm=15, diskThresh = 50, stype='pial', dtype='inflated', colbarTitle =' Beta weights',...
%     transparentPoint=2, colorscale = 'hot');

%save('decoder_Articulator_1000ms_20_fold_mixup_p1.mat',"modelweightnorm","cMatAll","accIter",'-v7.3')
%% Plot confusion matrix
cmatzero = zeros(size(cMatAll,2),size(cMatAll,3));
for iter = 1:size(cMatAll,1)
    cmatzero = cmatzero + squeeze(cMatAll(iter,:,:));
end
%cmatzero = decodeResultStruct.cmat;
cmatzero = cmatzero./size(cMatAll,1);

accTemp = mean(diag(cmatzero))
labels = {'a','ae','i','u','b','p','v','g','k'};        
%  labels = {'low','high','labials','dorsals'}      
 % labels = {'CVC','VCV'}
figure; 
imagesc(cmatzero);
colormap parula
set(gca,'XTick',[1:length(labels)]);
set(gca,'YTick',[1:length(labels)]);
set(gca,'xticklabels',labels)
set(gca,'yticklabels',labels)
cb = colorbar;
caxis([0 max(diag(cmatzero))])
ylabel(cb,'Accuracy (%)');
set(gca,'FontSize',15);
title(['Accuracy : ' num2str(accTemp*100) '%']);
formatTicks(gca)
%%
xlabelsAnat = {'Abs-base','z-score','Log-base','Abs-rel-base','Norm-base'};
figure; 
bar([27.2 34 27.9 29.5 28]);
yline(11.11, '--','chance','LineWidth',2, 'LaberhorizontalAlignment','left');
ylabel('P1 accuracy (%)');
formatTicks(gca)
set(gca,'xticklabels',xlabelsAnat)
%% Subject specific power
channelPower = [];
channelNameAll = [];
timeGamma = linspace(-1,1.5,250);
timeExtract = timeGamma>=-0.25&timeGamma<=0.25;
for iSubject = 1:length(ieegHGStruct)
    warning off;
    iSubject
    for iChan = 1:length(ieegHGStruct(iSubject).channelName)
        channelPower = [channelPower 10.*log10(mean2(squeeze(ieegHGStruct(iSubject).ieegHGNorm.data(iChan,:,timeExtract))))];
    end
    channelNameAll = [channelNameAll ieegHGStruct(iSubject).channelName];
end


%% Extract model weights

% modelweight = squeeze(nanmean(modelweightAll,1));
% modelweightshuffle = squeeze(nanmean(modelweightAllShuffle,1)); 
% modelweightnorm = squeeze(vecnorm(modelweightshuffle,2,1));

modelweight = squeeze(vecnorm(modelweightAll,2,2));
modelweightshuffle = squeeze(vecnorm(modelweightAllShuffle,2,2));
modelweightnorm = squeeze(mean(modelweight ,1));

modelweightChan = vecnorm(modelweightnorm',1)./size(modelweightnorm,2);

%%
load colors.mat
groupingIds =zeros(1,length(channelNameSelect));
hemisphere = 'rh';

bnLabels = extractBNlabels('ifg',hemisphere)
roiIds = extractRoiIds(Subject, channelNameSelect,bnLabels);
groupingIds(roiIds) = 1;


bnLabels = extractBNlabels('rmfg',hemisphere)
roiIds = extractRoiIds(Subject, channelNameSelect,bnLabels);
groupingIds(roiIds) = 2; 


bnLabels = extractBNlabels('cmfg',hemisphere)
roiIds = extractRoiIds(Subject, channelNameSelect,bnLabels)
groupingIds(roiIds) = 3;


bnLabels = extractBNlabels('smc',hemisphere)
roiIds = extractRoiIds(Subject, channelNameSelect,bnLabels);
groupingIds(roiIds) = 4;


bnLabels = extractBNlabels('ipc',hemisphere)
roiIds = extractRoiIds(Subject, channelNameSelect,bnLabels);
groupingIds(roiIds) = 5;


bnLabels = extractBNlabels('insula',hemisphere)
roiIds = extractRoiIds(Subject, channelNameSelect,bnLabels);
groupingIds(roiIds) = 6;


bnLabels = extractBNlabels('astg',hemisphere)
roiIds = extractRoiIds(Subject, channelNameSelect,bnLabels);
groupingIds(roiIds) = 7; 


bnLabels = extractBNlabels('pstg',hemisphere)
roiIds = extractRoiIds(Subject, channelNameSelect,bnLabels);
groupingIds(roiIds) = 8;


bnLabels = extractBNlabels('sts',hemisphere)
roiIds = extractRoiIds(Subject, channelNameSelect,bnLabels);
groupingIds(roiIds) = 9; 


bnLabels = extractBNlabels('heschl',hemisphere)
roiIds = extractRoiIds(Subject, channelNameSelect,bnLabels);
groupingIds(roiIds) = 10;




%% Plot model weights

figure;
hold on;
swarmchart(groupingIds(groupingIds>0),mean(modelweightnorm(:,groupingIds>0),1),10,'filled')
boxplot(mean(modelweightnorm(:,groupingIds>0),1),groupingIds(groupingIds>0),'symbol','','Colors','k','Whisker', inf);

set(gca,'xticklabels',{'IFG','rMFG','cMFG','SMC','IPC','Insula','aSTG','pSTG','STS','Heschl'})
ylabel('Beta weights')
formatTicks(gca)


 %%
plot_activation_density_geodesic_v2(channelNamePooled(groupingIds<9), modelweightChan(groupingIds<9)',...
    cLim = [0 1e-2],...
    gaussFwhm=12, diskThresh = 20, stype='inflated', ptype='inflated', colbarTitle=' Beta weights',...
    transparentPoint=quantile(modelweightChan(groupingIds<9)',0.75), colorscale = [ 1 1 1; colors(2,:)]);


plot_activation_timeSeries_bh_v2(channelNamePooled, modelweightnorm, tw = [-0.5 0.5],...
    cLim = [0 3e-2], movTitle = 'Beta_weight_phoneme_15_mm_gaussian',...
    frameRate=10, gaussFwhm=12, diskThresh = 20, stype='pial', ptype='inflated',...
    colbarTitle=' Beta weight', transparentPoint=5e-3, colorscale = [ 1 1 1; colors(2,:)]);
%% Subject specific decoding
decodeResultStruct = {};
parfor iSubject = 1:length(ieegHGStruct)
    warning off;
    iSubject
    dObj = decoderClass(10,[10:10:90],1);
    decodeResultStruct{iSubject} = dObj.baseClassify(ieegHGStruct(iSubject).ieegHGNorm,trialInfoStruct(iSubject).phonemeTrial.phonemeUnit(:,1)'  , d_time_window = [-0.5 0.5]);
end
%% 

for iSubject = 1:length(ieegHGStruct)
    chanNameSubject = ieegHGStruct(iSubject).channelName;
    roiIds = extractRoiIds(Subject, chanNameSubject,["temporal","sts" ]);
    if(sum(roiIds))
        figure;
        ieegGammaSubj = squeeze(mean(ieegHGStruct(iSubject).ieegHGNorm.data(roiIds,:,:),2));
        if(sum(roiIds)==1)
            ieegGammaSubj= ieegGammaSubj';
        end
        plot(timeGamma,ieegGammaSubj);
        xlabel('Time from response onset')
        title(chanNameSubject{1}(1:3))
    end
end


%% Trial subsampling


numTrials = 30:10:200;
numIterations = 20;

accTrialSamp = cell(length(numTrials), 1);

parfor iTrial = 1:length(numTrials)
    accTrialTemp = zeros(numIterations, 1);
    iSamp = 1;
    while(iSamp<=numIterations)
        try
            [ieegStructPooled, phonemeTrialPooled, channelNamePooled] = poolChannelWithSmoteTrial(ieegHGStruct, trialInfoStruct);
    
            phonemeUnit = phonemeTrialPooled.phonemeUnit(:, 1)';
            pass = 0
            while(pass==0)
                [trialSamp,trialIdx] = datasample(phonemeUnit,numTrials(iTrial),'Replace',false);
                [countTrial] = hist(trialSamp,unique(trialSamp)); %#ok<HIST>
        
                if (min(countTrial) < 3 || length(countTrial) < 9)
                    
                    disp('Repeating iteration')
                else
                    pass = 1;
                end
            end
            dObj = decoderClass(10,[80],1);
            decodeResultStruct = dObj.baseClassify(ieegStructPooled,phonemeUnit, d_time_window = [-0.5 0.5], selectTrial = trialIdx);
            accTrialTemp(iSamp) = decodeResultStruct.accPhoneme;
            iSamp = iSamp+1;
        catch
           disp('Repeating iteration')
        end
    end
    accTrialSamp{iTrial} = accTrialTemp;
end

accTrialSampMat = cell2mat(accTrialSamp');
numTrialsMat = repmat(numTrials,20,1);

accMod = fitlm(numTrialsMat(:),accTrialSampMat(:))

figure;
boxplot(accTrialSampMat,numTrials,'Colors','k')
hold on;
plot(predict(accMod,numTrials'))
xlabel('Trials')
ylabel('Decoding Accuracy')
yline(0.1111, '--','chance','LineWidth',1);
formatTicks(gca)
%% Linear modeling

betas_cvc = [];
Rsquared = [];

beta_p1 = [];
beta_p2 = [];
beta_p3 = [];
beta_syl = [];
mseLoss = [];
channelNameAll = [];

for iSubj = 1:length(ieegHGStruct)
    iSubj
    phonemeTrialSubj = trialInfoStruct(iSubj).phonemeTrial;
    cvcIds = find(phonemeTrialSubj.syllableUnit(:,1)'==2);
    
    p1=categorical(phonemeTrialSubj.phonemeUnit(:,1));
    p2=categorical(phonemeTrialSubj.phonemeUnit(:,2));
    p3=categorical(phonemeTrialSubj.phonemeUnit(:,3));
    syl = logical(phonemeTrialSubj.syllableUnit(:,1)-1);
    
    % p = 0:.5:1;
    % breaks = quantile(phonemeTrialPooled.phonotactic(cvcIds,6),p);
    % pfwd1 = ordinal(phonemeTrialPooled.phonotactic(cvcIds,6),{'Q1','Q2'},...
    %                    [],breaks);
    % 
    % breaks = quantile(phonemeTrialPooled.phonotactic(cvcIds,7),p);
    % pfwd2 = ordinal(phonemeTrialPooled.phonotactic(cvcIds,7),{'Q1','Q2'},...
    %                    [],breaks);
    % 
    % breaks = quantile(phonemeTrialPooled.phonotactic(cvcIds,8),p);
    % pbwd1 = ordinal(phonemeTrialPooled.phonotactic(cvcIds,8),{'Q1','Q2'},...
    %                    [],breaks);
    % 
    
    
    % breaks = quantile(phonemeTrialPooled.phonotactic(cvcIds,9),p);
    % pbwd2 = ordinal(phonemeTrialPooled.phonotactic(cvcIds,9),{'Q1','Q2'},...
    %                    [],breaks);
    
    % pfwd1 = phonemeTrialPooled.phonotactic(cvcIds,6);
    beta_p1_subj = [];
    beta_p2_subj = [];
    beta_p3_subj = [];
    beta_syl_subj = [];
    mseLoss_subj = [];
    ieegsubj = ieegHGStruct(iSubj).ieegHGNorm.data;
    dataSize = size(ieegsubj);
    for iChan = 1:dataSize(1)
        iChan
        
          parfor iTime = 1:dataSize(3) 
                
                hg = squeeze(ieegsubj(iChan,:,iTime))';
                ieegdata = table(hg,p1,p2,p3,syl);
                fit = fitrlinear(ieegdata,'hg~1+p1+p2+p3+syl','Lambda',0.001,...
                    'Regularization','ridge','KFold',10);
    %            fit = fitglm(ieegdata,'hg~p1+p2+p3+syl','Intercept',true,'Link','identity');
                betas_temp = zeros(1,29);
                for iFold = 1:length(fit.Trained)
                    betas_temp = betas_temp + (fit.Trained{iFold,1}.Beta)';
                end
                betas_chan = betas_temp./length(fit.Trained);
    %             length(betas_chan)
                beta_p1_subj(iChan,iTime) = norm(betas_chan(1:9));
                beta_p2_subj(iChan,iTime) = norm(betas_chan(10:18));
                beta_p3_subj(iChan,iTime) = norm(betas_chan(19:27));
                beta_syl_subj(iChan,iTime) = norm(betas_chan(28:29));
    
                mseLoss_subj(iChan,iTime) = kfoldLoss(fit);
                %betas(iChan,iTime,:) = fit.Beta;
    %             Rsquared(iChan,iTime,:) = fit.Rsquared.Adjusted;
            end
        
    end
    beta_p1 = cat(1,beta_p1, beta_p1_subj);
    beta_p2 = cat(1,beta_p2, beta_p2_subj);
    beta_p3 = cat(1,beta_p3, beta_p3_subj);
    beta_syl = cat(1,beta_syl, beta_syl_subj);
    mseLoss = cat(1,mseLoss, mseLoss_subj);
    channelNameAll = [channelNameAll ieegHGStruct(iSubj).channelName];
end

%% Visualize beta

selectIds = ismember(channelNameAll,elecNameProductionClean);
channelNameSelect = channelNameAll(selectIds);
beta_p1_select = beta_p1(:,selectIds,:);
beta_p2_select = beta_p2(:,selectIds,:);
beta_p3_select = beta_p3(:,selectIds,:);
beta_syl_select = beta_syl(:,selectIds,:);

beta_p1_shuff_select = beta_p1_shuff(:,selectIds,:);
beta_p2_shuff_select = beta_p2_shuff(:,selectIds,:);
beta_p3_shuff_select = beta_p3_shuff(:,selectIds,:);

colvals = lines(5)

timeRange = [-1 1.5];
figure;

hold on;
roiIds = extractRoiIds(Subject, channelNameSelect,["rh-parsopercular", "rh-parstriangular" ]);
sum(roiIds)
visTimePlot(timeRange,smoothdata(squeeze(abs(mean((beta_syl_select(:,roiIds,:)./2),1))),2,'gaussian',2),colval = colvals(1,:))
hold on;
visTimePlot(timeRange,smoothdata(squeeze(abs(mean((beta_p1_select(:,roiIds,:)./9),1))),2,'gaussian',2),colval = colvals(2,:))
hold on;
visTimePlot(timeRange,smoothdata(squeeze(abs(mean(abs(beta_p3_select(:,roiIds,:)./9),1))),2,'gaussian',2),colval = colvals(3,:))
xlabel('Time from response onset')
ylabel('Beta')

%%

selectIds = ismember(channelNameAll,elecNameProductionClean);

channelNameSelect = channelNameAll(selectIds);


yrange = [0.2 0.7]
data2plot = smoothdata(squeeze(abs(mean(beta_syl,1))),2,'gaussian',2);
data2plot = data2plot(selectIds,:);
groupingIds = zeros(1,length(channelNameAll));
timeRange = [-1 1.5];
figure;
[p,n]=numSubplots(8)
hold on;
subplot(p(1),p(2),1)
roiIds = extractRoiIds(Subject, channelNameSelect,["supramarginal", "inferiorparietal" ]);
visTimePlot(timeRange,data2plot(roiIds,:),colval = colors(2,:))
groupingIds(roiIds) = 2;
title('IPC')
xlabel('Time from response onset')
ylim(yrange)

subplot(p(1),p(2),2)
roiIds = extractRoiIds(Subject, channelNameSelect,["middlefrontal" ]);
visTimePlot(timeRange,data2plot(roiIds,:),colval = colors(3,:))
groupingIds(roiIds) = 3;
title('MFG')
xlabel('Time from response onset')
ylim(yrange)

subplot(p(1),p(2),3)
roiIds = extractRoiIds(Subject, channelNameSelect,["insula"]);
visTimePlot(timeRange,data2plot(roiIds,:),colval = colors(6,:))
groupingIds(roiIds) = 6;
title('Insula')
xlabel('Time from response onset')
ylim(yrange)

subplot(p(1),p(2),4)
roiIds = extractRoiIds(Subject, channelNameSelect,["superiorfrontal" ]);
visTimePlot(timeRange,data2plot(roiIds,:),colval = colors(7,:))
groupingIds(roiIds) = 7;
title('superiorfrontal')
xlabel('Time from response onset')
ylim(yrange)

subplot(p(1),p(2),5)
roiIds = extractRoiIds(Subject, channelNameSelect,["fusiform" ]);
visTimePlot(timeRange,data2plot(roiIds,:),colval = colors(8,:))
groupingIds(roiIds) = 8;
title('fusiform')
xlabel('Time from response onset')
ylim(yrange)

subplot(p(1),p(2),6)
roiIds = extractRoiIds(Subject, channelNameSelect,["central" ]);
visTimePlot(timeRange,data2plot(roiIds,:),colval = colors(5,:))
groupingIds(roiIds) = 5;
title('SMC')
xlabel('Time from response onset')
ylim(yrange)

subplot(p(1),p(2),7)
roiIds = extractRoiIds(Subject, channelNameSelect,["temporal", "sts"]);
sum(roiIds)
visTimePlot(timeRange,data2plot(roiIds,:),colval = colors(1,:))
groupingIds(roiIds) = 1;
title('MTG & STG')
xlabel('Time from response onset')
ylim(yrange)    

subplot(p(1),p(2),8)
roiIds = extractRoiIds(Subject, channelNameSelect,["rh-parsopercular", "rh-parstriangular" ]);
visTimePlot(timeRange,data2plot(roiIds,:),colval = colors(4,:))
groupingIds(roiIds) = 4;
title('IFG')
xlabel('Time from response onset')
ylim(yrange)
% visTimePlot(timeRange,modelweightnorm(groupingIds==0,:),colval = colors(9,:))
%  groupingIds(groupingIds==0) = 9;

% xlabel('Time from response onset')
% ylabel('Beta')
sgtitle('Beta-syllable')
%% Subject specific analysis
[subjProd,chanSubjProd] = extractSubjectPerChannel(elecNameProductionClean);
subjAnatStat = [];
hemisphere = 'lh';
planLabels = [extractBNlabels('ifg',hemisphere) extractBNlabels('rmfg',hemisphere) extractBNlabels('cmfg',hemisphere) ];
artLabels = [extractBNlabels('smc',hemisphere) extractBNlabels('ipc',hemisphere)];
monLabels = [extractBNlabels('astg',hemisphere) extractBNlabels('pstg',hemisphere) extractBNlabels('insula',hemisphere) extractBNlabels('heschl',hemisphere) extractBNlabels('sts',hemisphere)];
for iSubj = 1:length(subjProd)   
    roiIds = extractRoiIds(Subject, elecNameProductionClean(chanSubjProd==iSubj),planLabels);
    subjAnatStat(iSubj,1) = sum(roiIds);
     roiIds = extractRoiIds(Subject, elecNameProductionClean(chanSubjProd==iSubj),artLabels);
    subjAnatStat(iSubj,2) = sum(roiIds);
     roiIds = extractRoiIds(Subject, elecNameProductionClean(chanSubjProd==iSubj),monLabels);
    subjAnatStat(iSubj,3) = sum(roiIds);
end


%% NNMF factor - averaged data

ieegdatamean = squeeze(nanmean(ieegStructPooled.data,2));
for iFactor = 1:30
    iFactor
    for iRep = 1:20
        [W,H,D(iFactor,iRep)] = nnmf(ieegdatamean',iFactor);
    end
end


[H,W] = nnmf(ieegdatamean',4);
visTimePlot3_v2(timeEpoch,H');
[maxVal,maxId] = max(W);
cfg.elec_colors = lines(5);
plot_subjs_on_average_grouping(channelNamePooled', maxId', 'fsaverage', cfg);
plot_subjs_on_average_activation(channelNamePooled', W(3,:)', 'fsaverage', cfg);

for iFactor = 1:size(H,2)
    ieegStructPooled = ieegStructPooledAll;
    channelNamePooled = channelNamePooledAll;
    for iChan = 1:length(channelNamePooled)
        ieegStructPooled.data(iChan,:,:) = W(iFactor,iChan).*ieegStructPooled.data(iChan,:,:);
    end
    emptyChan = W(iFactor,:)==0;
    ieegStructPooled.data(emptyChan,:,:) = [];
    channelNamePooled(emptyChan) = [];
    save(fullfile(['pooledSubject_' Task.Name '_factor_' num2str(iFactor) ...
        '_stitch_prodelecs_data.mat']),...
        'ieegStructPooled','phonemeTrialPooled',...
        'channelNamePooled','timeEpoch','W','H');
end


%% NNMF factor - stretched data

ieegDataStretch = reshape(ieegStructPooled.data,[size(ieegStructPooled.data,1) size(ieegStructPooled.data,2)*size(ieegStructPooled.data,3)]);
for iFactor = 1:10
    iFactor
    for iRep = 1:20
        [W,H,DStretch(iFactor,iRep)] = nnmf(ieegDataStretch,iFactor,'algorithm','mult','replicates',100);
    end
end


[W,H] = nnmf(ieegDataStretch,4);
HPermute = reshape(H,[size(H,1) size(ieegStructPooled.data,2) size(ieegStructPooled.data,3)]);
figure; plot(squeeze(mean(HPermute,2))');
[maxVal,maxId] = max(W');
cfg.elec_colors = lines(4);
plot_subjs_on_average_grouping(channelNamePooled', maxId', 'fsaverage', cfg);

%% 3d tensor analysis

ieegHiGammaTensor = tensor(permute(ieegStructPooled.data,[1 3 2]));

R = 20;
n_fits = 15;
err = [];

factorVar= [];
similarFactor = [];
modelLoss = [];
OPTS.tol = 1e-5;
for r = 1:R
    
    for n = 1:n_fits
        % fit model
        est_factors = [];
        %ieegHiGammaPercProdTensor = tensor(permute(cat(3,ieegHiGammaPerc,ieegHiGammaProd) ,[1 3 2]));
        est_factors = cp_nmu(ieegHiGammaTensor,r);

        % store error
        err(r,n) = norm(full(est_factors) - ieegHiGammaTensor)/norm(ieegHiGammaTensor);
        factorVar(r,n) = 1-err(r,n).^2;
        similarFactor(r,n) = innerprod(full(est_factors),ieegHiGammaTensor)/(norm(full(est_factors))*norm(ieegHiGammaTensor));
%         trialFact = [];
%         for iCom = 1:r
%             trialFact(iCom,:) = est_factors.u{3}(:,iCom);
%         end
%         modelLossTemp = 0;
%         for iPhon =1:3
%             linearModel = crossval(fitcdiscr(trialFact',phonIndClass(:,iPhon)'));
%             modelLossTemp = modelLossTemp + kfoldLoss(linearModel);
%         end
%         modelLoss(r,n) = modelLossTemp/3;
        % visualize fit for first several fits
    %     if n < 4
    %         % score aligns the cp decompositions
    %         %[sc, est_factors] = score(est_factors, true_factors);
    %         
    %         % plot the estimated factors
    %         viz_ktensor_update(est_factors,chanMap,selectedChannels,timeGammaProd,phonClean(1,:),2);
    % 
    %     end
    end
    
end
%% cluster correction on encoding results
sigMatTmp2=zeros(size(r2,1),size(r2,3));
parfor iChan = 1:length(channelNameAll)
    iChan
    targetData = squeeze(r2(:,iChan,:));
    baseData = squeeze(r2_shuff(:,iChan,:));
    [zValsRawAct, pValsRaw, actClust] = timePermCluster(targetData, baseData, 1000, 1, 1.645);
    sigMatChan = zeros(1,size(targetData,2));
    for iClust=1:length(actClust.Size)
        if actClust.Size{iClust}>actClust.perm99
            sigMatChan(actClust.Start{iClust}: ...
                actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
                = 1;
        end
    end
    sigMatTmp2(iChan,:) = sigMatChan;

end


%% Extract whitematter electrodes
subjInfo = extractRoiInformation(Subject,voxRad=10);
whiteVoxelAll = [];
greyVoxelAll = [];
chanPowerWhite = [];
chanPowerGrey = [];
chanNameGrey = [];
for iSubj = 1:length(subjInfo)
    whiteIds = contains({subjInfo(iSubj).Trest{:,3}},["White","WM"]);
    %greyIds =  ~contains({subjInfo(iSubj).Trest{:,3}},["White","WM"]);
    prodIds = ismember(subjInfo(iSubj).Tname,channelNameAll);
    whiteVoxelAll = [whiteVoxelAll [subjInfo(iSubj).Trest{whiteIds&prodIds',4}]];
    if(~contains({subjInfo(iSubj).Trest{:,5}},["White","WM"]))
        greyVoxelAll = [greyVoxelAll [subjInfo(iSubj).Trest{whiteIds&prodIds',6}]];
    else
        greyVoxelAll = [greyVoxelAll [subjInfo(iSubj).Trest{whiteIds&prodIds',8}]];
    end
    chanNameWhiteSubj = subjInfo(iSubj).Tname(whiteIds&prodIds');
    %chanNameGreySubj = subjInfo(iSubj).Trest(greyIds',3);
    chanNameGrey = [chanNameGrey chanNameGreySubj'];
    chanPowerWhite = [chanPowerWhite channelPower(ismember(channelNameAll,chanNameWhiteSubj))];
    chanPowerGrey = [chanPowerGrey channelPower(ismember(channelNameAll,chanNameGreySubj))];
end
figure; scatter(whiteVoxelAll,greyVoxelAll,15,'filled');
xlabel('Fraction of white matter')
ylabel('Fraction of next maximum grey matter')

figure; scatter(greyVoxelAll,chanPowerWhite,15,'filled');
xlabel('Fraction of next maximum grey matter')
ylabel('HG-ESNR (dB)')

figure; histogram(whiteVoxelAll)
xlabel('Percentage of voxels')
ylabel('Number of electrodes')
%% Extract Grey matter electrodes
subjRoiInfo = extractRoiInformation(Subject,voxRad=10);
chanInfo = assignRoiInformation(subjRoiInfo, gmThresh=0.05);

chanInfo = [];
for iSubj = 1:length(subjInfo)
    greyIds =  ~contains({subjInfo(iSubj).Trest{:,1}},["White","hypointensities","nknown"]);
    TrestGM = subjInfo(iSubj).Trest(greyIds,:);
    TnameGM = subjInfo(iSubj).Tname(greyIds);
    for iChan = 1:size(TrestGM,1)
        % Iterating through grey matter channels
        iCol = 3;
        iColWhite = [];
        while(contains(TrestGM{iChan,iCol},["White","hypointensities","nknown"]))
            % checking for white-matter or unknown info and landing
            % on the first grey matter roi            
            if(contains(TrestGM{iChan,iCol},["White","hypointensities"]))
                % Marking the position of white matter roi, in case
                % there is no grey matter roi available
                iColWhite = [iColWhite iCol];
            end
            iCol = iCol + 2;
        end

        iColAssign = [];
        if(isempty(TrestGM{iChan,iCol}))
            disp('No Grey matter contact found')
            if(~isempty(iColWhite))
                disp('Assigning to the white matter contact location')
                iColAssign = iColWhite(1); 
            else
                iColAssign = 3;
            end
        else
            iColAssign = iCol;
        end
        chanInfo(iSubj).channelName{iChan} = TnameGM{iChan};
        chanInfo(iSubj).channelRoi{iChan} = TrestGM{iChan,iColAssign};
        chanInfo(iSubj).voxelPercent(iChan) = TrestGM{iChan,iColAssign+1};
    end    
end

vox3 = [];
vox10 = [];
sumIncorrect = 0;
chanNameIncorrect = [];
for iSubj = 1:44
    
    correctRoi = strcmp(chanInfo3(iSubj).channelRoi,chanInfo10(iSubj).channelRoi);
    sumIncorrect = sumIncorrect + sum(~correctRoi);
    chanNameIncorrect = [chanNameIncorrect chanInfo3(iSubj).channelName(~correctRoi)];
    vox3 = [vox3 chanInfo3(iSubj).voxelPercent(correctRoi)];
    vox10 = [vox10 chanInfo10(iSubj).voxelPercent(correctRoi)];
end

figure;
scatter(vox3,vox10,15,'filled');
%%

% Read an image
image = imread('path_to_your_image.jpg'); % Replace with your image path

% Convert image to grayscale if it is RGB
if size(image, 3) == 3
    grayImage = rgb2gray(image);
else
    grayImage = image;
end

% Detect edges using the Canny method
edges = edge(grayImage, 'Canny');

% Create a red border on the original image
borderedImage = image;
borderedImage(repmat(edges, [1, 1, size(image, 3)])) = 255;

% Display the original image and the image with borders
figure;
subplot(1,2,1);
imshow(image);
title('Original Image');

subplot(1,2,2);
imshow(borderedImage);
title('Image with Borders');
