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

selectRoi = '';
respTimeThresh = 0;
timeEpoch = [-0.5000    2.0000;   -0.5000    1.0000;   -1.0000    1.5000];

%trialInfoStruct = extractTrialInfo(Subject, remFastResponseTimeTrials=respTimeThresh);

ieegHGStructAuditory = extractHGDataWithROI(Subject,baseName = 'Start',...
    Epoch = 'Auditory', roi = selectRoi,Time= timeEpoch(1,:),respTimeThresh=respTimeThresh,...
    subsetElec=elecNameProductionClean,remWMchannels=true);
hgNormFactor = {ieegHGStructAuditory(:).normFactor};
ieegHGStructGo = extractHGDataWithROI(Subject,baseName = 'Start',...
    Epoch = 'Go', roi = selectRoi, Time=timeEpoch(2,:),respTimeThresh=respTimeThresh,...
    subsetElec=elecNameProductionClean,normFactor=hgNormFactor,remWMchannels=true);
ieegHGStructResponse = extractHGDataWithROI(Subject,baseName = 'Start',...
    Epoch = 'ResponseStart', roi = selectRoi, Time=timeEpoch(3,:),respTimeThresh=respTimeThresh,...
    subsetElec=elecNameProductionClean,normFactor=hgNormFactor,remWMchannels=true);

% Remove empty subjects
emptyIds = [];
ieegHGStruct = ieegHGStructAuditory;
trialInfoStruct = [];
for iSubject = 1:length(Subject)
    if(isempty(ieegHGStruct(iSubject).ieegHGNorm))
        emptyIds = [emptyIds iSubject];
    else
        ieegHGStruct(iSubject).ieegHGNorm.data = cat(3,ieegHGStructAuditory(iSubject).ieegHGNorm.data,...
            ieegHGStructGo(iSubject).ieegHGNorm.data,ieegHGStructResponse(iSubject).ieegHGNorm.data);
        ieegHGStruct(iSubject).ieegHGNorm.name = ieegHGStruct(iSubject).ieegHGNorm.name + '_Go' + '_ResponseOnset';
       trialInfoStruct(iSubject).phonemeTrial = phonemeSequenceTrialParser(ieegHGStruct(iSubject).trialInfo);   
    end

end

ieegHGStruct(emptyIds) = [];
trialInfoStruct(emptyIds) = [];

% Pooling across channels based on minimum trial matching
[ieegStructPooled,phonemeTrialPooled,channelNamePooled] = poolChannelWithMinTrial(ieegHGStruct,trialInfoStruct);
save(fullfile(['pooledSubject_' Task.Name '_' ...
    '_stitch_prodelecs_data.mat']),...
    'ieegStructPooled','phonemeTrialPooled',...
    'channelNamePooled','timeEpoch');

%% NNMF factor - averaged data

ieegdatamean = squeeze(mean(ieegStructPooled.data,2));
for iFactor = 1:10
    iFactor
    for iRep = 1:20
        [W,H,D(iFactor,iRep)] = nnmf(ieegdatamean,iFactor);
    end
end


[W,H] = nnmf(ieegdatamean,5);
figure; plot(H');
[maxVal,maxId] = max(W');
cfg.elec_colors = lines(5);
plot_subjs_on_average_grouping(channelNamePooled', maxId', 'fsaverage', cfg);

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

%% Extracting factors
R = 10
OPTS.tol = 1e-5;
ieegHiGammaTensor = tensor(permute(ieegStructPooled.data,[1 3 2]));
est_factors = cp_nmu(ieegHiGammaTensor,R);
elecFact = est_factors.u{1};
timeFact = est_factors.u{2};
trialFact = est_factors.u{3};

%% Visualizing factors

%viz_ktensor_freesurfer(est_factors,channelNamePooled',[0 2;-1.5 1; -1 1.5] ,phonemeTrialPooled,R);
% viz_ktensor_phonotactic(est_factors,channelNamePooled',[0 2;-1.5 1; -1 1.5] ,phonemeTrialPooled,1:5);
% viz_ktensor_phonotactic(est_factors,channelNamePooled',[0 2;-1.5 1; -1 1.5] ,phonemeTrialPooled,6:10);
viz_ktensor_phonotactic_v2(est_factors,channelNamePooled',[-0.5 2;-0.5 1; -1 1.5] ,phonemeTrialPooled,1:5);
viz_ktensor_phonotactic_v2(est_factors,channelNamePooled',[-0.5 2;-0.5 1; -1 1.5] ,phonemeTrialPooled,6:10);
%viz_ktensor_phonotactic_v3(est_factors,channelNamePooled',[0 2;-1.5 1; -1 1.5] ,phonemeTrialPooled,6:10);
%viz_ktensor_phonotactic(est_factors,channelNamePooled',[0 2;-1.5 1; -1 1.5] ,phonemeTrialPooled,11:15);

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
    phonemeTrialPooled.phonemeClass(:,1)',timeRes=timeRes,timeWin=timeWin,selectTrials= cvcIds);
% 2D temporalGeneralization            
disp('2D temporal generalization');
decodeTimeStruct2D = phonDecode.tempGenClassify2D(ieegStructPooled,...
    phonemeTrialPooled.phonemeClass(:,1)',timeRes=timeRes,timeWin=timeWin,selectTrials= cvcIds);



%save('pooledSubject_smc_time_gen_articulator_responseOnset.mat','decodeTimeStruct1D','decodeTimeStruct2D');


disp('Saving..');
if ~exist([DUKEDIR saveFolder]) %#ok<EXIST> 
    mkdir([DUKEDIR saveFolder])
end
save(fullfile([DUKEDIR saveFolder 'pooledSubject_' Task.Name '_' ...
    selectRoi '_' fieldEpoch '_Start_labial_vs_dorsal_resptime_0ms_decoded.mat']),...
    'decodeTimeStruct1D','decodeTimeStruct2D',...
    'channelNamePooled');
