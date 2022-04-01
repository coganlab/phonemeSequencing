addpath(genpath('E:\Box Sync\Box Sync\ECoG_Recon\matlab_code\'));
global DUKEDIR
DUKEDIR = 'E:\Box Sync\Box Sync\CoganLab\D_Data\Phoneme_Sequencing';
dLabels = dir(DUKEDIR);
dLabels = dLabels(3:end);
resFold = 'E:\inUnitData\phonemeSequencing\decoding\results\v8';
tw = [-2 2]; % epoch time window
etw = [-1 1]; % high gamma time window

prtw = [-0.5 0]; % preonset time window
pstw = [0.25 0.75]; % postonset time window
gammaF = [70 150]; % frequency in Hz
Task.Name = 'Phoneme_Sequencing';
Subject = popTaskSubjectData(Task);
%% Iterating through subjects
for iSubject = 1:length(Subject)
    % Update this code
    disp(['Loading Subject data:' dLabels(iSubject).name]);    
        d = []; ieegCarResp = []; ieegCarImpFilt = []; ieegGamma = []; ieegSplit = [];
        Experiment = Subject(iSubject).Experiment;
        fsD = Experiment.recording.sample_rate;
        Trials = Subject(iSubject).Trials;
        allChannels = ({Subject(iSubject).ChannelInfo.Name});
        badChannels = Subject(iSubject).badChannels;
        trialInfo = Subject(iSubject).trialInfo;
    disp('done')
    %% Recon Visualization
%     iSubject = 1;
%     dLabels(iSubject).name = 'D2';
%     cfg = [];
%     % cfg.alpha = 0.4;
%     % cfg.use_brainshifted = 1;
%     handles = plot_subj(dLabels(iSubject).name, cfg);
    disp('Extracting anatomical channels')
        channelName = {Subject(iSubject).ChannelInfo.Location};
        channelName(cellfun(@isempty,channelName)) = {'dummy'};
        sensorymotorChan = contains(channelName,'central');
        whiteChan = contains(channelName,'white');
        ifgChan = contains(channelName,'opercularis');
        frontalChan = contains(channelName,'front');
        temporalChan = contains(channelName,'temporal');
        anatChan = sensorymotorChan|ifgChan|frontalChan|temporalChan;
        disp(['Number of anatomical channels : ' num2str(sum(anatChan))]);
    disp('done')

%% Loading all the data
    
    
    disp('Loading IEEG data'); 
        [ieegSplit]=trialIEEG(Trials,1:length(allChannels),'ResponseStart',tw.*1000);
        ieegSplit = permute(ieegSplit,[2,1,3]);
        ieegBase = squeeze(trialIEEG(Trials,1:length(allChannels),'Start',tw.*1000));
        ieegBase = permute(ieegBase,[2,1,3]);
        respId = find(~[Trials.NoResponse]);
        ieegSplitResp = ieegSplit(:,respId,:);
        ieegBaseResp = ieegBase(:,respId,:);
    disp('done');
    
    

%% Common average referencing
    disp('Common average referencing');
        ieegCarBase = carFilterImpedance(ieegBaseResp,badChannels);
        ieegCarSplit = carFilterImpedance(ieegSplitResp,badChannels);

        goodChannels = setdiff(1:size(ieegSplitResp,1),badChannels);

        [~,goodtrials] = remove_bad_trials(ieegCarSplit,14);
        %goodTrialsCommon = extractCommonTrials(goodtrials(anatChan));
        ieegBaseClean = ieegCarBase(:,1:156,:);
        ieegCarClean = ieegCarSplit(:,1:156,:);
    disp('done');
 %% High Gamma Extraction 
    disp('Extracting High Gamma')
        fsDown =200;
        gInterval = [70 150];
        normFactor = [];
        ieegGammaBasedown = [];
        ieegGammaRespdown = [];
        disp('Extracting Baseline HG');
        for iTrial = 1:size(ieegBaseClean,2)
            
            [~,ieegGammaBasedown(:,iTrial,:)] = EcogExtractHighGammaTrial(double(squeeze(ieegBaseClean(:,iTrial,:))),fsD,fsDown,[gInterval(1) gInterval(2)],tw,prtw,[],[]);
            [~,ieegGammaRespdown(:,iTrial,:)] = EcogExtractHighGammaTrial(double(squeeze(ieegCarClean(:,iTrial,:))),fsD,fsDown,[gInterval(1) gInterval(2)],tw,[-0.25 0.25],[]);
        end
    %     
        ieegGammaRespPower = (squeeze(mean(log10(ieegGammaRespdown.^2),3)));
        ieegGammaBasePower = (squeeze(mean(log10(ieegGammaBasedown.^2),3)));
        pChan = [];

        for iChan = 1:size(ieegGammaBasedown,1)
            normFactor(iChan,:) = [mean2(squeeze(ieegGammaBasedown(iChan,:,:))) std2(squeeze(ieegGammaBasedown(iChan,:,:)))];
            pChan(iChan) = permtest_sk((ieegGammaRespPower(iChan,:)),(ieegGammaBasePower(iChan,:)),10000);
    %         ieegGammaPowerNorm(iChan) = 10.*log10(mean(ieegGammaRespPower(iChan,goodtrials{iChan})./mean(ieegGammaBasePower(iChan,goodtrials{iChan}))));
    %         evokeSnr(iChan) = esnr(squeeze(ieegGammaRespdown(iChan,goodtrials{iChan},:)),squeeze(ieegGammaBasedown(iChan,goodtrials{iChan},:)));
        end
         [p_fdr, pvalsMCleanProd] = fdr( pChan, 0.05);
        disp(['Number of significant HG channels : ' num2str(sum(pvalsMCleanProd))]);
        ieegGamma = [];
        disp('Normalizing HG');
        for iTrial = 1:size(ieegCarClean,2)
            

              [~,ieegTmp] = EcogExtractHighGammaTrial(double(squeeze(ieegCarClean(:,iTrial,:))),fsD,fsDown,[gInterval(1) gInterval(2)],tw,etw,squeeze(normFactor),2);
              ieegGamma(:,iTrial,:) = ieegTmp;    

        end
        ieegGamma = squeeze(ieegGamma);
        timeGamma = linspace(etw(1),etw(2),size(ieegGamma,3));
    disp('done');
%% Phoneme feature labeling
    disp('Labeling phoneme features');
        trialNames = [];
        trialInfoResp = trialInfo(respId);
        trialClean = trialInfoResp(1:156);
        phonCVClass = [];
        phonIndClass = [];

        for iTrial = 1:length(trialClean)
            phonPairParse = [];
            if(iscell(trialClean))
                trialNames = (trialClean{iTrial}.sound(1:end-4));
            else
                trialNames = (trialClean(iTrial).sound(1:end-4));
            end
            trialNames = strrep(trialNames,'ae','z');
            trialNames = num2cell(trialNames);
            trialNames = strrep(trialNames,'z','ae');
            for iPhon = 1:3
                 [~,phonCVClass(iTrial,iPhon),phonIndClass(iTrial,iPhon)] = phonemeEncoder(trialNames{iPhon});          
            end
        end
    disp('done');
%% HG power extraction & Phoneme decoding
    sigChannel = [find(anatChan) ];
    if(~isempty(sigChannel))
        accChan = [];
        phonError = [];
        disp('Spectrogram extraction');
            [specCarBase,pPerc]= getSpectrograms(ieegBaseClean,[],tw,etw,[1 200],[-0.5 0],[0.25 0.75],[70 150],fsD,1);
            [~,~,meanFreqChanOut] = specChanMap(specCarBase,[],allChannels,[],etw,[-0.5 -0],[0.25 0.75],[1 200],[70 150],[0.7 1.4],0,[]);
             [specCarResp,pProd]= getSpectrograms(ieegCarClean,[],tw,etw,[1 200],[-1.5 -1],[-0.25 0.25],[70 150],fsD,1);
        disp('done');
        specPower = [];
        CMatCat = [];  
        timeRange = -1:0.01:0.75;
        labels = {'a','ae','i','u','b','p','v','g','k'};
        accTime = zeros(3,length(timeRange));
        aucTime = [];
        accPhoneme = 0;
        accPhonemeUnbias = 0;
        accPhonemeChance = 0;
        if(sum(pvalsMCleanProd&anatChan)~=0)
            disp('Phoneme decoding...')
                for iPhon = 1:3
                    disp(['Position : ' num2str(iPhon)]);
                    CmatPhoneme = zeros(9,9);
                    for iTer = 1:3
                        [~,ytestAll,ypredAll,optimVarAll] = pcaLinearDecoderWrap(ieegGamma(pvalsMCleanProd&anatChan,:,:),phonIndClass(:,iPhon)',etw,[-0.5 0.5],80,20,0);
                       %[~,ytestAll,ypredAll,aucAll] = linearDecoder(ieegModel(pvalsMCleanProd&anatChan,:,:,:),phonIndClass,[0 1],[0 1],20,0);
                        CmatAll = confusionmat(ytestAll,ypredAll);
                        CmatPhoneme = CmatPhoneme + CmatAll;
                    end
                    CmatCatNorm = CmatPhoneme./sum(CmatPhoneme,2);
                    accPhonemeUnbias(iPhon) = trace(CmatPhoneme)/sum(CmatPhoneme(:));
                    accPhoneme(iPhon) = trace(CmatCatNorm)/size(CmatCatNorm,1);
                    [~,ytestAll,ypredAll] = pcaLinearDecoderWrap(ieegGamma(pvalsMCleanProd&anatChan,:,:),shuffle(phonIndClass(:,iPhon)'),etw,[-0.5 0.5],80,20,0);
                    cmatshuff = confusionmat(ytestAll,ypredAll);
                    cmatshuff = cmatshuff./sum(cmatshuff,2);
                    accPhonemeChance(iPhon) = trace(cmatshuff)/size(cmatshuff,1);       
                    [phonError(iPhon),cmatvect,phonemeDistVect] = phonemeDistanceError(CmatCatNorm,1:9);
                end
            disp('done')

            disp('Temporal Generalization');
                
                for iTime = 1:length(timeRange)
                    
                    timeSelect = timeGamma>=timeRange(iTime)&timeGamma<=timeRange(iTime)+0.25;
                    ieegGammaSelect =(ieegGamma(pvalsMCleanProd&anatChan,:,timeSelect));

                    for iPhone = 1:3

                        CmatCat = zeros(9,9);

                        for iTer = 1:3
                            [~,ytestAll,ypredAll] = pcaLinearDecoderWrap(ieegGammaSelect,phonIndClass(:,iPhone)',[0 1],[0 1],80,20,0);
                            %[~,ytestAll,ypredAll,aucAll,stmfTemplate] = stmfDecodeWrap(ieegModel,phonIndClass(:,iPhone)',[0 1],[0 1],20,0);
                            CmatAll = confusionmat(ytestAll,ypredAll);
    %                         acc = trace(CmatAll)/sum(sum(CmatAll));        
                            CmatCat = CmatCat+CmatAll;
                            %aucCat = aucCat + mean(aucAll,1);
                        end
                        CmatCatNorm = CmatCat./sum(CmatCat,2);
                        accAll = trace(CmatCatNorm)/size(CmatCatNorm,1);
                        accTime(iPhone,iTime) = accAll;
                        %aucTime(iPhone,iTime,:) = aucCat/iTer;
                    end

                end
              disp('done'); 

        end
            disp('Individual channel decoding');
                for iChan = 1:size(ieegGamma,1)
                    
                    gammaFlatBin = [(ieegGamma(iChan,:,:))];
                    phonClassBin = [phonIndClass(:,1)'  ];
                    Cmatiter = zeros(9,9);
                    for iTer = 1:3
                        [~,ytestAll,ypredAll] = pcaLinearDecoderWrap(gammaFlatBin,phonClassBin,etw,[-0.25 0.5],80,20,0);
                    %[~,ytestAll,ypredAll,aucAll] = linearDecoder(gammaFlatBin,phonClassBin,[0 1],[0 1],10,1);
                        CmatAll = confusionmat(ytestAll,ypredAll);
                        Cmatiter = Cmatiter + CmatAll;
                    end
                        CmatiterNorm = Cmatiter./sum(Cmatiter,2);
                        acc = trace(CmatiterNorm)/size(CmatiterNorm,1);
                        accChan(iChan) = acc;
                        CMatCat(iChan,:,:) = CmatiterNorm;

                        [~,specPowerTemp] = specChanMap(specCarResp,[],allChannels,[],etw,[-1.5 -1],[-0.25 0.25],[1 200],[70 150],[],iChan, meanFreqChanOut);
                        specPower = [specPower; specPowerTemp];
                end
                ieegGammaSig  = squeeze(ieegGamma(pvalsMCleanProd&anatChan,:,:));
            disp('done');        
          save(strcat(resFold,'\',dLabels(iSubject).name,'_phonDecodeMotorHGPackSigChannelCarAnatSignificant.mat'),'accChan','allChannels','CMatCat','anatChan','specPower','pvalsMCleanProd','accChan','accPhoneme','accPhonemeUnbias','accPhonemeChance','phonError','ieegGammaSig','phonIndClass','timeRange','accTime');
    end
close all;
end