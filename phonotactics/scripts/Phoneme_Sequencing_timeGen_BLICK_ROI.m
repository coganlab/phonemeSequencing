global TASK_DIR
%
global DUKEDIR
global BOX_DIR

%BOX_DIR='H:\Box Sync';
RECONDIR=[BOX_DIR '\ECoG_Recon'];

fDown = 200; %Downsampled Sampling Frequency
timeExtract = [-2 2];

Task=[];

Task.Name='Phoneme_Sequencing';
Subject = popTaskSubjectData(Task);

Task.Base.Name='Start';
Task.Base.Epoch='Start';
Task.Base.Time=[-0.5 0];


TASK_DIR=([BOX_DIR '/CoganLab/D_Data/' Task.Name]);
DUKEDIR=TASK_DIR;

% Task conditions
Task.Conds(1).Name='All';
Task.Conds(1).Field(1).Name='AuditorywDelay';
Task.Conds(1).Field(1).Epoch='Auditory';
Task.Conds(1).Field(1).Time=[0 1.5];
Task.Conds(1).Field(2).Name='DelaywGo';
Task.Conds(1).Field(2).Epoch='Go';
Task.Conds(1).Field(2).Time=[-1 1.5];
Task.Conds(1).Field(3).Name='Response';
Task.Conds(1).Field(3).Epoch='ResponseStart';
Task.Conds(1).Field(3).Time=[-1 1];

% Regression variable
% Use one of the following: 'POS1', 'POS1', etc.
phonRegress = 'BLICK';

% Decoder variable
numFold = 20;
varExplained = 80;
% varExplained = [10:10:90]; Picking desired variance through
% cross-validation; Computtionally expensive
timeWin = 0.25; % Window length for temporal generalization (s)
timeRes = 0.025; % Window hop for temporal generalization (s)

%SNList=1:29; %length(Subject);



% pick a label, generate channel (in goodChannel space) index and subject
% index (in SNList space)

%labelVals{1}='Plan_temp';
labelVals{1}='central';
labelVals{2}='postcentral';
labelVals{2}='temporal';
labelVals{3}='supram';
labelVals{4}='caudal';
%labelVals{6}='Operc';
%labelVals{5}='parietal';
%labelVals{6}='Triang';
%
%labelVals{8}='S_central';
%%
for iL=1:length(labelVals)
    SNList=setdiff(1:length(Subject),38:40);
    
    label=labelVals{iL};
    
    labelChan=[];
    labelChanName={};
    counter=0;
    for iSN=1:length(SNList)
        for iC=1:length(Subject(SNList(iSN)).goodChannels)
            if ~isempty(Subject(SNList(iSN)).ChannelInfo(Subject(SNList(iSN)).goodChannels(iC)).Location)
                if contains(Subject(SNList(iSN)).ChannelInfo(Subject(SNList(iSN)).goodChannels(iC)).Location,label,'IgnoreCase',1)
                    labelChan(counter+1,2)=iC;
                    labelChan(counter+1,1)=iSN;
                    labelChanName{counter+1}=Subject(SNList(iSN)).ChannelInfo(Subject(SNList(iSN)).goodChannels(iC)).Name;
                    counter=counter+1;
                end
            end
        end
    end
    
    load([TASK_DIR '/forKumar/prodDelayElecs3.mat']); %2 vs 3
    
    counter=0;
    sigDelayIdx=[];
    for iChan1=1:length(labelChanName)
        for iChan2=1:length(prodDelayElecs);
            if strcmp(labelChanName{iChan1},prodDelayElecs{iChan2})
                sigDelayIdx(counter+1)=iChan1;
                counter=counter+1;
            end
        end
    end
    
    % only significant channels!
    labelChan=labelChan(sigDelayIdx,:);
    labelChanName=labelChanName(sigDelayIdx);
    
    % create new subject list with only ROI subjects
    %SNList=unique(labelChan(:,1));
    %SNList=unique(labelChan(:,1));
    % for now, remove subjects with only 1 channel (Kumar needs to fix his
    % code)
%     counter=0;
%     SNList2=[];
%     for iSN=1:length(SNList);
%         ii=find(labelChan(:,1)==SNList(iSN));
%         if length(ii)>1
%             SNList2(counter+1)=SNList(iSN);
%             counter=counter+1;
%         end
%     end
%     
  %  SNList=SNList2;
    
    ieegFieldHGAll={};
    phonemeTrialAll={};
    
    %SNList=5;
    channelNamesAll=[];
    for iSN=1:length(SNList)
        SN=SNList(iSN);
        Subject(SN).Name
        Trials=Subject(SN).Trials;
        counterN=0;
        counterNR=0;
        noiseIdx=0;
        noResponseIdx=0;
        
        [ii jj]=find(labelChan(:,1)==SN);
        chanIdxAll=Subject(SN).goodChannels;
        chanIdx=Subject(SN).goodChannels(labelChan(ii,2));
        
        % prodDelayElecs is loaded from the mat file in 'forKumar'
        % Come up with automated way to parse necessary channels
        channelNames = {Subject(SN).ChannelInfo.Name};
        channelNames(cellfun(@isempty,channelNames)) = {'dummy'};
        channelNames = channelNames(chanIdx);
        channelNamesAll=cat(2,channelNamesAll,channelNames);
        
        % will sub-select for significance later
        
        %     [~,chan2select] = intersect(channelNames,prodDelayElecs);
        %     if(isempty(chan2select))
        %         disp('No channels found; Iterating next subject');
        %         continue;
        %         % Forces the iteration for next subject;
        %     end
        
        
        for iTrials=1:length(Trials)
            if Trials(iTrials).Noisy==1
                noiseIdx(counterN+1)=iTrials;
                counterN=counterN+1;
            end
            if Trials(iTrials).NoResponse==1
                noResponseIdx(counterNR+1)=iTrials;
                counterNR=counterNR+1;
            end
        end
        
        condIdx=ones(length(Subject(SN).Trials),1);
        
        baseEpoch=Task.Base.Epoch;
        baseTimeRange=Task.Base.Time;
       % Trials=Subject(SN).Trials(setdiff(1:length(Subject(SN).Trials),noiseIdx));
        
       RespTime=[];
       for iTrials=1:length(Trials);
           if ~isempty(Trials(iTrials).ResponseStart)
               RespTime(iTrials)=(Trials(iTrials).ResponseStart-Trials(iTrials).Go)./30000;
           else
               RespTime(iTrials)=0;
           end
       end
       iiResp=find(RespTime<0);
       
        
        
        
        ieegBase=trialIEEG(Trials,chanIdxAll,baseEpoch,timeExtract.*1000);
        ieegBase = permute(ieegBase,[2,1,3]);
        fs = Subject(SN).Experiment.recording.sample_rate;
        ieegBaseStruct = ieegStructClass(ieegBase, fs, timeExtract, [1 fs/2], baseEpoch);
        clear ieegBase;
        ieegBaseCAR=extractCar(ieegBaseStruct);
        [c,iA,iB]=intersect(chanIdx,chanIdxAll); % get it back to goodChannel Index  
        if(isempty(iB))
            disp('No channels found; Iterating next subject');
            continue;
            %Forces the iteration for next subject;
        end
        ieegBaseCAR.data=ieegBaseCAR.data(iB,:,:);
        clear ieegBaseStruct;      
        ieegBaseHG = extractHiGamma(ieegBaseCAR,fDown,baseTimeRange);
        normFactorBase = extractHGnormFactor(ieegBaseHG);
        for iC=1:length(Task.Conds)
            %trials2Select = setdiff(find(condIdx==iC),cat(2,noiseIdx,noResponseIdx));
            trials2Select = setdiff(find(condIdx==iC),cat(2,noiseIdx,noResponseIdx,iiResp));

            Trials=Subject(SN).Trials(trials2Select);
            trailInfo = Subject(SN).trialInfo(trials2Select);
            % Phoneme Sequence trial parsing
            phonemeTrial = phonemeSequenceTrialParser(trailInfo);
            
            % Iterting through field: 'Auditory', 'Delay', 'Production'
            for iF=1:length(Task.Conds(iC).Field)
                
                Epoch=Task.Conds(iC).Field(iF).Epoch;
                fieldTimeRange=Task.Conds(iC).Field(iF).Time;
                
                ieegField=trialIEEG(Trials,chanIdxAll,Epoch,timeExtract.*1000);
                if ndims(ieegField)==3
                    ieegField = permute(ieegField,[2,1,3]);
                else
                    ieegField=reshape(ieegField,1,size(ieegField,1),size(ieegField,2));
                end
                ieegFieldStruct = ieegStructClass(ieegField, fs, timeExtract, [1 fs/2], Epoch);
                clear ieegField
                % Common average referencing
                ieegFieldCAR = extractCar(ieegFieldStruct);
                ieegFieldCAR.data=ieegFieldCAR.data(iB,:,:);

                clear ieegFieldStruct;
                % Normalized High gamma extraction
                ieegFieldHG = extractHiGamma(ieegFieldCAR,fDown,fieldTimeRange,normFactorBase,2);
                
                ieegFieldHGAll{iSN}{iC}{iF}.ieegFieldHG=ieegFieldHG;
                phonemeTrialAll{iSN}{iC}{iF}.phonemeTrial=phonemeTrial;
            end
        end
    end
    ieegFieldHGAll = ieegFieldHGAll(~cellfun(@isempty, ieegFieldHGAll));
    phonemeTrialAll = phonemeTrialAll(~cellfun(@isempty, phonemeTrialAll));
    %%
    %for each condition and each field, create a minimum trial size feature matrix
    
    for iC=1:length(Task.Conds)
        for iF=1:length(Task.Conds(iC).Field)
            condMat=zeros(length(ieegFieldHGAll),126); % 126 specific triggers, but zeros in the middle
            for iSN=1:length(ieegFieldHGAll)
                for iTrial=1:length(phonemeTrialAll{iSN}{iC}{iF}.phonemeTrial.tokenIdentity)
                    condMat(iSN,phonemeTrialAll{iSN}{iC}{iF}.phonemeTrial.tokenIdentity(iTrial))=condMat(iSN,phonemeTrialAll{iSN}{iC}{iF}.phonemeTrial.tokenIdentity(iTrial))+1;
                end
            end
            
            % get minimum trial numbers for each token
            minCondMat=min(condMat);
            idxTokens=find(minCondMat>0);
            % now need to create feature vector that's common across all
            % subjects
            
            
            data=[];
            
            for iSN=1:length(ieegFieldHGAll)
                syllableUnit = [];
                phonemeUnit = [];
                phonemeClass = [];
                phonoTactic = [];
                tokenIdentity = [];
                dataTmp=[];
                
                for iT=1:length(idxTokens)
                    
                    ii=find(phonemeTrialAll{iSN}{iC}{iF}.phonemeTrial.tokenIdentity==idxTokens(iT));
                    sIdx=shuffle(1:length(ii));
                    sIdx=sIdx(1:minCondMat(idxTokens(iT)));
                    
                    
                    syllableUnit = cat(1,syllableUnit, phonemeTrialAll{iSN}{iC}{iF}.phonemeTrial.syllableUnit(ii(sIdx),:));
                    phonemeUnit = cat(1,phonemeUnit, phonemeTrialAll{iSN}{iC}{iF}.phonemeTrial.phonemeUnit(ii(sIdx),:));
                    phonemeClass = cat(1,phonemeClass, phonemeTrialAll{iSN}{iC}{iF}.phonemeTrial.phonemeClass(ii(sIdx),:));
                    phonoTactic = cat(1,phonoTactic, phonemeTrialAll{iSN}{iC}{iF}.phonemeTrial.phonoTactic(ii(sIdx),:));
                    tokenIdentity = cat(1,tokenIdentity, phonemeTrialAll{iSN}{iC}{iF}.phonemeTrial.tokenIdentity(ii(sIdx),:));
                    
                    dataTmp = cat(2,dataTmp, ieegFieldHGAll{iSN}{iC}{iF}.ieegFieldHG.data(:,ii(sIdx),:));
                end
                data=cat(1,data,dataTmp);
            end
            
            ieegFieldHG.data = data;
            ieegFieldHG.name=strcat(Task.Conds(iC).Field(iF).Name,'_CAR_High-Gamma-Normalized');
            ieegFieldHG.tw=Task.Conds(iC).Field(iF).Time;
            phonemeTrial.syllableUnit = syllableUnit;
            phonemeTrial.phonemeUnit = phonemeUnit;
            phonemeTrial.phonemeClass = phonemeClass;
            phonemeTrial.phonoTactic = phonoTactic;
            phonemeTrial.tokenIdentity = tokenIdentity;
            
            phonDecode = phonemeDecoderClass(numFold,varExplained);
            disp('1D temporal generalization');
            % 1D temporalGeneralization
            %decodeTimeStruct1D = tempGenRegress1D(phonDecode,ieegFieldHG,phonemeTrial,phonRegress,timeRes,timeWin,find(chan2select));
            decodeTimeStruct1D = tempGenRegress1D(phonDecode,ieegFieldHG,phonemeTrial,phonRegress,timeRes,timeWin,1:size(ieegFieldHG.data,1));
            
            % 2D temporalGeneralization
            disp('2D temporal generalization');
            %decodeTimeStruct2D = tempGenRegress2D(phonDecode,ieegFieldHG,phonemeTrial,phonRegress,timeRes,timeWin,find(chan2select));
            decodeTimeStruct2D = tempGenRegress2D(phonDecode,ieegFieldHG,phonemeTrial,phonRegress,timeRes,timeWin,1:size(ieegFieldHG.data,1));
            
            
            
            
            % decoding starts here
            
            % Initialize Decoder object
            
            channelsUsed =  channelNamesAll;
            
            %
            disp('Saving..');
            if ~exist([DUKEDIR '\TempDecode\tempGen\'])
                mkdir([DUKEDIR '\TempDecode\tempGen\'])
            end
            save([DUKEDIR '\TempDecode\tempGen\' label '_' Task.Name '_' ...
                Task.Conds(iC).Name '_' Task.Conds(iC).Field(iF).Name '_' Task.Base.Name '_' phonRegress 'decoded.mat'],'decodeTimeStruct1D','decodeTimeStruct2D','channelsUsed');
        end
    end
end

