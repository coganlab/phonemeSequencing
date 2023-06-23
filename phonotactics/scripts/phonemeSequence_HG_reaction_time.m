 global BOX_DIR
%BOX_DIR='C:\Users\gcoga\Box';
% BOX_DIR='H:\Box Sync';
%RECONDIR='C:\Users\gcoga\Box\ECoG_Recon';
saveFolder = '\TempDecode\temporal_nosig\';
BOX_DIR = 'C:\Users\sd355\Box'
RECONDIR='C:\Users\sd355\Box\ECoG_Recon';
fDown = 200; %Downsampled Sampling Frequency


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
Task.Conds(1).Field(1).Time=[-0.5 3];
Task.Conds(1).Field(2).Name='DelaywGowResponseTime';
Task.Conds(1).Field(2).Epoch='Go';
Task.Conds(1).Field(2).Time=[-1 1.5];
% Task.Conds(1).Field(2).Name='Response';
% Task.Conds(1).Field(2).Epoch='ResponseStart';
% Task.Conds(1).Field(2).Time=[-1 1.5];
timePad = 0.5;
%%
for iSubject=1:length(Subject)
     
    Subject(iSubject).Name
    Trials=Subject(iSubject).Trials;
    counterN=0;
    counterNR=0;
    noiseIdx=0;
    noResponseIdx=0;
    
    
    chanIdx = 1:length(Subject(iSubject).ChannelInfo);
    goodChan=Subject(iSubject).goodChannels;
    badChan = setdiff(chanIdx,goodChan);
    anatName = {Subject(iSubject).ChannelInfo.Location};
    noNameIds = (cellfun(@isempty,anatName))
    if(~isempty(find(noNameIds, 1)))
        save(fullfile(DUKEDIR,Subject(iSubject).Name,'noNameChannels.mat'),"noNameIds")
    end
    anatName = anatName(chanIdx);
    anatName(cellfun(@isempty,anatName)) = {'noName'};
    
%     sensorymotorChan = contains(anatName,'central');
%     whiteChan = contains(anatName,'white');
%     ifgChan = contains(anatName,'opercularis');
%     frontalChan = contains(anatName,'front');
%     temporalChan = contains(anatName,'temporal');
    anatChanId = ones(1,length(anatName));
    disp(['Number of anatomical channels : ' num2str(length(anatChanId))]);
    
    % prodDelayElecs is loaded from the mat file in 'forKumar'
    % Come up with automated way to parse necessary channels
    
    channelName = {Subject(iSubject).ChannelInfo.Name};
    channelName = channelName(chanIdx);
    channelName(cellfun(@isempty,channelName)) = {'noName'};
    
%     %muscleChannelsName = channelName(muscleChannels);
%     [~,delayChan] = intersect(channelName,elecNameProductionClean); 
   
%     chan2select = intersect(find(anatChanId),delayChan);
    chan2select = find(anatChanId);
    %chan2select = muscleChannels;
     if(isempty(chan2select))
        disp('No channels found; Iterating next subject');
        continue;
        % Forces the iteration for next subject;
    end
    chanNameSubject = channelName(chan2select);

%%
    RespTime=[];
    for iTrials=1:length(Trials)
        if Trials(iTrials).Noisy==1
            noiseIdx(counterN+1)=iTrials;
            counterN=counterN+1;
        end
        if Trials(iTrials).NoResponse==1
            noResponseIdx(counterNR+1)=iTrials;
            counterNR=counterNR+1;
        end

        if ~isempty(Trials(iTrials).ResponseStart)
           RespTime(iTrials)=(Trials(iTrials).ResponseStart-Trials(iTrials).Go)./30000;
       else
           RespTime(iTrials)=0;
       end
    end
    negResponseIdx=find(RespTime<0);
    
       

    condIdx=ones(length(Subject(iSubject).Trials),1);
    
    baseEpoch=Task.Base.Epoch;
    baseTimeRange=Task.Base.Time; 
    baseTimeExtract = [baseTimeRange(1)-timePad baseTimeRange(2)+timePad];
    Trials=Subject(iSubject).Trials(setdiff(1:length(Subject(iSubject).Trials),noiseIdx));
    
    ieegBase=trialIEEG(Trials,chanIdx,baseEpoch,baseTimeExtract.*1000);
    ieegBase = permute(ieegBase,[2,1,3]);
    fs = Subject(iSubject).Experiment.processing.ieeg.sample_rate;   
    ieegBaseStruct = ieegStructClass(ieegBase, fs, baseTimeExtract, [1 fs/2], baseEpoch);
    clear ieegBase;
    ieegBaseCAR=extractCar(ieegBaseStruct,badChan);
    clear ieegBaseStruct;
%     [NumTrials,goodtrials] = remove_bad_trials(ieegBaseCAR.data,10);
%     goodTrialsCommon = extractCommonTrials(goodtrials);
    ieegBaseHG = extractHiGamma(ieegBaseCAR,fDown,baseTimeRange);
    normFactorBase = extractHGnormFactor(ieegBaseHG);
    for iCond=1:length(Task.Conds)
        trials2Select = setdiff(find(condIdx==iCond),cat(2,noiseIdx,noResponseIdx,negResponseIdx));
        if(isempty(trials2Select))
            continue;
        end
        % Iterting through field: 'Auditory', 'Delay', 'Response'
        for iField=1:length(Task.Conds(iCond).Field)
            Trials = Subject(iSubject).Trials(trials2Select);
            Epoch=Task.Conds(iCond).Field(iField).Epoch;
            delTime = [];
            respTime = [];
            for iTrials = 1:length(Trials)
                delTime(iTrials) = (Trials(iTrials).Go-Trials(iTrials).Auditory)./30000;
                respTime(iTrials)=(Trials(iTrials).ResponseStart-Trials(iTrials).Go)./30000;
            end
            [respTimeSort,respSortId] = sort(respTime);
            [delTimeSort,delSortId] = sort(delTime);
            
            %kluge
            if(iField==1)
                timeSort = delTimeSort;
                sortId = delSortId;
            else
                timeSort = respTimeSort;
                sortId = respSortId;
            end

            fieldTimeRange=Task.Conds(iCond).Field(iField).Time;          
            fieldTimeExtract = [fieldTimeRange(1)-timePad fieldTimeRange(2)+timePad];            
            ieegField=trialIEEG(Trials,chanIdx,Epoch,fieldTimeExtract.*1000);
            ieegField = permute(ieegField,[2,1,3]);
            ieegFieldStruct = ieegStructClass(ieegField, fs, fieldTimeExtract, [1 fs/2], Epoch);
            clear ieegField
            % Common average referencing
            ieegFieldCAR = extractCar(ieegFieldStruct,badChan);
            clear ieegFieldStruct;
            ieegFieldHG = extractHiGamma(ieegFieldCAR,fDown,fieldTimeRange,normFactorBase,1);
            clear ieegFieldCar
            
            timeVect = linspace(ieegFieldHG.tw(1),ieegFieldHG.tw(2),size(ieegFieldCAR.data,3));
            totChanBlock=ceil(length(chan2select)./40);
        iChan2=0;
        for iB=0:totChanBlock-1;
            FigS=figure('Position', get(0, 'Screensize'));
            for iChan=1:min(40,length(chan2select)-iChan2);
                subplot(4,10,iChan);
                iChan2=iChan+iB*40;
                iChan2  
                
                imagesc(timeVect,[],squeeze(ieegFieldHG.data(iChan2,sortId,:)));
                caxis([0 1.5]);
                hold on;
                scatter(timeSort,1:length(timeSort),10,'k','filled');
                set(gca,'YTick',[]);
                %set(gca,'YTickLabels',round(waveSpecField.fscale(1:10:end)));
                set(gca,'FontSize',10);
                xlim(fieldTimeRange);
                title(chanNameSubject(iChan2))
                
            end;
            F=getframe(FigS);
            if ~exist([DUKEDIR '/Figs/response_time_hg/' Subject(iSubject).Name],'dir')
                mkdir([DUKEDIR '/Figs/response_time_hg/' Subject(iSubject).Name])
            end
            imwrite(F.cdata,[DUKEDIR '/Figs/response_time_hg/' Subject(iSubject).Name '/' Subject(iSubject).Name '_PhonemeSequence_' Task.Conds(iCond).Field(iField).Name '_HG_zscore_0_1-5_' num2str(iB+1) '.png'],'png');    
          %  imwrite(F.cdata,[DUKEDIR '/Figs/' Subject(SN).Name '/' Subject(SN).Name '_PhonemeSequence_Auditory_SpecGrams_0.7to1.4C_' num2str(iF+1) '.png'],'png');
            close
        end  
        
        
        end
        
    end
end
