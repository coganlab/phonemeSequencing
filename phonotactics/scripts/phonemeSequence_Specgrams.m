%duke;
global BOX_DIR
global RECONDIR
global TASK_DIR
global experiment
global DUKEDIR
BOX_DIR='C:\Users\gcoga\Box';
RECONDIR='C:\Users\gcoga\Box\ECoG_Recon';

Task = [];
Task.Name='Phoneme_Sequencing';

TASK_DIR=([BOX_DIR '/CoganLab/D_Data/' Task.Name]);
DUKEDIR=TASK_DIR;


% Populate Subject file
Subject = popTaskSubjectData(Task.Name);




DUKEDIR=['C:\Users\gcoga\Box\CoganLab\D_Data'];
DUKEDIR = [DUKEDIR '/Phoneme_Sequencing'];

SNList=1:length(Subject);
counterChan=0;
for iSN=1:length(SNList)
AnalParams=[];
CondParams=[];
SN = SNList(iSN);
Subject(SN).Experiment = loadExperiment(Subject(SN).Name);
%Subject(SN).Trials = dbTrials(Subject(SN).Name,Subject(SN).Day,'Speech_OvertMimeMove');
%Subject(SN).Trials = dbTrials(Subject(SN).Name,Subject(SN).Day,'Speech_CovertOvert'); % Test line
%Trials = Subject(SN).Trials;
experiment = Subject(SN).Experiment;
load([DUKEDIR '/' Subject(SN).Name '/' Subject(SN).Date '/mat/trialInfo.mat'])
load([DUKEDIR '/' Subject(SN).Name '/' Subject(SN).Date '/mat/Trials.mat']);
CondParams.Conds = [1:1];

for iTrials=1:length(Trials);
    Trials(iTrials).StartCode=1;Trials(iTrials).AuditoryCode=26;Trials(iTrials).GoCode=51;
    Trials(iTrials).Noisy=0;Trials(iTrials).NoResponse=0;
end

Subject(SN).Trials=Trials;

AnalParams.Channel=setdiff(Subject(SN).ChannelNums,Subject(SN).badChannels);
CondParams.Conds=[1:1];
%NumTrials = SelectChannels(Subject(SN).Trials, CondParams, AnalParams);
SelectedChannels=AnalParams.Channel; % really loose: accounts for practice trial confound
AnalParams.ReferenceChannels = SelectedChannels;
AnalParams.Channel = SelectedChannels;
AnalParams.TrialPooling = 1; %1;  %1; % used to be 1
AnalParams.dn=0.05;
AnalParams.Tapers = [.5,10];
AnalParams.fk = 200;
AnalParams.Reference = 'Grand average';% 'IndRef'; %'Grand average', 'Grand average induced'% induced' 'Single-ended','IndRef';%AnalParams.RefChans=subjRefChansInd(Subject(SN).Name);
AnalParams.ArtifactThreshold = 12; %8 %12;
srate=experiment.recording.sample_rate;
srate2=srate/4;
if srate<2048
    AnalParams.pad=2;
else
   AnalParams.pad=1;
end   
  
    CondParams.Conds=[1:1];
    CondParams.Field = 'Auditory';
    CondParams.bn = [-500,1000];
    for iCode = 1:length(CondParams.Conds)
        
        if isfield(CondParams,'Conds2')
            CondParams.Conds = CondParams.Conds2(iCode);
        else
            CondParams.Conds = iCode;
        end
        tic
        [Auditory_Spec{iCode}, Auditory_Data, Auditory_Trials{iCode}] = subjSpectrum(Subject(SN), CondParams, AnalParams);
        toc
        display(['Cond = ' num2str(iCode)])
    end

    CondParams.Conds=[1:1];
    CondParams.Field = 'Auditory';
    CondParams.bn = [1000,2000];
    for iCode = 1:length(CondParams.Conds)
        
        if isfield(CondParams,'Conds2')
            CondParams.Conds = CondParams.Conds2(iCode);
        else
            CondParams.Conds = iCode;
        end
        tic
        [Maint_Spec{iCode}, Maint_Data, Maint_Trials{iCode}] = subjSpectrum(Subject(SN), CondParams, AnalParams);
        toc
        display(['Cond = ' num2str(iCode)])
    end
    
      CondParams.Conds=[1:1];
CondParams.Field = 'Start';
CondParams.bn = [-500,500];
for iCode = 1:length(Auditory_Spec)
    if isfield(CondParams,'Conds2')
        CondParams.Conds = CondParams.Conds2(iCode);
    else
        CondParams.Conds = iCode;
    end
    tic
    [Start_Spec{iCode}, Start_Data, Start_Trials{iCode}] = subjSpectrum(Subject(SN), CondParams, AnalParams);
    toc
    display(['Cond = ' num2str(iCode)])
end  
    
    CondParams.Conds=[1:1];
CondParams.Field = 'Go';
CondParams.bn = [-500,1500];
for iCode = 1:length(Auditory_Spec)
    if isfield(CondParams,'Conds2')
        CondParams.Conds = CondParams.Conds2(iCode);
    else
        CondParams.Conds = iCode;
    end
    tic
    [Motor_Spec{iCode}, Motor_Data, Motor_Trials{iCode}] = subjSpectrum(Subject(SN), CondParams, AnalParams);
    toc
    display(['Cond = ' num2str(iCode)])
end
% %

   base=0;
    %base = zeros(1,size(Auditory_Spec{iCode}{iCh},2));
    for iCh = 1:length(Auditory_Spec{iCode})
        base=0;
        for iCode = 1:length(Auditory_Spec)
            %base = base + sq(Auditory_Spec{iCode}{iCh}(5,:)); % standard
            %   base= base+mean(sq(Auditory_Spec{iCode}{iCh}(1:10,:)),1); % used to be 1:9
            base= base+mean(sq(Start_Spec{iCode}{iCh}(1:10,:)),1); % used to be 1:9
            
            %base2(iCode,:)=std(sq(Auditory_Spec{iCode}{iCh}(1:6,:)),1); % std across time bins?
            
        end
        base = base./length(Auditory_Spec);
        for iCode = 1:length(Auditory_Spec)
            Auditory_nSpec(iCode,iCh,:,:) = Auditory_Spec{iCode}{iCh}(:,:)./base(ones(1,size(Auditory_Spec{iCode}{iCh},1)),:);
            Maint_nSpec(iCode,iCh,:,:) = Maint_Spec{iCode}{iCh}(:,:)./base(ones(1,size(Maint_Spec{iCode}{iCh},1)),:);
            Motor_nSpec(iCode,iCh,:,:) = Motor_Spec{iCode}{iCh}(:,:)./base(ones(1,size(Motor_Spec{iCode}{iCh},1)),:);
            Start_nSpec(iCode,iCh,:,:) = Start_Spec{iCode}{iCh}(:,:)./base(ones(1,size(Start_Spec{iCode}{iCh},1)),:);
            %  Motor_nSpec(iCode,iCh,:,:) = Motor_Spec{iCode}{iCh}(:,:)./base(ones(1,size(Motor_Spec{iCode}{iCh},1)),:);
            %  Start_nSpec(iCode,iCh,:,:) = Start_Spec{iCode}{iCh}(:,:)./base(ones(1,size(Start_Spec{iCode}{iCh},1)),:);
            % Del_nSpec(iCode,iCh,:,:) = Del_Spec{iCode}{iCh}(:,:)./base(ones(1,size(Del_Spec{iCode}{iCh},1)),:);
            %Auditory_nSpec(iCode,iCh,:,:)=(sq(Auditory_Spec{iCode}{iCh}(:,:))-base2)./(repmat(base,80,1)); % SD LINE
            
        end
        
    end

totChanBlock=ceil(length(AnalParams.Channel)./60);
iChan2=0;
for iF=0:totChanBlock-1;
    FigS=figure('Position', get(0, 'Screensize'));
    for iChan=1:min(60,length(AnalParams.Channel)-iChan2);
        subplot(6,10,iChan);
        iChan2=iChan+iF*60;
        tvimage(sq((Auditory_nSpec([1],iChan2,:,1:200))),'XRange',[-0.5,1]);
        title(experiment.channels(AnalParams.Channel(iChan2)).name);
        caxis([0.7 1.2]);
      %  caxis([0.7 1.4]);

    end;
    F=getframe(FigS);
    if ~exist([DUKEDIR '/Figs/' Subject(SN).Name],'dir')
        mkdir([DUKEDIR '/Figs/' Subject(SN).Name])
    end
    imwrite(F.cdata,[DUKEDIR '/Figs/' Subject(SN).Name '/' Subject(SN).Name '_PhonemeSequence_Auditory_SpecGrams_' num2str(iF+1) '.png'],'png');    
  %  imwrite(F.cdata,[DUKEDIR '/Figs/' Subject(SN).Name '/' Subject(SN).Name '_PhonemeSequence_Auditory_SpecGrams_0.7to1.4C_' num2str(iF+1) '.png'],'png');
    close
end  


totChanBlock=ceil(length(AnalParams.Channel)./60);
iChan2=0;
for iF=0:totChanBlock-1;
    FigS=figure('Position', get(0, 'Screensize'));
    for iChan=1:min(60,length(AnalParams.Channel)-iChan2);
        subplot(6,10,iChan);
        iChan2=iChan+iF*60;
        tvimage(sq((Maint_nSpec([1],iChan2,:,1:200))),'XRange',[0 1]);
        title(experiment.channels(AnalParams.Channel(iChan2)).name);
        caxis([0.7 1.2]);
     %   caxis([0.7 1.4]);

    end;
    F=getframe(FigS);
    if ~exist([DUKEDIR '/Figs/' Subject(SN).Name],'dir')
        mkdir([DUKEDIR '/Figs/' Subject(SN).Name])
    end
    imwrite(F.cdata,[DUKEDIR '/Figs/' Subject(SN).Name '/' Subject(SN).Name '_PhonemeSequence_Maintenance_SpecGrams_' num2str(iF+1) '.png'],'png');
  %  imwrite(F.cdata,[DUKEDIR '/Figs/' Subject(SN).Name '/' Subject(SN).Name '_PhonemeSequence_Maintenance_SpecGrams_0.7to1.4C_' num2str(iF+1) '.png'],'png');
    close
end  

totChanBlock=ceil(length(AnalParams.Channel)./60);
iChan2=0;
for iF=0:totChanBlock-1;
    FigS=figure('Position', get(0, 'Screensize'));
    for iChan=1:min(60,length(AnalParams.Channel)-iChan2);
        subplot(6,10,iChan);
        iChan2=iChan+iF*60;
        tvimage(sq((Motor_nSpec([1],iChan2,:,1:200))),'XRange',[-0.5 1.5]);
        title(experiment.channels(AnalParams.Channel(iChan2)).name);
        caxis([0.7 1.2]);
     %   caxis([0.7 1.4]);

    end;
    F=getframe(FigS);
    if ~exist([DUKEDIR '/Figs/' Subject(SN).Name],'dir')
        mkdir([DUKEDIR '/Figs/' Subject(SN).Name])
    end
    imwrite(F.cdata,[DUKEDIR '/Figs/' Subject(SN).Name '/' Subject(SN).Name '_PhonemeSequence_Motor_SpecGrams_' num2str(iF+1) '.png'],'png');
  %  imwrite(F.cdata,[DUKEDIR '/Figs/' Subject(SN).Name '/' Subject(SN).Name '_PhonemeSequence_Maintenance_SpecGrams_0.7to1.4C_' num2str(iF+1) '.png'],'png');
    close
end  

totChanBlock=ceil(length(AnalParams.Channel)./60);
iChan2=0;
for iF=0:totChanBlock-1;
    FigS=figure('Position', get(0, 'Screensize'));
    for iChan=1:min(60,length(AnalParams.Channel)-iChan2);
        subplot(6,10,iChan);
        iChan2=iChan+iF*60;
        tvimage(sq((Start_nSpec([1],iChan2,:,1:200))),'XRange',[-0.5 0.5]);
        title(experiment.channels(AnalParams.Channel(iChan2)).name);
        caxis([0.7 1.2]);
     %   caxis([0.7 1.4]);

    end;
    F=getframe(FigS);
    if ~exist([DUKEDIR '/Figs/' Subject(SN).Name],'dir')
        mkdir([DUKEDIR '/Figs/' Subject(SN).Name])
    end
    imwrite(F.cdata,[DUKEDIR '/Figs/' Subject(SN).Name '/' Subject(SN).Name '_PhonemeSequence_Start_SpecGrams_' num2str(iF+1) '.png'],'png');
  %  imwrite(F.cdata,[DUKEDIR '/Figs/' Subject(SN).Name '/' Subject(SN).Name '_PhonemeSequence_Maintenance_SpecGrams_0.7to1.4C_' num2str(iF+1) '.png'],'png');
    close
end  

totChanBlock=ceil(length(AnalParams.Channel)./60);
iChan2=0;
for iF=0:totChanBlock-1;
    FigS=figure('Position', get(0, 'Screensize'));
    for iChan=1:min(60,length(AnalParams.Channel)-iChan2);
        subplot(6,10,iChan);
        iChan2=iChan+iF*60;
        tvimage(sq((Auditory_nSpec([1],iChan2,:,1:200))),'XRange',[-0.5,1]);
        title(experiment.channels(AnalParams.Channel(iChan2)).name);
      %  caxis([0.7 1.2]);
        caxis([0.7 1.4]);

    end;
    F=getframe(FigS);
    if ~exist([DUKEDIR '/Figs/' Subject(SN).Name],'dir')
        mkdir([DUKEDIR '/Figs/' Subject(SN).Name])
    end
  %  imwrite(F.cdata,[DUKEDIR '/Figs/' Subject(SN).Name '/' Subject(SN).Name '_PhonemeSequence_Auditory_SpecGrams_' num2str(iF+1) '.png'],'png');    
    imwrite(F.cdata,[DUKEDIR '/Figs/' Subject(SN).Name '/' Subject(SN).Name '_PhonemeSequence_Auditory_SpecGrams_0.7to1.4C_' num2str(iF+1) '.png'],'png');
    close
end  


totChanBlock=ceil(length(AnalParams.Channel)./60);
iChan2=0;
for iF=0:totChanBlock-1;
    FigS=figure('Position', get(0, 'Screensize'));
    for iChan=1:min(60,length(AnalParams.Channel)-iChan2);
        subplot(6,10,iChan);
        iChan2=iChan+iF*60;
        tvimage(sq((Maint_nSpec([1],iChan2,:,1:200))),'XRange',[0 1]);
        title(experiment.channels(AnalParams.Channel(iChan2)).name);
      %  caxis([0.7 1.2]);
        caxis([0.7 1.4]);

    end;
    F=getframe(FigS);
    if ~exist([DUKEDIR '/Figs/' Subject(SN).Name],'dir')
        mkdir([DUKEDIR '/Figs/' Subject(SN).Name])
    end
  % imwrite(F.cdata,[DUKEDIR '/Figs/' Subject(SN).Name '/' Subject(SN).Name '_PhonemeSequence_Maintenance_SpecGrams_' num2str(iF+1) '.png'],'png');
    imwrite(F.cdata,[DUKEDIR '/Figs/' Subject(SN).Name '/' Subject(SN).Name '_PhonemeSequence_Maintenance_SpecGrams_0.7to1.4C_' num2str(iF+1) '.png'],'png');
    close
end 

totChanBlock=ceil(length(AnalParams.Channel)./60);
iChan2=0;
for iF=0:totChanBlock-1;
    FigS=figure('Position', get(0, 'Screensize'));
    for iChan=1:min(60,length(AnalParams.Channel)-iChan2);
        subplot(6,10,iChan);
        iChan2=iChan+iF*60;
        tvimage(sq((Motor_nSpec([1],iChan2,:,1:200))),'XRange',[-0.5 1.5]);
        title(experiment.channels(AnalParams.Channel(iChan2)).name);
      %  caxis([0.7 1.2]);
        caxis([0.7 1.4]);

    end;
    F=getframe(FigS);
    if ~exist([DUKEDIR '/Figs/' Subject(SN).Name],'dir')
        mkdir([DUKEDIR '/Figs/' Subject(SN).Name])
    end
  % imwrite(F.cdata,[DUKEDIR '/Figs/' Subject(SN).Name '/' Subject(SN).Name '_PhonemeSequence_Maintenance_SpecGrams_' num2str(iF+1) '.png'],'png');
    imwrite(F.cdata,[DUKEDIR '/Figs/' Subject(SN).Name '/' Subject(SN).Name '_PhonemeSequence_Motor_SpecGrams_0.7to1.4C_' num2str(iF+1) '.png'],'png');
    close
end

totChanBlock=ceil(length(AnalParams.Channel)./60);
iChan2=0;
for iF=0:totChanBlock-1;
    FigS=figure('Position', get(0, 'Screensize'));
    for iChan=1:min(60,length(AnalParams.Channel)-iChan2);
        subplot(6,10,iChan);
        iChan2=iChan+iF*60;
        tvimage(sq((Start_nSpec([1],iChan2,:,1:200))),'XRange',[-0.5 0.5]);
        title(experiment.channels(AnalParams.Channel(iChan2)).name);
      %  caxis([0.7 1.2]);
        caxis([0.7 1.4]);

    end;
    F=getframe(FigS);
    if ~exist([DUKEDIR '/Figs/' Subject(SN).Name],'dir')
        mkdir([DUKEDIR '/Figs/' Subject(SN).Name])
    end
  % imwrite(F.cdata,[DUKEDIR '/Figs/' Subject(SN).Name '/' Subject(SN).Name '_PhonemeSequence_Maintenance_SpecGrams_' num2str(iF+1) '.png'],'png');
    imwrite(F.cdata,[DUKEDIR '/Figs/' Subject(SN).Name '/' Subject(SN).Name '_PhonemeSequence_Start_SpecGrams_0.7to1.4C_' num2str(iF+1) '.png'],'png');
    close
end  

clear Auditory_nSpec;
clear Auditory_Spec;
clear Maint_nSpec;
clear Maint_Spec;
clear Start_Spec;
clear Motor_Spec;
clear Motor_nSpec;
clear Start_Spec;
clear Start_nSpec;
end

    
