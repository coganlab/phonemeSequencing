%% Phoneme decoder
colvals = lines(5)
scrsize = get(0, 'Screensize');
figure('Position', [scrsize(1) scrsize(2) scrsize(3) scrsize(4)/2]); 
for iPhon = 1:3
    visTimeGenAcc1D(decodeTimeStruct1D{iPhon},pVal2Cutoff=0.01,...
        axisLabel = 'Auditory',clowLimit = 0,timePad = 0.1,...
        maxVal = (iPhon-1)*0.015, chanceVal = 0.1111,clabel = 'Accuracy',colval=colvals(iPhon,:))
end

xline(0, '-','Auditory onset','LineWidth',2, 'Color','k');
xline(2, ':','Stitch','LineWidth',2, 'Color','k');
xline(2.5, '-','Go cue onset','LineWidth',2, 'Color','k');
xline(3, ':','Stitch','LineWidth',2, 'Color','k');
xline(4, '-','Response onset','LineWidth',2, 'Color','k');

%% Phoneme decoder - consonant
colvals = lines(5)
scrsize = get(0, 'Screensize');
figure('Position', [scrsize(1) scrsize(2) scrsize(3) scrsize(4)/2]); 
% for iPhon = 1:3
    visTimeGenAcc1D(decodeTimeStruct1D_cvc{1},pVal2Cutoff=0.01,...
        axisLabel = 'Auditory',clowLimit = 0,timePad = 0.1,...
        maxVal = 0, chanceVal = 0.2,clabel = 'Accuracy',colval=colvals(1,:))

    visTimeGenAcc1D(decodeTimeStruct1D_vcv{2},pVal2Cutoff=0.01,...
        axisLabel = 'Auditory',clowLimit = 0,timePad = 0.1,...
        maxVal = 0.015, chanceVal = 0.2,clabel = 'Accuracy',colval=colvals(2,:))

    visTimeGenAcc1D(decodeTimeStruct1D_cvc{3},pVal2Cutoff=0.01,...
        axisLabel = 'Auditory',clowLimit = 0,timePad = 0.1,...
        maxVal = 0.03, chanceVal = 0.2,clabel = 'Accuracy',colval=colvals(3,:))
% end

xline(0, '-','Auditory onset','LineWidth',2, 'Color','k');
xline(2, ':','Stitch','LineWidth',2, 'Color','k');
xline(2.5, '-','Go cue onset','LineWidth',2, 'Color','k');
xline(3, ':','Stitch','LineWidth',2, 'Color','k');
xline(4, '-','Response onset','LineWidth',2, 'Color','k');

%% Phoneme decoder - vowel
colvals = lines(5)
scrsize = get(0, 'Screensize');
figure('Position', [scrsize(1) scrsize(2) scrsize(3) scrsize(4)/2]); 
% for iPhon = 1:3
    visTimeGenAcc1D(decodeTimeStruct1D_vcv{1},pVal2Cutoff=0.01,...
        axisLabel = 'Auditory',clowLimit = 0,timePad = 0.1,...
        maxVal = 0, chanceVal = 0.25,clabel = 'Accuracy',colval=colvals(1,:))

    visTimeGenAcc1D(decodeTimeStruct1D_cvc{2},pVal2Cutoff=0.01,...
        axisLabel = 'Auditory',clowLimit = 0,timePad = 0.1,...
        maxVal = 0.015, chanceVal = 0.25,clabel = 'Accuracy',colval=colvals(2,:))

    visTimeGenAcc1D(decodeTimeStruct1D_vcv{3},pVal2Cutoff=0.01,...
        axisLabel = 'Auditory',clowLimit = 0,timePad = 0.1,...
        maxVal = 0.03, chanceVal = 0.25,clabel = 'Accuracy',colval=colvals(3,:))
% end

xline(0, '-','Auditory onset','LineWidth',2, 'Color','k');
xline(2, ':','Stitch','LineWidth',2, 'Color','k');
xline(2.5, '-','Go cue onset','LineWidth',2, 'Color','k');
xline(3, ':','Stitch','LineWidth',2, 'Color','k');
xline(4, '-','Response onset','LineWidth',2, 'Color','k');
%% Articulator decoder

figure; 
for iPhon = 1:3
    visTimeGenAcc1D(decodeTimeStruct1D{iPhon},pVal2Cutoff=0.01,...
        axisLabel = 'Auditory ',clowLimit = 0,timePad = 0.1,...
        maxVal = (iPhon-1)*0.015, chanceVal = 0.25,clabel = 'Accuracy',colval=colvals(iPhon,:))
end
xline(0, ':','Auditory onset','LineWidth',2, 'Color','k');
xline(2, '-','Stitch','LineWidth',2, 'Color','k');
xline(2.5, ':','Go cue onset','LineWidth',2, 'Color','k');
xline(3, '-','Stitch','LineWidth',2, 'Color','k');
xline(4, ':','Response onset','LineWidth',2, 'Color','k');

%% Syllable decoder

figure; 
    visTimeGenAcc1D(decodeTimeStruct1D,pVal2Cutoff=0.01,...
        axisLabel = 'Auditory',clowLimit = 0,timePad = 0.1,...
        maxVal = 0.4, chanceVal = 0.5,clabel = 'Accuracy',colval=colvals(1,:))
    xline(0, ':','Auditory onset','LineWidth',2, 'Color','k');
xline(2, '-','Stitch','LineWidth',2, 'Color','k');
xline(2.5, ':','Go cue onset','LineWidth',2, 'Color','k');
xline(3, '-','Stitch','LineWidth',2, 'Color','k');
xline(4, ':','Response onset','LineWidth',2, 'Color','k');

%% Roi analysis

cue = 'Auditory';
decoder = 'phoneme';
chanceVal = 0.1111;
load(['pooledSubject_Phoneme_Sequencing_temporal_' cue '__productionelecs_data_Start_' decoder '_decoded.mat'])
figure;
 visTimeGenAcc1D(decodeTimeStruct1D{1},pVal2Cutoff=0.01,...
        axisLabel = cue,clowLimit = 0,timePad = 0.1,...
        maxVal = 0, chanceVal = chanceVal,clabel = 'Accuracy',colval=colval(1,:));
 
 load(['pooledSubject_Phoneme_Sequencing_temporal_feedback_' cue '_productionelecs_data_Start_' decoder '_decoded.mat'])
 visTimeGenAcc1D(decodeTimeStruct1D{1},pVal2Cutoff=0.01,...
        axisLabel = cue,clowLimit = 0,timePad = 0.1,...
        maxVal = 0.02, chanceVal = chanceVal,clabel = 'Accuracy',colval=colval(2,:));
 
 load(['pooledSubject_Phoneme_Sequencing_parietal_' cue '__productionelecs_data_Start_' decoder '_decoded.mat'])
 visTimeGenAcc1D(decodeTimeStruct1D{1},pVal2Cutoff=0.01,...
        axisLabel = cue,clowLimit = 0,timePad = 0.1,...
        maxVal = 0.04, chanceVal = chanceVal,clabel = 'Accuracy',colval=colval(3,:));

 load(['pooledSubject_Phoneme_Sequencing_frontal_' cue '__productionelecs_data_Start_' decoder '_decoded.mat'])
 visTimeGenAcc1D(decodeTimeStruct1D{1},pVal2Cutoff=0.01,...
        axisLabel = cue,clowLimit = 0,timePad = 0.1,...
        maxVal = 0.06, chanceVal = chanceVal,clabel = 'Accuracy',colval=colval(4,:));

 load(['pooledSubject_Phoneme_Sequencing_sensorimotor_' cue '__productionelecs_data_Start_'  decoder '_decoded.mat'])
 visTimeGenAcc1D(decodeTimeStruct1D{1},pVal2Cutoff=0.01,...
        axisLabel = cue,clowLimit = 0,timePad = 0.1,...
        maxVal = 0.08, chanceVal = chanceVal,clabel = 'Accuracy',colval=colval(5,:));
 %% Stitched data analysis

 visTimeGenAcc2D(decodeTimeStruct2D{1},pVal2Cutoff=0.05, axisLabel = 'Auditory onset',...
     clowLimit = 1,timePad = 0.1,chanceVal = 0.1111,clabel='Accuracy/Chance')
 hold on;
 yline(0, '-','Auditory onset','LineWidth',2, 'Color','w');
 xline(0, '-','Auditory onset','LineWidth',2, 'Color','w');
 yline(2, ':','Stitch','LineWidth',2, 'Color','y');
 xline(2, ':','Stitch','LineWidth',2, 'Color','y');
 yline(2.5, '-','Go cue onset','LineWidth',2, 'Color','w');
 xline(2.5, '-','Go cue onset','LineWidth',2, 'Color','w');
 yline(3, ':','Stitch','LineWidth',2, 'Color','y');
 xline(3, ':','Stitch','LineWidth',2, 'Color','y');
 yline(4, '-','Response onset','LineWidth',2, 'Color','w');
 xline(4, '-','Response onset','LineWidth',2, 'Color','w');
 %% visualize HG
scrsize = get(0, 'Screensize');

%figure('Position', [scrsize(1) scrsize(2) scrsize(3) scrsize(4)/2]);

visTimePlot(timeEpoch,squeeze(mean(ieegStructPooled.data,2)),colval = colvalfeedback(6,:))

xline(0, '-','Auditory onset','LineWidth',2, 'Color','k');
xline(2, ':','Stitch','LineWidth',2, 'Color','k');
xline(2.5, '-','Go cue onset','LineWidth',2, 'Color','k');
xline(3, ':','Stitch','LineWidth',2, 'Color','k');
xline(4, '-','Response onset','LineWidth',2, 'Color','k');