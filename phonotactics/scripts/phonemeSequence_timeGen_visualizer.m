%% Phoneme decoder
colvals = lines(5)
scrsize = get(0, 'Screensize');
figure('Position', [scrsize(1) scrsize(2) scrsize(3) scrsize(4)/2]); 
%figure('Position', [scrsize(1) scrsize(2) scrsize(3) scrsize(4)/2]); 
timeEpoch =[-0.5000    1.5000;   -0.7500    0.5000;   -1.0000    1.5000];
   
for iPhon = 1:3
    if(iPhon==1)
     ax = visTimeGenAcc1DCluster_v2(decodeTimeStruct1D(iPhon,:),decodeTimeStruct1Dshuffle(iPhon,:),timeEpoch,pVal2Cutoff=1e-3,...
        axisLabel = 'Auditory',clowLimit = 0,timePad = 0.25,...
         maxVal = (iPhon-1)*0.015, chanceVal = 0.1111,clabel = 'Accuracy',colval=colvals(iPhon,:));
    else
        visTimeGenAcc1DCluster_v2(decodeTimeStruct1D(iPhon,:),decodeTimeStruct1Dshuffle(iPhon,:),timeEpoch,pVal2Cutoff=1e-3,...
        axisLabel = 'Auditory',clowLimit = 0,timePad =0.25,...
         maxVal = (iPhon-1)*0.015, chanceVal = 0.1111,clabel = 'Accuracy',colval=colvals(iPhon,:),tileaxis = ax)
    end

end

% xline(0, '-','Auditory onset','LineWidth',2, 'Color','k');
% xline(1.5, ':','Earliest Go time - Stitch','LineWidth',2, 'Color','k');
% %xline(2, ':','Stitch','LineWidth',2, 'Color','k');
% xline(2.25, '-','Go cue onset','LineWidth',2, 'Color','k');
% xline(2.75, ':','Stitch','LineWidth',2, 'Color','k');
% xline(3.75, '-','Response onset','LineWidth',2, 'Color','k');
% set(gca,'FontSize',15);

%% Phoneme decoder - consonant
colvals = lines(5)
scrsize = get(0, 'Screensize');
figure('Position', [scrsize(1) scrsize(2) scrsize(3) scrsize(4)/2]); 
% for iPhon = 1:3
    ax = visTimeGenAcc1DCluster_v2(decodeTimeStruct1D_cvc(1,:),decodeTimeStruct1D_cvc_shuffle(1,:),timeEpoch,pVal2Cutoff=5e-5,...
        axisLabel = 'Auditory',clowLimit = 0,timePad = 0.25,...
        maxVal = 0, chanceVal = 0.2,clabel = 'Accuracy',colval=colvals(1,:))

    visTimeGenAcc1DCluster_v2(decodeTimeStruct1D_vcv(2,:),decodeTimeStruct1D_vcv_shuffle(2,:),timeEpoch,pVal2Cutoff=5e-5,...
        axisLabel = 'Auditory',clowLimit = 0,timePad = 0.25,...
        maxVal = 0.015, chanceVal = 0.2,clabel = 'Accuracy',colval=colvals(2,:),tileaxis = ax)

    visTimeGenAcc1DCluster_v2(decodeTimeStruct1D_cvc(3,:),decodeTimeStruct1D_cvc_shuffle(3,:),timeEpoch,pVal2Cutoff=5e-5,...
        axisLabel = 'Auditory',clowLimit = 0,timePad = 0.25,...
        maxVal = 0.03, chanceVal = 0.2,clabel = 'Accuracy',colval=colvals(3,:),tileaxis = ax)
% end

% xline(0, '-','Auditory onset','LineWidth',2, 'Color','k');
% xline(1.5, ':','Earliest Go time - Stitch','LineWidth',2, 'Color','k');
% %xline(2, ':','Stitch','LineWidth',2, 'Color','k');
% xline(2.25, '-','Go cue onset','LineWidth',2, 'Color','k');
% xline(2.75, ':','Stitch','LineWidth',2, 'Color','k');
% xline(3.75, '-','Response onset','LineWidth',2, 'Color','k');
% set(gca,'FontSize',15);

%% Phoneme decoder - vowel
colvals = lines(5)
scrsize = get(0, 'Screensize');
figure('Position', [scrsize(1) scrsize(2) scrsize(3) scrsize(4)/2]); 
% for iPhon = 1:3
    ax = visTimeGenAcc1DCluster_v2(decodeTimeStruct1D_vcv(1,:),decodeTimeStruct1D_vcv_shuffle(1,:),timeEpoch,pVal2Cutoff=5e-4,...
        axisLabel = 'Auditory',clowLimit = 0,timePad = 0.25,...
        maxVal = 0, chanceVal = 0.25,clabel = 'Accuracy',colval=colvals(1,:))

    visTimeGenAcc1DCluster_v2(decodeTimeStruct1D_cvc(2,:),decodeTimeStruct1D_cvc_shuffle(2,:),timeEpoch,pVal2Cutoff=5e-4,...
        axisLabel = 'Auditory',clowLimit = 0,timePad = 0.25,...
        maxVal = 0.015, chanceVal = 0.25,clabel = 'Accuracy',colval=colvals(2,:),tileaxis = ax)

    visTimeGenAcc1DCluster_v2(decodeTimeStruct1D_cvc(3,:),decodeTimeStruct1D_cvc_shuffle(3,:),timeEpoch,pVal2Cutoff=5e-4,...
        axisLabel = 'Auditory',clowLimit = 0,timePad = 0.25,...
        maxVal = 0.03, chanceVal = 0.25,clabel = 'Accuracy',colval=colvals(3,:),tileaxis = ax)
% end

% xline(0, '-','Auditory onset','LineWidth',2, 'Color','k');
% xline(1.5, ':','Earliest Go time - Stitch','LineWidth',2, 'Color','k');
% %xline(2, ':','Stitch','LineWidth',2, 'Color','k');
% xline(2.25, '-','Go cue onset','LineWidth',2, 'Color','k');
% xline(2.75, ':','Stitch','LineWidth',2, 'Color','k');
% xline(3.75, '-','Response onset','LineWidth',2, 'Color','k');
% set(gca,'FontSize',15);
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

 visTimeGenAcc2D(decodeTimeStruct2D{3},pVal2Cutoff=0.05, axisLabel = 'Auditory onset',...
     clowLimit = 1,timePad = 0.1,chanceVal = 0.11111,clabel='Accuracy/Chance')
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

figure('Position', [scrsize(1) scrsize(2) scrsize(3) scrsize(4)/2]);

visTimePlot(timeEpoch,squeeze(mean(ieegStructPooled.data,2)),colval = colors(1,:))
xlabel('Time from auditory onset')
ylabel('z-score')
xline(0, '-','Auditory onset','LineWidth',2, 'Color','k');
xline(1.5, ':','Earliest Go time - Stitch','LineWidth',2, 'Color','k');
%xline(2, ':','Stitch','LineWidth',2, 'Color','k');
xline(2.25, '-','Go cue onset','LineWidth',2, 'Color','k');
xline(2.75, ':','Stitch','LineWidth',2, 'Color','k');
xline(3.75, '-','Response onset','LineWidth',2, 'Color','k');
set(gca,'FontSize',15);


%% visualize HG 2     4     5     6    15
scrsize = get(0, 'Screensize');

figure('Position', [scrsize(1) scrsize(2) scrsize(3) scrsize(4)/2]);

ax = visTimePlot3(timeEpoch,squeeze(mean(ieegStructPooled.data,2)),colval = colors(4,:))
hold on
visTimePlot3(timeEpoch,squeeze(mean(ieegStructPooled.data,2)),colval = colors(15,:),tileAxis = ax)
ylabel('HG (z-score)')


%% Phonotactic decoder
colvals = lines(5)
scrsize = get(0, 'Screensize');
figure('Position', [scrsize(1) scrsize(2) scrsize(3) scrsize(4)/2]); 
for iPhon = [1,2]
    visTimeGenAcc1D(decodeTimeStruct1D{iPhon},...
        pVal2Cutoff=0.05,...
        axisLabel = 'Auditory',clowLimit = 0,timePad = 0.2,...
        maxVal = (iPhon-1)*0.015, chanceVal = 0.25,clabel = 'Accuracy',colval=colvals(iPhon,:))
end

%     visTimeGenAcc1D(decodeTimeStruct1D_vcv{2},...
%         pVal2Cutoff=0.05,...
%         axisLabel = 'Auditory',clowLimit = 0,timePad = 0.2,...
%         maxVal = 0.015, chanceVal = 0.25,clabel = 'Accuracy',colval=colvals(2,:))


xline(0, '-','Auditory onset','LineWidth',2, 'Color','k');
xline(1.6, ':','Earliest Go time - Stitch','LineWidth',2, 'Color','k');
%xline(2, ':','Stitch','LineWidth',2, 'Color','k');
xline(2.35, '-','Go cue onset','LineWidth',2, 'Color','k');
xline(2.85, ':','Stitch','LineWidth',2, 'Color','k');
xline(3.85, '-','Response onset','LineWidth',2, 'Color','k');

%% Phonotactic decoder
colvals = lines(5)
scrsize = get(0, 'Screensize');
figure('Position', [scrsize(1) scrsize(2) scrsize(3) scrsize(4)/2]); 

    visTimeGenAcc1DCluster(decodeTimeStruct1D_cvc(5,:),decodeTimeStruct1D_cvc_shuffle(5,:),...
        pVal2Cutoff=0.05,...
        axisLabel = 'Auditory',clowLimit = 0,timePad = 0.2,...
        maxVal = 0, chanceVal = 0.25,clabel = 'Accuracy',colval=colvals(1,:))

    visTimeGenAcc1DCluster(decodeTimeStruct1D_vcv(5,:),decodeTimeStruct1D_vcv_shuffle(5,:),...
        pVal2Cutoff=0.05,...
        axisLabel = 'Auditory',clowLimit = 0,timePad = 0.2,...
        maxVal = 0.015, chanceVal = 0.25,clabel = 'Accuracy',colval=colvals(2,:))


xline(0, '-','Auditory onset','LineWidth',2, 'Color','k');
xline(1.5, ':','Earliest Go time - Stitch','LineWidth',2, 'Color','k');
%xline(2, ':','Stitch','LineWidth',2, 'Color','k');
xline(2.25, '-','Go cue onset','LineWidth',2, 'Color','k');
xline(2.75, ':','Stitch','LineWidth',2, 'Color','k');
xline(3.75, '-','Response onset','LineWidth',2, 'Color','k');
set(gca,'FontSize',15);