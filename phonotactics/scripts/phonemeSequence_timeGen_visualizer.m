%% Syllable decoder

figure;
    [ax,accResults] = visTimeGenAcc1DCluster(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,pVal2Cutoff=0.05,...
        axisLabel = 'Response',clowLimit = 0,timePad = 0.25,...
        maxVal = 0.31, chanceVal = 0.5,clabel = 'Accuracy',colval=colors(6,:),showShuffle=1,boxPlotPlace=0.95);
    ylim([0.2 1])
%% Phoneme decoder
colvals = lines(5)
scrsize = get(0, 'Screensize');
figure; 
%figure('Position', [scrsize(1) scrsize(2) scrsize(3) scrsize(4)/2]); 
%timeEpoch =[-0.5000    1.5000;   -0.7500    0.5000;   -1.0000    1.5000];
   
for iPhon = 1:3
    if(iPhon==1)
        [ax,accResultsTemp] = visTimeGenAcc1DCluster(decodeTimeStruct1D(iPhon,:),decodeTimeStruct1DSequenceshuffle(iPhon,:),pVal2Cutoff=0.05,...
            axisLabel = 'Response',clowLimit = 0,timePad = 0.25,...
             maxVal = (iPhon-1)*0.015, boxPlotPlace=(iPhon)*0.05+0.25, chanceVal = 0.1111,clabel = 'Accuracy',...
             colval=colvals(iPhon,:),showShuffle=1);

    else
       [axtemp,accResultsTemp] =  visTimeGenAcc1DCluster(decodeTimeStruct1D(iPhon,:),decodeTimeStruct1DSequenceshuffle(iPhon,:),pVal2Cutoff=0.05,...
            axisLabel = 'Response',clowLimit = 0,timePad =0.25,...
             maxVal = (iPhon-1)*0.015, boxPlotPlace=(iPhon)*0.05+0.25, chanceVal = 0.1111,clabel = 'Accuracy',...
             colval=colvals(iPhon,:),tileaxis = ax);
    end
    accResultsPhoneme{iPhon} = accResultsTemp;

end

% xline(0, '-','Auditory onset','LineWidth',2, 'Color','k');
% xline(1.5, ':','Earliest Go time - Stitch','LineWidth',2, 'Color','k');
% %xline(2, ':','Stitch','LineWidth',2, 'Color','k');
% xline(2.25, '-','Go cue onset','LineWidth',2, 'Color','k');
% xline(2.75, ':','Stitch','LineWidth',2, 'Color','k');
% xline(3.75, '-','Response onset','LineWidth',2, 'Color','k');
% set(gca,'FontSize',15);
%%

timeP21 = accResultsPhoneme{2}.timeRise(:)- accResultsPhoneme{1}.timeRise(:)';
timeP21Shuffle = accResultsPhoneme{2}.timeRiseShuffle(:)- accResultsPhoneme{1}.timeRiseShuffle(:)';
figure; histogram(timeP21Shuffle,50,'Normalization','probability','FaceColor',[0.5 0.5 0.5])
hold on; histogram(timeP21,10,'Normalization','probability','FaceColor',colvals(2,:),'EdgeColor',colvals(1,:),'LineWidth',2)
xlab('P2 - P1 (seconds)')
%[h,p,ci,stats] = vartest2(timeP21(:)' ,timeP21Shuffle(:)')
[h,p,ci,stats] = ttest2(timeP21(:)' ,timeP21Shuffle(:)',"Tail","right","Vartype","unequal")


timeP31 = accResultsPhoneme{3}.timeRise(:)- accResultsPhoneme{1}.timeRise(:)';
timeP31Shuffle = accResultsPhoneme{3}.timeRiseShuffle(:)- accResultsPhoneme{1}.timeRiseShuffle(:)';
% [h,p,ci,stats] = vartest2(timeP31(:)' ,timeP31Shuffle(:)')
[h,p,ci,stats] = ttest2(timeP31(:)' ,timeP31Shuffle(:)',"Tail","right","Vartype","unequal")
figure; histogram(timeP31Shuffle(:),50,'Normalization','probability','FaceColor',[0.5 0.5 0.5])
hold on; histogram(timeP31(:),10,'Normalization','probability','FaceColor',colvals(3,:),'EdgeColor',colvals(1,:),'LineWidth',2)
xlab('P3 - P1 (seconds)')

timeP32 = accResultsPhoneme{3}.timeRise(:)- accResultsPhoneme{2}.timeRise(:)';
timeP32Shuffle = accResultsPhoneme{3}.timeRiseShuffle(:)- accResultsPhoneme{2}.timeRiseShuffle(:)';
% [h,p,ci,stats] = vartest2(timeP32(:)' ,timeP32Shuffle(:)')
[h,p,ci,stats] = ttest2(timeP32(:)' ,timeP32Shuffle(:)',"Tail","right","Vartype","unequal")
figure; histogram(timeP32Shuffle(:),50,'Normalization','probability','FaceColor',[0.5 0.5 0.5])
hold on; histogram(timeP32(:),10,'Normalization','probability','FaceColor',colvals(3,:),'EdgeColor',colvals(2,:),'LineWidth',2)
xlab('P3 - P2 (seconds)')


%%

timeP21 = accResultsPhoneme{2}.timeCP1(:)- accResultsPhoneme{1}.timeCP1(:)';
timeP21Shuffle = accResultsPhoneme{2}.timeCP1Shuffle(:)- accResultsPhoneme{1}.timeCP1Shuffle(:)';
figure; histogram(timeP21Shuffle,50,'Normalization','probability','FaceColor',[0.5 0.5 0.5])
hold on; histogram(timeP21,10,'Normalization','probability','FaceColor',colvals(2,:),'EdgeColor',colvals(1,:),'LineWidth',2)
xlab('P2 - P1 (seconds)')
[h,p,ci,stats] = vartest2(timeP21(:)' ,timeP21Shuffle(:)','Tail','left')
%[h,p,ci,stats] = ttest2(timeP21(:)' ,timeP21Shuffle(:)')


timeP31 = accResultsPhoneme{3}.timeCP1(:)- accResultsPhoneme{1}.timeCP1(:)';
timeP31Shuffle = accResultsPhoneme{3}.timeCP1Shuffle(:)- accResultsPhoneme{1}.timeCP1Shuffle(:)';
[h,p,ci,stats] = vartest2(timeP31(:)' ,timeP31Shuffle(:)','Tail','left')
%[h,p,ci,stats] = ttest2(timeP31(:)' ,timeP31Shuffle(:)')

figure; histogram(timeP31Shuffle(:),50,'Normalization','probability','FaceColor',[0.5 0.5 0.5])
hold on; histogram(timeP31(:),10,'Normalization','probability','FaceColor',colvals(3,:),'EdgeColor',colvals(1,:),'LineWidth',2)
xlab('P3 - P1 (seconds)')

timeP32 = accResultsPhoneme{3}.timeCP1(:)- accResultsPhoneme{2}.timeCP1(:)';
timeP32Shuffle = accResultsPhoneme{3}.timeCP1Shuffle(:)- accResultsPhoneme{2}.timeCP1Shuffle(:)';
[h,p,ci,stats] = vartest2(timeP32(:)' ,timeP32Shuffle(:)','Tail','left')
figure; histogram(timeP32Shuffle(:),50,'Normalization','probability','FaceColor',[0.5 0.5 0.5])
hold on; histogram(timeP32(:),10,'Normalization','probability','FaceColor',colvals(3,:),'EdgeColor',colvals(2,:),'LineWidth',2)
xlab('P3 - P2 (seconds)')
%[h,p,ci,stats] = ttest2(timeP32(:)' ,timeP32Shuffle(:)')


%% Phonotactics regressor - cvc
colvals = lines(5)
scrsize = get(0, 'Screensize');

%figure('Position', [scrsize(1) scrsize(2) scrsize(3) scrsize(4)/2]); 
%timeEpoch =[-0.5000    1.5000;   -0.7500    0.5000;   -1.0000    1.5000];
   figure
for iPhon = [2,4]
   
    if(iPhon==2)
        [ax,accResultsTemp] = visTimeGenAcc1DCluster(decodeTimeStruct1D_cvc(iPhon,:),decodeTimeStruct1D_cvc_shuffle(iPhon,:),pVal2Cutoff=0.05,...
            axisLabel = 'Response',clowLimit = 0,timePad = 0.2,...
             maxVal = (iPhon-1)*0.015-0.1, boxPlotPlace=(iPhon)*0.05+0.3, chanceVal = 0,clabel = 'Coefficient of determination',...
             colval=colvals(iPhon,:),showShuffle=1);

    else
       [axtemp,accResultsTemp] =  visTimeGenAcc1DCluster(decodeTimeStruct1D_cvc(iPhon,:),decodeTimeStruct1D_cvc_shuffle(iPhon,:),pVal2Cutoff=0.05,...
            axisLabel = 'Response',clowLimit = 0,timePad =0.2,...
             maxVal = (iPhon-1)*0.015-0.1, boxPlotPlace=(iPhon)*0.05+0.3, chanceVal = 0,clabel = 'Coefficient of determination',...
             colval=colvals(iPhon,:),tileaxis=ax);
    end
    accResultsPhoneme{iPhon} = accResultsTemp;

end
ylim([-0.1 0.6])

%% Phonotactics regressor - vcv
colvals = lines(5)
scrsize = get(0, 'Screensize');

%figure('Position', [scrsize(1) scrsize(2) scrsize(3) scrsize(4)/2]); 
%timeEpoch =[-0.5000    1.5000;   -0.7500    0.5000;   -1.0000    1.5000];
   figure
for iPhon = [1,3]
   
    if(iPhon==1)
        [ax,accResultsTemp] = visTimeGenAcc1DCluster(decodeTimeStruct1D_vcv(iPhon,:),decodeTimeStruct1D_vcv_shuffle(iPhon,:),pVal2Cutoff=0.05,...
            axisLabel = 'Response',clowLimit = 0,timePad = 0.2,...
             maxVal = (iPhon-1)*0.015-0.1, boxPlotPlace=(iPhon)*0.05+0.5, chanceVal = 0,clabel = 'Coefficient of determination',...
             colval=colvals(iPhon,:),showShuffle=1);

    else
       [axtemp,accResultsTemp] =  visTimeGenAcc1DCluster(decodeTimeStruct1D_vcv(iPhon,:),decodeTimeStruct1D_vcv_shuffle(iPhon,:),pVal2Cutoff=0.05,...
            axisLabel = 'Response',clowLimit = 0,timePad =0.2,...
             maxVal = (iPhon-1)*0.015-0.1, boxPlotPlace=(iPhon)*0.05+0.5, chanceVal = 0,clabel = 'Coefficient of determination',...
             colval=colvals(iPhon,:),tileaxis=ax);
    end
    accResultsPhoneme{iPhon} = accResultsTemp;

end
ylim([-0.1 0.6])
%% Phoneme decoder - 1's consonant

hold on;
    visTimeGenAcc1DCluster(decodeTimeStruct1D_cvc(1,:),decodeTimeStruct1D_cvc_shuffle(1,:),pVal2Cutoff=5e-4,...
        axisLabel = 'Response',clowLimit = 0,timePad = 0.25,...
        maxVal = 0.065, chanceVal = 0.2,clabel = 'Accuracy',colval=colors(8,:))


    yline(0.2, '--','chance','LineWidth',1);
%% Phoneme decoder - consonant
colvals = lines(5)
scrsize = get(0, 'Screensize');
figure; 
% for iPhon = 1:3
    ax = visTimeGenAcc1DCluster(decodeTimeStruct1D_cvc(1,:),decodeTimeStruct1D_cvc_shuffle(1,:),pVal2Cutoff=1e-3,...
        axisLabel = 'Response',clowLimit = 0,timePad = 0.125,...
        maxVal = 0, chanceVal = 0.2,clabel = 'Accuracy',colval=colvals(1,:))

    visTimeGenAcc1DCluster(decodeTimeStruct1D_vcv(2,:),decodeTimeStruct1D_vcv_shuffle(2,:),pVal2Cutoff=1e-3,...
        axisLabel = 'Response',clowLimit = 0,timePad = 0.125,...
        maxVal = 0.015, chanceVal = 0.2,clabel = 'Accuracy',colval=colvals(2,:),tileaxis = ax)

    visTimeGenAcc1DCluster(decodeTimeStruct1D_cvc(3,:),decodeTimeStruct1D_cvc_shuffle(3,:),pVal2Cutoff=1e-3,...
        axisLabel = 'Response',clowLimit = 0,timePad = 0.125,...
        maxVal = 0.03, chanceVal = 0.2,clabel = 'Accuracy',colval=colvals(3,:),tileaxis = ax)
% end

%% Phoneme decoder - 1's vowel
hold on;

visTimeGenAcc1DCluster(decodeTimeStruct1D_vcv(1,:),decodeTimeStruct1D_vcv_shuffle(1,:),pVal2Cutoff=5e-4,...
    axisLabel = 'Response',clowLimit = 0,timePad = 0.25,...
    maxVal = 0.08, chanceVal = 0.2,clabel = 'Accuracy',colval=colors(8,:))
yline(0.25, '--','chance','LineWidth',1);
%% Phoneme decoder - vowel
colvals = lines(5)
scrsize = get(0, 'Screensize');
figure('Position', [scrsize(1) scrsize(2) scrsize(3) scrsize(4)/2]); 
% for iPhon = 1:3
    ax = visTimeGenAcc1DCluster(decodeTimeStruct1D_vcv(1,:),decodeTimeStruct1D_vcv_shuffle(1,:),pVal2Cutoff=1e-3,...
        axisLabel = 'Response',clowLimit = 0,timePad = 0.125,...
        maxVal = 0, chanceVal = 0.25,clabel = 'Accuracy',colval=colvals(1,:))

    visTimeGenAcc1DCluster(decodeTimeStruct1D_cvc(2,:),decodeTimeStruct1D_cvc_shuffle(2,:),pVal2Cutoff=1e-3,...
        axisLabel = 'Response',clowLimit = 0,timePad = 0.125,...
        maxVal = 0.015, chanceVal = 0.25,clabel = 'Accuracy',colval=colvals(2,:),tileaxis = ax)

    visTimeGenAcc1DCluster(decodeTimeStruct1D_vcv(3,:),decodeTimeStruct1D_vcv_shuffle(3,:),pVal2Cutoff=1e-3,...
        axisLabel = 'Response',clowLimit = 0,timePad = 0.125,...
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
    visTimeGenAcc1DCluster(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,pVal2Cutoff=0.01,...
        axisLabel = 'Response',clowLimit = 0,timePad = 0.25,...
        maxVal = 0.31, chanceVal = 0.5,clabel = 'Accuracy',colval=colors(4,:))
%     xline(0, ':','Auditory onset','LineWidth',2, 'Color','k');
% xline(2, '-','Stitch','LineWidth',2, 'Color','k');
% xline(2.5, ':','Go cue onset','LineWidth',2, 'Color','k');
% xline(3, '-','Stitch','LineWidth',2, 'Color','k');
% xline(4, ':','Response onset','LineWidth',2, 'Color','k');
%% Speech feature analysis
colvals = lines(8);
figure;
load('pooledSubject_Phoneme_Sequencing_Production_car_rh_z_score_prodelecs_data_v3_Start_syllable_decoded_1D.mat')
[ax,accResultsSyllable] = visTimeGenAcc1DCluster(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,pVal2Cutoff=0.05,...
        axisLabel = 'Response',clowLimit = 0,timePad = 0.2,boxPlotPlace=0.8,showPeaks = 5,...
        maxVal = -0.2, chanceVal = 0.5,clabel = 'Accuracy/Chance',colval=colvals(4,:),showAccperChance = 1,searchRange=[-1 0.5])
hold on;
load('pooledSubject_Phoneme_Sequencing_Production_car_rh_z_score_prodelecs_data_v3_Start_phoneme_decoded_1D.mat')
[ax,accResultsPhoneme] = visTimeGenAcc1DCluster(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,pVal2Cutoff=0.05,...
    axisLabel = 'Response',clowLimit = 0,timePad = 0.2,boxPlotPlace=1,showPeaks = 5,...
    maxVal = -0.1, chanceVal = 0.11111,clabel = 'Accuracy/Chance',colval=colvals(5,:),tileaxis = ax,showAccperChance = 1,searchRange=[-1 0.5])


% load('pooledSubject_Phoneme_Sequencing_Production_lh_car_z_score_prodelecs_data_v3_Start_consonant_decoded_1D.mat')
% visTimeGenAcc1DCluster(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,pVal2Cutoff=0.05,...
%     axisLabel = 'Response',clowLimit = 0,timePad = 0.2,boxPlotPlace=1,showPeaks = 5,...
%     maxVal = -0.5, chanceVal = 0.2,clabel = 'Accuracy/Chance',colval=colvals(5,:),tileaxis = ax,showAccperChance = 1,searchRange=[-1 0.5])
% hold on;
% load('pooledSubject_Phoneme_Sequencing_Production_lh_car_z_score_prodelecs_data_v3_Start_vowel_decoded_1D.mat')
% visTimeGenAcc1DCluster(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,pVal2Cutoff=0.05,...
%     axisLabel = 'Response',clowLimit = 0,timePad = 0.2,boxPlotPlace=1.2,showPeaks = 5,...
%     maxVal = -0.25, chanceVal = 0.25,clabel = 'Accuracy/Chance',colval=colvals(6,:),tileaxis = ax,showAccperChance = 1,searchRange=[-1 0.5])
figure;
histsyl = histfit(accResultsSyllable.timeCP1,5, 'kernel' );
histsyl(1).FaceColor = 'none';
histsyl(1).EdgeColor = 'none';
x_fit = histsyl(2).XData;
y_fit = histsyl(2).YData;
delete(histsyl(2));
hold on
patch([x_fit, fliplr(x_fit)], [zeros(1,length(y_fit)), fliplr(y_fit)], colvals(4,:), 'FaceAlpha', 0.3,'EdgeColor',colvals(4,:));
histsyl = histfit(accResultsPhoneme.timeCP1,5, 'kernel' );
histsyl(1).FaceColor = 'none';
histsyl(1).EdgeColor = 'none';
x_fit = histsyl(2).XData;
y_fit = histsyl(2).YData;
delete(histsyl(2));
hold on
patch([x_fit, fliplr(x_fit)], [zeros(1,length(y_fit)), fliplr(y_fit)], colvals(5,:), 'FaceAlpha', 0.3,'EdgeColor',colvals(5,:));

%% Roi analysis syllable

cue = 'Production';
decoder = 'syllable';
chanceVal = 0.5;
hemisphere = 'lh';
[p,n]=numSubplots(10);

figure;

% subplot(p(1),p(1),2)
%  load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_ifg_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
%  visTimeGenAcc1DCluster(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,pVal2Cutoff=0.05,...
%         axisLabel = 'Response',clowLimit = 0,timePad = 0.2,...
%         maxVal = 0.4, chanceVal = 0.5,clabel = 'Accuracy',colval=colors(7,:),showShuffle=1,boxPlotPlace=0.95,searchRange=[-1 1]);
%  ylim([0.3 1])
%  title('SMA')
 subplot(p(1),p(2),1)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_ifg_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
 visTimeGenAcc1DCluster(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,pVal2Cutoff=0.05,...
        axisLabel = 'Response',clowLimit = 0,timePad = 0.2,...
        maxVal = 0.4, chanceVal = 0.5,clabel = 'Accuracy',colval=colors(4,:),showShuffle=1,boxPlotPlace=0.95,searchRange=[-1 1]);
  ylim([0.3 1])
 title('IFG')
subplot(p(1),p(2),2)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_rmfg_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
 visTimeGenAcc1DCluster(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,pVal2Cutoff=0.05,...
        axisLabel = 'Response',clowLimit = 0,timePad = 0.2,...
        maxVal = 0.4, chanceVal = 0.5,clabel = 'Accuracy',colval=colors(3,:),showShuffle=1,boxPlotPlace=0.95,searchRange=[-1 1]);
 ylim([0.3 1])
 title('rMFG')
subplot(p(1),p(2),3)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_cmfg_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
 visTimeGenAcc1DCluster(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,pVal2Cutoff=0.05,...
        axisLabel = 'Response',clowLimit = 0,timePad = 0.2,...
        maxVal = 0.4, chanceVal = 0.5,clabel = 'Accuracy',colval=colors(10,:),showShuffle=1,boxPlotPlace=0.95,searchRange=[-1 1]);
 ylim([0.3 1])
 title('cMFG')
 subplot(p(1),p(2),4)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_smc_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
 visTimeGenAcc1DCluster(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,pVal2Cutoff=0.05,...
        axisLabel = 'Response',clowLimit = 0,timePad = 0.2,...
        maxVal = 0.4, chanceVal = 0.5,clabel = 'Accuracy',colval=colors(5,:),showShuffle=1,boxPlotPlace=0.95,searchRange=[-1 1]);
ylim([0.3 1])
title('SMC')
 subplot(p(1),p(2),5)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_ipc_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
 visTimeGenAcc1DCluster(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,pVal2Cutoff=0.05,...
        axisLabel = 'Response',clowLimit = 0,timePad = 0.2,...
        maxVal = 0.4, chanceVal = 0.5,clabel = 'Accuracy',colval=colors(2,:),showShuffle=1,boxPlotPlace=0.95,searchRange=[-1 1]);
ylim([0.3 1])
 title('IPC')
 subplot(p(1),p(2),6)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_insula_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
 visTimeGenAcc1DCluster(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,pVal2Cutoff=0.05,...
        axisLabel = 'Response',clowLimit = 0,timePad = 0.2,...
        maxVal = 0.4, chanceVal = 0.5,clabel = 'Accuracy',colval=colors(6,:),showShuffle=1,boxPlotPlace=0.95,searchRange=[-1 1]);
ylim([0.3 1])
 title('Insula')
 subplot(p(1),p(2),7)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_astg_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
 visTimeGenAcc1DCluster(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,pVal2Cutoff=0.05,...
        axisLabel = 'Response',clowLimit = 0,timePad = 0.2,...
        maxVal = 0.4, chanceVal = 0.5,clabel = 'Accuracy',colval=colors(9,:),showShuffle=1,boxPlotPlace=0.95,searchRange=[-1 1]);
    ylim([0.3 1])
 title('aSTG')
 subplot(p(1),p(2),8)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_pstg_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
 visTimeGenAcc1DCluster(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,pVal2Cutoff=0.05,...
        axisLabel = 'Response',clowLimit = 0,timePad = 0.2,...
        maxVal = 0.4, chanceVal = 0.5,clabel = 'Accuracy',colval=colors(1,:),showShuffle=1,boxPlotPlace=0.95,searchRange=[-1 1]);
ylim([0.3 1])
 title('pSTG')
  subplot(p(1),p(2),9)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_sts_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
 visTimeGenAcc1DCluster(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,pVal2Cutoff=0.05,...
        axisLabel = 'Response',clowLimit = 0,timePad = 0.2,...
        maxVal = 0.4, chanceVal = 0.5,clabel = 'Accuracy',colval=colors(8,:),showShuffle=1,boxPlotPlace=0.95,searchRange=[-1 1]);
 ylim([0.3 1])
 title('STS')
 subplot(p(1),p(2),10)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_heschl_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
 visTimeGenAcc1DCluster(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,pVal2Cutoff=0.05,...
        axisLabel = 'Response',clowLimit = 0,timePad = 0.2,...
        maxVal = 0.4, chanceVal = 0.5,clabel = 'Accuracy',colval=colors(8,:),showShuffle=1,boxPlotPlace=0.95,searchRange=[-1 1]);
 ylim([0.3 1])
 title('Heschl')
  %% Roi analysis phoneme

cue = 'Production';
decoder = 'phoneme';

hemisphere = 'rh';
[p,~]=numSubplots(10);

accResultsPhoneme =[];

figure;


 subplot(p(1),p(2),1)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_ifg_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
accResultsPhoneme{1} = plotPhonemeSequence(decodeTimeStruct1D,decodeTimeStruct1DSequenceshuffle)
ylim([-0.2 0.7]);title('IFG')
subplot(p(1),p(2),2)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_rmfg_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
accResultsPhoneme{2} = plotPhonemeSequence(decodeTimeStruct1D,decodeTimeStruct1DSequenceshuffle)
 ylim([-0.2 0.7])
 title('rMFG')
subplot(p(1),p(2),3)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_cmfg_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
 accResultsPhoneme{3} = plotPhonemeSequence(decodeTimeStruct1D,decodeTimeStruct1DSequenceshuffle)
 %plotPhonemeSyllableSequence(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,syllableType)
 ylim([-0.2 0.7])
 title('cMFG')
 subplot(p(1),p(2),4)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_smc_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
accResultsPhoneme{4} = plotPhonemeSequence(decodeTimeStruct1D,decodeTimeStruct1DSequenceshuffle)
 %plotPhonemeSyllableSequence(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,syllableType)
ylim([-0.2 0.7])
title('SMC')
 subplot(p(1),p(2),5)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_ipc_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
accResultsPhoneme{5} =  plotPhonemeSequence(decodeTimeStruct1D,decodeTimeStruct1DSequenceshuffle)
 %plotPhonemeSyllableSequence(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,syllableType)
 ylim([-0.2 0.7])
 title('IPC')
 subplot(p(1),p(2),6)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_insula_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
accResultsPhoneme{6} = plotPhonemeSequence(decodeTimeStruct1D,decodeTimeStruct1DSequenceshuffle)
 %plotPhonemeSyllableSequence(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,syllableType)
 ylim([-0.2 0.7])
 title('Insula')
 subplot(p(1),p(2),7)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_astg_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
accResultsPhoneme{7} = plotPhonemeSequence(decodeTimeStruct1D,decodeTimeStruct1DSequenceshuffle)
 %plotPhonemeSyllableSequence(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,syllableType)
 ylim([-0.2 0.7])
 title('aSTG')
 subplot(p(1),p(2),8)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_pstg_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
 accResultsPhoneme{8} = plotPhonemeSequence(decodeTimeStruct1D,decodeTimeStruct1DSequenceshuffle)
 %plotPhonemeSyllableSequence(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,syllableType)
ylim([-0.2 0.7])
 title('pSTG')
  subplot(p(1),p(2),9)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_sts_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
accResultsPhoneme{9} = plotPhonemeSequence(decodeTimeStruct1D,decodeTimeStruct1DSequenceshuffle)
 %plotPhonemeSyllableSequence(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,syllableType)
 ylim([-0.2 0.7])
 title('STS')
 subplot(p(1),p(2),10)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_heschl_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
accResultsPhoneme{10} = plotPhonemeSequence(decodeTimeStruct1D,decodeTimeStruct1DSequenceshuffle)
 ylim([-0.2 0.7])
 title('PAC')
 %% Roi analysis phoneme - group

cue = 'Production';
decoder = 'phoneme';

hemisphere = 'rh';
[p,~]=numSubplots(3);

accResultsPhoneme =[];

figure;


 subplot(p(1),p(2),1)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_group_1_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
accResultsPhoneme{1} = plotPhonemeSequence(decodeTimeStruct1D,decodeTimeStruct1DSequenceshuffle)
ylim([-0.2 0.7]);title('group 1')
subplot(p(1),p(2),2)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_group_2_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
accResultsPhoneme{2} = plotPhonemeSequence(decodeTimeStruct1D,decodeTimeStruct1DSequenceshuffle)
 ylim([-0.2 0.7])
 title('group 2')
subplot(p(1),p(2),3)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_group_3_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
 accResultsPhoneme{3} = plotPhonemeSequence(decodeTimeStruct1D,decodeTimeStruct1DSequenceshuffle)
 %plotPhonemeSyllableSequence(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,syllableType)
 ylim([-0.2 0.7])
 title('group 3')
 %% Roi analysis cvc/vcv

cue = 'Production';
decoder = 'cvc_class';

hemisphere = 'lh';
[p,n]=numSubplots(10);
syllableType = 2;


figure;

subplot(p(1),p(1),2)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_sma_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
 plotPhonemeSyllableSequence(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,syllableType)
 ylim([-1 2])
 title('SMA')
 subplot(p(1),p(2),2)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_ifg_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
plotPhonemeSyllableSequence(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,syllableType)
 ylim([-1 2]);title('IFG')
subplot(p(1),p(2),3)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_rmfg_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
 plotPhonemeSyllableSequence(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,syllableType)
 ylim([-1 2])
 title('rMFG')
subplot(p(1),p(2),4)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_cmfg_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
 plotPhonemeSyllableSequence(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,syllableType)
 ylim([-1 2])
 title('cMFG')
 subplot(p(1),p(2),5)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_smc_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
plotPhonemeSyllableSequence(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,syllableType)
 ylim([-1 2])
title('SMC')
 subplot(p(1),p(2),6)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_ipc_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
 plotPhonemeSyllableSequence(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,syllableType)
 ylim([-1 2])
 title('IPC')
 subplot(p(1),p(2),7)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_insula_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
plotPhonemeSyllableSequence(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,syllableType)
 ylim([-1 2])
 title('Insula')
 subplot(p(1),p(2),8)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_superiortemporal_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
plotPhonemeSyllableSequence(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,syllableType)
 ylim([-1 2])
 title('STG')
 subplot(p(1),p(2),9)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_middletemporal_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
 plotPhonemeSyllableSequence(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,syllableType)
 ylim([-1 2])
 title('MTG')
  subplot(p(1),p(2),10)
 load(['pooledSubject_Phoneme_Sequencing_' cue '_car_' hemisphere '_inferiortemporal_z_score_prodelecs_data_v3_Start_' decoder '_decoded_1D.mat'])
plotPhonemeSyllableSequence(decodeTimeStruct1D,decodeTimeStruct1Dshuffle,syllableType)
 ylim([-1 2])
 title('ITG')
 
  %% Roi analysis phoneme



hemisphere = 'rh';

[p,n]=numSubplots(10);

accResults = [];

figure;
 subplot(p(1),p(2),1)
[accResults{1}] = plotSpeechSequence(hemisphere,'ifg');
 title('IFG')
subplot(p(1),p(2),2)
[accResults{2}] =  plotSpeechSequence(hemisphere,'rmfg');
 title('rMFG')
subplot(p(1),p(2),3)
[accResults{3}] = plotSpeechSequence(hemisphere,'cmfg');
 title('cMFG')
 subplot(p(1),p(2),4)
[accResults{4}] =  plotSpeechSequence(hemisphere,'smc');
title('SMC')
 subplot(p(1),p(2),5)
[accResults{5}] = plotSpeechSequence(hemisphere,'ipc');
 title('IPC')
 subplot(p(1),p(2),6)
 [accResults{6}] = plotSpeechSequence(hemisphere,'insula');
 title('Insula')
 subplot(p(1),p(2),7)
 [accResults{7}] = plotSpeechSequence(hemisphere,'astg');
 title('aSTG')
 subplot(p(1),p(2),8)
 [accResults{8}] = plotSpeechSequence(hemisphere,'pstg');
 title('pSTG')
  subplot(p(1),p(2),9)
 [accResults{9}] = plotSpeechSequence(hemisphere,'sts');
 title('STS')
 subplot(p(1),p(2),10)
[accResults{10}] = plotSpeechSequence(hemisphere,'heschl');
 title('A1')
 %% 
hemisphere = 'lh';

[p,n]=numSubplots(3);

accResultsGroup = [];

figure;
 subplot(p(1),p(2),1)
[accResultsGroup{1}] = plotSpeechSequence(hemisphere,'group_1');
 title('Group 1')
subplot(p(1),p(2),2)
[accResultsGroup{2}] =  plotSpeechSequence(hemisphere,'group_2');
 title('Group 2')
subplot(p(1),p(2),3)
[accResultsGroup{3}] = plotSpeechSequence(hemisphere,'group_3');
 title('Group 3')
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
xlab('Time from auditory onset')
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
visTimePlot(timeEpoch,squeeze(mean(ieegStructPooled.data,2)),colval = colors(8,:), labels = 'Response Onset (s)' )


hold on
visTimePlot3(timeEpoch,squeeze(mean(ieegStructPooled.data,2)),colval = colors(15,:),tileAxis = ax)
ylabel('HG (z-score)')


% ieegMeanNorm = (max(ieegMean2plot'));
%%
ieegMeanNormTime = [];
ieegMeanNormShuffTime = [];
ieegMeanWin = extractTimeWindow(ieegMean2plot,20,2);
ieegMeanWinShuff = extractTimeWindow(ieegMeanShuff2plot,20,2);
timeGammaWin = linspace(-1,1.5,size(ieegMeanWin,3));
for iTime= 1:size(ieegMeanWin)
    ieegMeanNormTime(:,iTime) = vecnorm(squeeze(ieegMeanWin(iTime,:,:))',1)./size(ieegMeanWin,3);
    ieegMeanNormShuffTime(:,iTime) = vecnorm(squeeze(ieegMeanWinShuff(iTime,:,:))',1)./size(ieegMeanWinShuff,3);

end



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

%% 
xlabelsAnat = {'MTG&STG','SMC','MFG','IFG','IPC','Insula', 'SFG', 'Fusiform'};
figure; 
b = bar([21.6 30.7 16.3 11.14 18.97 16.18 15.3 18.72]);
yline(11.1111, '--','chance','LineWidth',2, 'LaberhorizontalAlignment','left');
cmap = colormap(colors)
for k = 1:length(xlabelsAnat)
    b.CData(k,:) = colors(k,:);
end
ylabel('decoding accuracy (%)');
formatTicks(gca)
set(gca,'xticklabels',xlabelsAnat)
xlim([0 9])

figure; 
bar([72.6 68.1 67.6 52.9 37.8 55.8 62.1 44 55.8]);
yline(25, '--','chance','LineWidth',2, 'LaberhorizontalAlignment','left');
ylabel('C1 accuracy (%)');
formatTicks(gca)
set(gca,'xticklabels',xlabelsAnat)




% figure; 
% bar([64 63.9 39.9 40.9 31.3 64.4 ]);
% yline(25, '--','chance','LineWidth',2, 'LaberhorizontalAlignment','left');
% ylabel('Articulator accuracy (%)');
% formatTicks(gca)
% set(gca,'xticklabels',xlabelsAnat)


figure; 
bar([90 88.7 85.7 73.5 60.7 76.3 79.8 65.3 73.8]);
yline(50, '--','chance','LineWidth',2, 'LaberhorizontalAlignment','left');
ylabel('Syllable accuracy (%)');
formatTicks(gca)
set(gca,'xticklabels',xlabelsAnat)

%%
channelInfoRoi = extractChannelLocation(Subject,channelNamePooledSelect(groupingIds==0))
figure; 
[uni,idx1,idx2] = unique({channelInfoRoi(:).Location},'sorted');
counts = accumarray( idx2, ones(size(idx2)) ) ;
[sortIdx1,sortIdIds] = sort(counts,'ascend');
barh(sortIdx1);
set(gca,'YTick',1:length(uni))  
yticklabels(uni(sortIdIds));


%%
figure;
x = rand(10,1)/20;
xg = randi(2,10,1)/40;
y = rand(10,1);
yg = randi(2,10,1)/2;
bv = boxplot(x,xg,'orientation','vertical');
xlim manual
hold on
bh = boxplot(y,yg,'orientation','horizontal');
hold off
xlim([0 2.5])
ylim([0 2.5])

%% 
timeGamma = linspace(-1.5,1.5,300);

for iSubj = 1:length(Subject)
    figure;
    plot(timeGamma,squeeze(mean(ieegHGStruct(iSubj).ieegHGNorm.data,2 )));
    xlab('Time from speech production onset')
    title(['temporal channels: ' ieegHGStruct(iSubj).channelName{1}(1:3)])
    %saveas(gcf,['temporal_channels_' ieegHGStruct(iSubj).channelName{1}(1:3) '.png'])
end

%%
timeGamma = linspace(-1.5,1.5,600);

for iSubj = 1:length(Subject)
    goTimeTrials = ([Subject(iSubj).Trials.Go] - [Subject(iSubj).Trials.Auditory])./30000;
   goTimeTrialInfo = [];
   for itrial = 1:208
       goTimeTrialInfo(itrial) = [Subject(iSubj).trialInfo{1,itrial}.goStart - Subject(iSubj).trialInfo{1,itrial}.audioStart];
   end
   figure;
  plot(goTimeTrials)
  hold on;
  plot(goTimeTrialInfo)
    xlab('Trials')
     ylabel('Go time ')
    title(Subject(iSubj).Name)
    
end
%% Visualize HG

%[ieegStructPooledBase,phonemeTrialPooled,channelNamePooled] = poolChannelWithMaxTrial(ieegHGStructBase,trialInfoStruct(maxTrial));    

 load('pooledSubject_Phoneme_Sequencing_Go_car_z_score_prodelecs_data_v3.mat')

ieegEpoch = [];
ieegEpochBga = [];
    channelNameEpoch = [];
    for iSubject = 1:length(ieegHGStruct)
        iSubject
       ieegStruct2process = ieegHGStruct(iSubject).ieegHGNorm;
%        
%        ieegStruct2process = applyNormalization(ieegStruct2process,ieegHGStruct(iSubject).normFactor,3);
        phonemeTrial = phonemeSequenceTrialParser(ieegHGStruct(iSubject).trialInfo);
       cvcIds = find(phonemeTrial.syllableUnit(:,1)'==2); 
       vcvIds = find(phonemeTrial.syllableUnit(:,1)'==1);
       % ieegEpoch = cat(1,ieegEpoch,squeeze(nanmean(ieegStruct2process.data,2)));
       ieegEpoch = cat(1,ieegEpoch,squeeze(nanmean(ieegStruct2process.data(:,vcvIds,:),2)));
%         ieegStruct2process = removeZScore(ieegHGStruct(iSubject).ieegHGNorm,ieegHGStruct(iSubject).normFactor);
%         ieegStructNorm = applyNormalization(ieegStruct2process,ieegHGStruct(iSubject).normFactor,3);
%         ieegEpochBga = cat(1,ieegEpochBga,squeeze(nanmean(ieegStructNorm.data,2)));
        channelNameEpoch = [channelNameEpoch ieegHGStruct(iSubject).channelName];
    end
tw = ieegHGStruct(iSubject).ieegHGNorm.tw;
timeGamma = linspace(tw(1),tw(2), size(ieegEpoch,2));
ieegEpoch = ieegEpoch(:,timeGamma>=-1);
timeGamma = timeGamma(timeGamma>=-1);
 timeRange = [-1 tw(2)];
    load('pooledSubject_Phoneme_Sequencing_start_car_z_score_prodelecs_data_v3.mat')
ieegBase = [];
ieegBaseBga = [];
channelNameBase = [];
for iSubject = 1:length(ieegHGStruct)
    iSubject
%     ieegStruct2process = removeZScore(ieegHGStruct(iSubject).ieegHGNorm,ieegHGStruct(iSubject).normFactor);
%     ieegStruct2process = applyNormalization(ieegStruct2process,ieegHGStruct(iSubject).normFactor,3);
%       
    ieegStruct2process = ieegHGStruct(iSubject).ieegHGNorm;
%     [NumTrials, goodtrials] = remove_bad_trials(ieegStruct2process.data, 0.1);
    ieegStructNew = extendTimeEpoch(ieegStruct2process,timeRange);
     
    ieegBase = cat(1,ieegBase,squeeze(nanmean(ieegStructNew.data,2)));

%     ieegStruct2process = removeZScore(ieegHGStruct(iSubject).ieegHGNorm,ieegHGStruct(iSubject).normFactor);
%     ieegStruct2process = applyNormalization(ieegStruct2process,ieegHGStruct(iSubject).normFactor,3);
%     ieegStructNew = extendTimeEpoch(ieegStruct2process,timeRange);
%     ieegBaseBga = cat(1,ieegBaseBga,squeeze(nanmean(ieegStructNew.data,2)));
    channelNameBase = [channelNameBase ieegHGStruct(iSubject).channelName];
end

chanSelect = ismember(channelNameEpoch,channelNameBase);
ieegEpoch = ieegEpoch(chanSelect,:);
%ieegBase = ieegBase(~chanSelect,:);
channelNameEpoch = channelNameEpoch(chanSelect);
%channelNameBase = channelNameBase(~chanSelect);



%%
load colors.mat
groupingIds =zeros(1,length(channelNameBase));
hemisphere = 'rh';
boxplace = 1.5;
ylineval = -0.1;
yrange = [-0.15 1.6];
boxStd = 0.1;
showBoxPlot = 0;
xlab = 'go onset'
    figure;
[p,n]=numSubplots(10);

hold on;

timeCentroid = [];
subplot(p(1),p(2),1)
bnLabels = extractBNlabels('ifg',hemisphere)
roiIds = extractRoiIds(Subject, channelNameEpoch,bnLabels);
roiColor = colors(4,:);
ieegRoiEpoch = ieegEpoch(roiIds,:);
ieegRoiBase = ieegBase(roiIds,:);

[ieegRoiEpochSubj,ieegRoiEpochSubjAll,ieegRoiEpochSubjElec] = extractMeanControlSubject(ieegRoiEpoch,channelNameEpoch(roiIds));
[ieegRoiBaseSubj,ieegRoiBaseSubjAll,ieegBaseEpochSubjElec] = extractMeanControlSubject(ieegRoiBase,channelNameEpoch(roiIds));

pmask = zeros(1,size(ieegRoiBase,2));
[zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegRoiEpochSubjAll,ieegRoiBaseSubjAll,1000,1,1.645);

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
visTimePlot(timeRange,ieegRoiEpochSubjAll,colval = roiColor, labels = xlab )
hold on;
visTimePlot(timeRange,ieegRoiBaseSubjAll,colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGamma(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',...
    MarkerEdgeColor=roiColor,MarkerFaceColor=roiColor);
groupingIds(roiIds) = 4;
title(['IFG ' num2str(sum(roiIds))])
ieegRoiEpochScaled = ieegRoiEpoch.*pmask;
[maxAcc,maxId] = max(ieegRoiEpochScaled');
timeMax = timeGamma(maxId);
timeMaxRoi{1} = timeMax;
maxId = [];
timeCentroidRoi = [];

for iElec = 1:size(ieegRoiEpochScaled,1)
    timeCentroidRoi = [timeCentroidRoi calculateTimeCentroid(timeGamma,ieegRoiEpochScaled(iElec,:))];
    maxPoints = findchangepts(ieegRoiEpochScaled(iElec,:),MaxNumChanges=4,Statistic="mean");
    if(~isempty(maxPoints))
        maxId = [maxId maxPoints(1)];
    end
end
timeCentroidRoi = timeCentroidRoi(timeCentroidRoi>-0.9);
timeCentroid{1} = timeCentroidRoi;
timeChangePoint{1} = timeGamma(maxId);
if(showBoxPlot&&sum(pmask)>0)
    boxchart(boxplace*ones(size(time_centroid)),...
        double(time_centroid),'orientation','horizontal','BoxFaceColor',...
        roiColor*0.75,'BoxWidth',boxStd, 'BoxLineColor',...
        roiColor,'MarkerStyle','none')
end
ylim(yrange)


subplot(p(1),p(2),2)
bnLabels = extractBNlabels('rmfg',hemisphere)
roiIds = extractRoiIds(Subject, channelNameEpoch,bnLabels);
roiColor = colors(3,:);
ieegRoiEpoch = ieegEpoch(roiIds,:);
ieegRoiBase = ieegBase(roiIds,:);

[ieegRoiEpochSubj,ieegRoiEpochSubjAll,ieegRoiEpochSubjElec] = extractMeanControlSubject(ieegRoiEpoch,channelNameEpoch(roiIds));
[ieegRoiBaseSubj,ieegRoiBaseSubjAll,ieegBaseEpochSubjElec] = extractMeanControlSubject(ieegRoiBase,channelNameEpoch(roiIds));

pmask = zeros(1,size(ieegRoiBase,2));
[zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegRoiEpochSubjAll,ieegRoiBaseSubjAll,1000,1,1.645);

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
visTimePlot(timeRange,ieegRoiEpochSubjAll,colval = roiColor, labels = xlab )
hold on;
visTimePlot(timeRange,ieegRoiBaseSubjAll,colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGamma(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',...
    MarkerEdgeColor=roiColor,MarkerFaceColor=roiColor);
groupingIds(roiIds) = 3;
title(['rMFG ' num2str(sum(roiIds))])
ieegRoiEpochScaled = ieegRoiEpoch.*pmask;
[maxAcc,maxId] = max(ieegRoiEpochScaled');
timeMax = timeGamma(maxId);
timeMaxRoi{2} = timeMax;
maxId = [];
timeCentroidRoi = [];
for iElec = 1:size(ieegRoiEpochScaled,1)
    timeCentroidRoi = [timeCentroidRoi calculateTimeCentroid(timeGamma,ieegRoiEpoch(iElec,:))];
    maxPoints = findchangepts(ieegRoiEpoch(iElec,:),MaxNumChanges=4,Statistic="mean");
    if(~isempty(maxPoints))
        maxId = [maxId maxPoints(1)];
    end
end
timeCentroidRoi = timeCentroidRoi(timeCentroidRoi>-0.9);
timeCentroid{2} = timeCentroidRoi;
timeChangePoint{2} = timeGamma(maxId);
if(showBoxPlot&&sum(pmask)>0)
    boxchart(boxplace*ones(size(time_centroid)),...
        double(time_centroid),'orientation','horizontal','BoxFaceColor',...
        roiColor*0.75,'BoxWidth',boxStd, 'BoxLineColor',...
        roiColor,'MarkerStyle','none')
end
ylim(yrange)



subplot(p(1),p(2),3)
bnLabels = extractBNlabels('cmfg',hemisphere)
roiIds = extractRoiIds(Subject, channelNameEpoch,bnLabels);
roiColor = colors(10,:);
ieegRoiEpoch = ieegEpoch(roiIds,:);
ieegRoiBase = ieegBase(roiIds,:);
 

  [ieegRoiEpochSubj,ieegRoiEpochSubjAll,ieegRoiEpochSubjElec] = extractMeanControlSubject(ieegRoiEpoch,channelNameEpoch(roiIds));
[ieegRoiBaseSubj,ieegRoiBaseSubjAll,ieegBaseEpochSubjElec] = extractMeanControlSubject(ieegRoiBase,channelNameEpoch(roiIds));

pmask = zeros(1,size(ieegRoiBase,2));
[zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegRoiEpochSubjAll,ieegRoiBaseSubjAll,1000,1,1.645);

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
visTimePlot(timeRange,ieegRoiEpochSubjAll,colval = roiColor, labels = xlab )
hold on;
visTimePlot(timeRange,ieegRoiBaseSubjAll,colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGamma(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',...
    MarkerEdgeColor=roiColor,MarkerFaceColor=roiColor);
groupingIds(roiIds) = 10;
title(['cMFG ' num2str(sum(roiIds))])
ieegRoiEpochScaled = ieegRoiEpoch.*pmask;
[maxAcc,maxId] = max(ieegRoiEpochScaled');
timeMax = timeGamma(maxId);
timeMaxRoi{3} = timeMax;
maxId = [];
timeCentroidRoi = [];
for iElec = 1:size(ieegRoiEpochScaled,1)
    timeCentroidRoi = [timeCentroidRoi calculateTimeCentroid(timeGamma,ieegRoiEpoch(iElec,:))];
    maxPoints = findchangepts(ieegRoiEpoch(iElec,:),MaxNumChanges=4,Statistic="mean");
    if(~isempty(maxPoints))
        maxId = [maxId maxPoints(1)];
    end
end
timeCentroidRoi = timeCentroidRoi(timeCentroidRoi>-0.9);
timeCentroid{3} = timeCentroidRoi;
timeChangePoint{3} = timeGamma(maxId);
if(showBoxPlot&&sum(pmask)>0)
    boxchart(boxplace*ones(size(time_centroid)),...
        double(time_centroid),'orientation','horizontal','BoxFaceColor',...
        roiColor*0.75,'BoxWidth',boxStd, 'BoxLineColor',...
        roiColor,'MarkerStyle','none')
end
ylim(yrange)

subplot(p(1),p(2),4)
bnLabels = extractBNlabels('smc',hemisphere)
roiIds = extractRoiIds(Subject, channelNameEpoch,bnLabels);
roiColor = colors(5,:);
ieegRoiEpoch = ieegEpoch(roiIds,:);
ieegRoiBase = ieegBase(roiIds,:);


  [ieegRoiEpochSubj,ieegRoiEpochSubjAll,ieegRoiEpochSubjElec] = extractMeanControlSubject(ieegRoiEpoch,channelNameEpoch(roiIds));
[ieegRoiBaseSubj,ieegRoiBaseSubjAll,ieegBaseEpochSubjElec] = extractMeanControlSubject(ieegRoiBase,channelNameEpoch(roiIds));

pmask = zeros(1,size(ieegRoiBase,2));
[zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegRoiEpochSubjAll,ieegRoiBaseSubjAll,1000,1,1.645);

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
visTimePlot(timeRange,ieegRoiEpochSubjAll,colval = roiColor, labels = xlab )
hold on;
visTimePlot(timeRange,ieegRoiBaseSubjAll,colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGamma(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',...
    MarkerEdgeColor=roiColor,MarkerFaceColor=roiColor);
groupingIds(roiIds) = 5;
title(['SMC ' num2str(sum(roiIds))])
ieegRoiEpochScaled = ieegRoiEpoch.*pmask; [maxAcc,maxId] = max(ieegRoiEpochScaled');
timeMax = timeGamma(maxId);
timeMaxRoi{4} = timeMax;
maxId = [];
timeCentroidRoi = [];
for iElec = 1:size(ieegRoiEpochScaled,1)
    timeCentroidRoi = [timeCentroidRoi calculateTimeCentroid(timeGamma,ieegRoiEpochScaled(iElec,:))];
    maxPoints = findchangepts(ieegRoiEpochScaled(iElec,:),MaxNumChanges=4,Statistic="mean");
    if(~isempty(maxPoints))
        maxId = [maxId maxPoints(1)];
    end
end
timeCentroidRoi = timeCentroidRoi(timeCentroidRoi>-0.9);
timeCentroid{4} = timeCentroidRoi;
timeChangePoint{4} = timeGamma(maxId);

if(showBoxPlot&&sum(pmask)>0)
    boxchart(boxplace*ones(size(time_centroid)),...
        double(time_centroid),'orientation','horizontal','BoxFaceColor',...
        roiColor*0.75,'BoxWidth',boxStd, 'BoxLineColor',...
        roiColor,'MarkerStyle','none')
end
ylim(yrange)


subplot(p(1),p(2),5)
bnLabels = extractBNlabels('ipc',hemisphere)
roiIds = extractRoiIds(Subject, channelNameEpoch,bnLabels);
roiColor = colors(2,:);
ieegRoiEpoch = ieegEpoch(roiIds,:);
ieegRoiBase = ieegBase(roiIds,:);
 

 [ieegRoiEpochSubj,ieegRoiEpochSubjAll,ieegRoiEpochSubjElec] = extractMeanControlSubject(ieegRoiEpoch,channelNameEpoch(roiIds));
[ieegRoiBaseSubj,ieegRoiBaseSubjAll,ieegBaseEpochSubjElec] = extractMeanControlSubject(ieegRoiBase,channelNameEpoch(roiIds));

pmask = zeros(1,size(ieegRoiBase,2));
[zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegRoiEpochSubjAll,ieegRoiBaseSubjAll,1000,1,1.645);

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
visTimePlot(timeRange,ieegRoiEpochSubjAll,colval = roiColor, labels = xlab )
hold on;
visTimePlot(timeRange,ieegRoiBaseSubjAll,colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGamma(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',...
    MarkerEdgeColor=roiColor,MarkerFaceColor=roiColor);
groupingIds(roiIds) = 2;
title(['IPC ' num2str(sum(roiIds))])
ieegRoiEpochScaled = ieegRoiEpoch.*pmask; [maxAcc,maxId] = max(ieegRoiEpochScaled');
timeMax = timeGamma(maxId);
timeMaxRoi{5} = timeMax;
maxId = [];
timeCentroidRoi = [];
for iElec = 1:size(ieegRoiEpochScaled,1)
    timeCentroidRoi = [timeCentroidRoi calculateTimeCentroid(timeGamma,ieegRoiEpochScaled(iElec,:))];
    maxPoints = findchangepts(ieegRoiEpochScaled(iElec,:),MaxNumChanges=4,Statistic="mean");
    if(~isempty(maxPoints))
        maxId = [maxId maxPoints(1)];
    end
end
timeCentroidRoi = timeCentroidRoi(timeCentroidRoi>-0.9);
timeCentroid{5} = timeCentroidRoi;
timeChangePoint{5} = timeGamma(maxId);
if(showBoxPlot&&sum(pmask)>0)
    boxchart(boxplace*ones(size(time_centroid)),...
        double(time_centroid),'orientation','horizontal','BoxFaceColor',...
        roiColor*0.75,'BoxWidth',boxStd, 'BoxLineColor',...
        roiColor,'MarkerStyle','none')
end
ylim(yrange)


subplot(p(1),p(2),6)
bnLabels = extractBNlabels('insula',hemisphere)
roiIds = extractRoiIds(Subject, channelNameEpoch,bnLabels);
roiColor = colors(6,:);
ieegRoiEpoch = ieegEpoch(roiIds,:);
ieegRoiBase = ieegBase(roiIds,:);
 

 [ieegRoiEpochSubj,ieegRoiEpochSubjAll,ieegRoiEpochSubjElec] = extractMeanControlSubject(ieegRoiEpoch,channelNameEpoch(roiIds));
[ieegRoiBaseSubj,ieegRoiBaseSubjAll,ieegBaseEpochSubjElec] = extractMeanControlSubject(ieegRoiBase,channelNameEpoch(roiIds));

pmask = zeros(1,size(ieegRoiBase,2));
[zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegRoiEpochSubjAll,ieegRoiBaseSubjAll,1000,1,1.645);

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
visTimePlot(timeRange,ieegRoiEpochSubjAll,colval = roiColor, labels = xlab )
hold on;
visTimePlot(timeRange,ieegRoiBaseSubjAll,colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGamma(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',...
    MarkerEdgeColor=roiColor,MarkerFaceColor=roiColor);
groupingIds(roiIds) = 6;
title(['Insula ' num2str(sum(roiIds))])
ieegRoiEpochScaled = ieegRoiEpoch.*pmask; [maxAcc,maxId] = max(ieegRoiEpochScaled');
timeMax = timeGamma(maxId);
timeMaxRoi{6} = timeMax;
maxId = [];
timeCentroidRoi = [];
for iElec = 1:size(ieegRoiEpochScaled,1)
    timeCentroidRoi = [timeCentroidRoi calculateTimeCentroid(timeGamma,ieegRoiEpochScaled(iElec,:))];
    maxPoints = findchangepts(ieegRoiEpochScaled(iElec,:),MaxNumChanges=4,Statistic="mean");
    if(~isempty(maxPoints))
        maxId = [maxId maxPoints(1)];
    end
end
timeCentroidRoi = timeCentroidRoi(timeCentroidRoi>-0.9);
timeCentroid{6} = timeCentroidRoi;
timeChangePoint{6} = timeGamma(maxId);
if(showBoxPlot&&sum(pmask)>0)
    boxchart(boxplace*ones(size(time_centroid)),...
        double(time_centroid),'orientation','horizontal','BoxFaceColor',...
        roiColor*0.75,'BoxWidth',boxStd, 'BoxLineColor',...
        roiColor,'MarkerStyle','none')
end
ylim(yrange)


subplot(p(1),p(2),7)
bnLabels = extractBNlabels('astg',hemisphere)
roiIds = extractRoiIds(Subject, channelNameEpoch,bnLabels);
roiColor = colors(9,:);
ieegRoiEpoch = ieegEpoch(roiIds,:);
ieegRoiBase = ieegBase(roiIds,:);
 

[ieegRoiEpochSubj,ieegRoiEpochSubjAll,ieegRoiEpochSubjElec] = extractMeanControlSubject(ieegRoiEpoch,channelNameEpoch(roiIds));
[ieegRoiBaseSubj,ieegRoiBaseSubjAll,ieegBaseEpochSubjElec] = extractMeanControlSubject(ieegRoiBase,channelNameEpoch(roiIds));

pmask = zeros(1,size(ieegRoiBase,2));
[zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegRoiEpochSubjAll,ieegRoiBaseSubjAll,1000,1,1.645);

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
visTimePlot(timeRange,ieegRoiEpochSubjAll,colval = roiColor, labels = xlab )
hold on;
visTimePlot(timeRange,ieegRoiBaseSubjAll,colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGamma(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',...
    MarkerEdgeColor=roiColor,MarkerFaceColor=roiColor);
groupingIds(roiIds) = 9;
title(['aSTG ' num2str(sum(roiIds))])
ieegRoiEpochScaled = ieegRoiEpoch.*pmask; [maxAcc,maxId] = max(ieegRoiEpochScaled');
timeMax = timeGamma(maxId);
timeMaxRoi{7} = timeMax;
maxId = [];
timeCentroidRoi = [];
for iElec = 1:size(ieegRoiEpochScaled,1)
    timeCentroidRoi = [timeCentroidRoi calculateTimeCentroid(timeGamma,ieegRoiEpochScaled(iElec,:))];
    maxPoints = findchangepts(ieegRoiEpochScaled(iElec,:),MaxNumChanges=4,Statistic="mean");
    if(~isempty(maxPoints))
        maxId = [maxId maxPoints(1)];
    end
end
timeCentroidRoi = timeCentroidRoi(timeCentroidRoi>-0.9);
timeCentroid{7} = timeCentroidRoi;
timeChangePoint{7} = timeGamma(maxId);
if(showBoxPlot&&sum(pmask)>0)
    boxchart(boxplace*ones(size(time_centroid)),...
        double(time_centroid),'orientation','horizontal','BoxFaceColor',...
        roiColor*0.75,'BoxWidth',boxStd, 'BoxLineColor',...
        roiColor,'MarkerStyle','none')
end
ylim(yrange)


subplot(p(1),p(2),8)
bnLabels = extractBNlabels('pstg',hemisphere)
roiIds = extractRoiIds(Subject, channelNameEpoch,bnLabels);
roiColor = colors(1,:);
ieegRoiEpoch = ieegEpoch(roiIds,:);
ieegRoiBase = ieegBase(roiIds,:);


[ieegRoiEpochSubj,ieegRoiEpochSubjAll,ieegRoiEpochSubjElec] = extractMeanControlSubject(ieegRoiEpoch,channelNameEpoch(roiIds));
[ieegRoiBaseSubj,ieegRoiBaseSubjAll,ieegBaseEpochSubjElec] = extractMeanControlSubject(ieegRoiBase,channelNameEpoch(roiIds));

pmask = zeros(1,size(ieegRoiBase,2));
[zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegRoiEpochSubjAll,ieegRoiBaseSubjAll,1000,1,1.645);

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
visTimePlot(timeRange,ieegRoiEpochSubjAll,colval = roiColor, labels = xlab )
hold on;
visTimePlot(timeRange,ieegRoiBaseSubjAll,colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGamma(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',...
    MarkerEdgeColor=roiColor,MarkerFaceColor=roiColor);
groupingIds(roiIds) = 1;
title(['pSTG ' num2str(sum(roiIds))])
ieegRoiEpochScaled = ieegRoiEpoch.*pmask; [maxAcc,maxId] = max(ieegRoiEpochScaled');
timeMax = timeGamma(maxId);
timeMaxRoi{8} = timeMax;
maxId = [];
timeCentroidRoi = [];
for iElec = 1:size(ieegRoiEpochScaled,1)
    timeCentroidRoi = [timeCentroidRoi calculateTimeCentroid(timeGamma,ieegRoiEpochScaled(iElec,:))];
    maxPoints = findchangepts(ieegRoiEpochScaled(iElec,:),MaxNumChanges=4,Statistic="mean");
    if(~isempty(maxPoints))
        maxId = [maxId maxPoints(1)];
    end
end
timeCentroidRoi = timeCentroidRoi(timeCentroidRoi>-0.9);
timeCentroid{8} = timeCentroidRoi;
timeChangePoint{8} = timeGamma(maxId);
if(showBoxPlot&&sum(pmask)>0)
    boxchart(boxplace*ones(size(time_centroid)),...
        double(time_centroid),'orientation','horizontal','BoxFaceColor',...
        roiColor*0.75,'BoxWidth',boxStd, 'BoxLineColor',...
        roiColor,'MarkerStyle','none')
end
ylim(yrange)


subplot(p(1),p(2),9)
bnLabels = extractBNlabels('sts',hemisphere)
roiIds = extractRoiIds(Subject, channelNameEpoch,bnLabels);
roiColor = colors(8,:);
ieegRoiEpoch = ieegEpoch(roiIds,:);
ieegRoiBase = ieegBase(roiIds,:);
 [ieegRoiEpochSubj,ieegRoiEpochSubjAll,ieegRoiEpochSubjElec] = extractMeanControlSubject(ieegRoiEpoch,channelNameEpoch(roiIds));
[ieegRoiBaseSubj,ieegRoiBaseSubjAll,ieegBaseEpochSubjElec] = extractMeanControlSubject(ieegRoiBase,channelNameEpoch(roiIds));

pmask = zeros(1,size(ieegRoiBase,2));
[zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegRoiEpochSubjAll,ieegRoiBaseSubjAll,1000,1,1.645);

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
visTimePlot(timeRange,ieegRoiEpochSubjAll,colval = roiColor, labels = xlab )
hold on;
visTimePlot(timeRange,ieegRoiBaseSubjAll,colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGamma(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',...
    MarkerEdgeColor=roiColor,MarkerFaceColor=roiColor);
groupingIds(roiIds) = 8;
title(['STS ' num2str(sum(roiIds))])
ieegRoiEpochScaled = ieegRoiEpoch.*pmask; [maxAcc,maxId] = max(ieegRoiEpochScaled');
timeMax = timeGamma(maxId);
timeMaxRoi{9} = timeMax;
maxId = [];
timeCentroidRoi = [];
for iElec = 1:size(ieegRoiEpochScaled,1)
    timeCentroidRoi = [timeCentroidRoi calculateTimeCentroid(timeGamma,ieegRoiEpochScaled(iElec,:))];
    maxPoints = findchangepts(ieegRoiEpochScaled(iElec,:),MaxNumChanges=4,Statistic="mean");
    if(~isempty(maxPoints))
        maxId = [maxId maxPoints(1)];
    end
end
timeCentroidRoi = timeCentroidRoi(timeCentroidRoi>-0.9);
timeCentroid{9} = timeCentroidRoi;
timeChangePoint{9} = timeGamma(maxId);
if(showBoxPlot&&sum(pmask)>0)
    boxchart(boxplace*ones(size(time_centroid)),...
        double(time_centroid),'orientation','horizontal','BoxFaceColor',...
        roiColor*0.75,'BoxWidth',boxStd, 'BoxLineColor',...
        roiColor,'MarkerStyle','none')
end
ylim(yrange)


subplot(p(1),p(2),10)
bnLabels = extractBNlabels('heschl',hemisphere)
roiIds = extractRoiIds(Subject, channelNameEpoch,bnLabels);
roiColor = colors(7,:);
ieegRoiEpoch = ieegEpoch(roiIds,:);
ieegRoiBase = ieegBase(roiIds,:);
 [ieegRoiEpochSubj,ieegRoiEpochSubjAll,ieegRoiEpochSubjElec] = extractMeanControlSubject(ieegRoiEpoch,channelNameEpoch(roiIds));
[ieegRoiBaseSubj,ieegRoiBaseSubjAll,ieegBaseEpochSubjElec] = extractMeanControlSubject(ieegRoiBase,channelNameEpoch(roiIds));

pmask = zeros(1,size(ieegRoiBase,2));
[zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegRoiEpochSubjAll,ieegRoiBaseSubjAll,1000,1,1.645);

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
visTimePlot(timeRange,ieegRoiEpochSubjAll,colval = roiColor, labels = xlab )
hold on;
visTimePlot(timeRange,ieegRoiBaseSubjAll,colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGamma(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',...
    MarkerEdgeColor=roiColor,MarkerFaceColor=roiColor);
groupingIds(roiIds) = 7;
title(['PAC ' num2str(sum(roiIds))])
ieegRoiEpochScaled = ieegRoiEpoch.*pmask; [maxAcc,maxId] = max(ieegRoiEpochScaled');
timeMax = timeGamma(maxId);
timeMaxRoi{10} = timeMax;
maxId = [];
timeCentroidRoi = [];
for iElec = 1:size(ieegRoiEpochScaled,1)
    timeCentroidRoi = [timeCentroidRoi calculateTimeCentroid(timeGamma,ieegRoiEpochScaled(iElec,:))];
    maxPoints = findchangepts(ieegRoiEpochScaled(iElec,:),MaxNumChanges=4,Statistic="mean");
    if(~isempty(maxPoints))
        maxId = [maxId maxPoints(1)];
    end
end
timeCentroidRoi = timeCentroidRoi(timeCentroidRoi>-0.9);
timeCentroid{10} = timeCentroidRoi;
timeChangePoint{10} = timeGamma(maxId);
if(showBoxPlot&&sum(pmask)>0)
    boxchart(boxplace*ones(size(time_centroid)),...
        double(time_centroid),'orientation','horizontal','BoxFaceColor',...
        roiColor*0.75,'BoxWidth',boxStd, 'BoxLineColor',...
        roiColor,'MarkerStyle','none')
end
ylim(yrange)


sgtitle('')


%% Heschl regions
hemisphere = 'lh';
figure;
bnLabels = extractBNlabels('triangular',hemisphere)
roiIds = extractRoiIds(Subject, channelNameEpoch,bnLabels);
% roiIds = find(roiIds);
% roiIds = roiIds(1:floor(length(roiIds)/2));
%roiIds = roiIds(floor(length(roiIds)/2)+1:end);
roiColor = colors(11,:);
ieegRoiEpoch = ieegEpoch(roiIds,:);
ieegRoiBase = ieegBase(roiIds,:);
% ieegRoiEpochNoNorm = ieegEpochBga(roiIds,:);
% ieegRoiBaseNoNorm = ieegBaseBga(roiIds,:);
[ieegRoiEpochSubj,ieegRoiEpochSubjAll,ieegRoiEpochSubjElec] = extractMeanControlSubject(ieegRoiEpoch,channelNameEpoch(roiIds));
[ieegRoiBaseSubj,ieegRoiBaseSubjAll,ieegBaseEpochSubjElec] = extractMeanControlSubject(ieegRoiBase,channelNameEpoch(roiIds));
figure;
pmask = zeros(1,size(ieegRoiBase,2));
[zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegRoiEpochSubjAll,ieegRoiBaseSubjAll,1000,1,1.645);

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
visTimePlot(timeRange,ieegRoiEpochSubjAll,colval = roiColor, labels = xlab )
hold on;
visTimePlot(timeRange,ieegRoiBaseSubjAll,colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGamma(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',...
    MarkerEdgeColor=roiColor,MarkerFaceColor=roiColor);

title(['Pars Triangularis ' num2str(sum(roiIds))])
% [p] = numSubplots(length(ieegBaseEpochSubjElec))
% 
% subjNames = extractSubjectPerChannel(channelNameEpoch(roiIds));
% for iSubj = 1:length(ieegBaseEpochSubjElec)
%    subplot(p(1),p(2),iSubj)
%     
%      [zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegRoiEpochSubjElec{iSubj},ieegBaseEpochSubjElec{iSubj},1000,1,1.645);
%     pmask = zeros(1,size(ieegRoiBase,2));
%     
%     for iClust=1:length(actClust.Size)
%         if actClust.Size{iClust}>actClust.perm95
%             pmask(actClust.Start{iClust}: ...
%                 actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
%                 = 1;
%         end
%     end
%     if(size(ieegRoiEpochSubjElec{iSubj},1)==1)
%         plot(timeGamma, ieegRoiEpochSubjElec{iSubj})
%         formatTicks(gca)
%     else
%     visTimePlot(timeRange,ieegRoiEpochSubjElec{iSubj},colval = roiColor, labels = xlab )
%     hold on;
%     visTimePlot(timeRange,ieegBaseEpochSubjElec{iSubj},colval = [0.5 0.5 0.5], labels = xlab)
%     scatter(timeGamma(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',...
%         MarkerEdgeColor=roiColor,MarkerFaceColor=roiColor);
%     end
%     groupingIds(roiIds) = 7;
%     title(['Heschl  ' num2str(sum(roiIds)) ' subject ' subjNames{iSubj}])
% end
%%
for iSubj = 1:length(ieegHGStruct)
    for iChan = 1:size(ieegHGStruct(iSubj).ieegHGNorm.data,1)
        figure; 
        imagesc(timeGamma,[],squeeze(ieegHGStruct(iSubj).ieegHGNorm.data(iChan,:,:))); 
        caxis([0 2]);
    end
end
%%
load colors.mat
groupingIds =zeros(1,length(channelNameBase));
hemisphere = 'lh';
boxplace = 1;
ylineval = -0.1;
yrange = [-0.15 1.4];
boxStd = 0.1;
showBoxPlot = 1;


figure;
hold on;



bnLabels = extractBNlabels('ifg',hemisphere)
roiIds = extractRoiIds(Subject, channelNameEpoch,bnLabels);
roiColor(1,:) = colors(4,:);
ieegRoiEpoch = ieegEpoch(roiIds,:);
ieegRoiBase = ieegBase(roiIds,:);
groupingIds(roiIds) = 4;
visTimePlot(timeRange,ieegRoiEpoch,colval = roiColor(1,:), labels = xlab    )

bnLabels = extractBNlabels('rmfg',hemisphere)
roiIds = extractRoiIds(Subject, channelNameEpoch,bnLabels);
roiColor(2,:) = colors(3,:);
ieegRoiEpoch = ieegEpoch(roiIds,:);
ieegRoiBase = ieegBase(roiIds,:);
groupingIds(roiIds) = 3; 
visTimePlot(timeRange,ieegRoiEpoch,colval = roiColor(2,:), labels = xlab )

bnLabels = extractBNlabels('cmfg',hemisphere)
roiIds = extractRoiIds(Subject, channelNameEpoch,bnLabels);
roiColor(3,:) = colors(10,:);
ieegRoiEpoch = ieegEpoch(roiIds,:);
ieegRoiBase = ieegBase(roiIds,:);
groupingIds(roiIds) = 10;
visTimePlot(timeRange,ieegRoiEpoch,colval = roiColor(3,:), labels = xlab )

figure;

bnLabels = extractBNlabels('smc',hemisphere)
roiIds = extractRoiIds(Subject, channelNameEpoch,bnLabels);
roiColor(4,:) = colors(5,:);
ieegRoiEpoch = ieegEpoch(roiIds,:);
ieegRoiBase = ieegBase(roiIds,:);
groupingIds(roiIds) = 5;
visTimePlot(timeRange,ieegRoiEpoch,colval = roiColor(4,:), labels = xlab )

bnLabels = extractBNlabels('ipc',hemisphere)
roiIds = extractRoiIds(Subject, channelNameEpoch,bnLabels);
roiColor(5,:) = colors(2,:);
ieegRoiEpoch = ieegEpoch(roiIds,:);
ieegRoiBase = ieegBase(roiIds,:);
groupingIds(roiIds) = 2;
visTimePlot(timeRange,ieegRoiEpoch,colval = roiColor(5,:), labels = xlab )
hold on;
visTimePlot(timeRange,ieegRoiBase,colval = [0.5 0.5 0.5], labels = xlab)

figure;

bnLabels = extractBNlabels('insula',hemisphere)
roiIds = extractRoiIds(Subject, channelNameEpoch,bnLabels);
roiColor(6,:) = colors(6,:);
ieegRoiEpoch = ieegEpoch(roiIds,:);
ieegRoiBase = ieegBase(roiIds,:);
groupingIds(roiIds) = 6;
visTimePlot(timeRange,ieegRoiEpoch,colval = roiColor(6,:), labels = xlab )

bnLabels = extractBNlabels('astg',hemisphere)
roiIds = extractRoiIds(Subject, channelNameEpoch,bnLabels);
roiColor(7,:) = colors(9,:);
ieegRoiEpoch = ieegEpoch(roiIds,:);
ieegRoiBase = ieegBase(roiIds,:);
groupingIds(roiIds) = 9; 
visTimePlot(timeRange,ieegRoiEpoch,colval = roiColor(7,:), labels = xlab )

bnLabels = extractBNlabels('pstg',hemisphere)
roiIds = extractRoiIds(Subject, channelNameEpoch,bnLabels);
roiColor(8,:) = colors(1,:);
ieegRoiEpoch = ieegEpoch(roiIds,:);
ieegRoiBase = ieegBase(roiIds,:);
groupingIds(roiIds) = 1;
visTimePlot(timeRange,ieegRoiEpoch,colval = roiColor(8,:), labels = xlab )

bnLabels = extractBNlabels('sts',hemisphere)
roiIds = extractRoiIds(Subject, channelNameEpoch,bnLabels);
roiColor(9,:) = colors(8,:);
ieegRoiEpoch = ieegEpoch(roiIds,:);
ieegRoiBase = ieegBase(roiIds,:);
groupingIds(roiIds) = 8; 
visTimePlot(timeRange,ieegRoiEpoch,colval = roiColor(9,:), labels = xlab )

bnLabels = extractBNlabels('heschl',hemisphere)
roiIds = extractRoiIds(Subject, channelNameEpoch,bnLabels);
ieegRoiEpoch = ieegEpoch(roiIds,:);
ieegRoiBase = ieegBase(roiIds,:);



roiColor(10,:) = colors(7,:); 
visTimePlot(timeRange,ieegRoiEpoch,colval = roiColor(10,:), labels = xlab )
groupingIds(roiIds) = 7;


% visTimePlot(timeRange,modelweightnorm(groupingIds==0,:),colval = colors(9,:))
%  groupingIds(groupingIds==0) = 9;

% xlabel('Time from response onset')
% ylabel('Beta')

%%

y = num2cell(1:numel(timeMaxRoi));
x = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], timeMaxRoi, y, 'UniformOutput', 0); % adding labels to the cells
X = vertcat(x{:});
w95 = 0.7193;
for i = 1:10
    medX = median(x{i});
end

figure;
boxplot(X(:,1), X(:,2),'symbol','','Colors','k','Whisker', w95)

%%

% Assuming timeMaxRoi is a cell array where each cell contains numeric data for each group
y = num2cell(1:numel(timeMaxRoi));
x = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], timeMaxRoi, y, 'UniformOutput', 0); % Add labels to the cells
X = vertcat(x{:});

w95 = 0.7193;

% Calculate the median for each group
medX = cellfun(@median, timeMaxRoi);

% Sort medians and get the sort order
[~, sortIdx] = sort(medX);

% Create a new mapping for the group labels based on the sorted order
% This step ensures the group labels are sorted according to the median values
sortedLabels = arrayfun(@(x) find(sortIdx == x), X(:,2), 'UniformOutput', false);
sortedLabelsMat = cell2mat(sortedLabels);

figure;
% Use the sorted labels for plotting
boxplot(X(:,1), sortedLabelsMat, 'symbol', '', 'Colors', 'k', 'Whisker', w95)

% Optionally, set custom x-tick labels to reflect the original group IDs in their sorted order
set(gca, 'XTick', 1:numel(timeMaxRoi_lh), 'XTickLabel', labels(sortIdx));

%% ROI histogram plot
labels = {'IFG','cMFG','rMFG','SMC','IPC','Insula','aSTG','pSTG','STS','PAC'};
roiColors = [colors(4,:); colors(3,:); colors(10,:); colors(5,:); colors(2,:);...
    colors(6,:); colors(9,:); colors(1,:); colors(8,:); colors(7,:);];
medX = cellfun(@median, timeCentroid);

% Sort medians and get the sort order
[~, sortIdx] = sort(medX);
    
timeMaxRoi2plot = timeCentroid; 
labels2plot = labels;

plotRoiHistograms(timeMaxRoi2plot,labels2plot,roiColors)

%% Group histogram plot

timeGroup{1} = [timeCentroid{1:3}];
timeGroup{2} = [timeCentroid{4:5}];
timeGroup{3} = [timeCentroid{6:10}];

labels2plot = {'planning','articulation','monitoring'};
 plotRoiHistograms(timeGroup,labels2plot,lines(3))


 [h,p,ci,stats] = ttest2( timeGroup{2},timeGroup{1})
 [h,p,ci,stats] = ttest2( timeGroup{3},timeGroup{2})
 [h,p,ci,stats] = ttest2( timeGroup{3},timeGroup{1})
%% plot ROI - syllable phoneme


timeCpSpeech = [];
for iRoi = 1:10
    if(~isfield(accResults{iRoi}{1},'timeCP1'))
        timeCpSpeech{iRoi,1} = nan;
    else
        timeCpSpeech{iRoi,1} = rmoutliers(accResults{iRoi}{1}.timeCP1,'percentiles', [5 95]) ;
    end

    if(~isfield(accResults{iRoi}{2},'timeCP1'))
        timeCpSpeech{iRoi,2} = nan;
    else
        timeCpSpeech{iRoi,2} = rmoutliers(accResults{iRoi}{2}.timeCP1,'percentiles', [5 95]) ;
    end
    
end
labels = {'IFG','rMFG','cMFG','SMC','IPC','Insula','aSTG','pSTG','STS','PAC'};
colvals = lines(8);
% plotRoiHistograms(timeCpSyllable,labels,repmat(colvals(4,:),10,1))
% plotRoiHistograms(timeCpPhoneme,labels,repmat(colvals(5,:),10,1))
colvalSpeech = [];
colvalSpeech(:,1,:) = repmat(colvals(4,:),10,1);
colvalSpeech(:,2,:) = repmat(colvals(5,:),10,1);
%  plotRoiHistogramsSpeech(timeCpSpeech([1 3],:), labels([1 3]), colvalSpeech([1 3],:,:))
%   plotRoiHistogramsSpeech(timeCpSpeech(4:5,:), labels(4:5), colvalSpeech(4:5,:,:))
   plotRoiHistogramsSpeech(timeCpSpeech(1:10,:), labels(1:10), colvalSpeech(1:10,:,:))
plotRoiHistogramsSpeechBoxPlot(timeCpSpeech(1:10,:), labels(1:10), colvalSpeech(1:10,:,:))

   %% plot ROI group- syllable phoneme 


timeCpSpeech = [];
for iRoi = 1:3
    if(~isfield(accResultsGroup{iRoi}{1},'timeCP1'))
        timeCpSpeech{iRoi,1} = nan;
    else
        timeCpSpeech{iRoi,1} = rmoutliers(accResultsGroup{iRoi}{1}.timeCP1,'percentiles', [1 99]) ;
        timeCpSpeech{iRoi,1} = timeCpSpeech{iRoi,1}(timeCpSpeech{iRoi,1}<=0.5);
    end

    if(~isfield(accResultsGroup{iRoi}{2},'timeCP1'))
        timeCpSpeech{iRoi,2} = nan;
    else
        timeCpSpeech{iRoi,2} = rmoutliers(accResultsGroup{iRoi}{2}.timeCP1,'percentiles', [1 99]) ;
        timeCpSpeech{iRoi,2} = timeCpSpeech{iRoi,2}(timeCpSpeech{iRoi,2}<=0.5);
    end
    
end
labels = {'Planning','Articulation','Monitoring'};
colvals = lines(8);
% plotRoiHistograms(timeCpSyllable,labels,repmat(colvals(4,:),10,1))
% plotRoiHistograms(timeCpPhoneme,labels,repmat(colvals(5,:),10,1))
colvalSpeech = [];
colvalSpeech(:,1,:) = repmat(colvals(4,:),3,1);
colvalSpeech(:,2,:) = repmat(colvals(5,:),3,1);
%  plotRoiHistogramsSpeech(timeCpSpeech([1 3],:), labels([1 3]), colvalSpeech([1 3],:,:))
%   plotRoiHistogramsSpeech(timeCpSpeech(4:5,:), labels(4:5), colvalSpeech(4:5,:,:))
   plotRoiHistogramsSpeech(timeCpSpeech(1:3,:), labels(1:3), colvalSpeech(1:3,:,:))

%    [p,h] = permutationTest(timeCpSpeech{1,2}, timeCpSpeech{1,1},1000,'sidedness','larger' )
%    
% 
%    [p,h] = ranksum(timeCpSpeech{3,2}, timeCpSpeech{3,1},'Tail','right' )
%    
   

timeCpSpeechDiff1 = timeCpSpeech{1,2} - timeCpSpeech{1,1}';
timeCpSpeechDiff2 = timeCpSpeech{2,2} - timeCpSpeech{2,1}';
timeCpSpeechDiff3 = timeCpSpeech{3,2} - timeCpSpeech{3,1}';


[h,p,ci,stats] = ttest2(timeCpSpeechDiff1(:)' ,timeCpSpeechDiff2(:)' )
[h,p,ci,stats] = ttest2(timeCpSpeechDiff2(:)' ,timeCpSpeechDiff3(:)')
[h,p,ci,stats] = ttest2(timeCpSpeechDiff1(:)' ,timeCpSpeechDiff3(:)')

figure;
mean1 = nanmean(timeCpSpeechDiff1(:)');
mean2 = nanmean(timeCpSpeechDiff2(:)');
mean3 = nanmean(timeCpSpeechDiff3(:)');

std1 = nanstd(timeCpSpeechDiff1(:)')./sqrt(sum(~isnan(timeCpSpeechDiff1(:)')));
std2 = nanstd(timeCpSpeechDiff2(:)')./sqrt(sum(~isnan(timeCpSpeechDiff2(:)')));
std3 = nanstd(timeCpSpeechDiff3(:)')./sqrt(sum(~isnan(timeCpSpeechDiff3(:)')));

x = 1:3;
data = [mean1 mean2 mean3]';
errhigh = [std1 std2 std3];
errlow  = [std1 std2 std3];
bar(x,data)                

hold on

er = errorbar(x,data,errlow,errhigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';

hold off
formatTicks(gca)  

   %% plot ROI - phoneme sequence


timeCpPhoneme = [];
for iRoi = 1:10
    if(~isfield(accResultsPhoneme{iRoi}{1},'timeCP1'))
        timeCpPhoneme{iRoi,1} = nan;
    else
        timeCpPhoneme{iRoi,1} = (accResultsPhoneme{iRoi}{1}.timeCP1) ;
    end

    if(~isfield(accResultsPhoneme{iRoi}{2},'timeCP1'))
        timeCpPhoneme{iRoi,2} = nan;
    else
        timeCpPhoneme{iRoi,2} = (accResultsPhoneme{iRoi}{2}.timeCP1) ;
    end

    if(~isfield(accResultsPhoneme{iRoi}{3},'timeCP1'))
        timeCpPhoneme{iRoi,3} = nan;
    else
        timeCpPhoneme{iRoi,3} = rmoutliers(accResultsPhoneme{iRoi}{3}.timeCP1,'percentiles', [1 99]) ;
    end
    
end
labels = {'IFG','rMFG','cMFG','SMC','IPC','Insula','aSTG','pSTG','STS','PAC'};
%labels = {'Planning','Articulation','Monitoring'};
colvals = lines(8);
% plotRoiHistograms(timeCpSyllable,labels,repmat(colvals(4,:),10,1))
% plotRoiHistograms(timeCpPhoneme,labels,repmat(colvals(5,:),10,1))
colvalSpeech = [];
colvalSpeech(:,1,:) = repmat(colvals(1,:),10,1);
colvalSpeech(:,2,:) = repmat(colvals(2,:),10,1);
colvalSpeech(:,3,:) = repmat(colvals(3,:),10,1);
 %plotRoiHistogramsSpeech(timeCpSpeech([1 3],:), labels([1 3]), colvalSpeech([1 3],:,:))
%  plotRoiHistogramsSpeech(timeCpPhoneme(2:3,:), labels(2:3), colvalSpeech(2:3,:,:))
   plotRoiHistogramsSpeech(timeCpPhoneme(1:10,:), labels(1:10), colvalSpeech(1:10,:,:))

   
   %% plot ROI - phoneme sequence


timeCpPhoneme = [];
for iRoi = 1:3
    if(~isfield(accResultsPhoneme{iRoi}{1},'timeCP1'))
        timeCpPhoneme{iRoi,1} = nan;
    else
        timeCpPhoneme{iRoi,1} = rmoutliers(accResultsPhoneme{iRoi}{1}.timeCP1,'percentiles', [1 99]) ;
    end

    if(~isfield(accResultsPhoneme{iRoi}{2},'timeCP1'))
        timeCpPhoneme{iRoi,2} = nan;
    else
        timeCpPhoneme{iRoi,2} = rmoutliers(accResultsPhoneme{iRoi}{2}.timeCP1,'percentiles', [1 99]) ;
    end

    if(~isfield(accResultsPhoneme{iRoi}{3},'timeCP1'))
        timeCpPhoneme{iRoi,3} = nan;
    else
        timeCpPhoneme{iRoi,3} = rmoutliers(accResultsPhoneme{iRoi}{3}.timeCP1,'percentiles', [1 99]) ;
    end
    
end
%labels = {'IFG','rMFG','cMFG','SMC','IPC','Insula','aSTG','pSTG','STS','PAC'};
labels = {'Planning','Articulation','Monitoring'};
colvals = lines(8);
% plotRoiHistograms(timeCpSyllable,labels,repmat(colvals(4,:),10,1))
% plotRoiHistograms(timeCpPhoneme,labels,repmat(colvals(5,:),10,1))
colvalSpeech = [];
colvalSpeech(:,1,:) = repmat(colvals(1,:),3,1);
colvalSpeech(:,2,:) = repmat(colvals(2,:),3,1);
colvalSpeech(:,3,:) = repmat(colvals(3,:),3,1);
 %plotRoiHistogramsSpeech(timeCpSpeech([1 3],:), labels([1 3]), colvalSpeech([1 3],:,:))
  plotRoiHistogramsSpeech(timeCpPhoneme(2:3,:), labels(2:3), colvalSpeech(2:3,:,:))
   %plotRoiHistogramsSpeech(timeCpPhoneme(1:10,:), labels(1:10), colvalSpeech(1:10,:,:))
[p,h] = ranksum( timeCpPhoneme{2,2},timeCpPhoneme{2,1},'Tail','right')
[p,h] = ranksum( timeCpPhoneme{2,3},timeCpPhoneme{2,2},'Tail','right')
[p,h] = ranksum( timeCpPhoneme{2,3},timeCpPhoneme{2,1},'Tail','right')


[p,h] = ranksum( timeCpPhoneme{3,2},timeCpPhoneme{3,1},'Tail','right')
[p,h] = ranksum( timeCpPhoneme{3,3},timeCpPhoneme{3,2},'Tail','right')
[p,h] = ranksum( timeCpPhoneme{3,3},timeCpPhoneme{3,1},'Tail','right')



[h,p,ci,stats] = ttest2( timeCpPhoneme{2,2},timeCpPhoneme{2,1},'Tail','right')
[h,p,ci,stats] = ttest2( timeCpPhoneme{2,3},timeCpPhoneme{2,2},'Tail','right')
[h,p,ci,stats] = ttest2( timeCpPhoneme{2,3},timeCpPhoneme{2,1},'Tail','right')


[h,p,ci,stats] = ttest2( timeCpPhoneme{3,2},timeCpPhoneme{3,1},'Tail','right')
[h,p,ci,stats] = ttest2( timeCpPhoneme{3,3},timeCpPhoneme{3,2},'Tail','right')
[h,p,ci,stats] = ttest2( timeCpPhoneme{3,3},timeCpPhoneme{3,1},'Tail','right')
%% Group histogram plot

timeGroup{1} = [timeSyllableCentroid{1:3}];
timeGroup{2} = [timeSyllableCentroid{4:5}];
timeGroup{3} = [timeSyllableCentroid{6:10}];

labels2plot = {'pre-frontal','sensorimotor','temporal-monitoring'};
 plotRoiHistograms(timeGroup,labels2plot,lines(3))

%% 

% Example cell array with three groups
data = timeCentroid;

% Convert cell array data to a format suitable for ANOVA
group = [];
values = [];

for i = 1:numel(data)
    group = [group; i*ones(size(data{i},2), 1)]; % Assign group numbers
    values = [values; data{i}']; % Concatenate all values into a single array
end
group3 = group;
group3(group<=3) = 1;
group3(group>3&group<=5) = 2;
group3(group>5) = 3;

% Perform one-way ANOVA
[p, tbl, stats] = anova1(values, group3, 'off');  % 'off' turns off the ANOVA table display

% Display results
disp('ANOVA Table:');
disp(tbl);
fprintf('p-value = %.4f\n', p);

% If you want to investigate significant differences further with post-hoc tests
if p < 0.05
    figure;
    [c,m,h,gnames] = multcompare(stats);
    title('Multiple Comparisons of Means');
end

%%

load('pooledSubject_Phoneme_Sequencing_Production_car_z_score_prodelecs_data_v3.mat')

ieegStruct = ieegHGStruct; 
for iSubject = 1:length(ieegHGStruct)
    iSubject
     ieegStruct2process = removeZScore(ieegStruct(iSubject).ieegHGNorm,ieegStruct(iSubject).normFactor);
     ieegStruct2process = applyNormalization(ieegStruct2process,ieegStruct(iSubject).normFactor,3);
     ieegHGStruct(iSubject).ieegHGNorm = ieegStruct2process;
end


save(fullfile(['pooledSubject_' Task.Name '_Production_car_bga' ...
'_prodelecs_data_v3.mat']),...
'ieegHGStruct','trialInfoStruct','maxTrial',...
'timeEpoch','-v7.3');
%       

%%
% Generate some random data
data = randn(1000, 1);
    
% Create histogram and fit a normal distribution
figure;
h = histfit(data, 50, 'normal'); % 50 bins and normal distribution

% Make the histogram bars invisible
h(1).FaceColor = 'none'; % Set face color to none to hide the bars
h(1).EdgeColor = 'none';

% Get the x and y data from the fitted curve
x_fit = h(2).XData;
y_fit = h(2).YData;

% Remove the line plot used in histfit for fitting
delete(h(2));

% Use a patch to create a filled curve for the fit
hold on;
patch([x_fit, fliplr(x_fit)], [zeros(1,length(y_fit)), fliplr(y_fit)], 'r', 'FaceAlpha', 0.3);

% Labels and titles
xlabel('Data values');
ylabel('Frequency');
title('Fitted Distribution Without Histogram Bars');

% Show grid
grid on;
hold off;
%% Subject specific analysis - power traces
subjNames = {'D35','D39','D41','D49','D59','D64','D68','D79','D81','D88','D93','D96'};
%subjNames = extractSubjectPerChannel(channelNameEpoch);    
for iSubj = 1:length(subjNames)
    iSubj
    channelIds = extractChannelPerSubject(channelNameEpoch,subjNames{iSubj});
    channelNameEpoch(channelIds)

    ieegEpochSubj = ieegEpoch(channelIds,:);
    ieegBaseSubj = ieegBase(channelIds,:);
    figure;
    bnLabels = [extractBNlabels('ifg',hemisphere) extractBNlabels('rmfg',hemisphere) extractBNlabels('cmfg',hemisphere) ];
    roiIds = extractRoiIds(Subject, channelNameEpoch(channelIds),bnLabels);    
    if(sum(roiIds)>1)
    ieegRoiEpoch = ieegEpochSubj(roiIds,:);
    ieegRoiBase = ieegBaseSubj(roiIds,:);
    
    pmask = zeros(1,size(ieegRoiBase,2));
    [zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegRoiEpoch,ieegRoiBase,1000,1,1.645);
    
    for iClust=1:length(actClust.Size)
        if actClust.Size{iClust}>actClust.perm95
            pmask(actClust.Start{iClust}: ...
                actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
                = 1;
        end
    end    
    ieegRoiEpochScaled = ieegRoiEpoch.*pmask;

    timeCentroidRoi = [];
    for iElec = 1:size(ieegRoiEpochScaled,1)
        timeCentroidRoi = [timeCentroidRoi calculateTimeCentroid(timeGamma,ieegRoiEpochScaled(iElec,:))];   
        
    end
    timeCentroidRoi = timeCentroidRoi(timeCentroidRoi>-0.9);
    timeCentroid{iSubj,1} = timeCentroidRoi;
    
    visTimePlot(timeRange,ieegRoiEpoch,colval = [0.05 0.03 0.53], labels = xlab )
    hold on;
    yline(0, '-','','LineWidth',1, 'Color','k');

    scatter(timeGamma(logical(pmask)),-0.1.*ones(1,sum(pmask)),'filled',...
    MarkerEdgeColor=[0.05 0.03 0.53],MarkerFaceColor=[0.05 0.03 0.53]);
    end
    bnLabels = [extractBNlabels('smc',hemisphere) extractBNlabels('ipc',hemisphere)];
    roiIds = extractRoiIds(Subject, channelNameEpoch(channelIds),bnLabels);    
    if(sum(roiIds)>1)
    ieegRoiEpoch = ieegEpochSubj(roiIds,:);
    ieegRoiBase = ieegBaseSubj(roiIds,:);
    
    pmask = zeros(1,size(ieegRoiBase,2));
    [zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegRoiEpoch,ieegRoiBase,1000,1,1.645);
    
    for iClust=1:length(actClust.Size)
        if actClust.Size{iClust}>actClust.perm95
            pmask(actClust.Start{iClust}: ...
                actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
                = 1;
        end
    end    
    ieegRoiEpochScaled = ieegRoiEpoch.*pmask;
    timeCentroidRoi = [];
    for iElec = 1:size(ieegRoiEpochScaled,1)
        timeCentroidRoi = [timeCentroidRoi calculateTimeCentroid(timeGamma,ieegRoiEpochScaled(iElec,:))];   
        
    end
    timeCentroidRoi = timeCentroidRoi(timeCentroidRoi>-0.9);
    timeCentroid{iSubj,2} = timeCentroidRoi;
    visTimePlot(timeRange,ieegRoiEpoch,colval = [0.8 0.28 0.47], labels = xlab )
    hold on;
    yline(0, '-','','LineWidth',1, 'Color','k');

    scatter(timeGamma(logical(pmask)),-0.2.*ones(1,sum(pmask)),'filled',...
    MarkerEdgeColor=[0.8 0.28 0.47],MarkerFaceColor=[0.8 0.28 0.47]);
    end

    bnLabels = [extractBNlabels('astg',hemisphere) extractBNlabels('pstg',hemisphere) extractBNlabels('insula',hemisphere) extractBNlabels('heschl',hemisphere) extractBNlabels('sts',hemisphere)];
    roiIds = extractRoiIds(Subject, channelNameEpoch(channelIds),bnLabels);    
    if(sum(roiIds)>1)
    ieegRoiEpoch = ieegEpochSubj(roiIds,:);
    ieegRoiBase = ieegBaseSubj(roiIds,:);
    
    pmask = zeros(1,size(ieegRoiBase,2));
    [zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegRoiEpoch,ieegRoiBase,1000,1,1.645);
    
    for iClust=1:length(actClust.Size)
        if actClust.Size{iClust}>actClust.perm95
            pmask(actClust.Start{iClust}: ...
                actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
                = 1;
        end
    end    
    ieegRoiEpochScaled = ieegRoiEpoch.*pmask;
    timeCentroidRoi = [];
    for iElec = 1:size(ieegRoiEpochScaled,1)
        timeCentroidRoi = [timeCentroidRoi calculateTimeCentroid(timeGamma,ieegRoiEpochScaled(iElec,:))];   
        
    end
    timeCentroidRoi = timeCentroidRoi(timeCentroidRoi>-0.9);
    timeCentroid{iSubj,3} = timeCentroidRoi;
    visTimePlot(timeRange,ieegRoiEpoch,colval = [0.13 0.57 0.55], labels = xlab )
    hold on;
    yline(0, '-','','LineWidth',1, 'Color','k');

    scatter(timeGamma(logical(pmask)),-0.3.*ones(1,sum(pmask)),'filled',...
    MarkerEdgeColor=[0.13 0.57 0.55],MarkerFaceColor=[0.13 0.57 0.55]);
    end
    title(subjNames{iSubj})
end
timeCentMedian = []
timeCentExecute = [];
for iSubj = 1:length(subjNames)
%     figure;
%     labels2plot = {'planning','articulation','monitoring'};
%     plotRoiHistograms( timeCentroid(iSubj,:),labels2plot,[0.05 0.03 0.53; 0.8 0.28 0.47; 0.13 0.57 0.55])
%     sgtitle(subjNames{iSubj})
    timeCentMedian(iSubj,1) = mean(timeCentroid{iSubj,1});
    timeCentMedian(iSubj,2) = mean(timeCentroid{iSubj,2});
    timeCentMedian(iSubj,3) = mean(timeCentroid{iSubj,3});
    if(timeCentMedian(iSubj,2)<timeCentMedian(iSubj,3)||isnan(timeCentMedian(iSubj,3)))
        timeCentExecute(iSubj,2) = timeCentMedian(iSubj,2);
    end
    if(timeCentMedian(iSubj,2)>timeCentMedian(iSubj,3)||isnan(timeCentMedian(iSubj,2)))
        timeCentExecute(iSubj,2) = timeCentMedian(iSubj,3);
    end
end

figure;
scatter(timeCentMedian(:,1)', 1:12);
hold on;
scatter(timeCentMedian(:,2)', 1:12);
scatter(timeCentMedian(:,3)', 1:47);

figure;
scatter(timeCentMedian(:,1)', 1:12);
hold on;
scatter(timeCentExecute(:,2)', 1:12);


%% Subject specific analysis - power traces
%subjNames = {'D35','D39','D41','D49','D59','D64','D68','D79','D81','D88','D93','D96'};
subjNames = extractSubjectPerChannel(channelNameEpoch);    
for iSubj = 1:length(subjNames)
    iSubj
    channelIds = extractChannelPerSubject(channelNameEpoch,subjNames{iSubj});
    channelNameEpoch(channelIds)

    ieegEpochSubj = ieegEpoch(channelIds,:);
    ieegBaseSubj = ieegBase(channelIds,:);
    figure;
    bnLabels = [extractBNlabels('ifg',hemisphere) extractBNlabels('rmfg',hemisphere) extractBNlabels('cmfg',hemisphere) ];
    roiIds = extractRoiIds(Subject, channelNameEpoch(channelIds),bnLabels);    
    if(sum(roiIds)>1)
    ieegRoiEpoch = ieegEpochSubj(roiIds,:);
    ieegRoiBase = ieegBaseSubj(roiIds,:);
    
    pmask = zeros(1,size(ieegRoiBase,2));
    [zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegRoiEpoch,ieegRoiBase,1000,1,1.645);
    
    for iClust=1:length(actClust.Size)
        if actClust.Size{iClust}>actClust.perm95
            pmask(actClust.Start{iClust}: ...
                actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
                = 1;
        end
    end    
    ieegRoiEpochScaled = ieegRoiEpoch.*pmask;

    timeCentroidRoi = [];
    for iElec = 1:size(ieegRoiEpochScaled,1)
        timeCentroidRoi = [timeCentroidRoi calculateTimeCentroid(timeGamma,ieegRoiEpochScaled(iElec,:))];   
        
    end
    timeCentroidRoi = timeCentroidRoi(timeCentroidRoi>-0.9);
    timeCentroid{iSubj,1} = timeCentroidRoi;
    
    visTimePlot(timeRange,ieegRoiEpoch,colval = [0.05 0.03 0.53], labels = xlab )
    hold on;
    yline(0, '-','','LineWidth',1, 'Color','k');

    scatter(timeGamma(logical(pmask)),-0.1.*ones(1,sum(pmask)),'filled',...
    MarkerEdgeColor=[0.05 0.03 0.53],MarkerFaceColor=[0.05 0.03 0.53]);
    end
    bnLabels = [extractBNlabels('smc',hemisphere) extractBNlabels('ipc',hemisphere) extractBNlabels('astg',hemisphere) extractBNlabels('pstg',hemisphere) extractBNlabels('insula',hemisphere) extractBNlabels('heschl',hemisphere) extractBNlabels('sts',hemisphere)];
    roiIds = extractRoiIds(Subject, channelNameEpoch(channelIds),bnLabels);    
    if(sum(roiIds)>1)
    ieegRoiEpoch = ieegEpochSubj(roiIds,:);
    ieegRoiBase = ieegBaseSubj(roiIds,:);
    
    pmask = zeros(1,size(ieegRoiBase,2));
    [zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegRoiEpoch,ieegRoiBase,1000,1,1.645);
    
    for iClust=1:length(actClust.Size)
        if actClust.Size{iClust}>actClust.perm95
            pmask(actClust.Start{iClust}: ...
                actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
                = 1;
        end
    end    
    ieegRoiEpochScaled = ieegRoiEpoch.*pmask;
    timeCentroidRoi = [];
    for iElec = 1:size(ieegRoiEpochScaled,1)
        timeCentroidRoi = [timeCentroidRoi calculateTimeCentroid(timeGamma,ieegRoiEpochScaled(iElec,:))];   
        
    end
    timeCentroidRoi = timeCentroidRoi(timeCentroidRoi>-0.9);
    timeCentroid{iSubj,2} = timeCentroidRoi;
    visTimePlot(timeRange,ieegRoiEpoch,colval = [0.8 0.28 0.47], labels = xlab )
    hold on;
    yline(0, '-','','LineWidth',1, 'Color','k');

    scatter(timeGamma(logical(pmask)),-0.2.*ones(1,sum(pmask)),'filled',...
    MarkerEdgeColor=[0.8 0.28 0.47],MarkerFaceColor=[0.8 0.28 0.47]);
    end
% 
%     bnLabels = [extractBNlabels('astg',hemisphere) extractBNlabels('pstg',hemisphere) extractBNlabels('insula',hemisphere) extractBNlabels('heschl',hemisphere) extractBNlabels('sts',hemisphere)];
%     roiIds = extractRoiIds(Subject, channelNameEpoch(channelIds),bnLabels);    
%     if(sum(roiIds)>1)
%     ieegRoiEpoch = ieegEpochSubj(roiIds,:);
%     ieegRoiBase = ieegBaseSubj(roiIds,:);
%     
%     pmask = zeros(1,size(ieegRoiBase,2));
%     [zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegRoiEpoch,ieegRoiBase,1000,1,1.645);
%     
%     for iClust=1:length(actClust.Size)
%         if actClust.Size{iClust}>actClust.perm95
%             pmask(actClust.Start{iClust}: ...
%                 actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
%                 = 1;
%         end
%     end    
%     ieegRoiEpochScaled = ieegRoiEpoch.*pmask;
%     timeCentroidRoi = [];
%     for iElec = 1:size(ieegRoiEpochScaled,1)
%         timeCentroidRoi = [timeCentroidRoi calculateTimeCentroid(timeGamma,ieegRoiEpochScaled(iElec,:))];   
%         
%     end
%     timeCentroidRoi = timeCentroidRoi(timeCentroidRoi>-0.9);
%     timeCentroid{iSubj,3} = timeCentroidRoi;
%     visTimePlot(timeRange,ieegRoiEpoch,colval = [0.13 0.57 0.55], labels = xlab )
%     hold on;
%     yline(0, '-','','LineWidth',1, 'Color','k');
% 
%     scatter(timeGamma(logical(pmask)),-0.3.*ones(1,sum(pmask)),'filled',...
%     MarkerEdgeColor=[0.13 0.57 0.55],MarkerFaceColor=[0.13 0.57 0.55]);
%     end
%     title(subjNames{iSubj})
end
timeCentMedian = []
for iSubj = 1:length(subjNames)
%     figure;
%     labels2plot = {'planning','articulation','monitoring'};
%     plotRoiHistograms( timeCentroid(iSubj,:),labels2plot,[0.05 0.03 0.53; 0.8 0.28 0.47; 0.13 0.57 0.55])
%     sgtitle(subjNames{iSubj})
    timeCentMedian(iSubj,1) = mean(timeCentroid{iSubj,1});
    timeCentMedian(iSubj,2) = mean(timeCentroid{iSubj,2});
   % timeCentMedian(iSubj,3) = mean(timeCentroid{iSubj,3});

end

figure;
scatter(timeCentMedian(:,1)', 1:47);
hold on;
scatter(timeCentMedian(:,2)', 1:47);
%scatter(timeCentMedian(:,3)', 1:47);

%% Subject specific analysis - syllable

subjNames = {'D35','D39','D41','D49','D59','D64','D68','D79','D81','D88','D93','D96'};

for iSubj = 1:length(subjNames)
    iSubj
    channelIds = extractChannelPerSubject(channelNamePooled,subjNames{iSubj});
    channelNamePooled(channelIds)

    ieegMeanSubj = ieegMeanMod(:,channelIds,:);
    ieegShuffSubj = ieegShuffMod(:,channelIds,:);


    bnLabels = [extractBNlabels('ifg',hemisphere) extractBNlabels('rmfg',hemisphere) extractBNlabels('cmfg',hemisphere) ];
    roiIds = extractRoiIds(Subject, channelNamePooled(channelIds),bnLabels);

    
    
    trials2select = randperm(size(ieegMeanMod,1),groupNum);
    ieegMeanRoi = squeeze(mean(mean(ieegMeanSubj(trials2select,roiIds,:),2),1))';
    
    ieegShuffleRoi = [];
    for iTer = 1:1000
        trials2select = randperm(size(ieegMeanMod,1),groupNum);
        ieegShuffleRoi(iTer, :) = squeeze(mean(mean(ieegShuffSubj(trials2select,roiIds,:),2),1))';       
    end
    
     [zValsRawAct, pValsRaw, actClust] = timePermClusterAfterPerm(ieegMeanRoi,ieegShuffleRoi,1,1.645);
    pmask = zeros(1,size(ieegMeanRoi,2));
    % [~, ~, actClust]=timePermClusterZOptim(acctimeall,acctimeshuffle,10000,1,norminv(1-pVal2Cutoff));
    
    for iClust=1:length(actClust.Size)
        if actClust.Size{iClust}>actClust.perm95
            pmask(actClust.Start{iClust}: ...
                actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
                = 1;
        end
    end

    figure;
    visTimePlot(timeRange,squeeze((mean(ieegMeanSubj(:,roiIds,:),1))),colval = colors(1,:), labels = xlab )
    hold on;
    visTimePlot(timeRange,squeeze((mean(ieegShuffSubj(:,roiIds,:),1))),colval = [0.5 0.5 0.5], labels = xlab)
    scatter(timeGammaMod(logical(pmask)),-0.005.*ones(1,sum(pmask)),'filled',MarkerEdgeColor=colors(1,:),MarkerFaceColor=colors(1,:));


    bnLabels = [extractBNlabels('smc',hemisphere) extractBNlabels('ipc',hemisphere)];
    roiIds = extractRoiIds(Subject, channelNamePooled(channelIds),bnLabels);

    
    
    trials2select = randperm(size(ieegMeanMod,1),groupNum);
    ieegMeanRoi = squeeze(mean(mean(ieegMeanSubj(trials2select,roiIds,:),2),1))';
    
    ieegShuffleRoi = [];
    for iTer = 1:1000
        trials2select = randperm(size(ieegMeanMod,1),groupNum);
        ieegShuffleRoi(iTer, :) = squeeze(mean(mean(ieegShuffSubj(trials2select,roiIds,:),2),1))';       
    end
    
     [zValsRawAct, pValsRaw, actClust] = timePermClusterAfterPerm(ieegMeanRoi,ieegShuffleRoi,1,1.645);
    pmask = zeros(1,size(ieegMeanRoi,2));
    % [~, ~, actClust]=timePermClusterZOptim(acctimeall,acctimeshuffle,10000,1,norminv(1-pVal2Cutoff));
    
    for iClust=1:length(actClust.Size)
        if actClust.Size{iClust}>actClust.perm95
            pmask(actClust.Start{iClust}: ...
                actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
                = 1;
        end
    end

    
    visTimePlot(timeRange,squeeze((mean(ieegMeanSubj(:,roiIds,:),1))),colval = colors(2,:), labels = xlab )
    hold on;
    visTimePlot(timeRange,squeeze((mean(ieegShuffSubj(:,roiIds,:),1))),colval = [0.5 0.5 0.5], labels = xlab)
    scatter(timeGammaMod(logical(pmask)),-0.01.*ones(1,sum(pmask)),'filled',MarkerEdgeColor=colors(2,:),MarkerFaceColor=colors(2,:));



    bnLabels = [extractBNlabels('astg',hemisphere) extractBNlabels('pstg',hemisphere) extractBNlabels('insula',hemisphere) extractBNlabels('heschl',hemisphere) extractBNlabels('sts',hemisphere)];
    roiIds = extractRoiIds(Subject, channelNamePooled(channelIds),bnLabels);

    
    
    trials2select = randperm(size(ieegMeanMod,1),groupNum);
    ieegMeanRoi = squeeze(mean(mean(ieegMeanSubj(trials2select,roiIds,:),2),1))';
    
    ieegShuffleRoi = [];
    for iTer = 1:1000
        trials2select = randperm(size(ieegMeanMod,1),groupNum);
        ieegShuffleRoi(iTer, :) = squeeze(mean(mean(ieegShuffSubj(trials2select,roiIds,:),2),1))';       
    end
    
     [zValsRawAct, pValsRaw, actClust] = timePermClusterAfterPerm(ieegMeanRoi,ieegShuffleRoi,1,1.645);
    pmask = zeros(1,size(ieegMeanRoi,2));
    % [~, ~, actClust]=timePermClusterZOptim(acctimeall,acctimeshuffle,10000,1,norminv(1-pVal2Cutoff));
    
    for iClust=1:length(actClust.Size)
        if actClust.Size{iClust}>actClust.perm95
            pmask(actClust.Start{iClust}: ...
                actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
                = 1;
        end
    end

    
    visTimePlot(timeRange,squeeze((mean(ieegMeanSubj(:,roiIds,:),1))),colval = colors(3,:), labels = xlab )
    hold on;
    visTimePlot(timeRange,squeeze((mean(ieegShuffSubj(:,roiIds,:),1))),colval = [0.5 0.5 0.5], labels = xlab)
    scatter(timeGammaMod(logical(pmask)),-0.015.*ones(1,sum(pmask)),'filled',MarkerEdgeColor=colors(3,:),MarkerFaceColor=colors(3,:));

    title(subjNames{iSubj})

end