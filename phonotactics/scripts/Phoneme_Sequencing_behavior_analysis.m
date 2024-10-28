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

% fieldEpoch = 'Auditory';
% selectRoi = {'lh-middletemporal','lh-superiortemporal','lh-inferiortemporal', 'lh-bankssts'};
% selectRoi = {'lh-supramarginal','lh-inferiorparietal'};
% selectRoi = {'lh-parsopercula','lh-parstriangular'};
% selectRoi = {'lh-rostralmiddlefrontal','lh-caudalmiddlefrontal'};
% selectRoi = {'lh-precentral','lh-postcentral'};
% selectRoi = {'lh-insula'};
% selectRoi = {'lh-superiorfrontal'};
% selectRoi = {'lh-fusiform'};
%selectRoi = {'supramarginal','inferiorparietal','temporal','sts','opercular','triangular','middlefrontal','central','insula','superiorfrontal','fusiform'};
selectRoi = {};
respTimeThresh = 0.01;
timeEpoch = [-0.5  0];

%trialInfoStruct = extractTrialInfo(Subject, remFastResponseTimeTrials=respTimeThresh);

% elecs2remove = elecNameFeedBack_000;
% elecs2remove =  elecs2remove(~cellfun('isempty',elecs2remove));    
    
% elecNameProductionInfo = extractChannelLocation(Subject,elecNameProductionClean);
% elecNameProductionWM = elecNameProductionClean(contains({elecNameProductionInfo.Location},["White"]));
% roiIds = extractRoiIds(Subject, elecNameProductionClean,["White","hypointen","known"]);
subSelectElecs = elecNameProductionClean;

ieegHGStruct = extractRawDataWithROI(Subject, Epoch = 'Start', roi = selectRoi, Time=timeEpoch,remFastResponseTimeTrials=respTimeThresh,...
    subsetElec=subSelectElecs,remWMchannels=true);

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


%  save(fullfile(['pooledSubject_' Task.Name '_Auditory_offset_car_z_score' ...
%     '_prodelecs_data_v1.mat']),...
%     'ieegHGStruct','trialInfoStruct',...
%     'timeEpoch','-v7.3');
%% Behavior analysis
% Load the struct from a .mat file

% Get the number of fields in the struct (assuming each field is a separate dataset)
[m,n] = numSubplots(numel(ieegHGStruct));

% Create a figure for subplots
figure;
% Set the figure size
% set(gcf, 'Position', [100, 100, 1200, 300]); % Adjust width and height as needed
tiledLayout = tiledlayout(m(1), m(2), 'TileSpacing', 'compact', 'Padding', 'compact');

% Loop through each field in the struct
for i = 1:length(ieegHGStruct)
    
    % Create subplot
 nexttile;
    % Combine data1 and data2 for boxplot and create a grouping variable
    combinedData = [ieegHGStruct(i).responseDuration];
   
    groupVariable = [trialInfoStruct(i).phonemeTrial.syllableUnit(:,1)];
     [p, h] = ranksum(combinedData(groupVariable==1), combinedData(groupVariable==2));
    % Create boxplot
    boxplot(combinedData, groupVariable, 'symbol', '');   
    
%     ylim([0 1.5])
%     xlim([0 3])
    % Display significance on the plot
     % Determine the level of significance and display it
     set(gca,'XTick',1:2,'XTickLabel',{'VCV','CVC'})
    if h == 1 % Test is significant
        
        
        % Draw horizontal bar
        %line([1, 2], [ypos, ypos], 'Color', 'k', 'LineWidth', 1.5); 

        % Add significance markers
        
        H=sigstar({{'VCV','CVC'}},p);
    end
    hold off;
    
    formatTicks(gca)
   
end

% Adjust subplot spacing if necessary
% sgtitle('Boxplots Grouped by Categorical Variable'); % Use sgtitle for a common title if using MATLAB R2018b or later
%% Linear mixed effect modeling

reactionTime = [];
reactionDuration = [];
subIds = [];
syllableIds = [];
for i = 1:length(ieegHGStruct)
    
   
    % Combine data1 and data2 for boxplot and create a grouping variable
    reactionTime = [reactionTime; ieegHGStruct(i).responseTime'];
    reactionDuration = [reactionDuration; ieegHGStruct(i).responseDuration'];
    subIds = [subIds; repmat(i,length(ieegHGStruct(i).responseTime),1)];
    syllableIds = [syllableIds; categorical(trialInfoStruct(i).phonemeTrial.syllableUnit(:,2))];
end
subIds = categorical(subIds);

behaviorTable = table(reactionTime,reactionDuration,subIds,syllableIds);
mdl_RT=fitlme(behaviorTable,'reactionTime~syllableIds+(syllableIds|subIds)')
mdl_RD=fitlme(behaviorTable,'reactionDuration~syllableIds+(syllableIds|subIds)')
mdl_RT_RD_control = fitlme(behaviorTable,'reactionTime~syllableIds+(1|subIds)+(reactionDuration|subIds)')
mdl_RT_RD_control_v2 = fitlme(behaviorTable,'reactionTime~syllableIds*reactionDuration+(1|subIds)')
mdl_RT_RD_control_v3 = fitlme(behaviorTable,'reactionTime~syllableIds+(reactionDuration*syllableIds|subIds)')
% mdl_cvc = fitlme(behaviorTable,'syllableIds~reactionTime*reactionDuration+(1|subIds)')
beta = fixedEffects(mdl_RT)
[~,~,STATS] = randomEffects(mdl_RT) % Compute the random-effects statistics (STATS)
beta_rt_shuff = [];
beta_rd_shuff = [];
beta_rt_rd_shuff = [];
for iTer = 1:1000
    iTer
    reactionTimeShuffle = shuffle(reactionTime);
    reactionDurationShuffle = shuffle(reactionDuration);
    behaviorTableShuffle = table(reactionTimeShuffle,reactionDurationShuffle,subIds,syllableIds);
    behaviorTableShuffle2 = table(reactionTimeShuffle,reactionDuration,subIds,syllableIds);
    mdl_RT_shuff=fitlme(behaviorTableShuffle,'reactionTimeShuffle~1+syllableIds+(1|subIds)');
    mdl_RD_shuff=fitlme(behaviorTableShuffle,'reactionDurationShuffle~1+syllableIds+(1|subIds)');
    mdl_RT_RD_control_shuff = fitlme(behaviorTableShuffle2,'reactionTimeShuffle~syllableIds*reactionDuration+(1|subIds)');
    beta_rt_shuff(iTer,:) = fixedEffects(mdl_RT_shuff);
    beta_rd_shuff(iTer,:) = fixedEffects(mdl_RD_shuff);
    beta_rt_rd_shuff(iTer,:) = fixedEffects(mdl_RT_RD_control_shuff);
end

%%
 histHand(:,i) = histfit(timeMaxRoi2plot{i},50, 'kernel' );
    x = [xl(1),histHand(2,i).XData,xl([2,1])]; 
    y = [0,histHand(2,i).YData,0,0]; 
    fillHand = fill(subHand(i),x,y,roiColor(i,:),'FaceAlpha',0.4,'EdgeColor','k','LineWidth',1);
%%
figure; histogram(beta_rt_shuff(:,2),50)
hold on; xline(mdl_RT.Coefficients.Estimate(2), '-','Beta actual fit','LineWidth',2, 'Color','g');
xlabel('Beta for syllable coding');
title('Mixed effects modeling: Syllable vs. Response time')

figure; histogram(beta_rd_shuff(:,2),50)
hold on; xline(mdl_RD.Coefficients.Estimate(2), '-','Beta actual fit','LineWidth',2, 'Color','g');
xlabel('Beta for syllable coding');
title('Mixed effects modeling: Syllable vs. Response Duration')

figure; histogram(beta_rt_rd_shuff(:,4),50)
hold on; xline(mdl_RT_RD_control_v2.Coefficients.Estimate(4), '-','Beta actual fit','LineWidth',2, 'Color','g');
xlabel('Beta for syllable coding');
title('Syllable vs. Response time (controlled for Duration) ')
%% CVC and VCV specific stats

reactionTimeCVC = reactionTime(double(syllableIds)==1);
reactionDurationCVC = reactionDuration(double(syllableIds)==1);

reactionTimeVCV = reactionTime(double(syllableIds)==2);
reactionDurationVCV = reactionDuration(double(syllableIds)==2);
% 
% histsyl = histfit(reactionTimeCVC,50, 'kernel' );
% histVcv =  histfit(reactionTimeVCV,50, 'kernel' );
% col = lines(2);
figure;
histsyl = histfit(reactionTimeCVC,20, 'kernel' );
histsyl(1).FaceColor = 'none';
histsyl(1).EdgeColor = 'none';
x_fit = histsyl(2).XData;
y_fit = histsyl(2).YData;
delete(histsyl(2));
hold on
patch([x_fit, fliplr(x_fit)], [zeros(1,length(y_fit)), fliplr(y_fit)], col(1,:), 'FaceAlpha', 0.3,'EdgeColor',col(1,:));
histsyl = histfit(reactionTimeVCV,20, 'kernel' );
histsyl(1).FaceColor = 'none';
histsyl(1).EdgeColor = 'none';
x_fit = histsyl(2).XData;
y_fit = histsyl(2).YData;
delete(histsyl(2));
hold on
patch([x_fit, fliplr(x_fit)], [zeros(1,length(y_fit)), fliplr(y_fit)], col(2,:), 'FaceAlpha', 0.3,'EdgeColor',col(2,:));

xlabel('Response duration (s)')
ylabel('Trials')

%%
subIdsCVC = subIds(double(syllableIds)==1);
subIdsVCV = subIds(double(syllableIds)==2);

behaviorDataCVC = table(reactionTimeCVC,reactionDurationCVC,subIdsCVC);
behaviorDataVCV = table(reactionTimeVCV,reactionDurationVCV,subIdsVCV);

mdl_CVC=fitlme(behaviorDataCVC,'reactionTimeCVC~reactionDurationCVC+(1|subIdsCVC)')
mdl_VCV=fitlme(behaviorDataVCV,'reactionTimeVCV~reactionDurationVCV+(1|subIdsVCV)')


beta_cvc_shuff = [];
beta_vcv_shuff = [];
for iTer = 1:1000
    iTer
    reactionDurationCVCShuffle = shuffle(reactionDurationCVC);
    reactionDurationVCVShuffle = shuffle(reactionDurationVCV);
    behaviorDataCVCShuffle = table(reactionDurationCVCShuffle,reactionTimeCVC,subIdsCVC);
    behaviorDataVCVShuffle = table(reactionDurationVCVShuffle,reactionTimeVCV,subIdsVCV);
    
    mdl_CVC_shuffle=fitlme(behaviorDataCVCShuffle,'reactionDurationCVCShuffle~reactionTimeCVC+(1|subIdsCVC)');
    mdl_VCV_shuffle=fitlme(behaviorDataVCVShuffle,'reactionDurationVCVShuffle~reactionTimeVCV+(1|subIdsVCV)');
    beta_cvc_shuff(iTer,:) = fixedEffects(mdl_CVC_shuffle);
    beta_vcv_shuff(iTer,:) = fixedEffects(mdl_VCV_shuffle);

end

[rho,pval] = corr(reactionTime(double(syllableIds)==1),reactionDuration(double(syllableIds)==1))

