accGrid = chanLab;
%%
accStereo = [];
%%
accStereo = [accStereo chanLab];
%%
figure; distributionPlot({accStereo.*100,accGrid.*100, accMicro.*100},'colormap',copper,'showMM',6,'xNames',{'SEEG','macro-ECoG','\muECoG'},'yLabel','Accuracy (%)');
 sigstar({[1,3],[2,3]},[2.4e-40,1.37e-4]);
 %%
 maxGrid = []; 
 maxStereo = [];
 maxMicro = [];
 
 meanGrid = []; 
 meanStereo = [];
 meanMicro = [];
sampInd = 2;
 for iSamp = 1: 10000
     gridSamp = randsample(accGrid,sampInd,false);
     stereoSamp = randsample(accStereo,sampInd,false);
     microSamp = randsample(accMicro,sampInd,false);
     
     x = gridSamp; SEM = std(x)/sqrt(length(x));               % Standard Error
    ts = tinv([0.025  0.975],length(x)-1);      % T-Score
    maxGrid(iSamp,:) = mean(x) + ts*SEM;
     meanGrid(iSamp) = mean(x);
    
     x = stereoSamp; SEM = std(x)/sqrt(length(x));               % Standard Error
    ts = tinv([0.025  0.975],length(x)-1);      % T-Score
    maxStereo(iSamp,:) = mean(x) + ts*SEM;
     meanStereo(iSamp) = mean(x);
    
     x = microSamp; SEM = std(x)/sqrt(length(x));               % Standard Error
    ts = tinv([0.025  0.975],length(x)-1);      % T-Score
    maxMicro(iSamp,:) = mean(x) + ts*SEM;
     meanMicro(iSamp) = mean(x);
     
 end
 
mStereo =  mean(maxStereo)*100
mGrid =  mean(maxGrid)*100
mMicro =  mean(maxMicro)*100
 

x = meanStereo; SEM = std(x)/sqrt(length(x));               % Standard Error
    ts = tinv([0.025  0.975],length(x)-1);      % T-Score
    ciStereo = mean(x) + ts*SEM;
     
    x = meanGrid; SEM = std(x)/sqrt(length(x));               % Standard Error
    ts = tinv([0.025  0.975],length(x)-1);      % T-Score
    ciGrid = mean(x) + ts*SEM;
    
    x = meanMicro; SEM = std(x)/sqrt(length(x));               % Standard Error
    ts = tinv([0.025  0.975],length(x)-1);      % T-Score
    ciMicro = mean(x) + ts*SEM;
    
    
    x = accStereo; SEM = std(x)/sqrt(length(x));               % Standard Error
    ts = tinv([0.025  0.975],length(x)-1);      % T-Score
    ogStereo = mean(x) + ts*SEM
     
    x = accGrid; SEM = std(x)/sqrt(length(x));               % Standard Error
    ts = tinv([0.025  0.975],length(x)-1);      % T-Score
    ogGrid = mean(x) + ts*SEM
    
    x = accMicro; SEM = std(x)/sqrt(length(x));               % Standard Error
    ts = tinv([0.025  0.975],length(x)-1);      % T-Score
    ogMicro = mean(x) + ts*SEM
 %%
figure; distributionPlot({accStereo'.*100,accGrid'.*100, accMicro'.*100},'cidata',[ogStereo; ogGrid; ogMicro].*100,'colormap',copper,'showMM',7,'xNames',{'SEEG','macro-ECoG','\muECoG'},'yLabel','Accuracy (%)');
 sigstar({[1,3],[2,3]},[2.4e-40,1.37e-4]);
 %%
 accAll = [accStereo'.*100;accGrid'.*100; accMicro'.*100];
 g1 = repmat({'SEEG'},length(accStereo),1);
g2 = repmat({'macro-ECoG'},length(accGrid),1);
g3 = repmat({'\muECoG'},length(accMicro),1);
g = [g1; g2; g3];
 figure; boxplot(accAll,g,'whisker', 0.7193);
 ylabel('Accuracy (%)')
 set(gca,'FontSize',20);
 %% Power & Accuracy Compilation
 %accStandard = [accStandard accChan];
 gammaPowerStandard = [gammaPowerStandard ieegGammaPowerNorm(motorChan) ieegGammaPowerNorm(sensoryChan) ieegGammaPowerNorm(ifgChan)];
 
 %% Accuracy comparison
accAll = [accStandard'.*100; accChan(pvalsMCleanProdResp)'.*100];
g1 = repmat({'Standard ECoG'},length(accStandard),1);
g2 = repmat({'?ECoG'},length(find(pvalsMCleanProdResp)),1);
g = [g1; g2];
 figure; h = boxplot(accAll,g,'symbol','');
 set(h,{'linew'},{2})
 ylabel('Accuracy (%)')
 set(gca,'FontSize',15);
pAccuracy = ranksum(accStandard,accChan)
sigstar({[1,2]},[pAccuracy]);
 %%
figure;
 accMean = [mean(gammaPowerStandard) mean(ieegGammaPowerNorm(pvalsMCleanProdResp))];
 accErr = [std(gammaPowerStandard)./sqrt(length(gammaPowerStandard)) std(ieegGammaPowerNorm(pvalsMCleanProdResp))./sqrt(length(find(pvalsMCleanProdResp))) ];
 %accHigh = [std(sigPowerStereoPerc,0.75) quantile(sigPowerGridPerc,0.75) quantile(sigPowerMicroPerc',0.75)];
 elecNames = {'Standard ECoG','uECoG'};
 bar(1:2,accMean)                
pPower = ranksum(gammaPowerStandard, evokeSnr(pvalsMCleanProdResp));
hold on
 ylabel('ESNR (dB)')
er = errorbar(1:2,accMean,accErr);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 2;
set(gca, 'XTickLabel',elecNames)
 set(gca,'FontSize',15);
sigstar({[1,2]},[pPower]);

%% spec Power compile

rmsECoG = [rmsECoG; (specPower(anatChan & pvalsMCleanProd))]
%%

rmsMicro = sigPowerMicro;
rmsAll = [20.*log10(rmsECoG') 20.*log10(rmsMicro)];
g1 = repmat({'Standard IEEG'},length(rmsECoG),1);
g2 = repmat({'µECoG'},length(rmsMicro),1);
g = [g1; g2];
 figure; h = boxplot(rmsAll,g, 'symbol','');
 set(h,{'linew'},{2})
 pGamma = ranksum(20.*log10(rmsECoG ),20.*log10(rmsMicro))
sigstar({[1,2]},[pGamma]);
ylabel('High Gamma Power (dB)')
set(gca,'FontSize',15);
axis square
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         );
%%
figure;
bar(1,[56.7 56.6 30.2 50.6]);
hold on;
bar(2,[43.77 33.32 26.04 22.31 21.92 20.88 20.24]);
yline(11.11, '--','chance','LineWidth',2);
ylabel('P1 accuracy (%)');
set(gca,'FontSize',15);
axis square

 
 