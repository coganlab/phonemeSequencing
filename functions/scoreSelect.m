function [lossVect,aucVect,C] = scoreSelect(sig2analyzeAllFeature,labels,numDim)
    cvp = cvpartition(labels,'KFold',20,'stratify',true);
   % linearTemplate = templateDiscriminant('DiscrimType','linear');
    labUnique = unique(labels);
    aucVect = []; scoreF = [];
    %[~,scoreall,~,~,~] = pca(sig2analyzeAllFeature,'Centered',false);
    for nCv = 1:cvp.NumTestSets
      
        train = cvp.training(nCv);
        test = cvp.test(nCv);
        featureTrain = sig2analyzeAllFeature(train,:);
        featureTest = sig2analyzeAllFeature(test,:);
        meanTrain = mean(featureTrain,1);
        stdTrain = std(featureTrain,0,1);
        featureTrainNorm = (featureTrain - meanTrain)./stdTrain;
        [coeffTrain,scoreTrain] = pca(featureTrainNorm,'Centered',false);
        featureTestNorm = (featureTest - meanTrain)./stdTrain;
        scoreTest = featureTestNorm*coeffTrain;
        for v = 1:numDim          
            scoreTrainGrid = scoreTrain(:,1:v);
            scoreTestGrid = scoreTest(:,1:v);
            
            linearModel = fitcdiscr((scoreTrainGrid),labels(train),'DiscrimType','linear','CrossVal','off');            
            %linearModel = fitcecoc((scoreTrainGrid),labels(:,train),'Learners','discriminant','Coding','onevsall','CrossVal','off'); 
            [yhat,yscore] = predict(linearModel,scoreTestGrid);
            lossVect(nCv,v) = loss(linearModel,scoreTestGrid,labels(:,test));
%             for t = 1:length(labUnique)
%             [~,~,~,aucVect(nCv,v,t)] = perfcurve(labels(test),yscore(:,t),labUnique(t));
%             end
            C{nCv,v} = confusionmat(labels(test),yhat);
        end
    end
end