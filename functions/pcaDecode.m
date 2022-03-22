function [lossMod,Cmat,scoreTrainMode,scoreTestMode] = pcaDecode(sigTrain,sigTest,YTrain,YTest,nModes)
        meanTrain = mean(sigTrain,1);
        stdTrain = std(sigTrain,0,1);
        sigTrainNorm = (sigTrain - meanTrain)./stdTrain;
        [coeffTrain,scoreTrain] = pca(sigTrainNorm,'Centered',false);
        sigTestNorm = (sigTest - meanTrain)./stdTrain;
        scoreTest = sigTestNorm*coeffTrain;
        
        scoreTrainMode = scoreTrain(:,1:nModes);
        scoreTestMode = scoreTest(:,1:nModes);
       linearModel = fitcdiscr((scoreTrainMode),YTrain,'DiscrimType','linear','CrossVal','off');      
       % linearModel = fitcecoc((scoreTrainGrid),YTrain,'Learners','linear','CrossVal','off');
       %linearModel = fitcecoc((scoreTrainGrid),YTrain,'Learners','discriminant','Coding','onevsall','CrossVal','off'); 
           
        [yhat,yscore] = predict(linearModel,scoreTestMode);
        lossMod = loss(linearModel,scoreTestMode,YTest);
        Cmat =  confusionmat(YTest,yhat);
end