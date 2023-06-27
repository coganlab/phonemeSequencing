function  decodeTimeGenPhonotactic2D(matfilename,numCore)
%DECODETIMEGENBLICK Summary of this function goes here
%   Detailed explanation goes here
% warning('off','all')
addpath(genpath(fullfile(pwd,'scripts')))
load(fullfile(pwd,'workspace','data',matfilename))
parpool(str2num(numCore))
%addpath(genpath([pwd '/workspace/data/']))
% addpath(genpath([pwd '/scripts/']))



numFold = 5; % K-Fold cross-validation
varExplained = 90; % PCA variance
decodeTimeStruct2D = [];



numIter = 1; % Number of decoder iterations for stable results
timeWin = 0.2; % Window length for temporal generalization (s)
timeRes1D = 0.02; % Window hop for temporal generalization (s)

% save([pwd '/workspace/results/' matfilename(1:length(matfilename)-4) '_Start_phoneme_decoded_1D.mat'],...
%     'channelNamePooled');
phonDecode = decoderClass(numFold,varExplained,numIter);
cvcIds = find(phonemeTrialPooled.syllableUnit(:,1)'==2);
vcvIds = find(phonemeTrialPooled.syllableUnit(:,1)'==1);

phontactid = [6:9];
ieegStruct = ieegStructPooled;
ieegStruct.data = ieegStruct.data(:,cvcIds,:);
for iPhon = 1:length(phontactid)
    disp('2D temporal generalization CVC');
   
    phontactfeature = -log(phonemeTrialPooled.phonotactic(cvcIds,phontactid(iPhon))');
    phontactstats = quantile(phontactfeature,[0,0.25,0.5,0.75,1]);
    phontactcat= ones(1,length(phontactfeature));
    phontactcat(phontactfeature>phontactstats(2) & phontactfeature<=phontactstats(3)) = 2;
    phontactcat(phontactfeature>phontactstats(3) & phontactfeature<=phontactstats(4)) = 3;
    phontactcat(phontactfeature>phontactstats(4) & phontactfeature<=phontactstats(5)) = 4;
    decodeTimeStruct2D_cvc{iPhon} = phonDecode.tempGenClassify2D(ieegStruct,...
            phontactcat,timeRes = timeRes1D, timeWin = timeWin);
end
save([pwd '/workspace/results/' matfilename(1:length(matfilename)-4) '_Start_phonotactic_cvc_decoded_2D.mat'],...   
     'decodeTimeStruct2D_cvc','channelNamePooled','-v7.3');
clear decodeTimeStruct2D_cvc;
clear ieegStruct;


ieegStruct = ieegStructPooled;
ieegStruct.data = ieegStruct.data(:,vcvIds,:);
for iPhon = 1:length(phontactid)
    disp('2D temporal generalization VCV');
   
    phontactfeature = -log(phonemeTrialPooled.phonotactic(vcvIds,phontactid(iPhon))');
    phontactstats = quantile(phontactfeature,[0,0.25,0.5,0.75,1]);
    phontactcat= ones(1,length(phontactfeature));
    phontactcat(phontactfeature>phontactstats(2) & phontactfeature<=phontactstats(3)) = 2;
    phontactcat(phontactfeature>phontactstats(3) & phontactfeature<=phontactstats(4)) = 3;
    phontactcat(phontactfeature>phontactstats(4) & phontactfeature<=phontactstats(5)) = 4;
    decodeTimeStruct2D_vcv{iPhon} = phonDecode.tempGenClassify2D(ieegStruct,...
            phontactcat,timeRes = timeRes1D, timeWin = timeWin);
end
save([pwd '/workspace/results/' matfilename(1:length(matfilename)-4) '_Start_phonotactic_vcv_decoded_2D.mat'],...   
     'decodeTimeStruct2D_vcv', 'channelNamePooled','-v7.3');
clear decodeTimeStruct2D_vcv;



end

