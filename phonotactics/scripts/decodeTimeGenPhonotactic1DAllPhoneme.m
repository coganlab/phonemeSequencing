function  decodeTimeGenPhonotactic1DAllPhoneme(matfilename,numCore)
%DECODETIMEGENBLICK Summary of this function goes here
%   Detailed explanation goes here
% warning('off','all')
addpath(genpath(fullfile(pwd,'scripts')))
load(fullfile(pwd,'workspace','data',matfilename))
parpool(str2num(numCore))
%addpath(genpath([pwd '/workspace/data/']))
% addpath(genpath([pwd '/scripts/']))



numFold = 20; % K-Fold cross-validation
varExplained = 80; % PCA variance
decodeTimeStruct1D = [];

decodeTimeStruct2D = [];

numIter = 1; % Number of decoder iterations for stable results
timeWin = 0.2; % Window length for temporal generalization (s)
timeRes1D = 0.005; % Window hop for temporal generalization (s)

save([pwd '/workspace/results/' matfilename(1:length(matfilename)-4) '_Start_phonotactic_decoded_1D.mat'],...
    'channelNamePooled');
phonDecode = decoderClass(numFold,varExplained,numIter);
cvcIds = find(phonemeTrialPooled.syllableUnit(:,1)'==2);
vcvIds = find(phonemeTrialPooled.syllableUnit(:,1)'==1);

% for iPhon = 1:3
%     phontactfeature = -log(phonemeTrialPooled.phonotactic(:,ipos)');
%     disp('1D temporal generalization');
%     % 1D temporalGeneralization            
%     decodeTimeStruct1D{iPhon} = phonDecode.tempGenClassify1D(ieegStructPooled,...
%         phonemeTrialPooled.phonemeUnit(:,iPhon)',timeRes = timeRes1D, timeWin = timeWin);      
% end
% save([pwd '/workspace/results/' matfilename(1:length(matfilename)-4) '_Start_phoneme_decoded_1D.mat'],...   
%      'decodeTimeStruct1D','channelNamePooled','-v7.3');
% clear decodeTimeStruct1D;
phontactid = [6:9];
% ieegStruct = ieegStructPooled;
% ieegStruct.data = ieegStruct.data(:,cvcIds,:);
for iPhon = 1:length(phontactid)
    disp('1D temporal generalization CVC');
   
    phontactfeature = (phonemeTrialPooled.phonotactic(:,phontactid(iPhon))');
    phontactstats = quantile(phontactfeature,[0,0.25,0.5,0.75,1]);
    phontactcat= ones(1,length(phontactfeature));
    phontactcat(phontactfeature>phontactstats(2) & phontactfeature<=phontactstats(3)) = 2;
    phontactcat(phontactfeature>phontactstats(3) & phontactfeature<=phontactstats(4)) = 3;
    phontactcat(phontactfeature>phontactstats(4) & phontactfeature<=phontactstats(5)) = 4;
    decodeTimeStruct1D{iPhon} = phonDecode.tempGenClassify1D(ieegStructPooled,...
            phontactcat,timeRes = timeRes1D, timeWin = timeWin);
    save([pwd '/workspace/results/' matfilename(1:length(matfilename)-4) '_Start_phonotactic_decoded_1D.mat'],...   
     'decodeTimeStruct1D','channelNamePooled','-v7.3');
end

clear decodeTimeStruct1D;
clear ieegStruct;

% 
% ieegStruct = ieegStructPooled;
% ieegStruct.data = ieegStruct.data(:,vcvIds,:);
% for iPhon = 1:length(phontactid)
%     disp('1D temporal generalization VCV');
%    
%     phontactfeature = (phonemeTrialPooled.phonotactic(vcvIds,phontactid(iPhon))');
%     phontactstats = quantile(phontactfeature,[0,0.25,0.5,0.75,1]);
%     phontactcat= ones(1,length(phontactfeature));
%     phontactcat(phontactfeature>phontactstats(2) & phontactfeature<=phontactstats(3)) = 2;
%     phontactcat(phontactfeature>phontactstats(3) & phontactfeature<=phontactstats(4)) = 3;
%     phontactcat(phontactfeature>phontactstats(4) & phontactfeature<=phontactstats(5)) = 4;
%     decodeTimeStruct1D_vcv{iPhon} = phonDecode.tempGenClassify1D(ieegStruct,...
%             phontactcat,timeRes = timeRes1D, timeWin = timeWin);
%     save([pwd '/workspace/results/' matfilename(1:length(matfilename)-4) '_Start_phonotactic_vcv_decoded_1D.mat'],...   
%      'decodeTimeStruct1D_vcv', 'channelNamePooled','-v7.3');
% end
% 
% clear decodeTimeStruct1D_vcv;

% for iPhon = 1:3    
%     % 2D temporalGeneralization            
%     disp('2D temporal generalization');
%     decodeTimeStruct2D{iPhon} = phonDecode.tempGenClassify2D(ieegStructPooled,...
%         phonemeTrialPooled.phonemeUnit(:,iPhon)');
% end
% 
% save([pwd '/workspace/results/' matfilename(1:length(matfilename)-4) '_Start_phoneme_decoded_2D.mat'],...   
%     'decodeTimeStruct2D',  'channelNamePooled','-v7.3');



end

