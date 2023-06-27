function  decodeTimeGenPhoneme2D(matfilename,numCore)
%DECODETIMEGENBLICK Summary of this function goes here
%   Detailed explanation goes here
% warning('off','all')
addpath(genpath(fullfile(pwd,'scripts')))
load(fullfile(pwd,'workspace','data',matfilename))
parpool(str2num(numCore))
%addpath(genpath([pwd '/workspace/data/']))
% addpath(genpath([pwd '/scripts/']))



numFold = 10; % K-Fold cross-validation
varExplained = 80; % PCA variance
decodeTimeStruct2D = [];



numIter = 1; % Number of decoder iterations for stable results
timeWin = 0.2; % Window length for temporal generalization (s)
timeRes1D = 0.01; % Window hop for temporal generalization (s)

% save([pwd '/workspace/results/' matfilename(1:length(matfilename)-4) '_Start_phoneme_decoded_1D.mat'],...
%     'channelNamePooled');
phonDecode = decoderClass(numFold,varExplained,numIter);
cvcIds = find(phonemeTrialPooled.syllableUnit(:,1)'==2);
vcvIds = find(phonemeTrialPooled.syllableUnit(:,1)'==1);

% for iPhon = 1:3
%     disp('1D temporal generalization');
%     % 1D temporalGeneralization            
%     decodeTimeStruct1D{iPhon} = phonDecode.tempGenClassify1D(ieegStructPooled,...
%         phonemeTrialPooled.phonemeUnit(:,iPhon)',timeRes = timeRes1D, timeWin = timeWin);      
% end
% save([pwd '/workspace/results/' matfilename(1:length(matfilename)-4) '_Start_phoneme_decoded_1D.mat'],...   
%      'decodeTimeStruct1D','channelNamePooled');
% clear decodeTimeStruct1D;
% for iPhon = 1:3
% decodeTimeStruct1D_cvc{iPhon} = phonDecode.tempGenClassify1D(ieegStructPooled,...
%         phonemeTrialPooled.phonemeUnit(:,iPhon)',timeRes = timeRes1D, timeWin = timeWin, selectTrials = cvcIds);
% end
% save([pwd '/workspace/results/' matfilename(1:length(matfilename)-4) '_Start_phoneme_cvc_decoded_1D.mat'],...   
%      'decodeTimeStruct1D_cvc','channelNamePooled');
% clear decodeTimeStruct1D_cvc;
% for iPhon = 1:3
% decodeTimeStruct1D_vcv{iPhon} = phonDecode.tempGenClassify1D(ieegStructPooled,...
%     phonemeTrialPooled.phonemeUnit(:,iPhon)',timeRes = timeRes1D, timeWin = timeWin, selectTrials = vcvIds);
% end
% save([pwd '/workspace/results/' matfilename(1:length(matfilename)-4) '_Start_phoneme_vcv_decoded_1D.mat'],...   
%      'decodeTimeStruct1D_vcv', 'channelNamePooled');
% clear decodeTimeStruct1D_vcv;

for iPhon = 1:3    
    % 2D temporalGeneralization            
    disp('2D temporal generalization');
    decodeTimeStruct2D{iPhon} = phonDecode.tempGenClassify2D(ieegStructPooled,...
        phonemeTrialPooled.phonemeUnit(:,iPhon)');
end

save([pwd '/workspace/results/' matfilename(1:length(matfilename)-4) '_Start_phoneme_decoded_2D.mat'],...   
    'decodeTimeStruct2D',  'channelNamePooled','-v7.3');



end

