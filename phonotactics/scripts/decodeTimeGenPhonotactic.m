function  decodeTimeGenPhonotactic(matfilename)
%DECODETIMEGENBLICK Summary of this function goes here
%   Detailed explanation goes here
addpath(genpath([pwd '/data/']))
addpath(genpath([pwd '/scripts/']))

load(matfilename);

numFold = 10; % K-Fold cross-validation
varExplained = 80; % PCA variance
decodeTimeStruct1D_cvc = [];
decodeTimeStruct1D_vcv = [];
decodeTimeStruct2D_cvc = [];
decodeTimeStruct2D_vcv = [];
decodeTimeStruct1D = [];
decodeTimeStruct2D = [];
numIter = 1; % Number of decoder iterations for stable results
timeWin = 0.2; % Window length for temporal generalization (s)
timeRes = 0.02; % Window hop for temporal generalization (s)
cvcIds = find(phonemeTrialPooled.syllableUnit(:,1)'==2);
vcvIds = find(phonemeTrialPooled.syllableUnit(:,1)'==1);
save([pwd '/results/' matfilename(1:length(matfilename)-4) '_Start_phonotactic_decoded.mat'],...
    'channelNamePooled');

for ipos = 1:9
    % Initialize Decoder object
    phonDecode = decoderClass(numFold,varExplained,numIter);
    disp('CVC - 1D temporal generalization');
    phontactfeature = -log(phonemeTrialPooled.phonotactic(:,ipos)');

    % cvc - 1D temporalGeneralization            
    decodeTimeStruct1D_cvc{ipos} = phonDecode.tempGenRegress1D(ieegStructPooled,...
        phontactfeature,  selectTrials= cvcIds);
    % cvc - 2D temporalGeneralization            
    disp('CVC - 2D temporal generalization');
    decodeTimeStruct2D_cvc{ipos} = phonDecode.tempGenRegress2D(ieegStructPooled,...
        phontactfeature, selectTrials = cvcIds);
    
    disp('VCV - 1D temporal generalization');
    % 1D temporalGeneralization            
    decodeTimeStruct1D_vcv{ipos} = phonDecode.tempGenRegress1D(ieegStructPooled,...
        -log(1+phonemeTrialPooled.phonotactic(:,ipos)'),selectTrials = vcvIds);
    % 2D temporalGeneralization            
    disp('VCV - 2D temporal generalization');
    decodeTimeStruct2D_vcv{ipos} = phonDecode.tempGenRegress2D(ieegStructPooled,...
        phontactfeature,selectTrials = vcvIds);
    
    
    disp('1D temporal generalization');
    % 1D temporalGeneralization            
    decodeTimeStruct1D{ipos} = phonDecode.tempGenRegress1D(ieegStructPooled,...
        phontactfeature);
    % 2D temporalGeneralization            
    disp('2D temporal generalization');
    decodeTimeStruct2D{ipos} = phonDecode.tempGenRegress2D(ieegStructPooled,...
        phontactfeature);
end

save([pwd '/results/' matfilename(1:length(matfilename)-4) '_Start_log_phonotactic_decoded.mat'],...
    'decodeTimeStruct1D_cvc','decodeTimeStruct2D_cvc',...
    'decodeTimeStruct1D_vcv','decodeTimeStruct2D_vcv',...
    'decodeTimeStruct2D', 'decodeTimeStruct1D', 'channelNamePooled');



end

