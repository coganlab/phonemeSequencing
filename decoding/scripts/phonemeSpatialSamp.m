subjectId = 'S26';
load([subjectId '-dataPack.mat'])
nonSigChannel = [];

nSamp = 1.5;
firstElec = [0.025 0.05];
totalDist = 120;



for iWay = 2:9
    accDist = [];
    numSamp = [];
    labels2choose = nchoosek(1:9,iWay);
    for iDist = 1:totalDist
        elecSampDensity = firstElec+(iDist-1)*firstElec(1);
        iDist
        for iSamp = 1:round(nSamp*(totalDist/iDist))
            
           
            elecSampClose = elecSampDensity;
             nElec = round((elecSampClose(1) + (elecSampClose(2)-elecSampClose(1)) .* rand(1,1))*length(selectedChannels));
             %elecPtIds = round(poissonDisc([8,16],poisSpace,nElec));
             elecPtIds = ceil(poissonDisc2(size(chanMap),nElec));
                
    
            elecPt = [];
            for iElec = 1:size(elecPtIds,1)
                elecPt(iElec) = chanMap(elecPtIds(iElec,1),elecPtIds(iElec,2));
            end
            elecPt = elecPt(~isnan(elecPt)); 
            elecPtcm = ismember(selectedChannels,elecPt);
            accTrial = zeros(1,size(labels2choose,1));
            parfor iTrial = 1:size(labels2choose,1)
                rTrials= ismember(phonemeUnits,labels2choose(iTrial,:));
                dObj = decoderClass(10,80);
                decodeResultStruct = dObj.baseClassify(ieegHGRespNorm,phonemeUnits, [-0.5 0.5], find(elecPtcm),rTrials);
                accTrial(iTrial) = decodeResultStruct.accPhoneme;
            end
            mean(accTrial)
            accDist = [accDist mean(accTrial)];
            
            numSamp = [numSamp sum(elecPtcm)];
        end
        sum(elecPtcm)
        save([subjectId '_PhonemeDecodeHighResFirstPhoneme_' num2str(iWay) '_way_v2.mat'],'accDist', 'numSamp');   
    end
end
