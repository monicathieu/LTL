function gpar = LT_GroupParams(subjArray, analysis)
%gpar = LT_GroupParams(subjArray)

[~,substrArray] = LT_SA;
subsToInclude = substrArray(subjArray);
gpar.subjArray = subjArray;
gpar.thisAnalysis = analysis;

gpar.exptdir = '/Volumes/Tyler_Drive1/LTLab/data/fmri';
gpar.modelTemplate = '/Volumes/Tyler_Drive1/LTLab/scripts/analysis/univariate/GroupTemplate.mat';
gpar.mask = '';
gpar.constat = 'T';
gpar.stat = 't1';


gpar.conTemplate = fullfile(gpar.exptdir, substrArray{1}, 'analysis', gpar.thisAnalysis, 'SPM');
ldTemp = load(gpar.conTemplate);
gpar.SPMcons = ldTemp.SPM.xCon;
% gpar.nCovs = 0;
%     gpar.covVec = [];
%     gpar.covName = [];
for conNum = 1:length(gpar.SPMcons)
    thisConName = gpar.SPMcons(conNum).name
    gpar.cons{conNum}.dir = {fullfile(gpar.exptdir,'group_analyses/fmri', gpar.thisAnalysis,thisConName)};
    gpar.cons{conNum}.name = thisConName;
    
    sCI = 0;
    for sNum = 1:length(subsToInclude)
        
        thisAnalysisDir = fullfile(gpar.exptdir, subsToInclude{sNum}, 'analysis', gpar.thisAnalysis);
        thisConStruct = load(fullfile(thisAnalysisDir, 'conStruct.mat'));
        thisConIdx = find(strcmp(thisConStruct.conStruct.con_name, thisConName));
        
        if ~isempty(thisConIdx)
            sCI = sCI + 1;
            gpar.cons{conNum}.scans{sCI} = fullfile(thisAnalysisDir, thisConStruct.conStruct.con_fileName{thisConIdx});
        end
    end
    gpar.cons{conNum}.groupContrast = {[1] [-1]};
end
end
    

