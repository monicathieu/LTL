function LT_GroupRunUnivariate(gpar, analysis, flags)
%LT_GroupRunUnivariate(gpar, flags)

if ~exist('flags','var')
    flags='mec'
    fprintf('***WARNING: Assuming default flags for group univariate analyses***')

if ~isstruct(gpar)
    gpar = LT_GroupParams(gpar, analysis);
end

if ismember('m',flags)
    LT_GroupModelSpec(gpar, analysis);
end
if ismember('e', flags)
    for cNum = 1:length(gpar.cons)
        LT_GroupModelEst(gpar, cNum, analysis);
    end
end

if ismember('c', flags)
    for cNum = 1:length(gpar.cons)
        LT_GroupSetContrasts(gpar, cNum, analysis)
    end
end
end