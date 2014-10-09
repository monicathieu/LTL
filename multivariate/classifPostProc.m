function [res] = classifPostProc(res, saveit, plotit )
if nargin < 3
    plotit = 1;
end
if nargin < 2
    saveit = 0
end

S = res.subj{1}.penalty.nVox.weights.S;
trainConds = {}; testConds = {};
for i = 1:length(S.condsTrain)
    trainConds{i} = S.condsTrain{i}{1}
end
for i = 1:length(S.condsTest)
    testConds{i} = S.condsTest{i}{1}
end
    
nSubs  = length(res.subj);
for i =1:nSubs
var = [res.subj{i}.penalty.nVox.weights.iter{1}.iterations.perfmet];
res.conf(:,:,i) = confusion([var.guesses],[var.desireds]);
end
res.meanConf = mean(res.conf,3);
if plotit
    plotConfmat(res.meanConf,S.saveName,trainConds,testConds);
end
end