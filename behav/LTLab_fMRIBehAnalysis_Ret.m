function [res psyphys idx gRes] = LTLab_fMRIBehAnalysis_Ret(par)
% [res psyphys idx gRes] = LTLab_fMRIBehAnalysis_Retrieval(par) <par> =
% subject parameters setup by LT_Param

if ~isstruct(par) % if it is par_params struct
    par = LT_Params(par ,'lab',1);
end

%specify variables for analysis
EXCLUDEHIGHRTS = 1; %based on LT_Params
FIXDOUBLERESPS = 1; %what to do: get rid of those that are not repeats (i.e. keep 22, get rid of 23) 

res = []; psyphys = []; idx = []; gRes = [];

%cd to subjects behavioral directory and find test logiles
cd(par.behavdir)
logFnames = dir(fullfile(par.logdir, '*test*.mat'));
fileNames = {logFnames.name};
dotFiles = cell2mat(cellfun(@(x) x(1)=='.', fileNames, 'UniformOutput', false));
fileNames(find(dotFiles)) = [];%remove hidden files that are prepended with dots.
%concatenate into one useful file
trialData = combineDataFile(fileNames, par.logdir);
numTrials = length(trialData);
trialRunNum = [];
for runNum = 1:par.lab.numRuns
    trialRunNum = [trialRunNum; zeros(par.lab.trialsPerRun(runNum),1)+runNum];
end

buttonMapLeftSource = {'' 'OBJ' 'FACE' 'SCENE' '' '' '' 'NEW' 'OLD'};
buttonMapRightSource = {'' '' 'NEW' 'OLD' '' '' 'OBJ' 'FACE' 'SCENE'};
buttonRemap = {'OBJ' 'FACE' 'SCENE' 'OLD' 'NEW'};

%% extract stim presentation and behavioral variables
%first, handle responses. fixing the doubles and converting numbers to
%names
resp = cat(2,trialData.resp);
resp = str2double(resp);
if FIXDOUBLERESPS == 1
    doubleresps = resp;
    doubleresps(isnan(resp))=-1;
    doubleresps=(mod(doubleresps,11)==0);
    resp(doubleresps) = floor(resp(doubleresps)/10);
end
resp(resp>10 | resp < 1) = nan;

% transform the responses to [1 2 3 4 5] = [object face scene old new]. to
% do this, we have to modify resps from each hand separately
if strcmp(par.lab.sourceHand, 'L')
    leftresp = resp(:,1:240);
    rightresp = resp(:,241:end);
else
    rightresp = resp(:,1:240);
    leftresp = resp(:,241:end);
end
leftresp = leftresp - 1;
leftresp(leftresp > 3 & leftresp < 7) = 0;
leftresp(leftresp == 8) = 4;
leftresp(leftresp == 7) = 5;
rightresp(rightresp<10) = rightresp(rightresp<10) - 6;
rightresp(rightresp == -2) = 4;
rightresp(rightresp == -3) = 5;
if strcmp(par.lab.sourceHand, 'L') %recombine hands
    resp = cat(2,leftresp,rightresp);
else
    resp = cat(2,rightresp,leftresp);
end

%now set up rts, excluding high ones if desired
idx.rt = cat(2,trialData.respRT)';
if EXCLUDEHIGHRTS
    idx.excludedRts = idx.rt > par.lab.rtMax;
    idx.rt(idx.excludedRts) = nan;
end

valid_trials = ~idx.excludedRts .* ismember(resp,[1 2 3 4 5])';

oldNew = cat(1,trialData.target);
sourcestr = cat(1, trialData.picCat);
idx.runNum = trialRunNum;
idx.trialTime0 = cat(2,trialData.onset)';
%stim attributes
idx.old = ismember(oldNew, {'OLD','Y'});
idx.obj = strcmp('OBJ', sourcestr) .* idx.old;
idx.face = strcmp('FACE', sourcestr) .* idx.old;
idx.scene = strcmp('SCENE', sourcestr) .* idx.old;
idx.new = strcmp('NEW', oldNew);

%if the target is coded as NEW, then make the phase 0
idx.phase = cat(1,trialData.stInterv);
idx.phase = idx.phase .* idx.old;

%resp attributes
idx.validTrials = valid_trials;
idx.respObj = (resp' == 1);
idx.respFace = (resp' == 2);
idx.respScene = (resp' == 3);
idx.respAllSrc = (resp' < 4 & resp' > 0);
idx.respOld = (resp' == 4);
idx.respNew = (resp' == 5);

%keep track of onsets and junk responses
idx.allTrialOns_fixation = idx.trialTime0 - 10 + ((idx.runNum-1) .* par.lab.numvols(idx.runNum)' * par.TR);%ownership-1 .* numvols * TR + stim_on'
idx.allTrialOns_word = idx.allTrialOns_fixation + 1.5;
idx.allTrialOns = idx.allTrialOns_word;
idx.allTrialOns_prompt = idx.allTrialOns_word + 1.5;
idx.junk = ~(valid_trials);

%given these images and responses, what are the outcomes?
idx.OBJobj = idx.obj .* idx.respObj;
idx.FACEface = idx.face .* idx.respFace;
idx.SCENEscene = idx.scene .* idx.respScene;
idx.srcHit = idx.OBJobj + idx.FACEface + idx.SCENEscene;
idx.OBJsrcInc = idx.obj .* (idx.respFace | idx.respScene);
idx.FACEsrcInc = idx.face .* (idx.respObj | idx.respScene);
idx.SCENEsrcInc = idx.scene .* (idx.respObj | idx.respFace);
idx.hitSrcMiss = idx.OBJsrcInc + idx.FACEsrcInc + idx.SCENEsrcInc;
idx.OBJnoSrc = idx.obj .* idx.respOld;
idx.FACEnoSrc = idx.face .* idx.respOld;
idx.SCENEnoSrc = idx.scene .* idx.respOld;
idx.itemHit = idx.OBJnoSrc + idx.FACEnoSrc + idx.SCENEnoSrc;
idx.allHit = idx.srcHit + idx.hitSrcMiss + idx.itemHit;
idx.miss = idx.old .* idx.respNew;
idx.srcFA = idx.respAllSrc .* idx.new;
idx.noSrcFA = idx.respOld .* idx.new;
idx.allFA = idx.srcFA + idx.noSrcFA;
idx.CR = idx.new .* idx.respNew;

%check for some possible errors
assert(isequal(idx.old, (idx.obj + idx.face + idx.scene)));
assert(isequal(idx.old, logical(idx.phase)));
assert(isequal(idx.allHit, (idx.respObj + idx.respScene + idx.respFace + idx.respOld) .* idx.old));


%% results sections
res.pRespOld = sum(idx.respOld)/(sum(idx.respOld)+sum(idx.respNew));
%accuracy results
res.pHit = sum(idx.allHit)/sum(idx.old);
res.pSrcCorHit = sum(idx.srcHit)/sum(idx.old);
res.pSrcIncHit = sum(idx.hitSrcMiss)/sum(idx.old);
res.pNoSrcHit = sum(idx.itemHit)/sum(idx.old);
res.pSrcCorOfHit = sum(idx.srcHit)/sum(idx.allHit); 
res.pFA = sum(idx.allFA)/sum(idx.new);
res.pSrcFA = sum(idx.srcFA)/sum(idx.new);
res.pNoSrcFA = sum(idx.noSrcFA)/sum(idx.new);
%d' results
res.alldprime = norminv(res.pHit) - norminv(res.pFA);
% res.srcCordprime = res.pSrcCorHit - res.pSrcFA; res.srcIncdprime =
% res.pSrcIncHit - res.pSrcFA; res.noSrcdprime = res.pNoSrcHit -
% res.pNoSrcFA;
%rt results
res.rtHit = nanmean(idx.rt(find(idx.allHit)),1);
res.rtsrcHit = nanmean(idx.rt(find(idx.srcHit)),1);
res.rthitSrcMiss = nanmean(idx.rt(find(idx.hitSrcMiss)),1);
res.rtitemHit = nanmean(idx.rt(find(idx.itemHit)),1);
res.rtCR = nanmean(idx.rt(find(idx.CR)),1);
res.rtMiss = nanmean(idx.rt(find(idx.CR)),1);
res.rtAllFA = nanmean(idx.rt(find(idx.allFA)),1);


%% phase psychophysics
%hits by phase
i=0;
for phaseNum = [1,3,4]
    i = i+1;
    thisPhase = (idx.phase==phaseNum);
    %1 hits by phase
    psyphys.hitByPhz.mean(i) = sum(idx.allHit .* thisPhase)/sum(thisPhase);
    psyphys.hitByPhz.SE(i) = sqrt((psyphys.hitByPhz.mean(i) * (1-(psyphys.hitByPhz.mean(i))))/sum(thisPhase));
    psyphys.hitByPhz.ticks.x = 0:4;
    psyphys.hitByPhz.ticks.y = [0 1];
    psyphys.hitByPhz.xlabel = 'Phase ( 1 2 4)';
    psyphys.hitByPhz.ylabel = 'Percent Hits';
    psyphys.hitByPhz.xVals = 1:3;
    psyphys.hitByPhz.marker = 'o-';
    %2 hits source correct by phase
    psyphys.srcHitByPhz.mean(i) = sum(idx.srcHit .* thisPhase)/sum(thisPhase);
    psyphys.srcHitByPhz.SE(i) = sqrt((psyphys.srcHitByPhz.mean(i) * (1-(psyphys.srcHitByPhz.mean(i))))/sum(thisPhase));
    psyphys.srcHitByPhz.ticks.x = 0:4;
    psyphys.srcHitByPhz.ticks.y = [0 1];
    psyphys.srcHitByPhz.xlabel = 'Phase ( 1 2 4)';
    psyphys.srcHitByPhz.ylabel = 'Percent Hits (Source Correct)';
    psyphys.srcHitByPhz.xVals = 1:3;
    psyphys.srcHitByPhz.marker = 'o-';
    %3 percent of hits that are source correct by phase
    psyphys.pSrcCorOfHitByPhz.mean(i) = sum(idx.srcHit .* thisPhase)/sum(idx.allHit .* thisPhase);
    psyphys.pSrcCorOfHitByPhz.SE(i) = sqrt((psyphys.pSrcCorOfHitByPhz.mean(i) * (1-(psyphys.pSrcCorOfHitByPhz.mean(i))))/sum(thisPhase));
    psyphys.pSrcCorOfHitByPhz.ticks.x = 0:4;
    psyphys.pSrcCorOfHitByPhz.ticks.y = [0 1];
    psyphys.pSrcCorOfHitByPhz.xlabel = 'Phase ( 1 2 4)';
    psyphys.pSrcCorOfHitByPhz.ylabel = '(#Source Correct)/(#Item Correct)';
    psyphys.pSrcCorOfHitByPhz.xVals = 1:3;
    psyphys.pSrcCorOfHitByPhz.marker = 'o-';
    %4 hits source incorrect by phase
    psyphys.hitSrcMissByPhz.mean(i) = sum(idx.hitSrcMiss .* thisPhase)/sum(thisPhase);
    psyphys.hitSrcMissByPhz.SE(i) = sqrt((psyphys.hitSrcMissByPhz.mean(i) * (1-(psyphys.hitSrcMissByPhz.mean(i))))/sum(thisPhase));
    psyphys.hitSrcMissByPhz.ticks.x = 0:4;
    psyphys.hitSrcMissByPhz.ticks.y = [0 1];
    psyphys.hitSrcMissByPhz.xlabel = 'Phase ( 1 2 4)';
    psyphys.hitSrcMissByPhz.ylabel = 'Percent Hits (Source Incorrect)';
    psyphys.hitSrcMissByPhz.xVals = 1:3;
    psyphys.hitSrcMissByPhz.marker = 'o-';
    %5 hits source incorrect by phase
    psyphys.itemHitByPhz.mean(i) = sum(idx.itemHit .* thisPhase)/sum(thisPhase);
    psyphys.itemHitByPhz.SE(i) = sqrt((psyphys.itemHitByPhz.mean(i) * (1-(psyphys.itemHitByPhz.mean(i))))/sum(thisPhase));
    psyphys.itemHitByPhz.ticks.x = 0:4;
    psyphys.itemHitByPhz.ticks.y = [0 1];
    psyphys.itemHitByPhz.xlabel = 'Phase ( 1 2 4)';
    psyphys.itemHitByPhz.ylabel = 'Percent Hits (No Source)';
    psyphys.itemHitByPhz.xVals = 1:3;
    psyphys.itemHitByPhz.marker = 'o-';
    %6 collapsed d'
    psyphys.allDprimeByPhz.mean(i) = norminv(psyphys.hitByPhz.mean(i)) - norminv(res.pFA);
    psyphys.allDprimeByPhz.SE(i) = 0;
    psyphys.allDprimeByPhz.ticks.x = 0:4;
    psyphys.allDprimeByPhz.ticks.y = [0 1.5];
    psyphys.allDprimeByPhz.xlabel = 'Phase ( 1 2 4)';
    psyphys.allDprimeByPhz.ylabel = 'd'' overall';
    psyphys.allDprimeByPhz.xVals = 1:3;
    psyphys.allDprimeByPhz.marker = 'o-';
    %7 RT Hit
    psyphys.rtHitByPhz.mean(i) = nanmean(idx.rt(find(idx.allHit .* thisPhase)),1);
    psyphys.rtHitByPhz.SE(i) = nanste(idx.rt(find(idx.allHit .* thisPhase)),1);
    psyphys.rtHitByPhz.ticks.x = 0:4;
    psyphys.rtHitByPhz.ticks.y = [1 5];
    psyphys.rtHitByPhz.xlabel = 'Phase (1 2 4)';
    psyphys.rtHitByPhz.ylabel = 'RT (all hits)';
    psyphys.rtHitByPhz.xVals = 1:3;
    psyphys.rtHitByPhz.marker = 'o-';
    %8 RT Hit (Source Correct)
    psyphys.rtsrcHitByPhz.mean(i) = nanmean(idx.rt(find(idx.srcHit .* thisPhase)),1);
    psyphys.rtsrcHitByPhz.SE(i) = nanste(idx.rt(find(idx.srcHit .* thisPhase)),1);
    psyphys.rtsrcHitByPhz.ticks.x = 0:4;
    psyphys.rtsrcHitByPhz.ticks.y = [1 5];
    psyphys.rtsrcHitByPhz.xlabel = 'Phase (1 2 4)';
    psyphys.rtsrcHitByPhz.ylabel = 'RT (Source Correct)';
    psyphys.rtsrcHitByPhz.xVals = 1:3;
    psyphys.rtsrcHitByPhz.marker = 'o-';
    %9 RT Hit (Source Incorrect)
    psyphys.rthitSrcMissByPhz.mean(i) = nanmean(idx.rt(find(idx.hitSrcMiss .* thisPhase)),1);
    psyphys.rthitSrcMissByPhz.SE(i) = nanste(idx.rt(find(idx.hitSrcMiss .* thisPhase)),1);
    psyphys.rthitSrcMissByPhz.ticks.x = 0:4;
    psyphys.rthitSrcMissByPhz.ticks.y = [1 5];
    psyphys.rthitSrcMissByPhz.xlabel = 'Phase (1 2 4)';
    psyphys.rthitSrcMissByPhz.ylabel = 'RT (Source Incorrect)';
    psyphys.rthitSrcMissByPhz.xVals = 1:3;
    psyphys.rthitSrcMissByPhz.marker = 'o-';
    %10 RT Hit (No Source)
    psyphys.rtitemHitByPhz.mean(i) = nanmean(idx.rt(find(idx.itemHit .* thisPhase)),1);
    psyphys.rtitemHitByPhz.SE(i) = nanste(idx.rt(find(idx.itemHit .* thisPhase)),1);
    psyphys.rtitemHitByPhz.ticks.x = 0:4;
    psyphys.rtitemHitByPhz.ticks.y = [1 5];
    psyphys.rtitemHitByPhz.xlabel = 'Phase (1 2 4)';
    psyphys.rtitemHitByPhz.ylabel = 'RT (No Source)';
    psyphys.rtitemHitByPhz.xVals = 1:3;
    psyphys.rtitemHitByPhz.marker = 'o-';
    
    
    
end
%RT by phase; no accounting for accuracy
i=0;
for phaseNum = [0,1,3,4]
    i = i+1;
    thisPhase = (idx.phase==phaseNum);
    psyphys.rtByPhz.mean(i) = nanmean(idx.rt(find(thisPhase),1));
    psyphys.rtByPhz.SE(i) = nanste(idx.rt(find(thisPhase),1));
    psyphys.rtByPhz.ticks.x = 0:4;
    psyphys.rtByPhz.ticks.y = [1 5];
    psyphys.rtByPhz.xlabel = 'Phase (0 1 2 4) note: 0 is new images';
    psyphys.rtByPhz.ylabel = 'RT (correct and incorrect trials)';
    psyphys.rtByPhz.xVals = 0:3;
    psyphys.rtByPhz.marker = 'o-';
    
    
end
