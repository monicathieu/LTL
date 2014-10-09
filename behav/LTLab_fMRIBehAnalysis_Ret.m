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
    doubleresps = resp
    doubleresps(isnan(resp))=-1
    doubleresps=(mod(doubleresps,11)==0)
    resp(doubleresps) = floor(resp(doubleresps)/10);
end
resp(resp>10 | resp < 1) = nan

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
idx.time0 = cat(2,trialData.onset)';
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
idx.allTrialOns = idx.time0 - 10 + ((idx.runNum-1) .* par.lab.numvols(idx.runNum)' * par.TR);%ownership-1 .* numvols * TR + stim_on'
idx.junk = ~(valid_trials);

%given these images and responses, what are the outcomes?
idx.OBJobj = idx.obj .* idx.respObj;
idx.FACEface = idx.face .* idx.respFace;
idx.SCENEscene = idx.scene .* idx.respScene;
idx.hitSrcCor = idx.OBJobj + idx.FACEface + idx.SCENEscene;
idx.OBJsrcInc = idx.obj .* (idx.respFace | idx.respScene);
idx.FACEsrcInc = idx.face .* (idx.respObj | idx.respScene);
idx.SCENEsrcInc = idx.scene .* (idx.respObj | idx.respFace);
idx.hitSrcInc = idx.OBJsrcInc + idx.FACEsrcInc + idx.SCENEsrcInc;
idx.OBJnoSrc = idx.obj .* idx.respOld;
idx.FACEnoSrc = idx.face .* idx.respOld;
idx.SCENEnoSrc = idx.scene .* idx.respOld;
idx.hitNoSrc = idx.OBJnoSrc + idx.FACEnoSrc + idx.SCENEnoSrc;
idx.allHit = idx.hitSrcCor + idx.hitSrcInc + idx.hitNoSrc;
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
res.pSrcCorHit = sum(idx.hitSrcCor)/sum(idx.old);
res.pSrcIncHit = sum(idx.hitSrcInc)/sum(idx.old);
res.pNoSrcHit = sum(idx.hitNoSrc)/sum(idx.old);
res.pSrcCorOfHit = sum(idx.hitSrcCor)/sum(idx.allHit); 
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
res.rtHitSrcCor = nanmean(idx.rt(find(idx.hitSrcCor)),1);
res.rtHitSrcInc = nanmean(idx.rt(find(idx.hitSrcInc)),1);
res.rtHitNoSrc = nanmean(idx.rt(find(idx.hitNoSrc)),1);
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
    psyphys.hitSrcCorByPhz.mean(i) = sum(idx.hitSrcCor .* thisPhase)/sum(thisPhase);
    psyphys.hitSrcCorByPhz.SE(i) = sqrt((psyphys.hitSrcCorByPhz.mean(i) * (1-(psyphys.hitSrcCorByPhz.mean(i))))/sum(thisPhase));
    psyphys.hitSrcCorByPhz.ticks.x = 0:4;
    psyphys.hitSrcCorByPhz.ticks.y = [0 1];
    psyphys.hitSrcCorByPhz.xlabel = 'Phase ( 1 2 4)';
    psyphys.hitSrcCorByPhz.ylabel = 'Percent Hits (Source Correct)';
    psyphys.hitSrcCorByPhz.xVals = 1:3;
    psyphys.hitSrcCorByPhz.marker = 'o-';
    %3 percent of hits that are source correct by phase
    psyphys.pSrcCorOfHitByPhz.mean(i) = sum(idx.hitSrcCor .* thisPhase)/sum(idx.allHit .* thisPhase);
    psyphys.pSrcCorOfHitByPhz.SE(i) = sqrt((psyphys.pSrcCorOfHitByPhz.mean(i) * (1-(psyphys.pSrcCorOfHitByPhz.mean(i))))/sum(thisPhase));
    psyphys.pSrcCorOfHitByPhz.ticks.x = 0:4;
    psyphys.pSrcCorOfHitByPhz.ticks.y = [0 1];
    psyphys.pSrcCorOfHitByPhz.xlabel = 'Phase ( 1 2 4)';
    psyphys.pSrcCorOfHitByPhz.ylabel = '(#Source Correct)/(#Item Correct)';
    psyphys.pSrcCorOfHitByPhz.xVals = 1:3;
    psyphys.pSrcCorOfHitByPhz.marker = 'o-';
    %4 hits source incorrect by phase
    psyphys.hitSrcIncByPhz.mean(i) = sum(idx.hitSrcInc .* thisPhase)/sum(thisPhase);
    psyphys.hitSrcIncByPhz.SE(i) = sqrt((psyphys.hitSrcIncByPhz.mean(i) * (1-(psyphys.hitSrcIncByPhz.mean(i))))/sum(thisPhase));
    psyphys.hitSrcIncByPhz.ticks.x = 0:4;
    psyphys.hitSrcIncByPhz.ticks.y = [0 1];
    psyphys.hitSrcIncByPhz.xlabel = 'Phase ( 1 2 4)';
    psyphys.hitSrcIncByPhz.ylabel = 'Percent Hits (Source Incorrect)';
    psyphys.hitSrcIncByPhz.xVals = 1:3;
    psyphys.hitSrcIncByPhz.marker = 'o-';
    %5 hits source incorrect by phase
    psyphys.hitNoSrcByPhz.mean(i) = sum(idx.hitNoSrc .* thisPhase)/sum(thisPhase);
    psyphys.hitNoSrcByPhz.SE(i) = sqrt((psyphys.hitNoSrcByPhz.mean(i) * (1-(psyphys.hitNoSrcByPhz.mean(i))))/sum(thisPhase));
    psyphys.hitNoSrcByPhz.ticks.x = 0:4;
    psyphys.hitNoSrcByPhz.ticks.y = [0 1];
    psyphys.hitNoSrcByPhz.xlabel = 'Phase ( 1 2 4)';
    psyphys.hitNoSrcByPhz.ylabel = 'Percent Hits (No Source)';
    psyphys.hitNoSrcByPhz.xVals = 1:3;
    psyphys.hitNoSrcByPhz.marker = 'o-';
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
    psyphys.rtHitSrcCorByPhz.mean(i) = nanmean(idx.rt(find(idx.hitSrcCor .* thisPhase)),1);
    psyphys.rtHitSrcCorByPhz.SE(i) = nanste(idx.rt(find(idx.hitSrcCor .* thisPhase)),1);
    psyphys.rtHitSrcCorByPhz.ticks.x = 0:4;
    psyphys.rtHitSrcCorByPhz.ticks.y = [1 5];
    psyphys.rtHitSrcCorByPhz.xlabel = 'Phase (1 2 4)';
    psyphys.rtHitSrcCorByPhz.ylabel = 'RT (Source Correct)';
    psyphys.rtHitSrcCorByPhz.xVals = 1:3;
    psyphys.rtHitSrcCorByPhz.marker = 'o-';
    %9 RT Hit (Source Incorrect)
    psyphys.rtHitSrcIncByPhz.mean(i) = nanmean(idx.rt(find(idx.hitSrcInc .* thisPhase)),1);
    psyphys.rtHitSrcIncByPhz.SE(i) = nanste(idx.rt(find(idx.hitSrcInc .* thisPhase)),1);
    psyphys.rtHitSrcIncByPhz.ticks.x = 0:4;
    psyphys.rtHitSrcIncByPhz.ticks.y = [1 5];
    psyphys.rtHitSrcIncByPhz.xlabel = 'Phase (1 2 4)';
    psyphys.rtHitSrcIncByPhz.ylabel = 'RT (Source Incorrect)';
    psyphys.rtHitSrcIncByPhz.xVals = 1:3;
    psyphys.rtHitSrcIncByPhz.marker = 'o-';
    %10 RT Hit (No Source)
    psyphys.rtHitNoSrcByPhz.mean(i) = nanmean(idx.rt(find(idx.hitNoSrc .* thisPhase)),1);
    psyphys.rtHitNoSrcByPhz.SE(i) = nanste(idx.rt(find(idx.hitNoSrc .* thisPhase)),1);
    psyphys.rtHitNoSrcByPhz.ticks.x = 0:4;
    psyphys.rtHitNoSrcByPhz.ticks.y = [1 5];
    psyphys.rtHitNoSrcByPhz.xlabel = 'Phase (1 2 4)';
    psyphys.rtHitNoSrcByPhz.ylabel = 'RT (No Source)';
    psyphys.rtHitNoSrcByPhz.xVals = 1:3;
    psyphys.rtHitNoSrcByPhz.marker = 'o-';
    
    
    
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
