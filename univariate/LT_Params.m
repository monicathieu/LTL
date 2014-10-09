function par = LT_Params(subNo, task, loadScans)
% sets up parameters for batching of LT (lab and camera) experiment.
% <subNo>: a vlue in the range of 1:16, or the subId (e.g. 'LTL001')
% <task>: 'lab' or 'cam'
% <load scans>:
% 1 - load all the scans (this is necessary for analyses of fmri data).
% 0 - do not load scans (saves time)

% this script is adapted by TBM from AG's PM_Params

%% basic task information
par.subNo = subNo;
par.task = task;

%% convert input
if (nargin > 2)
    par.loadScans = loadScans;
else
    par.loadScans = 1; % by deault, load scans
end

par.str.lab = sprintf('LTL%03d',par.subNo);
par.str.cam = sprintf('LTC%03d',par.subNo);

%% subject-specific information

% each subject has slightly different acquisition information.  Here is
% where we store that information.

%par.lab.sourceHand = hand used to indicate source memory on each run
%par.cam.rememberHand = hand used on first run to indicate remembered image
%par.flagIt = should we flag this subject for some special reason?
par.cam.numvols = [280 280 280 280 280 280 280 280];
par.cam.scan0 = [5 5 5 5 5 5 5 5];
par.cam.trialsPerRun = [45 45 45 45 45 45 45 45];

par.lab.numvols = [245 245 245 245 245 245 245 245];
par.lab.scan0 = [5 5 5 5 5 5 5 5];
par.lab.trialsPerRun = [60 60 60 60 60 60 60 60];
par.lab.numRuns = 8;
par.TR = 2;

par.lab.leftButtonMap = {'' 'OBJ' 'FACE' 'SCENE' '' '' '' 'NEW' 'OLD'};
par.lab.rightButtonMap = {'' '' 'NEW' 'OLD' '' '' 'OBJ' 'FACE' 'SCENE'};
par.lab.buttonRemap = {'OBJ' 'FACE' 'SCENE' 'OLD' 'NEW'};

par.cam.leftButtonMap = {'HiR', 'MedR', 'HiF', 'MedF', '', 'HiNew', 'MedNew', 'Uns', 'K'};
par.cam.rightButtonMap = {'HiNew', 'MedNew', 'Uns', 'K', '', 'HiR', 'MedR', 'HiF', 'MedF'};

if mod(par.subNo,2)
    par.lab.sourceHand = 'L';
    par.cam.leftHandRec = [1 1 1 1 0 0 0 0];
else
    par.lab.sourceHand = 'R';
    par.cam.leftHandRec = [0 0 0 0 1 1 1 1];
end

switch par.subNo
    case 1
        par.lab.sourceHand = 'L'; %hand used on first run to indicate source memory
        par.cam.rememberHand = 'L'; %hand used on first run to indicate remembered image
        par.cam.leftHandRec = [1 1 1 1 0 0 0 ];
        par.cam.trialsPerRun = [ 45 45 45 45 45 45 45];
        par.cam.numvols = [ 280 280 280 280 280 280 280];
        par.cam.scan0 = [ 5 5 5 5 5 5 5];
        %par.flagIt = should we flag this subject for some special reason?
        %NOTE LTL: phase 1 lab stimlist repeated at encoding
        %NOTE LTC: no run 8 because we accidentally represented a run here
    case 2 
        par.lab.sourceHand = 'R'; %hand used on first run to indicate source memory
        par.cam.rememberHand = 'R'; %hand used on first run to indicate remembered image
        %par.flagIt = should we flag this subject for some special reason?
        
        %NOTE LTL: lots of missed responses in the 1st and 4th runs.
    case 3
        par.lab.sourceHand = 'L'; %hand used on first run to indicate source memory
        par.cam.rememberHand = 'L'; %hand used on first run to indicate remembered image
        %par.flagIt = should we flag this subject for some special reason?
    case 4
        par.lab.sourceHand = 'R'; %hand used on first run to indicate source memory
        par.cam.rememberHand = 'R'; %hand used on first run to indicate remembered image
        %par.flagIt = should we flag this subject for some special reason?
        
        %NOTE: check that this is the correct sourceHand
    case 5
        par.lab.sourceHand = 'L'; %hand used on first run to indicate source memory
        par.cam.rememberHand = 'L'; %hand used on first run to indicate remembered image
        %par.flagIt = should we flag this subject for some special reason?
    case 6
        par.lab.sourceHand = 'R'; %hand used on first run to indicate source memory
        par.cam.rememberHand = 'R'; %hand used on first run to indicate remembered image
        %par.flagIt = should we flag this subject for some special reason?
    case 7
        par.lab.sourceHand = 'L'; %hand used on first run to indicate source memory
        par.cam.rememberHand = 'L'; %hand used on first run to indicate remembered image
        %par.flagIt = should we flag this subject for some special reason?
    case 8
        par.lab.sourceHand = 'R'; %hand used on first run to indicate source memory
        par.cam.rememberHand = 'R'; %hand used on first run to indicate remembered image
        %par.flagIt = should we flag this subject for some special reason?
        
        %NOTE LTL: Lots of missed responses in the first run (nearly all
        %trials).
    case 9
        par.lab.sourceHand = 'L'; %hand used on first run to indicate source memory
        par.cam.rememberHand = 'L'; %hand used on first run to indicate remembered image
        %par.flagIt = should we flag this subject for some special reason?
        par.cam.leftHandRec = [1 1 1 0 0 0 0];
        par.cam.trialsPerRun = [ 45 45 45 45 45 45 45];
        par.cam.numvols = [ 280 280 280 280 280 280 280];
        par.cam.scan0 = [ 5 5 5 5 5 5 5];
        
        %NOTE: LTC. problem with the first run prescription. discarded
    case 10
        par.lab.sourceHand = 'R'; %hand used on first run to indicate source memory
        par.cam.rememberHand = 'R'; %hand used on first run to indicate remembered image
        %par.flagIt = should we flag this subject for some special reason?
    case 11
        par.lab.sourceHand = 'L'; %hand used on first run to indicate source memory
        par.cam.rememberHand = 'L'; %hand used on first run to indicate remembered image
        %par.flagIt = should we flag this subject for some special reason?
    case 12
        par.lab.sourceHand = 'R'; %hand used on first run to indicate source memory
        par.cam.rememberHand = 'R'; %hand used on first run to indicate remembered image
        %par.flagIt = should we flag this subject for some special reason?
    case 13
        par.lab.sourceHand = 'L'; %hand used on first run to indicate source memory
        par.cam.rememberHand = 'L'; %hand used on first run to indicate remembered image
        par.lab.trialsPerRun = [60 60 60 60 60 60 60 61];
        %par.flagIt = should we flag this subject for some special reason?
        
        %NOTE LTL: Lots of missed responses in first 3 runs
    case 14
        par.lab.sourceHand = 'R'; %hand used on first run to indicate source memory
        par.cam.rememberHand = 'R'; %hand used on first run to indicate remembered image
        %par.flagIt = should we flag this subject for some special reason?
    case 15
        par.lab.sourceHand = 'L'; %hand used on first run to indicate source memory
        par.cam.rememberHand = 'L'; %hand used on first run to indicate remembered image
        %par.flagIt = should we flag this subject for some special reason?
    case 16
        par.lab.sourceHand = 'R'; %hand used on first run to indicate source memory
        par.cam.rememberHand = 'R'; %hand used on first run to indicate remembered image\
        par.cam.leftHandRec = [0 0 0 1 1 1 1 0 ];
        par.cam.trialsPerRun = [16 45 45 45 45 45 45 20];
        par.cam.numvols = [106 280 280 280 280 280 280 125];
        %par.flagIt = should we flag this subject for some special reason?
    otherwise
        error('unrecognized subject');
end

par.substr = par.str.(par.task);
subject = par.str.(task);
par.lab.numscans = sum(par.lab.numvols);
par.lab.numRuns = length(par.lab.trialsPerRun);
par.cam.numRuns = length(par.cam.trialsPerRun);
par.cam.numscans = sum(par.cam.numvols);

%% select which scans to include
% fill in at a later date
par.lab.rtMin = 0;
par.lab.rtMax = 4;
par.cam.rtMin = 0;
par.cam.rtMax = 4;

%% select which scans to include
par.scansSelect.cam = 1:length(par.cam.trialsPerRun);
par.scansSelect.lab = 1:length(par.lab.trialsPerRun);

%% directory information

if strcmp(par.task,'lab')
    studyname = 'LTLab';
elseif strcmp(par.task,'cam')
    studyname = 'LTCam'
end

par.exptdir = ['/Volumes/Tyler_Drive1/' studyname];
par.scriptsdir = fullfile(par.exptdir,'scripts');
par.fmridir = fullfile(par.exptdir,'fmri_data');
par.subdir = fullfile(par.fmridir, subject);
par.funcdir = fullfile(par.subdir,'functional');
par.behavdir = fullfile(par.subdir, 'behav');
par.logdir = fullfile(par.subdir,'logfiles');
par.analysisdir = fullfile(par.subdir, 'analysis');
par.mvpadir = fullfile(par.analysisdir, [task 'MVPA']);


%univariate SPM setup
par.timing.fmri_t = 36;
par.timing.fmri_t0 = 35;
par.timing.units = 'secs';
par.volt = 1;
par.bases.hrf.derivs = [1 1];
par.cvi = 'AR(1)'; %note that this actually gets changed to AR(0.2) in spm_fmri_spm_ui.  
par.global = 'None';
par.mask = '';
par.constat = 'T';
par.sess.condstr = 'ons.mat';
par.sess.regstr = 'regs.mat';
par.sess.hpf = 128;

warning('off');
if ~isempty(par.substr) && par.loadScans
    [par.swascanfiles.(par.task) par.swascanfilesByRun.(par.task)] = findScans(par, 'swavol*.nii',3,par.funcdir);
else
end
warning('on');
end

function [scanfiles scansByRun] = findScans(par, matchstring, scanDim, scandir)
% input:
% <par> - the relevant par structure
% <matchstring> - a string representing a class of image data of interest, e.g. 'Rrun*.nii'
% <scanDim> - how many dimensions is the data? (e.g. 3)
% <scanDir> - where the functional data are located

% returns: 
% <scanfiles> - lists all the images of interest, for the specified
% <scansByRun> - lists all scan files of interest by run
scanfiles = [];
scansByRun = [];

clear sf_h sf_h2

if scanDim == 3
    for rnum=1:length(par.scansSelect.(par.task))
        runDir = fullfile(scandir, sprintf('test%02d',rnum));
        dFN = dir(fullfile(runDir, matchstring));
        scansByRun{rnum} = {dFN.name}';
        scansByRun{rnum} = strcat([runDir '/'],scansByRun{rnum});
        scanfiles = cat(1,scanfiles, scansByRun{rnum});
        
    end
elseif scanDim == 4
    error('we can''t handle 4d nifitis yet');
end
end

    
