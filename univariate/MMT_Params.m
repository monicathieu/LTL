function par = MMT_Params(subNo, task, loadScans)
% sets up parameters for batching of LT (lab and camera) experiment.
% <subNo>: a value in the range of 2:41
% <task>: 'study' or 'loc'
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

par.substr = sprintf('MMT%03d',par.subNo);


%% subject-specific information

% each subject has slightly different acquisition information.  Here is
% where we store that information.

%par.lab.sourceHand = hand used to indicate source memory on each run
%par.cam.rememberHand = hand used on first run to indicate remembered image
%par.flagIt = should we flag this subject for some special reason?
par.study.numvols = [225 225 225 225 225 225];
par.study.scan0 = [5 5 5 5 5 5]; % number of the first scan image
par.study.trialsPerRun = [55 55 55 55 55 55];

par.loc.numvols = [117 117];
par.loc.scan0 = [5 5];
par.loc.trialsPerRun = [60 60 60 60 60 60 60 60]; % ??? check this
par.loc.numRuns = 2;
par.TR = 2;

switch par.subNo % special cases for pesky subjs
    case 4 % stopped in middle of run 1, restarted in the middle of run 1, weirdly partitioned vols for trial 1
    case 21 % removed and replaced in scanner between runs 3 and 4, two inplanes
    case 32 % removed and replaced in scanner between runs 1 and 2, two inplanes
    otherwise
         % error('unrecognized subject');
end

% assigned par.substr above; we use same naming convention for both study
% and localizer
subject = par.str.(task);
par.study.numscans = sum(par.study.numvols);
par.study.numRuns = length(par.study.trialsPerRun);
par.loc.numRuns = length(par.loc.trialsPerRun);
par.loc.numscans = sum(par.loc.numvols);

%% select which scans to include
% this one is scrubbing for bad RTs
% par.lab.rtMin = 0;
% par.lab.rtMax = 4;
% par.cam.rtMin = 0;
% par.cam.rtMax = 4;

%% select which scans to include
par.scansSelect.cam = 1:length(par.cam.trialsPerRun);
par.scansSelect.lab = 1:length(par.lab.trialsPerRun);

%% directory information

par.exptdir = '/biac4/wagner/wagner/MMT_fMRI';
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

    
