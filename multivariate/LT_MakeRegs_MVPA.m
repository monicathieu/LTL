function LT_MakeRegs_MVPA(par, task, saveit)
%LTCam_MakeRegs(par, saveit) creates onsets and regresors (boxcar and motion)
%for use in SPM's GLM.
%inputs:
%par = subject-specific parameters created by LTCam_Params
%saveit = 1 if we'd like to save mat files for onsets and boxcar + motion
%regressors(default is 1)
STIMDUR = 0;

if nargin<3, saveit=1; end

%recursive call if <par> is array
if length(par) > 1
    for i =1:length(par)
        LT_MakeRegs_MVPA(par(i), task, saveit)
    end
    return
end
      
if ~isstruct(par)
    par = LT_Params(par, task,1);
end


%get <idx> for appropriate task
switch par.task
    case 'lab'
        [~, ~, idx] = LTLab_fMRIBehAnalysis_Ret(par);
    case 'cam'
        [~, ~, idx] = LTCam_fMRIBehAnalysis_Ret(par);
end

idxFnames = fieldnames(idx);

onsets = cell(size(fieldnames(idx)));
names = cell(size(fieldnames(idx)));
durations = cell(size(fieldnames(idx)));

for xx = 1:length(idxFnames)
    fName = idxFnames{xx};
    idx.thisOns = idx.(fName);
    onsets{xx} = idx.allTrialOns(find(idx.thisOns));
    names{xx} = fName;
    durations{xx} = STIMDUR;
end

mvpadir_exists = exist(par.mvpadir,'dir');
if ~mvpadir_exists, mkdir(par.mvpadir); end

if saveit, save(fullfile(par.mvpadir, 'mvpa_ons.mat'),'onsets','durations','names'); end