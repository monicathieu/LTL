function ltMakeRegsMvpa(par, task, saveit)
%LTL_MakeRegs(par, saveit) creates onsets and regresors (boxcar and motion)
%for use in SPM's GLM.
%inputs:
%par = subject-specific parameters created by LTL_Params
%saveit = 1 if we'd like to save mat files for onsets and boxcar + motion
%regressors(default is 1)
%

analysis = [task 'MVPA'];

if nargin<3
    saveit=1;
end

if length(par) > 1
    for i =1:length(par)
        ltMakeRegsMvpa(par(i), task)
    end
    
    return
end
        
if ~isstruct(par)
    par = LT_Params(par, task,1);
end
cd (par.behavdir);

thisAnalysisDir = fullfile(par.analysisdir, analysis);

stimDur = 0;

switch par.task
    case 'lab'
        [~, ~, idx] = LTLab_fMRIBehAnalysis_Ret(par);
    case 'cam'
        [~, ~, idx] = LTCam_fMRIBehAnalysis_Ret(par);
    case 'loc'
        [~, ~, idx] = LTLoc_fMRIBehAnalysis_Ret(par);
end


i=0;
idxFnames = fieldnames(idx);
for xx  = 1:length(idxFnames)
    i = i+1;
    fName = idxFnames{xx};
    idx.thisOns = idx.(fName);
    stimOnsets{i} = idx.allTrialOns(find(idx.thisOns));
    stimNames{i} = fName;
    stimDurations{i} = stimDur;
end

if sum(idx.junk>0)
    i= i+1;
    stimOnsets{i} = idx.allTrialOns(find(idx.junk));
    stimNames{i} = 'junk';
    stimDurations{i} = 0;
end

onsets = stimOnsets;
names = stimNames;
durations = stimDurations;

if ~exist(thisAnalysisDir)
    mkdir (thisAnalysisDir);
end

sessReg = zeros(sum(par.(task).numvols),length(par.(task).numvols) -1);
runStart = 0;
runEnd = 0;
for i =1:(length(par.(task).numvols) -1)
    runStart = runEnd + 1;
    runEnd = runStart+par.(task).numvols(i)-1;
    sessReg(runStart:runEnd,i) = ones(par.(task).numvols(i),1);
end

motRegs = [];
for runNum =1:length(par.(par.task).numvols)
    
    fid = fopen(fullfile(par.funcdir, sprintf('test%02d/rp_avol0005.txt',runNum)));
    motRegs = vertcat(motRegs,cell2mat(textscan(fid,'%n%n%n%n%n%n')));
end
fclose('all');

if strcmp(par.task,'cam')
    switch par.subNo
        case 1
            motRegs = motRegs(1:size(sessReg,1),1:size(sessReg,2));
        case 16
            display('You''re doing anything to fix motion for s16, is this okay?')
    end
end
R = horzcat(sessReg,motRegs);


cd(thisAnalysisDir);
if saveit
    save mvpa_ons.mat onsets durations names;
    
%     save regs.mat R
end