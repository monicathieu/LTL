function LT_MakeRegs(par, task, analysis, saveit)
%LT_MakeRegs(par, saveit) creates onsets and regresors (boxcar and motion)
%for use in SPM's GLM.
%inputs:
%par = subject-specific parameters created by LT_Params
%saveit = 1 if we'd like to save mat files for onsets and boxcar + motion
%regressors(default is 1)
%
%5/15/14 - currently only works for 'analysisByPerf_noPhase'. Will change
%so that it can make regressors for several different models.

STIMDUR = 0;

if nargin<4
    saveit=1;
end

if ~isstruct(par)
    par = LT_Params(par, task,1);
end

par.stimDur = STIMDUR;

thisAnalysisDir = fullfile(par.analysisdir, analysis);

switch par.task
    case 'lab'
        [onsets, names, durations] = makeLabOns(par,analysis);
    case 'cam'
        error('whoops! don''t know how to do this yet');
end

R = makeMotRegs(par);

if ~exist(thisAnalysisDir)
    mkdir (thisAnalysisDir);
end

cd(thisAnalysisDir);
if saveit
    save ons.mat onsets durations names;
    save regs.mat R
end
end



function R = makeMotRegs(par)

sessReg = zeros(sum(par.(par.task).numvols),length(par.(par.task).numvols) -1);
for i =1:(length(par.(par.task).numvols) -1)
    sessReg((i-1)*par.(par.task).numvols(i)+1:i*par.(par.task).numvols(i),i) = ones(par.(par.task).numvols(i),1);
end

motRegs = [];
for runNum =1:length(par.(par.task).numvols)
    fid = fopen(fullfile(par.funcdir, sprintf('test%02d/rp_avol0005.txt',runNum)));
    motRegs = vertcat(motRegs,cell2mat(textscan(fid,'%n%n%n%n%n%n')));
end

R = horzcat(sessReg,motRegs);
end