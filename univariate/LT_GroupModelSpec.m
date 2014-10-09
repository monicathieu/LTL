function LT_GroupModelSpec(gpar, analysis)
%function LT_GroupModelSpec(gpar)
origdir = pwd;
if ~isstruct(gpar)
    gpar = LT_GroupParams(gpar, analysis);
end

clear jobs
load(char(gpar.modelTemplate));

yTemplate = load(gpar.conTemplate);
yCons = yTemplate.SPM.xCon;

for cNum= 1:length(yCons);
    
    jobs{cNum}.stats{1}.factorial_design.des.t1.scans = [];
    jobs{cNum}.stats{1}.factorial_design.dir = gpar.cons{cNum}.dir;
    
    if ~exist(jobs{cNum}.stats{1}.factorial_design.dir{1});
        mkdir(jobs{cNum}.stats{1}.factorial_design.dir{1});
    end
    
    %jobs{c}.stats{1}.factorial_design.dir = {fullfile(expt_dir, 'group_analyses', tasks{t}, yCons(c).name)};
    %jobs{c}.stats{1}.factorial_design.dir = {fullfile(gpar.expt_dir, 'group_analyses', gpar.tasks{t}, gpar.task{t}.SPMcons(c).name)};
    
    if isfield(gpar,'exMask') & ~isempty(gpar.exMask)
        jobs{cNum}.stats{1}.factorial_design.masking.em{1} = gpar.exMask;
    end
    
    if isfield(gpar, 'covVec') & ~isempty(gpar.covVec)
        jobs{cNum}.stats{1}.factorial_design.cov(1).cname = gpar.covName{cNum};
        jobs{cNum}.stats{1}.factorial_design.cov(1).c = gpar.covVec{cNum};
        jobs{cNum}.stats{1}.factorial_design.cov(1).iCC = 1;
        jobs{cNum}.stats{1}.factorial_design.cov(1).iCFI = 1;
    end
    
    for s = 1:length(gpar.cons{cNum}.scans)
        
        
        %jobs{c}.stats{1}.factorial_design.des.t1.scans{s} = fullfile(gpar.expt_dir, gpar.subjArray{s}, gpar.conGroup{t}, ['con_' prepend(num2str(c), 4) '.img']);
        jobs{cNum}.stats{1}.factorial_design.des.t1.scans{s} = gpar.cons{cNum}.scans{s};
    end
    
end
spm_jobman('run',jobs)
end
