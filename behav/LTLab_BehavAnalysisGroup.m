function [resG, psyG, idxG] = LTLab_BehavAnalysisGroup(sarr)
%[M, groupPsy, idxG] = LTLab_BehavAnalysisGroup(sarr) aggregates behavioral
%data for the Long Term Lab Experiment
% inputs:
%   <sarr> = array of subjects to be included in group analysis
% outputs:
% <resG> = group aggregate of subject-specific results struct
% <psyG> = group aggregate of subject-specific psyphys struct
% <idxG> = group aggregate of subject-specific idx struct 

resG = []; psyG = []; idxG = []; 
for snum = 1:length(sarr)

    par = LT_Params(sarr(snum), 'lab', 0);
    
    [res, psyphys, idxG.sub{snum}] = LTLab_fMRIBehAnalysis_Ret(par);
    idxG.sub{snum}.subNo = par.subNo;
    
    fnRet = fieldnames(res);
    
    for f = 1:length(fnRet)
        resG.res.vals(snum, f) = res.(fnRet{f});
        resG.res.names = fnRet';
    end
    
    psyG(snum).dat = psyphys;
%     psyG(snum).goodSub = par.goodSub;
end

    
%group plot of psychophysics data    
LTLab_GroupPsyPhysPlot(psyG);
    


end

