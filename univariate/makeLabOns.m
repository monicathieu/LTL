function [onsets, names, durations] = makeLabOns(par, analysis)

[~, ~, idx] = LTLab_fMRIBehAnalysis_Ret(par);

perfSet = {'srcHit', 'hitSrcMiss' 'itemHit', 'miss', 'CR', 'srcFA', 'noSrcFA'};
srcSet = {'obj', 'face', 'scene', 'new'};
respSet = {'respObj', 'respFace', 'respScene', 'respOld' ,'respNew'};

switch analysis
    case 'analysisByPerf_noPhase'
        i = 0;
        for p = 1:length(perfSet)
            for srcNum = 1:length(srcSet)
                for respNum = 1:length(respSet)
                    idx.thisOns = idx.(perfSet{p}) .* idx.(srcSet{srcNum}) .*  idx.(respSet{respNum});
                    
                    if ~isempty(idx.allTrialOns(find(idx.thisOns)))
                        i = i+1;
                        fName = sprintf('%s_%s_%s', perfSet{p}, srcSet{srcNum}, respSet{respNum});
                        onsets{i} = idx.allTrialOns(find(idx.thisOns));
                        names{i} = fName;
                        durations{i} = par.stimDur;
                    end
                end
            end
        end
        
    case 'ByPerf_noPhase_simple'
        i = 0;
        for p = 1:length(perfSet)
            
            idx.thisOns = idx.(perfSet{p});
            if ~isempty(idx.allTrialOns(find(idx.thisOns)))
                i = i+1;
                fName = sprintf('%s', perfSet{p});
                onsets{i} = idx.allTrialOns(find(idx.thisOns));
                names{i} = fName;
                durations{i} = par.stimDur;
            end
        end
        
        
    case 'analysisByPerfAndPhase'
        i = 0;
        phaseName = {'', 'Phase1_','Phase2_','Phase3_','Phase4_'}
        for p = 1:length(perfSet)
            for phase = 0:4
                idx.thisOns = (idx.phase == phase) .* idx.(perfSet{p});
                
                if ~isempty(idx.allTrialOns(find(idx.thisOns)))
                    i = i+1;
                    fName = sprintf('%s%s', phaseName{phase+1}, perfSet{p});
                    onsets{i} = idx.allTrialOns(find(idx.thisOns));
                    names{i} = fName;
                    durations{i} = par.stimDur;
                end
            end
        end
end

if sum(idx.junk>0)
    i= i+1;
    onsets{i} = idx.allTrialOns(find(idx.junk));
    names{i} = 'junk';
    durations{i} = 0;
end
end
