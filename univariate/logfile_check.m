function logfile_check(S)
if nargin  == 0, S = [1:16];end
FDIR = '/Volumes/Tyler_Drive1/LTLab/fmri_data/';
RES1DIR_FMT = 'LTL%03d/univariateQA/results01_short/LTL%03d_onsets_model01_short.mat';
RES2DIR_FMT = 'LTL%03d/univariateQA/results02_0s/LTL%03d_onsets_model02_0s.mat';
EMPTYONSET = zeros(480,1);

for s = S
    s_id = sprintf('LTL%03d',s);
    
    s_res1dir = sprintf(RES1DIR_FMT,s,s);
    %load handmade logfile
    lf = fullfile(FDIR,s_res1dir);
    res1 = load(lf);   
    
    
    s_res2dir = sprintf(RES2DIR_FMT,s,s);
    %load handmade logfile
    lf = fullfile(FDIR,s_res2dir);
    res2 = load(lf);
      
    %run LTLab_fMRIBehAnalysis_Ret
    [~, ~, idx] = LTLab_fMRIBehAnalysis_Ret(s);
    
    subplot(1,3,1); 
    srcmsocounts{1} = res1.onsets{8};
    srcmsocounts{2} = get_cond_ons(idx, idx.SCENEsrcInc .* idx.respObj)
    srcmsocounts
    show_onset_diffs(srcmsocounts);
    subplot(1,3,2);
    counts{1} = res2.onsets{12};
    counts{2} = get_cond_ons(idx,idx.srcFA);
    counts
    show_onset_diffs(counts);
        subplot(1,3,3);

    crcounts{1} = res1.onsets{16};
    crcounts{2} = get_cond_ons(idx, idx.CR)
    show_onset_diffs(crcounts);
colorbar;
%     
     keyboard;
    %clf;
    
end    
    
end

function o = get_cond_ons(idx,x) 
o=sort(idx.allTrialOns(find(x))');
end

function show_onset_diffs(counts)
counts{1} = get_logical(counts{1});
counts{2} = get_logical(counts{2});
imagesc((counts{1}-counts{2}));
end

function o = get_logical(a)
o = zeros(480,1);
o(find(a)) = 1;
end
