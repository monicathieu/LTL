function LT_CheckDesignLev1(task, analysis, subarray)
if nargin < 3
    subarray = 1:16
end

for i = subarray
    par = LTC_Params(i, task, 0);
    cd([par.analysisdir '/' analysis]);

load SPM
figure(1)
imagesc(SPM.xX.X)
colormap('gray')
caxis([-.1 .5])
title(sprintf('Design Matrix for LTL%03d', i))
spm_fMRI_design_show(SPM,1)
pause
figure(1)
clf
end