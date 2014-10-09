function LT_CheckDesignLev1(task, analysis)
for i = 1:16
    par = LTC_Params(i, task, 0);
    cd([par.analysisdir '/' analysis]);

load SPM
imagesc(SPM.xX.X)
colormap('gray')
caxis([-.1 .5])
title(sprintf('Design Matrix for LTL%03d', i))
pause
end