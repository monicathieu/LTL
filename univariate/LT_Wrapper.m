function LT_Wrapper(subpar, flags)
if isstruct(subpar)
    par = subpar;
else
    par = LT_Params(subpar);
end

if ismember('r', flags); LT_MakeRegs(par);end
