function [sa SA] = LT_SA()

sa.sa_all = 1:16;
for i = sa.sa_all
    SA{i} = sprintf('LTL%03d',i);
end
end