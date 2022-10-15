%% Find the line

%%
clc
clear

%%
T = [18.53 22.47 26.87 30.05 35.87 38.56 41.50 46.26];
dist = 1.375:1:8.375;

[a] = polyfit(dist,T,1);

T_0 = polyval(a,0);