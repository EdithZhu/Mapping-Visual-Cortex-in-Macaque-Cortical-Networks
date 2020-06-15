close all;
clear all;
clc;

load mac95.mat; % bd network
adj = matrix;
cluster = clustering_coef_bd(adj);
avg_quality = sum(quality_array)/(loop+1)
[h,pvalue,ci]=ttest(quality_array,avg_quality)
[hh,pp] = lillietest(quality_array)
