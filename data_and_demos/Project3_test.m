close all;
clear all;
clc;

load cat.mat; % wd network
%****************************mac95.mat*********************
step_end = 3;
step = 0.1
adj = CIJctx;
adj(adj>0)=1;
writetoPAJ(adj, 'cat' , 1); %get data for a edges.csv