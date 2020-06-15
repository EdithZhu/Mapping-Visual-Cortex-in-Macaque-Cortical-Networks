close all;
clear all;
clc;

load mac95.mat; % bd network
load fve30.mat;
adj = matrix;
%adj = CIJ;

disp('---------cluster--------')
cluster = clustering_coef_bd(adj);
avg_cluster = sum(cluster)/length(adj)
[h,pvalue,ci]=ttest(cluster,avg_cluster)
%[hh,pp] = lillietest(cluster)
disp('---------degree--------')
[id,od,degree] = degrees_dir(adj);%Degree distribution
degree = degree /2;
avg_degree = sum(degree)/length(adj)
[h,pvalue,ci]=ttest(degree,avg_degree)
%[hh,pp] = lillietest(degree)
disp('---------path length--------')
dismatrix = distance_bin(adj);%distance matrix
[lambda,dv,efficiency,ecc,radius,diameter] = charpath1(dismatrix,0,0);%path length
char_path_length = lambda %characteristic path length
[h,pvalue,ci]=ttest(dv,char_path_length)
%[hh,pp] = lillietest(dv)
disp('---------distance--------')
[modu_class,modu_value] = modularity_dir(adj,1.4);     % get community assignments
partitions = modu_class;
locations = position;
[ comm_ave_pairwise_spatial_dist_array,store_dist ] = comm_ave_pairwise_spatial_dist1(partitions,locations);
edge_num = length(adj)*(length(adj)+1)/2;
distance = store_dist;
%avg_distance = sum(distance)/edge_num
avg_distance = mean(distance)
[h,pvalue,ci]=ttest(distance,avg_distance)
%[hh,pp] = lillietest(distance)