close all;
clear all;
clc;

load mac95.mat; % bd network
adj = matrix;
% random model
%ITER=100;
%[random_matrix,eff]=randmio_dir_connected(adj, ITER);

gamma_dis = 1:0.1:2;
class_matrix = zeros(94,11);
index = 1;
for gamma = 1:0.1:2
i = 1;
[class_plot,value_plot] = modularity_dir(adj,gamma);     % get community assignments
%[class_plot_ran,value_plot_ran] = modularity_dir(random_matrix,gamma);     % get community assignments
while i < 95
class_matrix(i,index) = class_plot(i);
i = i + 1;
end
index = index + 1;
end

figure(1);
x = 1:0.1:2;
y = 1:94;
imagesc(x,y,class_matrix);