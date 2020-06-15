close all;
clear all;
clc;

load celegans131.mat; % bd network
load celegans277.mat; % bd network
load mac95.mat; % bd network
load cat.mat; % wd network
adj1 = CIJctx;
adj1(adj1>0)=1;
%************************mac95.mat****************************
AIJ = matrix;                                % load your adjacency matrix
COOR = position;                               % load 3D coordinates for each node
[x,y,z] = adjacency_plot_und(AIJ,COOR);  % call function
%************************celegans131.mat****************************
AIJ1 = celegans131matrix;                                % load your adjacency matrix
COOR1 = celegans131positions;                               % load 3D coordinates for each node
[x1,y1] = adjacency_plot_und(AIJ1,COOR1);  % call function
%************************celegans277.mat****************************
AIJ2 = celegans277matrix;                                % load your adjacency matrix
COOR2 = celegans277positions;                               % load 3D coordinates for each node
[x2,y2] = adjacency_plot_und(AIJ2,COOR2);  % call function
figure(1);
subplot(1,3,1);plot3(x,y,z,'b-','MarkerSize',10);grid;
xlabel('x axis');ylabel('y axis');zlabel('z axis');title('3D macaque cortical connectivity network(mm)');legend('connectivity');
subplot(1,3,2);plot(x,y,'g-','MarkerSize',10);grid;
xlabel('x axis');ylabel('y axis');title('2D C. elegans local network of 131 frontal neurons(mm)');legend('connectivity');
subplot(1,3,3);plot(x,y,'r-','MarkerSize',10);grid;
xlabel('x axis');ylabel('y axis');title('2D C. elegans global network of 277 neurons(mm)');legend('connectivity');
% hold on;
% for k=1:94
% text(x(k),y(k),labels(k));
% end

%***************************************************************************************
adj = matrix;%mac95
COOR_new = COOR;%mac95
% adj = celegans131matrix;%celegans131
% COOR_new = celegans131positions;%celegans131
% adj = celegans277matrix;%celegans277
% COOR_new = celegans277positions;%celegans277
[modu_class,modu_value] = modularity_dir(adj,1.4);     % get community assignments
[X,Y,INDSORT] = grid_communities(modu_class); % call function
figure(2);
imagesc(adj(INDSORT,INDSORT));           % plot ordered adjacency matrix
hold on;                                 % hold on to overlay community visualization
plot(X,Y,'w','linewidth',2);             % plot community boundaries
title('Observed Reordered');

%*******************************************************************************************
% adj = matrix;
% n = length(adj);
% M  = 1:n;                   % initial community affiliations
% %Q0 = -1; Q1 = 0;            % initialize modularity values
% %while Q1-Q0>1e-5;           % while modularity increases
% for loop = 1:100
% %Q0 = Q1;                % perform community detection
% %[M, Q1] = community_louvain(adj, [], M);
% [M, Q1] = community_louvain(adj,1.4, M);
% end
% [A,B,INDEX] = grid_communities(M); % call function
% figure(4);
% imagesc(adj(INDEX,INDEX));           % plot reordered adjacency matrix
% colorbar;
% hold on;                                 % hold on to overlay community visualization
% plot(A,B,'w','linewidth',2);             % plot community boundaries
% xlabel('Sequence of Reordered Nodes');ylabel('Sequence of Reordered Nodes');
% title('reordered Network');

partitions = modu_class;
locations = COOR_new;
disp('community average pairwise spatial distance:'); 
disp(comm_ave_pairwise_spatial_dist(partitions,locations));
%disp('community laterality:');
%disp(comm_laterality(partitions,locations));
disp('community radius:'); 
disp(comm_radius(partitions,locations));
disp('community spatial diameter:'); 
disp(comm_spatial_diameter(partitions,locations));
disp('community spatial extent:'); 
disp(comm_spatial_extent(partitions,locations));

%********************************************************************************************
step_end = 3;
step = 0.1;
gamma_dis = 0:step:step_end;
counter = zeros(step_end/step +1,1);
corevalue = zeros(step_end/step +1,1);
error1 = zeros(step_end/step +1,1);
error2 = zeros(step_end/step +1,1);
k=1;
for gamma = 0:step:step_end
[Core,coreval] = core_periphery_dir(adj,gamma);
counter(k) = sum(Core(:)==1);
corevalue(k) = coreval;
error1(k) = 0.025 * counter(k);
error2(k) = 0.025 * corevalue(k);
k = k+1;
end

counter1 = zeros(step_end/step +1,1);
corevalue1 = zeros(step_end/step +1,1);
error11 = zeros(step_end/step +1,1);
error12 = zeros(step_end/step +1,1);
k1=1;
for gamma1 = 0:step:step_end
[Core1,coreval1] = core_periphery_dir(adj1,gamma1);
counter1(k1) = sum(Core1(:)==1);
corevalue1(k1) = coreval1;
error11(k1) = 0.025 * counter1(k1);
error12(k1) = 0.025 * corevalue1(k1);
k1 = k1+1;
end
figure(3);
subplot(1,2,1); errorbar(gamma_dis,counter,error1,'-k','linewidth',2);grid;hold on;
errorbar(gamma_dis,counter1,error11,'--b','linewidth',2);
xlabel({'gamma','(a)'});ylabel('# of core');title('number of detected cores');legend('mac94','cat52');
subplot(1,2,2); errorbar(gamma_dis,corevalue,error2,'-k','linewidth',2);grid;hold on;
errorbar(gamma_dis,corevalue1,error12,'--b','linewidth',2);
xlabel({'gamma','(b)'});ylabel('maximized coreness statistic');title('coreness statistic');legend('mac94','cat52');

[Core2,corevalue2] = core_periphery_dir(adj,1.4);
[A,B,index] = grid_communities(Core2); % call function
% Core22 = reshape(Core2,[94,1]);
disp(sum(Core1(:)==1));
figure(4);
imagesc(adj(index,index));           % plot ordered adjacency matrix
hold on;                                 % hold on to overlay community visualization
plot(A,B,'r','linewidth',2);             % plot community boundaries