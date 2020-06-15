close all;
clear all;
clc;

load mac95.mat; % bd network
load fve30.mat; % bd network
load cat.mat; % bd network
load celegans131.mat; % bd network
load celegans277.mat; % bd network

ran_matrix = matrix;
% random model
ITER=100;
[random_matrix,eff]=randmio_dir_connected(ran_matrix, ITER);
%****************************mac95.mat*********************
step_end = 3;
step = 0.1;
adj = matrix;
%visual cortex in mac95
index = [31,37,42,43,44,45,46,47,48,56,62,65,67,68,69,70,71,76,77,78,79,83,84,85,86,87,88,91,94];
index_store = zeros(29,1);
index_class = zeros(step_end/step +1,1);
error_class = zeros(step_end/step +1,1);
error1_class = zeros(step_end/step +1,1);
error_sum = zeros(step_end/step +1,1);
gamma_dis = 0:step:step_end;
class = zeros(step_end/step +1,1);
value = zeros(step_end/step +1,1);
size = zeros(step_end/step +1,1);
error1 = zeros(step_end/step +1,1);
error2 = zeros(step_end/step +1,1);
error3 = zeros(step_end/step +1,1);
errora = zeros(step_end/step +1,1);
errorb = zeros(step_end/step +1,1);
errorc = zeros(step_end/step +1,1);
errord = zeros(step_end/step +1,1);
optimal = zeros(step_end/step +1,1);
k=1;
for gamma = 0:step:step_end
[class_plot,value_plot] = modularity_dir(adj,gamma);     % get community assignments
class(k) = max(class_plot);
value(k) = value_plot;
size(k) = 94/max(class_plot);
error1(k) = 0.025 * class(k);
error2(k) = 0.025 * value(k);
error3(k) = 0.025 * size(k);
for index_index = 1:29
    index_store(index_index) = class_plot(index(index_index));%save the # of community of visual neurons in mac94
    index_index = index_index + 1;
end
result30 = unique(index_store);%communnity number
count30 = hist(index_store,unique(index_store));% # of nodes in one module
[index_class(k),ind] = max(count30);
errord(k) = 0.025 * index_class(k);
%result94 = unique(class_plot);%communnity number
%count94 = hist(class_plot,unique(class_plot));
indexdev = find(class_plot == result30(ind));
%error_class(k) = ( length(indexdev) - max(count30) )/30;
error_class(k) =  length(indexdev) - 29;
error1_class(k) =  max(count30) - 29;
error_sum(k) = abs(error_class(k)) + abs(error1_class(k));
errora(k) = 0.025 * error_class(k);
errorb(k) = 0.025 * error1_class(k);
errorc(k) = 0.025 * error_sum(k);
disp(max(count30));
k = k+1;
end
%**********************************fve30.mat****************************
adj1 = CIJ;
gamma_dis1 = 0:step:step_end;
class1 = zeros(step_end/step +1,1);
value1 = zeros(step_end/step +1,1);
size1 = zeros(step_end/step +1,1);
error11 = zeros(step_end/step +1,1);
error12 = zeros(step_end/step +1,1);
error13 = zeros(step_end/step +1,1);
k1=1;
for gamma1= 0:step:step_end
[class_plot1,value_plot1] = modularity_dir(adj1,gamma1);     % get community assignments
class1(k1) = max(class_plot1);
value1(k1) = value_plot1;
size1(k1) = 30/max(class_plot1);
error11(k1) = 0.025 * class1(k1);
error12(k1) = 0.025 * value1(k1);
error13(k1) = 0.025 * size1(k1);
k1 = k1+1;
end

%**********************************cat.mat****************************
adj2 = random_matrix;
gamma_dis2 = 0:step:step_end;
class2 = zeros(step_end/step +1,1);
value2 = zeros(step_end/step +1,1);
size2 = zeros(step_end/step +1,1);
error21 = zeros(step_end/step +1,1);
error22 = zeros(step_end/step +1,1);
error23 = zeros(step_end/step +1,1);
k2=1;
for gamma2= 0:step:step_end
[class_plot2,value_plot2] = modularity_dir(adj2,gamma2);     % get community assignments
class2(k2) = max(class_plot2);
value2(k2) = value_plot2;
size2(k2) = 94/max(class_plot2);
error21(k2) = 0.025 * class2(k2);
error22(k2) = 0.025 * value2(k2);
error23(k2) = 0.025 * size2(k2);
k2 = k2+1;
end

%**********************************celegans131.mat****************************
adj3 = celegans131matrix;
gamma_dis3 = 0:step:step_end;
class3 = zeros(step_end/step +1,1);
value3 = zeros(step_end/step +1,1);
size3 = zeros(step_end/step +1,1);
error31 = zeros(step_end/step +1,1);
error32 = zeros(step_end/step +1,1);
error33 = zeros(step_end/step +1,1);
k3=1;
for gamma3= 0:step:step_end
[class_plot3,value_plot3] = modularity_dir(adj3,gamma3);     % get community assignments
class3(k3) = max(class_plot3);
value3(k3) = value_plot3;
size3(k3) = 131/max(class_plot3);
error31(k3) = 0.025 * class3(k3);
error32(k3) = 0.025 * value3(k3);
error33(k3) = 0.025 * size3(k3);
k3 = k3+1;
end

%**********************************celegans277.mat****************************
adj4 = celegans277matrix;
gamma_dis4 = 0:step:step_end;
class4 = zeros(step_end/step +1,1);
value4 = zeros(step_end/step +1,1);
size4 = zeros(step_end/step +1,1);
error41 = zeros(step_end/step +1,1);
error42 = zeros(step_end/step +1,1);
error43 = zeros(step_end/step +1,1);
k4=1;
for gamma4= 0:step:step_end
[class_plot4,value_plot4] = modularity_dir(adj4,gamma4);     % get community assignments
class4(k4) = max(class_plot4);
value4(k4) = value_plot4;
size4(k4) = 277/max(class_plot4);
error41(k4) = 0.025 * class4(k4);
error42(k4) = 0.025 * value4(k4);
error43(k4) = 0.025 * size4(k4);
k4 = k4+1;
end

figure(1);
subplot(1,3,1); errorbar(gamma_dis,class,error1,'-k','linewidth',2);grid;hold on;
errorbar(gamma_dis1,class1,error11,'--b','linewidth',2);grid;hold on;
errorbar(gamma_dis2,class2,error21,'-.m','linewidth',2);grid;hold on;
errorbar(gamma_dis3,class3,error31,':c','linewidth',2);grid;hold on;
errorbar(gamma_dis4,class4,error41,'-y','linewidth',2);
xlabel({'gamma','(a)'});ylabel('# of communities');title('number detected communities');legend('mac94','fve30','mac94(random model)','celegans131','celegans277');

subplot(1,3,2); errorbar(gamma_dis,value,error2,'-k','linewidth',2);grid;hold on;
errorbar(gamma_dis1,value1,error12,'--b','linewidth',2);grid;hold on;
errorbar(gamma_dis2,value2,error22,'-.m','linewidth',2);grid;hold on;
errorbar(gamma_dis3,value3,error32,':c','linewidth',2);grid;hold on;
errorbar(gamma_dis4,value4,error42,'-y','linewidth',2);
xlabel({'gamma','(b)'});ylabel('Modularity(Q)');title('Louvain maximized modularity');legend('mac94','fve30','mac94(random model)','celegans131','celegans277');

subplot(1,3,3); errorbar(gamma_dis,size,error3,'-k','linewidth',2);grid;hold on;
errorbar(gamma_dis1,size1,error13,'--b','linewidth',2);grid;hold on;
errorbar(gamma_dis2,size2,error23,'-.m','linewidth',2);grid;hold on;
errorbar(gamma_dis3,size3,error33,':c','linewidth',2);grid;hold on;
errorbar(gamma_dis4,size4,error43,'-y','linewidth',2);
xlabel({'gamma','(c)'});ylabel('# of nodes');title('Average community size');legend('mac94','fve30','mac94(random model)','celegans131','celegans277');

figure(2);
subplot(1,2,1);errorbar(gamma_dis,index_class,errord,'-k','linewidth',2);grid;
xlabel('gamma');ylabel('# of visual areas in one module');title('# of areas');legend('# of areas');
subplot(1,2,2);errorbar(gamma_dis,error_class,errora,'-b','linewidth',2);grid;
hold on;errorbar(gamma_dis,error1_class,errorb,'--m','linewidth',2);
hold on;errorbar(gamma_dis,error_sum,errorc,'--*k','linewidth',3);
xlabel('gamma');ylabel('errors in visual module');title('Distribution of errors');legend('global error','local error','sum error');
%***************************************************************************************
[modu_class,modu_value] = modularity_dir(adj,1.4);     % get community assignments
[X,Y,INDSORT] = grid_communities(modu_class); % call function
figure(3);
imagesc(adj(INDSORT,INDSORT));           % plot ordered adjacency matrix
colorbar;
hold on;                                 % hold on to overlay community visualization
plot(X,Y,'w','linewidth',2);             % plot community boundaries
xlabel({'Sequence of Reordered Nodes','(b)'},'FontSize',12);ylabel('Sequence of Reordered Nodes','FontSize',12);
title('Observed Reordered','FontSize',12);
%**************************************************************************************
% Iterative community finetuning.
% adj is the input connection matrix.
%n  = size(adj,1);             % number of nodes
n = length(adj);
M  = 1:n;                   % initial community affiliations
Q0 = -1; Q1 = 0;            % initialize modularity values
while Q1-Q0>1e-5;           % while modularity increases
Q0 = Q1;                % perform community detection
%[M, Q1] = community_louvain(adj, [], M);
[M, Q1] = community_louvain(adj,1.4, M);
end
[A,B,INDEX] = grid_communities(M); % call function
figure(4);
imagesc(adj(INDEX,INDEX));           % plot reordered adjacency matrix
colorbar;
hold on;                                 % hold on to overlay community visualization
plot(A,B,'w','linewidth',2);             % plot community boundaries
xlabel('Sequence of Reordered Nodes');ylabel('Sequence of Reordered Nodes');
title('reordered Network');
%******************************************************************************************
writetoPAJ(adj, 'mac94' , 1); %get data for a edges.csv

%[C,q] = core_periphery_dir(W,gamm,C0)
% [Core,corevalue] = core_periphery_dir(adj,1.5)% core periphery
% [A,B,index] = grid_communities(Core); % call function
% figure(2);
% imagesc(adj(index,index));           % plot ordered adjacency matrix
% hold on;                                 % hold on to overlay community visualization
% plot(A,B,'r','linewidth',2);             % plot community boundaries
%******************************************************************************************
% Iterative community finetuning.
% adj is the input connection matrix.
%n  = size(adj,1);             % number of nodes
figure(5);
subplot(1,2,1);
imagesc(adj);           % plot reordered adjacency matrix
colorbar;
xlabel({'Sequence of Original Nodes','(a)'});ylabel('Sequence of Original Nodes');
title('Original Network');
subplot(1,2,2);
imagesc(adj(INDSORT,INDSORT));           % plot ordered adjacency matrix
colorbar;
hold on;                                 % hold on to overlay community visualization
plot(X,Y,'w','linewidth',2);             % plot community boundaries
xlabel({'Sequence of Reordered Nodes','(b)'});ylabel('Sequence of Reordered Nodes');
title('Observed Reordered');