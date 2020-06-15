close all;
clear all;
clc;

load fve32.mat; % bd network
node = 32;
% lattice model
ITER=50;
dismatrix = distance_bin(CIJ);%distance matrix
[Rlatt,Rrp,ind_rp,eff] = latmio_dir_connected(CIJ,ITER,dismatrix);
disp(Rlatt);%latticized network
%disp(eff); %number of actual rewirings carried out

adj = Rlatt;

cluster = clustering_coef_bd(adj);%clustering coefficient
%disp(['Clustering coefficient:' cluster]);
disp(sum(cluster)/30);
group_num = 6;
[nnn, xout] = hist(cluster,group_num);
%disp(sum(sum(dismatrix))/node);% average path length
nnn = nnn / sum(nnn);
figure;
bar(xout, nnn);
xlabel('clustering coef c');ylabel('probability p(c)');title('lattice clustering coef distribution');

gloeffi =efficiency_bin(adj); %global efficiency
loeffi =efficiency_bin(adj,1); %local efficiency
disp(gloeffi);
disp(sum(loeffi)/node);

nodecentr=betweenness_bin(adj); %node betweenness centrality
disp(sum(nodecentr)/node);
edgecentr = edge_betweenness_bin(adj); %edge betweenness centrality
%[EBC BC] = edge_betweenness_bin(adj);
disp(sum(sum(edgecentr))/(node*node));

dismatrix = distance_bin(adj);%distance matrix
%disp(dismatrix);
[lambda,efficiency,ecc,radius,diameter] = charpath(dismatrix);%path length
disp(lambda); %characteristic path length
group_num = 6;
[nn, xout] = hist(dismatrix,group_num);
disp(sum(sum(dismatrix))/node);% average path length
nn = nn / sum(nn);
figure;
bar(xout, nn);
xlabel('shortest path length l');ylabel('probability p(l)');title('lattice path length distribution');

[id,od,deg] = degrees_dir(adj);%Degree distribution
%disp(deg);
group_num = 6;
[n, xout] = hist(deg,group_num);
disp(sum(deg)/node);% average degree
n = n / sum(n);
figure;
bar(xout, n);
xlabel('degree k');ylabel('probability p(k)');title('lattice node degree(indegree & outdegree) distribution');