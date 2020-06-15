clustering = [];
k = 0;
while k<10
ITER = 500;
[Random,eff] = randmio_dir_connected(adj, ITER); % random model

clustering = [clustering Random];
k = k + 1;
end

C_R_mean = mean(clustering);