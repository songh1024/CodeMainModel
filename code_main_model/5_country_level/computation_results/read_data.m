clear all;
%% Parameters
fileID = fopen('parameters_matlab.txt','r');
parameters_matlab=fscanf(fileID,'%s',[1, inf]);
read_parameter(parameters_matlab);

%% read bank
fileID = fopen('bank.txt','r');
bank=fscanf(fileID,'%f,',[1, inf]);
bank=read1(bank, nmkt);
bank=bank';

%% read d
fileID = fopen('d.txt','r');
d=fscanf(fileID,'%f,',[1, inf]);
d=read2(d, period_trans, nmkt);

%% read bank_opt
fileID = fopen('bank_opt.txt','r');
bank_opt=fscanf(fileID,'%f,',[1, inf]);
bank_opt=read1(bank_opt, nmkt);

%% read d_opt
fileID = fopen('d_opt.txt','r');
d_opt=fscanf(fileID,'%f,',[1, inf]);
d_opt=read2(d_opt, period_trans, nmkt);

%% read n_opt
fileID = fopen('n_opt.txt','r');
n_opt=fscanf(fileID,'%f,',[1, inf]);
n_opt=read2(n_opt, nyear, nseg);

%% read n_cum_opt
fileID = fopen('n_cum_opt.txt','r');
n_cum_opt=fscanf(fileID,'%f,',[1, inf]);
n_cum_opt=read2(n_cum_opt, nyear, nseg);


save('data.mat');
fclose('all');
