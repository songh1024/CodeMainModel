clear all;
%% Parameters
fileID = fopen('parameters_matlab.txt','r');
parameters_matlab=fscanf(fileID,'%s',[1, inf]);
read_parameter(parameters_matlab);

%% read Pi
fileID = fopen('Pi.txt','r');
Pi=fscanf(fileID,'%f,',[1, inf]);
Pi=read1(Pi, nmkt);

%% read Z
fileID = fopen('Z.txt','r');
Z=fscanf(fileID,'%f,',[1, inf]);
Z=read1(Z, nmkt);

%% read d
fileID = fopen('d.txt','r');
d=fscanf(fileID,'%f,',[1, inf]);
d=read2(d, period_trans, nmkt);

%% read zeta
fileID = fopen('zeta.txt','r');
zeta=fscanf(fileID,'%f,',[1, inf]);
zeta=read2(zeta, period_trans, nmkt);

%% read psi
fileID = fopen('psi.txt','r');
psi=fscanf(fileID,'%f,',[1, inf]);
psi=read2(psi, period_trans, nmkt);

%% read gridz
fileID = fopen('gridz.txt','r');
gridz=fscanf(fileID,'%f,',[1, inf]);
gridz=read1(gridz, nz);

%% read probz
fileID = fopen('probz.txt','r');
probz=fscanf(fileID,'%f,',[1, inf]);
probz=read1(probz, nz);

%% read trans_z
fileID = fopen('trans_z.txt','r');
trans_z=fscanf(fileID,'%f,',[1, inf]);
trans_z=read2(trans_z, nz, nz);

%% read gridb
fileID = fopen('gridb.txt','r');
gridb=fscanf(fileID,'%f,',[1, inf]);
gridb=read1(gridb, nb);

%% read gridtheta
fileID = fopen('gridtheta.txt','r');
gridtheta=fscanf(fileID,'%f,',[1, inf]);
gridtheta=read1(gridtheta, ntheta);

%% read gridb_pdf
fileID = fopen('gridb_pdf.txt','r');
gridb_pdf=fscanf(fileID,'%f,',[1, inf]);
gridb_pdf=read1(gridb_pdf, nb_pdf);

%% read gridtheta_pdf
fileID = fopen('gridtheta_pdf.txt','r');
gridtheta_pdf=fscanf(fileID,'%f,',[1, inf]);
gridtheta_pdf=read1(gridtheta_pdf, ntheta_pdf);

%% read income_entre
fileID = fopen('income_entre.txt','r');
income_entre=fscanf(fileID,'%f,',[1, inf]);
income_entre=read3(income_entre, period_trans-1, nmkt, nincome_pdf);

%% read income_worker
fileID = fopen('income_worker.txt','r');
income_worker=fscanf(fileID,'%f,',[1, inf]);
income_worker=read3(income_worker, period_trans-1, nmkt, nincome_pdf);

%% read grid_income
fileID = fopen('grid_income.txt','r');
grid_income=fscanf(fileID,'%f,',[1, inf]);

%% read frac_entre
fileID = fopen('frac_entre.txt','r');
frac_entre=fscanf(fileID,'%f,',[1, inf]);
frac_entre=read2(frac_entre, period_trans-1, nmkt);

%% read frac_credit
fileID = fopen('frac_credit.txt','r');
frac_credit=fscanf(fileID,'%f,',[1, inf]);
frac_credit=read2(frac_credit, period_trans-1, nmkt);

%% read labor_demand_stat
fileID = fopen('labor_demand_stat.txt','r');
labor_demand_stat=fscanf(fileID,'%f,',[1, inf]);
labor_demand_stat=read2(labor_demand_stat, period_trans-1, nmkt);

save('data.mat');
fclose('all');