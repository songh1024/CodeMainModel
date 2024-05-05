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
d=read1(d, nmkt);

%% read zeta
fileID = fopen('zeta.txt','r');
zeta=fscanf(fileID,'%f,',[1, inf]);
zeta=read1(zeta, nmkt);

%% read psi
fileID = fopen('psi.txt','r');
psi=fscanf(fileID,'%f,',[1, inf]);
psi=read1(psi, nmkt);

%% read gridz
fileID = fopen('gridz.txt','r');
gridz=fscanf(fileID,'%f,',[1, inf]);
gridz=read1(gridz, nz);

%% read probz
fileID = fopen('probz.txt','r');
probz=fscanf(fileID,'%f,',[1, inf]);
probz=read1(probz, zn);

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

%% read sum_pdf_B
fileID = fopen('sum_pdf_B.txt','r');
sum_pdf_B=fscanf(fileID,'%f,',[1, inf]);
sum_pdf_B=read1(sum_pdf_B, nmkt);

%% read error_pdf_B
fileID = fopen('error_pdf_B.txt','r');
error_pdf_B=fscanf(fileID,'%f,',[1, inf]);
error_pdf_B=read1(error_pdf_B, nmkt);

%% read w_iterate
fileID = fopen('w_iterate.txt','r');
w_iterate=fscanf(fileID,'%f,',[1, inf]);
w_iterate=read3(w_iterate, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);

%% read r_iterate
fileID = fopen('r_iterate.txt','r');
r_iterate=fscanf(fileID,'%f,',[1, inf]);
r_iterate=read1(r_iterate, maxit_bisec_ss_r);

%% read labor_demand_iterate
fileID = fopen('labor_demand_iterate.txt','r');
labor_demand_iterate=fscanf(fileID,'%f,',[1, inf]);
labor_demand_iterate=read3(labor_demand_iterate, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);

%% read labor_supply_iterate
fileID = fopen('labor_supply_iterate.txt','r');
labor_supply_iterate=fscanf(fileID,'%f,',[1, inf]);
labor_supply_iterate=read3(labor_supply_iterate, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);

%% read capital_demand_iterate
fileID = fopen('capital_demand_iterate.txt','r');
capital_demand_iterate=fscanf(fileID,'%f,',[1, inf]);
capital_demand_iterate=read3(capital_demand_iterate, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);

%% read capital_supply_iterate
fileID = fopen('capital_supply_iterate.txt','r');
capital_supply_iterate=fscanf(fileID,'%f,',[1, inf]);
capital_supply_iterate=read3(capital_supply_iterate, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);

%% read GDP_output
fileID = fopen('GDP_output.txt','r');
GDP_output=fscanf(fileID,'%f,',[1, inf]);
GDP_output=read3(GDP_output, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);

%% read pdf_A
fileID = fopen('pdf_A.txt','r');
pdf_A=fscanf(fileID,'%f,',[1, inf]);
pdf_A=read4(pdf_A, nmkt, nz, nb_pdf, ntheta_pdf);

%% read pdf_B
fileID = fopen('pdf_B.txt','r');
pdf_B=fscanf(fileID,'%f,',[1, inf]);
pdf_B=read4(pdf_B, nmkt, nz, nb_pdf, ntheta_pdf);

%% read mig
fileID = fopen('mig.txt','r');
mig=fscanf(fileID,'%f,',[1, inf]);
mig=read5(mig, period, nmkt, nz, nb_pdf, ntheta_pdf);

%% read V
fileID = fopen('V.txt','r');
V=fscanf(fileID,'%f,',[1, inf]);
V=read5(V, maxit_DP, nmkt, nz, nb, ntheta);

%% read gc
fileID = fopen('gc.txt','r');
gc=fscanf(fileID,'%f,',[1, inf]);
gc=read4(gc, nmkt, nz, nb, ntheta);

%% read gtheta
fileID = fopen('gtheta.txt','r');
gtheta=fscanf(fileID,'%f,',[1, inf]);
gtheta=read4(gtheta, nmkt, nz, nb, ntheta);

%% read gb
fileID = fopen('gb.txt','r');
gb=fscanf(fileID,'%f,',[1, inf]);
gb=read4(gb, nmkt, nz, nb, ntheta);

%% read adjust
fileID = fopen('adjust.txt','r');
adjust=fscanf(fileID,'%f,',[1, inf]);
adjust=read4(adjust, nmkt, nz, nb, ntheta);

%% read wealth
fileID = fopen('wealth.txt','r');
wealth=fscanf(fileID,'%f,',[1, inf]);
wealth=read4(wealth, nmkt, nz, nb, ntheta);

%% read gk
fileID = fopen('gk.txt','r');
gk=fscanf(fileID,'%f,',[1, inf]);
gk=read4(gk, nmkt, nz, nb, ntheta);

%% read gl
fileID = fopen('gl.txt','r');
gl=fscanf(fileID,'%f,',[1, inf]);
gl=read4(gl, nmkt, nz, nb, ntheta);


save('data.mat');
fclose('all');

