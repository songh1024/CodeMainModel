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

%% read sum_pdf_B_ini
fileID = fopen('sum_pdf_B_ini.txt','r');
sum_pdf_B_ini=fscanf(fileID,'%f,',[1, inf]);
sum_pdf_B_ini=read1(sum_pdf_B_ini, nmkt);

%% read error_pdf_B_ini
fileID = fopen('error_pdf_B_ini.txt','r');
error_pdf_B_ini=fscanf(fileID,'%f,',[1, inf]);
error_pdf_B_ini=read1(error_pdf_B_ini, nmkt);

%% read w_iterate_ini
fileID = fopen('w_iterate_ini.txt','r');
w_iterate_ini=fscanf(fileID,'%f,',[1, inf]);
w_iterate_ini=read3(w_iterate_ini, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);

%% read r_iterate_ini
fileID = fopen('r_iterate_ini.txt','r');
r_iterate_ini=fscanf(fileID,'%f,',[1, inf]);
r_iterate_ini=read1(r_iterate_ini, maxit_bisec_ss_r);

%% read labor_demand_iterate_ini
fileID = fopen('labor_demand_iterate_ini.txt','r');
labor_demand_iterate_ini=fscanf(fileID,'%f,',[1, inf]);
labor_demand_iterate_ini=read3(labor_demand_iterate_ini, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);

%% read labor_supply_iterate_ini
fileID = fopen('labor_supply_iterate_ini.txt','r');
labor_supply_iterate_ini=fscanf(fileID,'%f,',[1, inf]);
labor_supply_iterate_ini=read3(labor_supply_iterate_ini, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);

%% read capital_demand_iterate_ini
fileID = fopen('capital_demand_iterate_ini.txt','r');
capital_demand_iterate_ini=fscanf(fileID,'%f,',[1, inf]);
capital_demand_iterate_ini=read3(capital_demand_iterate_ini, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);

%% read capital_supply_iterate_ini
fileID = fopen('capital_supply_iterate_ini.txt','r');
capital_supply_iterate_ini=fscanf(fileID,'%f,',[1, inf]);
capital_supply_iterate_ini=read3(capital_supply_iterate_ini, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);

%% read GDP_output_ini
fileID = fopen('GDP_output_ini.txt','r');
GDP_output_ini=fscanf(fileID,'%f,',[1, inf]);
GDP_output_ini=read3(GDP_output_ini, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);

%% read pdf_A_ini
fileID = fopen('pdf_A_ini.txt','r');
pdf_A_ini=fscanf(fileID,'%f,',[1, inf]);
pdf_A_ini=read4(pdf_A_ini, nmkt, nz, nb_pdf, ntheta_pdf);

%% read pdf_B_ini
fileID = fopen('pdf_B_ini.txt','r');
pdf_B_ini=fscanf(fileID,'%f,',[1, inf]);
pdf_B_ini=read4(pdf_B_ini, nmkt, nz, nb_pdf, ntheta_pdf);

%% read mig_ini
fileID = fopen('mig_ini.txt','r');
mig_ini=fscanf(fileID,'%f,',[1, inf]);
mig_ini=read5(mig_ini, period, nmkt, nz, nb_pdf, ntheta_pdf);

%% read V_noad_ini
fileID = fopen('V_noad_ini.txt','r');
V_noad_ini=fscanf(fileID,'%f,',[1, inf]);
V_noad_ini=read4(V_noad_ini, nmkt, nz, nb, ntheta);

%% read V_ad_ini
fileID = fopen('V_ad_ini.txt','r');
V_ad_ini=fscanf(fileID,'%f,',[1, inf]);
V_ad_ini=read4(V_ad_ini, nmkt, nz, nb, ntheta);

%% read V_ini
fileID = fopen('V_ini.txt','r');
V_ini=fscanf(fileID,'%f,',[1, inf]);
V_ini=read5(V_ini, maxit_DP, nmkt, nz, nb, ntheta);

%% read gc_ini
fileID = fopen('gc_ini.txt','r');
gc_ini=fscanf(fileID,'%f,',[1, inf]);
gc_ini=read4(gc_ini, nmkt, nz, nb, ntheta);

%% read gtheta_ini
fileID = fopen('gtheta_ini.txt','r');
gtheta_ini=fscanf(fileID,'%f,',[1, inf]);
gtheta_ini=read4(gtheta_ini, nmkt, nz, nb, ntheta);

%% read gb_ini
fileID = fopen('gb_ini.txt','r');
gb_ini=fscanf(fileID,'%f,',[1, inf]);
gb_ini=read4(gb_ini, nmkt, nz, nb, ntheta);

%% read adjust_noad_ini
fileID = fopen('adjust_noad_ini.txt','r');
adjust_noad_ini=fscanf(fileID,'%f,',[1, inf]);
adjust_noad_ini=read4(adjust_noad_ini, nmkt, nz, nb, ntheta);

%% read adjust_ini
fileID = fopen('adjust_ini.txt','r');
adjust_ini=fscanf(fileID,'%f,',[1, inf]);
adjust_ini=read4(adjust_ini, nmkt, nz, nb, ntheta);

%% read wealth_ini
fileID = fopen('wealth_ini.txt','r');
wealth_ini=fscanf(fileID,'%f,',[1, inf]);
wealth_ini=read4(wealth_ini, nmkt, nz, nb, ntheta);

%% read gk_ini
fileID = fopen('gk_ini.txt','r');
gk_ini=fscanf(fileID,'%f,',[1, inf]);
gk_ini=read4(gk_ini, nmkt, nz, nb, ntheta);

%% read gl_ini
fileID = fopen('gl_ini.txt','r');
gl_ini=fscanf(fileID,'%f,',[1, inf]);
gl_ini=read4(gl_ini, nmkt, nz, nb, ntheta);


save('data.mat');
fclose('all');







%% read sum_pdf_B_new
fileID = fopen('sum_pdf_B_new.txt','r');
sum_pdf_B_new=fscanf(fileID,'%f,',[1, inf]);
sum_pdf_B_new=read1(sum_pdf_B_new, nmkt);

%% read error_pdf_B_new
fileID = fopen('error_pdf_B_new.txt','r');
error_pdf_B_new=fscanf(fileID,'%f,',[1, inf]);
error_pdf_B_new=read1(error_pdf_B_new, nmkt);

%% read w_iterate_new
fileID = fopen('w_iterate_new.txt','r');
w_iterate_new=fscanf(fileID,'%f,',[1, inf]);
w_iterate_new=read3(w_iterate_new, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);

%% read r_iterate_new
fileID = fopen('r_iterate_new.txt','r');
r_iterate_new=fscanf(fileID,'%f,',[1, inf]);
r_iterate_new=read1(r_iterate_new, maxit_bisec_ss_r);

%% read labor_demand_iterate_new
fileID = fopen('labor_demand_iterate_new.txt','r');
labor_demand_iterate_new=fscanf(fileID,'%f,',[1, inf]);
labor_demand_iterate_new=read3(labor_demand_iterate_new, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);

%% read labor_supply_iterate_new
fileID = fopen('labor_supply_iterate_new.txt','r');
labor_supply_iterate_new=fscanf(fileID,'%f,',[1, inf]);
labor_supply_iterate_new=read3(labor_supply_iterate_new, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);

%% read capital_demand_iterate_new
fileID = fopen('capital_demand_iterate_new.txt','r');
capital_demand_iterate_new=fscanf(fileID,'%f,',[1, inf]);
capital_demand_iterate_new=read3(capital_demand_iterate_new, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);

%% read capital_supply_iterate_new
fileID = fopen('capital_supply_iterate_new.txt','r');
capital_supply_iterate_new=fscanf(fileID,'%f,',[1, inf]);
capital_supply_iterate_new=read3(capital_supply_iterate_new, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);

%% read GDP_output_new
fileID = fopen('GDP_output_new.txt','r');
GDP_output_new=fscanf(fileID,'%f,',[1, inf]);
GDP_output_new=read3(GDP_output_new, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);

%% read pdf_A_new
fileID = fopen('pdf_A_new.txt','r');
pdf_A_new=fscanf(fileID,'%f,',[1, inf]);
pdf_A_new=read4(pdf_A_new, nmkt, nz, nb_pdf, ntheta_pdf);

%% read pdf_B_new
fileID = fopen('pdf_B_new.txt','r');
pdf_B_new=fscanf(fileID,'%f,',[1, inf]);
pdf_B_new=read4(pdf_B_new, nmkt, nz, nb_pdf, ntheta_pdf);

%% read mig_new
fileID = fopen('mig_new.txt','r');
mig_new=fscanf(fileID,'%f,',[1, inf]);
mig_new=read5(mig_new, period, nmkt, nz, nb_pdf, ntheta_pdf);

%% read V_new
fileID = fopen('V_new.txt','r');
V_new=fscanf(fileID,'%f,',[1, inf]);
V_new=read5(V_new, maxit_DP, nmkt, nz, nb, ntheta);

%% read gc_new
fileID = fopen('gc_new.txt','r');
gc_new=fscanf(fileID,'%f,',[1, inf]);
gc_new=read4(gc_new, nmkt, nz, nb, ntheta);

%% read gtheta_new
fileID = fopen('gtheta_new.txt','r');
gtheta_new=fscanf(fileID,'%f,',[1, inf]);
gtheta_new=read4(gtheta_new, nmkt, nz, nb, ntheta);

%% read gb_new
fileID = fopen('gb_new.txt','r');
gb_new=fscanf(fileID,'%f,',[1, inf]);
gb_new=read4(gb_new, nmkt, nz, nb, ntheta);

%% read adjust_new
fileID = fopen('adjust_new.txt','r');
adjust_new=fscanf(fileID,'%f,',[1, inf]);
adjust_new=read4(adjust_new, nmkt, nz, nb, ntheta);

%% read wealth_new
fileID = fopen('wealth_new.txt','r');
wealth_new=fscanf(fileID,'%f,',[1, inf]);
wealth_new=read4(wealth_new, nmkt, nz, nb, ntheta);

%% read gk_new
fileID = fopen('gk_new.txt','r');
gk_new=fscanf(fileID,'%f,',[1, inf]);
gk_new=read4(gk_new, nmkt, nz, nb, ntheta);

%% read gl_new
fileID = fopen('gl_new.txt','r');
gl_new=fscanf(fileID,'%f,',[1, inf]);
gl_new=read4(gl_new, nmkt, nz, nb, ntheta);


save('data.mat');
fclose('all');




%% read w_guess
fileID = fopen('w_guess.txt','r');
w_guess=fscanf(fileID,'%f,',[1, inf]);
w_guess=read4(w_guess, maxit_trans_r, maxit_trans_w, period_trans, nmkt);

%% read w_implied
fileID = fopen('w_implied.txt','r');
w_implied=fscanf(fileID,'%f,',[1, inf]);
w_implied=read4(w_implied, maxit_trans_r, maxit_trans_w, period_trans, nmkt);

%% read labor_demand
fileID = fopen('labor_demand.txt','r');
labor_demand=fscanf(fileID,'%f,',[1, inf]);
labor_demand=read4(labor_demand, maxit_trans_r, maxit_trans_w, period_trans, nmkt);

%% read labor_supply
fileID = fopen('labor_supply.txt','r');
labor_supply=fscanf(fileID,'%f,',[1, inf]);
labor_supply=read4(labor_supply, maxit_trans_r, maxit_trans_w, period_trans, nmkt);

%% read r_guess
fileID = fopen('r_guess.txt','r');
r_guess=fscanf(fileID,'%f,',[1, inf]);
r_guess=read2(r_guess, maxit_trans_r, period_trans);

%% read r_implied
fileID = fopen('r_implied.txt','r');
r_implied=fscanf(fileID,'%f,',[1, inf]);
r_implied=read2(r_implied, maxit_trans_r, period_trans);

%% read capital_demand
fileID = fopen('capital_demand.txt','r');
capital_demand=fscanf(fileID,'%f,',[1, inf]);
capital_demand=read4(capital_demand, maxit_trans_r, maxit_trans_w, period_trans, nmkt);

%% read capital_supply
fileID = fopen('capital_supply.txt','r');
capital_supply=fscanf(fileID,'%f,',[1, inf]);
capital_supply=read4(capital_supply, maxit_trans_r, maxit_trans_w, period_trans, nmkt);

%% read GDP_output_trans
fileID = fopen('GDP_output_trans.txt','r');
GDP_output_trans=fscanf(fileID,'%f,',[1, inf]);
GDP_output_trans=read4(GDP_output_trans, maxit_trans_r, maxit_trans_w, period_trans, nmkt);

%% read pdf_A_forward
fileID = fopen('pdf_A_forward.txt','r');
pdf_A_forward=fscanf(fileID,'%f,',[1, inf]);
pdf_A_forward=read5(pdf_A_forward, period_trans, nmkt, nz, nb_pdf, ntheta_pdf);

%% read pdf_B_forward
fileID = fopen('pdf_B_forward.txt','r');
pdf_B_forward=fscanf(fileID,'%f,',[1, inf]);
pdf_B_forward=read5(pdf_B_forward, period_trans, nmkt, nz, nb_pdf, ntheta_pdf);

%% read V_backward
fileID = fopen('V_backward.txt','r');
V_backward=fscanf(fileID,'%f,',[1, inf]);
V_backward=read5(V_backward, period_trans, nmkt, nz, nb, ntheta);

%% read gc_backward
fileID = fopen('gc_backward.txt','r');
gc_backward=fscanf(fileID,'%f,',[1, inf]);
gc_backward=read5(gc_backward, period_trans, nmkt, nz, nb, ntheta);

%% read gtheta_backward
fileID = fopen('gtheta_backward.txt','r');
gtheta_backward=fscanf(fileID,'%f,',[1, inf]);
gtheta_backward=read5(gtheta_backward, period_trans, nmkt, nz, nb, ntheta);

%% read gb_backward
fileID = fopen('gb_backward.txt','r');
gb_backward=fscanf(fileID,'%f,',[1, inf]);
gb_backward=read5(gb_backward, period_trans, nmkt, nz, nb, ntheta);

%% read adjust_backward
fileID = fopen('adjust_backward.txt','r');
adjust_backward=fscanf(fileID,'%f,',[1, inf]);
adjust_backward=read5(adjust_backward, period_trans, nmkt, nz, nb, ntheta);

%% read wealth_backward
fileID = fopen('wealth_backward.txt','r');
wealth_backward=fscanf(fileID,'%f,',[1, inf]);
wealth_backward=read5(wealth_backward, period_trans, nmkt, nz, nb, ntheta);

%% read gk_backward
fileID = fopen('gk_backward.txt','r');
gk_backward=fscanf(fileID,'%f,',[1, inf]);
gk_backward=read5(gk_backward, period_trans, nmkt, nz, nb, ntheta);

%% read gl_backward
fileID = fopen('gl_backward.txt','r');
gl_backward=fscanf(fileID,'%f,',[1, inf]);
gl_backward=read5(gl_backward, period_trans, nmkt, nz, nb, ntheta);

%% read mig_trans
fileID = fopen('mig_trans.txt','r');
mig_trans=fscanf(fileID,'%f,',[1, inf]);
mig_trans=read6(mig_trans, period_trans, nmkt, nz, nb_pdf, ntheta_pdf,nmkt);


save('data.mat');
fclose('all');


%% read income_entre
fileID = fopen('income_entre.txt','r');
income_entre=fscanf(fileID,'%f,',[1, inf]);
income_entre=read3(income_entre, period_trans-1, nmkt, nincome_pdf);

%% read income_worker
fileID = fopen('income_worker.txt','r');
income_worker=fscanf(fileID,'%f,',[1, inf]);
income_worker=read3(income_worker, period_trans-1, nmkt, nincome_pdf);

%% read noninterest_income_entre
fileID = fopen('noninterest_income_entre.txt','r');
noninterest_income_entre=fscanf(fileID,'%f,',[1, inf]);
noninterest_income_entre=read3(noninterest_income_entre, period_trans-1, nmkt, nincome_pdf);

%% read noninterest_income_worker
fileID = fopen('noninterest_income_worker.txt','r');
noninterest_income_worker=fscanf(fileID,'%f,',[1, inf]);
noninterest_income_worker=read3(noninterest_income_worker, period_trans-1, nmkt, nincome_pdf);

%% read wealth_entre
fileID = fopen('wealth_entre.txt','r');
wealth_entre=fscanf(fileID,'%f,',[1, inf]);
wealth_entre=read3(wealth_entre, period_trans-1, nmkt, nb_pdf);

%% read wealth_worker
fileID = fopen('wealth_worker.txt','r');
wealth_worker=fscanf(fileID,'%f,',[1, inf]);
wealth_worker=read3(wealth_worker, period_trans-1, nmkt, nb_pdf);

%% read grid_labor
fileID = fopen('grid_labor.txt','r');
grid_labor=fscanf(fileID,'%f,',[1, inf]);

%% read labor_dist
fileID = fopen('labor_dist.txt','r');
labor_dist=fscanf(fileID,'%f,',[1, inf]);

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

%% read L
fileID = fopen('L.txt','r');
L=fscanf(fileID,'%f,',[1, inf]);
L=read2(L, period_trans-1, nmkt);

%% read K
fileID = fopen('K.txt','r');
K=fscanf(fileID,'%f,',[1, inf]);
K=read2(K, period_trans-1, nmkt);

%% read Y
fileID = fopen('Y.txt','r');
Y=fscanf(fileID,'%f,',[1, inf]);
Y=read2(Y, period_trans-1, nmkt);

%% read TFP
fileID = fopen('TFP.txt','r');
TFP=fscanf(fileID,'%f,',[1, inf]);
TFP=read2(TFP, period_trans-1, nmkt);

%% read GDP
fileID = fopen('GDP.txt','r');
GDP=fscanf(fileID,'%f,',[1, inf]);
GDP=read2(GDP, period_trans-1, nmkt);

%% read capital_demand_stat
fileID = fopen('capital_demand_stat.txt','r');
capital_demand_stat=fscanf(fileID,'%f,',[1, inf]);
capital_demand_stat=read2(capital_demand_stat, period_trans-1, nmkt);

%% read capital_supply_stat
fileID = fopen('capital_supply_stat.txt','r');
capital_supply_stat=fscanf(fileID,'%f,',[1, inf]);
capital_supply_stat=read2(capital_supply_stat, period_trans-1, nmkt);

%% read labor_demand_stat
fileID = fopen('labor_demand_stat.txt','r');
labor_demand_stat=fscanf(fileID,'%f,',[1, inf]);
labor_demand_stat=read2(labor_demand_stat, period_trans-1, nmkt);

%% read labor_supply_stat
fileID = fopen('labor_supply_stat.txt','r');
labor_supply_stat=fscanf(fileID,'%f,',[1, inf]);
labor_supply_stat=read2(labor_supply_stat, period_trans-1, nmkt);

%% read cash
fileID = fopen('cash.txt','r');
cash=fscanf(fileID,'%f,',[1, inf]);
cash=read2(cash, period_trans-1, nmkt);

%% read deposit
fileID = fopen('deposit.txt','r');
deposit=fscanf(fileID,'%f,',[1, inf]);
deposit=read2(deposit, period_trans-1, nmkt);

%% read credit
fileID = fopen('credit.txt','r');
credit=fscanf(fileID,'%f,',[1, inf]);
credit=read2(credit, period_trans-1, nmkt);


save('data.mat');
fclose('all');