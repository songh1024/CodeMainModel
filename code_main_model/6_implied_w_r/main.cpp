#include "parameters.h"
#include <stdio.h>  /* defines FILENAME_MAX */
#include <direct.h>  //windows
//#include <unistd.h>  //linux
//#include <chrono> 

#define GetCurrentDir _getcwd
//using namespace std::chrono; 
using namespace std;

double kappa;
double eta;

int main()
{
	char cCurrentPath[FILENAME_MAX];
	GetCurrentDir(cCurrentPath, sizeof(cCurrentPath));
	std::string currentpath = std::string(cCurrentPath);
	std::string inputpath = std::string(cCurrentPath);
	currentpath = currentpath.append("/computation_results/"); //linux
	inputpath = inputpath.append("/input/"); //linux
	std::string filename;

	// Initialize some commonly used variables in iteration, simulation and GDP accounting
	// iteration
	double* theta_candidate = ini_matrix1(coefficient);
	linspace(theta_candidate, mintheta, maxtheta, coefficient);

	// simulation
	int*** b_AtoB_M_index = ini_int_matrix3(nz, nb_pdf, ntheta_pdf);
	int*** theta_AtoB_M_index = ini_int_matrix3(nz, nb_pdf, ntheta_pdf);
	double***** V_AtoB = ini_matrix5(nmkt, nz, nb_pdf, ntheta_pdf, nmkt);
	double***** step_AtoB_prob = ini_matrix5(nmkt, nz, nb_pdf, ntheta_pdf, nmkt);
	int**** b_BtoA_next_index = ini_int_matrix4(nmkt, nz, nb_pdf, ntheta_pdf);
	int**** theta_BtoA_next_index = ini_int_matrix4(nmkt, nz, nb_pdf, ntheta_pdf);
	double** diff = ini_matrix2(maxit_DP, nmkt);

	// demand_supply
	double*** temp_fin = ini_matrix3(nz, nb, ntheta);
	double*** temp_cap_supply = ini_matrix3(nz, nb, ntheta);
	double*** temp_adjust = ini_matrix3(nz, nb, ntheta);
	double*** l = ini_matrix3(nz, nb, ntheta);

	// prices
	double* w = ini_matrix1(nmkt);
	double r = 0.0;


	// Initialize Markets
	double* Pi = ini_matrix1(nmkt);	
	filename=inputpath+("Pi.dat");	
	read_dat1(Pi, nmkt, filename);
	double Pi_sum = sum(nmkt, Pi);
	for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
		Pi[i_mkt] = Pi[i_mkt] / Pi_sum;

	double* Z = ini_matrix1(nmkt);
	filename = inputpath + ("Z.dat");
	read_dat1(Z, nmkt, filename);

	double** d = ini_matrix2(period_trans, nmkt);
	filename = inputpath + ("d.dat");
	read_dat2(d, nyears, nmkt, filename);
	for (int t = nyears; t < period_trans; ++t)
		for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
			d[t][i_mkt] = d[nyears - 1][i_mkt];

	double** zeta=ini_matrix2(period_trans, nmkt);
	double** psi= ini_matrix2(period_trans, nmkt);
	for (int t = 0; t < period_trans; ++t)
		for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
		{
			psi[t][i_mkt] = conspsi + d[t][i_mkt] * thetapsi;
			zeta[t][i_mkt] = conszeta + d[t][i_mkt] * thetazeta;
		}

	filename = currentpath + ("parameters.txt");
	write_parameters(filename);
	filename = currentpath + ("parameters_matlab.txt");
	write_parameters_matlab(filename);
	filename = currentpath + ("Pi.txt");
	write1(filename, Pi, nmkt);
	filename = currentpath + ("Z.txt");
	write1(filename, Z, nmkt);
	filename = currentpath + ("d.txt");
	write2(filename, d, period_trans, nmkt);
	filename=currentpath+("zeta.txt");
	write2(filename, zeta, period_trans, nmkt);
	filename=currentpath+("psi.txt");
	write2(filename, psi, period_trans, nmkt);

	// Initialization for Value Function Iteration
	double* gridz = ini_matrix1(nz); // grids of individual talent
	double* probz = ini_matrix1(nz);
	double** trans_z = ini_matrix2(nz,nz);
	z_generator(gridz, probz, trans_z);
	//for (int i = 0; i < nz; ++i)
		//cout << i << "    " << gridz[i] << "    " << probz[i] << endl;

	double* gridb=ini_matrix1(nb);
	linspace(gridb,minb,maxb,nb);
	double* gridtheta=ini_matrix1(ntheta);
	linspace(gridtheta, mintheta, maxtheta, ntheta);

	// Initialization for Simulation
	double* gridb_pdf = ini_matrix1(nb_pdf);
	linspace(gridb_pdf, minb, maxb, nb_pdf);
	double* gridtheta_pdf = ini_matrix1(ntheta_pdf);
	linspace(gridtheta_pdf, mintheta, maxtheta, ntheta_pdf);

	filename = currentpath + ("gridz.txt");
	write1(filename, gridz, nz);
	filename = currentpath + ("probz.txt");
	write1(filename, probz, nz);
	filename = currentpath + ("trans_z.txt");
	write2(filename, trans_z, nz, nz);
	filename = currentpath + ("gridb.txt");
	write1(filename, gridb, nb);
	filename = currentpath + ("gridtheta.txt");
	write1(filename, gridtheta, ntheta);

	filename = currentpath + ("gridb_pdf.txt");
	write1(filename, gridb_pdf, nb_pdf);
	filename = currentpath + ("gridtheta_pdf.txt");
	write1(filename, gridtheta_pdf, ntheta_pdf);

	// Initialization for Initial Steady-State
	kappa = kappa_ini; // When solving the initial steady state, we don't allow migration
	eta = eta_ini;

	double**** V_ad_ini=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gc_ad_ini=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gtheta_ad_ini=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gb_ad_ini=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gk_ad_ini=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gl_ad_ini=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** adjust_ad_ini=ini_matrix4(nmkt, nz, nb, ntheta); // 1 worker, 2 noborrow entrepreneur, 4 borrow entrepreneur
	double**** wealth_ad_ini=ini_matrix4(nmkt, nz, nb, ntheta);

	double**** V_noad_ini=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gc_noad_ini=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gtheta_noad_ini=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gb_noad_ini=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gk_noad_ini=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gl_noad_ini=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** adjust_noad_ini=ini_matrix4(nmkt, nz, nb, ntheta); // 1 worker, 2 noborrow entrepreneur, 4 borrow entrepreneur
	double**** wealth_noad_ini=ini_matrix4(nmkt, nz, nb, ntheta);

	double***** V_ini=ini_matrix5(maxit_DP, nmkt, nz, nb, ntheta);
	double**** gc_ini=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gtheta_ini=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gb_ini=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gk_ini=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gl_ini=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** adjust_ini=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** wealth_ini=ini_matrix4(nmkt, nz, nb, ntheta);

	double***** mig_ini= ini_matrix5(period, nmkt, nz, nb_pdf, ntheta_pdf);
	double***** pdf_A_ini=ini_matrix5(period, nmkt, nz, nb_pdf, ntheta_pdf);
	double***** pdf_B_ini=ini_matrix5(period, nmkt, nz, nb_pdf, ntheta_pdf);
	double** sum_pdf_B_ini=ini_matrix2(period, nmkt);
	double** error_pdf_B_ini=ini_matrix2(period, nmkt);
	
	double*** w_iterate_ini=ini_matrix3(maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	double* r_iterate_ini=ini_matrix1(maxit_bisec_ss_r);
	double*** labor_demand_iterate_ini=ini_matrix3(maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	double*** labor_supply_iterate_ini=ini_matrix3(maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	double*** capital_demand_iterate_ini=ini_matrix3(maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	double*** capital_supply_iterate_ini=ini_matrix3(maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	double*** GDP_output_ini=ini_matrix3(maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);

	GE(Pi, Z, psi[0], zeta[0], w, &r,

			 V_ad_ini, gc_ad_ini, gtheta_ad_ini, gb_ad_ini, gk_ad_ini, gl_ad_ini, adjust_ad_ini, wealth_ad_ini, 
			 V_noad_ini, gc_noad_ini, gtheta_noad_ini, gb_noad_ini, gk_noad_ini, gl_noad_ini, adjust_noad_ini, wealth_noad_ini, 
			 V_ini, gc_ini, gtheta_ini, gb_ini, gk_ini, gl_ini, adjust_ini, wealth_ini, 
			 gridb, gridtheta, gridz, probz, trans_z,
			 gridb_pdf, gridtheta_pdf,

			 pdf_A_ini, pdf_B_ini,  
			 w_iterate_ini, r_iterate_ini,
			 labor_demand_iterate_ini, labor_supply_iterate_ini, capital_demand_iterate_ini, capital_supply_iterate_ini, GDP_output_ini,
			 mig_ini,

			 sum_pdf_B_ini, error_pdf_B_ini, 
			 
			 diff, theta_candidate, b_AtoB_M_index, theta_AtoB_M_index, V_AtoB, step_AtoB_prob, b_BtoA_next_index, theta_BtoA_next_index, temp_fin, temp_cap_supply, l, currentpath);
	
	for (int i_mkt=0; i_mkt <nmkt; ++i_mkt)
		w_iterate_ini[maxit_bisec_ss_r-1][maxit_bisec_ss_w-1][i_mkt]=w[i_mkt];
	r_iterate_ini[maxit_bisec_ss_r-1]=r;


	filename=currentpath+("sum_pdf_B_ini.txt");
	write2(filename, sum_pdf_B_ini, period, nmkt);
	filename=currentpath+("error_pdf_B_ini.txt");
	write2(filename, error_pdf_B_ini, period, nmkt);

	filename=currentpath+("w_iterate_ini.txt");
	write3(filename, w_iterate_ini, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	filename=currentpath+("r_iterate_ini.txt");
	write1(filename, r_iterate_ini, maxit_bisec_ss_r);

	filename=currentpath+("labor_demand_iterate_ini.txt");
	write3(filename, labor_demand_iterate_ini, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	filename=currentpath+("labor_supply_iterate_ini.txt");
	write3(filename, labor_supply_iterate_ini, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	filename=currentpath+("capital_demand_iterate_ini.txt");
	write3(filename, capital_demand_iterate_ini, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	filename=currentpath+("capital_supply_iterate_ini.txt");
	write3(filename, capital_supply_iterate_ini, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	filename=currentpath+("GDP_output_ini.txt");
	write3(filename, GDP_output_ini, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);

	filename=currentpath+("pdf_A_ini.txt");
	write4(filename, pdf_A_ini[period-2], nmkt, nz, nb_pdf, ntheta_pdf);
	//write5(filename, pdf_A, period, nmkt, nz, nb_pdf, ntheta_pdf);
	filename=currentpath+("pdf_B_ini.txt");
	write4(filename, pdf_B_ini[period-2], nmkt, nz, nb_pdf, ntheta_pdf);
	//write5(filename, pdf_B, period, nmkt, nz, nb_pdf, ntheta_pdf);
	filename = currentpath + ("mig_ini.txt");
	write5(filename, mig_ini, period, nmkt, nz, nb_pdf, ntheta_pdf);

	filename = currentpath + ("V_ad_ini.txt");
	write4(filename, V_ad_ini, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gc_ad_ini.txt");
	write4(filename, gc_ad_ini, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gtheta_ad_ini.txt");
	write4(filename, gtheta_ad_ini, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gb_ad_ini.txt");
	write4(filename, gb_ad_ini, nmkt, nz, nb, ntheta);
	filename = currentpath + ("adjust_ad_ini.txt");
	write4(filename, adjust_ad_ini, nmkt, nz, nb, ntheta);
	filename = currentpath + ("wealth_ad_ini.txt");
	write4(filename, wealth_ad_ini, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gk_ad_ini.txt");
	write4(filename, gk_ad_ini, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gl_ad_ini.txt");
	write4(filename, gl_ad_ini, nmkt, nz, nb, ntheta);

	filename=currentpath+("V_noad_ini.txt");
	write4(filename, V_noad_ini, nmkt, nz, nb, ntheta);
	filename=currentpath+("gc_noad_ini.txt");
	write4(filename, gc_noad_ini, nmkt, nz, nb, ntheta);
	filename=currentpath+("gtheta_noad_ini.txt");
	write4(filename, gtheta_noad_ini, nmkt, nz, nb, ntheta);
	filename=currentpath+("gb_noad_ini.txt");
	write4(filename, gb_noad_ini, nmkt, nz, nb, ntheta);
	filename=currentpath+("adjust_noad_ini.txt");
	write4(filename, adjust_noad_ini, nmkt, nz, nb, ntheta);
	filename=currentpath+("wealth_noad_ini.txt");
	write4(filename, wealth_noad_ini, nmkt, nz, nb, ntheta);
	filename=currentpath+("gk_noad_ini.txt");
	write4(filename, gk_noad_ini, nmkt, nz, nb, ntheta);
	filename=currentpath+("gl_noad_ini.txt");
	write4(filename, gl_noad_ini, nmkt, nz, nb, ntheta);

	filename=currentpath+("V_ini.txt");
	write5(filename, V_ini, maxit_DP, nmkt, nz, nb, ntheta);
	filename=currentpath+("gc_ini.txt");
	write4(filename, gc_ini, nmkt, nz, nb, ntheta);
	filename=currentpath+("gtheta_ini.txt");
	write4(filename, gtheta_ini, nmkt, nz, nb, ntheta);
	filename=currentpath+("gb_ini.txt");
	write4(filename, gb_ini, nmkt, nz, nb, ntheta);
	filename=currentpath+("adjust_ini.txt");
	write4(filename, adjust_ini, nmkt, nz, nb, ntheta);
	filename=currentpath+("wealth_ini.txt");
	write4(filename, wealth_ini, nmkt, nz, nb, ntheta);
	filename=currentpath+("gk_ini.txt");
	write4(filename, gk_ini, nmkt, nz, nb, ntheta);
	filename=currentpath+("gl_ini.txt");
	write4(filename, gl_ini, nmkt, nz, nb, ntheta);









	// Initialization for New Steady-State
	kappa = kappa_new;
	eta = eta_new;

	double**** V_ad_new = ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gc_ad_new = ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gtheta_ad_new = ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gb_ad_new = ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gk_ad_new = ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gl_ad_new = ini_matrix4(nmkt, nz, nb, ntheta);
	double**** adjust_ad_new = ini_matrix4(nmkt, nz, nb, ntheta); // 1 worker, 2 noborrow entrepreneur, 4 borrow entrepreneur
	double**** wealth_ad_new = ini_matrix4(nmkt, nz, nb, ntheta);

	double**** V_noad_new = ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gc_noad_new = ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gtheta_noad_new = ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gb_noad_new = ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gk_noad_new = ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gl_noad_new = ini_matrix4(nmkt, nz, nb, ntheta);
	double**** adjust_noad_new = ini_matrix4(nmkt, nz, nb, ntheta); // 1 worker, 2 noborrow entrepreneur, 4 borrow entrepreneur
	double**** wealth_noad_new = ini_matrix4(nmkt, nz, nb, ntheta);

	double***** V_new = ini_matrix5(maxit_DP, nmkt, nz, nb, ntheta);
	double**** gc_new = ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gtheta_new = ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gb_new = ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gk_new = ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gl_new = ini_matrix4(nmkt, nz, nb, ntheta);
	double**** adjust_new = ini_matrix4(nmkt, nz, nb, ntheta);
	double**** wealth_new = ini_matrix4(nmkt, nz, nb, ntheta);

	double***** mig_new = ini_matrix5(period, nmkt, nz, nb_pdf, ntheta_pdf);
	double***** pdf_A_new = ini_matrix5(period, nmkt, nz, nb_pdf, ntheta_pdf);
	double***** pdf_B_new = ini_matrix5(period, nmkt, nz, nb_pdf, ntheta_pdf);
	double** sum_pdf_B_new = ini_matrix2(period, nmkt);
	double** error_pdf_B_new = ini_matrix2(period, nmkt);

	double*** w_iterate_new = ini_matrix3(maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	double* r_iterate_new = ini_matrix1(maxit_bisec_ss_r);
	double*** labor_demand_iterate_new = ini_matrix3(maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	double*** labor_supply_iterate_new = ini_matrix3(maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	double*** capital_demand_iterate_new = ini_matrix3(maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	double*** capital_supply_iterate_new = ini_matrix3(maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	double*** GDP_output_new = ini_matrix3(maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);

	GE(Pi, Z, psi[period_trans-1], zeta[period_trans-1], w, &r,

		V_ad_new, gc_ad_new, gtheta_ad_new, gb_ad_new, gk_ad_new, gl_ad_new, adjust_ad_new, wealth_ad_new,
		V_noad_new, gc_noad_new, gtheta_noad_new, gb_noad_new, gk_noad_new, gl_noad_new, adjust_noad_new, wealth_noad_new,
		V_new, gc_new, gtheta_new, gb_new, gk_new, gl_new, adjust_new, wealth_new,
		gridb, gridtheta, gridz, probz, trans_z,
		gridb_pdf, gridtheta_pdf,

		pdf_A_new, pdf_B_new,
		w_iterate_new, r_iterate_new,
		labor_demand_iterate_new, labor_supply_iterate_new, capital_demand_iterate_new, capital_supply_iterate_new, GDP_output_new,
		mig_new,

		sum_pdf_B_new, error_pdf_B_new,

		diff, theta_candidate, b_AtoB_M_index, theta_AtoB_M_index, V_AtoB, step_AtoB_prob, b_BtoA_next_index, theta_BtoA_next_index, temp_fin, temp_cap_supply, l, currentpath);

	for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
		w_iterate_new[maxit_bisec_ss_r - 1][maxit_bisec_ss_w - 1][i_mkt] = w[i_mkt];
	r_iterate_new[maxit_bisec_ss_r - 1] = r;


	filename = currentpath + ("sum_pdf_B_new.txt");
	write2(filename, sum_pdf_B_new, period, nmkt);
	filename = currentpath + ("error_pdf_B_new.txt");
	write2(filename, error_pdf_B_new, period, nmkt);

	filename = currentpath + ("w_iterate_new.txt");
	write3(filename, w_iterate_new, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	filename = currentpath + ("r_iterate_new.txt");
	write1(filename, r_iterate_new, maxit_bisec_ss_r);

	filename = currentpath + ("labor_demand_iterate_new.txt");
	write3(filename, labor_demand_iterate_new, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	filename = currentpath + ("labor_supply_iterate_new.txt");
	write3(filename, labor_supply_iterate_new, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	filename = currentpath + ("capital_demand_iterate_new.txt");
	write3(filename, capital_demand_iterate_new, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	filename = currentpath + ("capital_supply_iterate_new.txt");
	write3(filename, capital_supply_iterate_new, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	filename = currentpath + ("GDP_output_new.txt");
	write3(filename, GDP_output_new, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);

	filename = currentpath + ("pdf_A_new.txt");
	write4(filename, pdf_A_new[period - 2], nmkt, nz, nb_pdf, ntheta_pdf);
	//write5(filename, pdf_A, period, nmkt, nz, nb_pdf, ntheta_pdf);
	filename = currentpath + ("pdf_B_new.txt");
	write4(filename, pdf_B_new[period - 2], nmkt, nz, nb_pdf, ntheta_pdf);
	//write5(filename, pdf_B, period, nmkt, nz, nb_pdf, ntheta_pdf);
	filename = currentpath + ("mig_new.txt");
	write5(filename, mig_new, period, nmkt, nz, nb_pdf, ntheta_pdf);

	filename = currentpath + ("V_ad_new.txt");
	write4(filename, V_ad_new, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gc_ad_new.txt");
	write4(filename, gc_ad_new, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gtheta_ad_new.txt");
	write4(filename, gtheta_ad_new, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gb_ad_new.txt");
	write4(filename, gb_ad_new, nmkt, nz, nb, ntheta);
	filename = currentpath + ("adjust_ad_new.txt");
	write4(filename, adjust_ad_new, nmkt, nz, nb, ntheta);
	filename = currentpath + ("wealth_ad_new.txt");
	write4(filename, wealth_ad_new, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gk_ad_new.txt");
	write4(filename, gk_ad_new, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gl_ad_new.txt");
	write4(filename, gl_ad_new, nmkt, nz, nb, ntheta);

	filename = currentpath + ("V_noad_new.txt");
	write4(filename, V_noad_new, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gc_noad_new.txt");
	write4(filename, gc_noad_new, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gtheta_noad_new.txt");
	write4(filename, gtheta_noad_new, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gb_noad_new.txt");
	write4(filename, gb_noad_new, nmkt, nz, nb, ntheta);
	filename = currentpath + ("adjust_noad_new.txt");
	write4(filename, adjust_noad_new, nmkt, nz, nb, ntheta);
	filename = currentpath + ("wealth_noad_new.txt");
	write4(filename, wealth_noad_new, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gk_noad_new.txt");
	write4(filename, gk_noad_new, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gl_noad_new.txt");
	write4(filename, gl_noad_new, nmkt, nz, nb, ntheta);

	filename = currentpath + ("V_new.txt");
	write5(filename, V_new, maxit_DP, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gc_new.txt");
	write4(filename, gc_new, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gtheta_new.txt");
	write4(filename, gtheta_new, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gb_new.txt");
	write4(filename, gb_new, nmkt, nz, nb, ntheta);
	filename = currentpath + ("adjust_new.txt");
	write4(filename, adjust_new, nmkt, nz, nb, ntheta);
	filename = currentpath + ("wealth_new.txt");
	write4(filename, wealth_new, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gk_new.txt");
	write4(filename, gk_new, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gl_new.txt");
	write4(filename, gl_new, nmkt, nz, nb, ntheta);
	









	// Initialization for Transitional Dynamics
	double****** mig_trans = ini_matrix6(period_trans, nmkt, nz, nb_pdf, ntheta_pdf, nmkt);
	double** sum_pdf_A_trans = ini_matrix2(period_trans, nmkt);
	double** error_pdf_A_trans = ini_matrix2(period_trans, nmkt);
	double** sum_pdf_B_trans = ini_matrix2(period_trans, nmkt);
	double** error_pdf_B_trans = ini_matrix2(period_trans, nmkt);
	double***** pdf_A_forward = ini_matrix5(period_trans, nmkt, nz, nb_pdf, ntheta_pdf); // when period_trans=0, the initial dist should be the same as the initial steady state dist
	double***** pdf_B_forward = ini_matrix5(period_trans, nmkt, nz, nb_pdf, ntheta_pdf);

	double**** w_guess = ini_matrix4(maxit_trans_r, maxit_trans_w, period_trans, nmkt);
	double**** w_implied = ini_matrix4(maxit_trans_r, maxit_trans_w, period_trans, nmkt);
	double** r_guess = ini_matrix2(maxit_trans_r, period_trans);
	double** r_implied = ini_matrix2(maxit_trans_r, period_trans);

	double w_temp[period_trans];
	for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
	{
		linspace(w_temp, w_iterate_ini[maxit_bisec_ss_r - 1][maxit_bisec_ss_w - 1][i_mkt], w_iterate_new[maxit_bisec_ss_r - 1][maxit_bisec_ss_w - 1][i_mkt], period_trans);
		for (int i_t = 0; i_t < period_trans; ++i_t)
			w_guess[0][0][i_t][i_mkt] = w_temp[i_t];
	}
	linspace(r_guess[0], r_iterate_ini[maxit_bisec_ss_r - 1], r_iterate_new[maxit_bisec_ss_r - 1], period_trans);

	double**** labor_demand = ini_matrix4(maxit_trans_r, maxit_trans_w, period_trans, nmkt);
	double**** labor_supply = ini_matrix4(maxit_trans_r, maxit_trans_w, period_trans, nmkt);
	double**** capital_demand = ini_matrix4(maxit_trans_r, maxit_trans_w, period_trans, nmkt);
	double**** capital_supply = ini_matrix4(maxit_trans_r, maxit_trans_w, period_trans, nmkt);
	double**** GDP_output_trans = ini_matrix4(maxit_trans_r, maxit_trans_w, period_trans, nmkt);

	// Initialization for Backward Induction
	double***** V_ad_backward = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gc_ad_backward = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gtheta_ad_backward = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gb_ad_backward = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gk_ad_backward = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gl_ad_backward = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** adjust_ad_backward = ini_matrix5(period_trans, nmkt, nz, nb, ntheta); // 1 worker, 2 noborrow entrepreneur, 5 borrow entrepreneur
	double***** wealth_ad_backward = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);

	double***** V_noad_backward = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gc_noad_backward = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gtheta_noad_backward = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gb_noad_backward = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gk_noad_backward = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gl_noad_backward = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** adjust_noad_backward = ini_matrix5(period_trans, nmkt, nz, nb, ntheta); // 1 worker, 2 noborrow entrepreneur, 5 borrow entrepreneur
	double***** wealth_noad_backward = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);

	double***** V_backward = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gc_backward = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gtheta_backward = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gb_backward = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gk_backward = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gl_backward = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** adjust_backward = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** wealth_backward = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);

	// Initialization for Update
	double***** V_ad_update = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gc_ad_update = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gtheta_ad_update = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gb_ad_update = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gk_ad_update = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gl_ad_update = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** adjust_ad_update = ini_matrix5(period_trans, nmkt, nz, nb, ntheta); // 1 worker, 2 noborrow entrepreneur, 5 borrow entrepreneur
	double***** wealth_ad_update = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);

	double***** V_noad_update = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gc_noad_update = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gtheta_noad_update = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gb_noad_update = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gk_noad_update = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gl_noad_update = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** adjust_noad_update = ini_matrix5(period_trans, nmkt, nz, nb, ntheta); // 1 worker, 2 noborrow entrepreneur, 5 borrow entrepreneur
	double***** wealth_noad_update = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);

	double***** V_update = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gc_update = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gtheta_update = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gb_update = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gk_update = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** gl_update = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** adjust_update = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);
	double***** wealth_update = ini_matrix5(period_trans, nmkt, nz, nb, ntheta);

	for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
		for (int i_z = 0; i_z < nz; ++i_z)
			for (int i_b = 0; i_b < nb; ++i_b)
				for (int i_theta = 0; i_theta < ntheta; ++i_theta)
				{
					V_backward[period_trans - 1][i_mkt][i_z][i_b][i_theta] = V_new[0][i_mkt][i_z][i_b][i_theta]; // The last value function iteration is stored at the beginning
					V_update[period_trans - 1][i_mkt][i_z][i_b][i_theta] = V_new[0][i_mkt][i_z][i_b][i_theta];;
				}

	// Start solving transition
	double capital_demand_sum[period_trans];
	double capital_supply_sum[period_trans];
	double GDP_output_trans_sum[period_trans];

	double temp_diff = 0.0;
	int temp = 0;
	int i_w_last;
	int i_r_last;
	for (int i_r = 0; i_r < maxit_trans_r; ++i_r)
	{
		i_r_last = i_r;
		temp = (i_r < maxit_trans_r - 1) ? maxit_trans_w_mid : maxit_trans_w;
		cout << "i_r: " << i_r << endl;
		for (int t = 0; t < period_trans; ++t)
			cout << t << ": " << r_guess[i_r][t] << endl;

		for (int i_w = 0; i_w < temp; ++i_w)
		{
			filename = currentpath + ("progress.txt");	write_progress(filename, i_r, i_w);
			i_w_last = i_w;

			for (int i_t = period_trans - 2; i_t >= 0; --i_t)
				for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
					optimal_decision(Z[i_mkt],
						adjust_ad_backward[i_t][i_mkt], wealth_ad_backward[i_t][i_mkt], gk_ad_backward[i_t][i_mkt], gl_ad_backward[i_t][i_mkt],
						adjust_noad_backward[i_t][i_mkt], wealth_noad_backward[i_t][i_mkt], gk_noad_backward[i_t][i_mkt], gl_noad_backward[i_t][i_mkt],
						gridb, gridtheta, gridz, probz, trans_z, w_guess[i_r][i_w][i_t][i_mkt], r_guess[i_r][i_t], psi[i_t][i_mkt], zeta[i_t][i_mkt]);

			for (int i_t = period_trans - 2; i_t >= 0; --i_t) //calculate value functions along transtion path using backward induction
				for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
					temp_diff = backward(i_t, i_mkt,
						theta_candidate,
						V_ad_backward[i_t][i_mkt], gc_ad_backward[i_t][i_mkt], gtheta_ad_backward[i_t][i_mkt], gb_ad_backward[i_t][i_mkt], gk_ad_backward[i_t][i_mkt], gl_ad_backward[i_t][i_mkt], adjust_ad_backward[i_t][i_mkt], wealth_ad_backward[i_t][i_mkt],
						V_noad_backward[i_t][i_mkt], gc_noad_backward[i_t][i_mkt], gtheta_noad_backward[i_t][i_mkt], gb_noad_backward[i_t][i_mkt], gk_noad_backward[i_t][i_mkt], gl_noad_backward[i_t][i_mkt], adjust_noad_backward[i_t][i_mkt], wealth_noad_backward[i_t][i_mkt],
						V_backward[i_t][i_mkt], V_backward[i_t + 1], gc_backward[i_t][i_mkt], gtheta_backward[i_t][i_mkt], gb_backward[i_t][i_mkt], gk_backward[i_t][i_mkt], gl_backward[i_t][i_mkt], adjust_backward[i_t][i_mkt], wealth_backward[i_t][i_mkt],
						gridb, gridtheta, gridz, probz, trans_z, r_guess[i_r][i_t]);

			for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
				for (int i_z = 0; i_z < nz; ++i_z)
					for (int i_b = 0; i_b < nb_pdf; ++i_b)
						for (int i_theta = 0; i_theta < ntheta_pdf; ++i_theta)
						{
							for (int i_t = 0; i_t < period_trans; ++i_t)
							{
								pdf_A_forward[i_t][i_mkt][i_z][i_b][i_theta] = 0.0;
								pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta] = 0.0;
							}
							pdf_A_forward[0][i_mkt][i_z][i_b][i_theta] = pdf_A_ini[period - 2][i_mkt][i_z][i_b][i_theta];
							pdf_B_forward[0][i_mkt][i_z][i_b][i_theta] = pdf_A_ini[period - 2][i_mkt][i_z][i_b][i_theta];
						}
			forward_simulation(mig_trans,
				pdf_A_forward, pdf_B_forward,
				V_backward, gtheta_backward, gb_backward, gridb, gridtheta, gridb_pdf, gridtheta_pdf, gridz, probz, trans_z,
				b_AtoB_M_index, theta_AtoB_M_index, V_AtoB, step_AtoB_prob, b_BtoA_next_index, theta_BtoA_next_index);
			error_simulation_trans(sum_pdf_A_trans, error_pdf_A_trans, pdf_A_forward);
			error_simulation_trans(sum_pdf_B_trans, error_pdf_B_trans, pdf_B_forward);

			for (int i_t = period_trans - 1; i_t >= 0; --i_t)
				for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
				{
					labor_demand_supply(pdf_B_forward[i_t][i_mkt], gl_backward[i_t][i_mkt], adjust_backward[i_t][i_mkt], gridb, gridtheta, gridb_pdf, gridtheta_pdf,
						&labor_demand[i_r][i_w][i_t][i_mkt], &labor_supply[i_r][i_w][i_t][i_mkt],
						l);
					GDP(Z[i_mkt],
						pdf_B_forward[i_t][i_mkt], gk_backward[i_t][i_mkt], gl_backward[i_t][i_mkt], adjust_backward[i_t][i_mkt], gridb, gridtheta, gridb_pdf, gridtheta_pdf, gridz, probz, trans_z,
						&GDP_output_trans[i_r][i_w][i_t][i_mkt], psi[i_t][i_mkt],
						temp_fin);
					capital_demand_supply(pdf_B_forward[i_t][i_mkt], gk_backward[i_t][i_mkt], adjust_backward[i_t][i_mkt], gridb, gridtheta, gridb_pdf, gridtheta_pdf,
						&capital_demand[i_r][i_w][i_t][i_mkt], &capital_supply[i_r][i_w][i_t][i_mkt], psi[i_t][i_mkt],
						temp_cap_supply);
					capital_demand_sum[i_t] = sum(nmkt, capital_demand[i_r][i_w][i_t]);
					capital_supply_sum[i_t] = sum(nmkt, capital_supply[i_r][i_w][i_t]);
					GDP_output_trans_sum[i_t] = sum(nmkt, GDP_output_trans[i_r][i_w][i_t]);
				}
			update_w(Z,
				sum_pdf_B_trans,
				theta_candidate,
				V_ad_update, gc_ad_update, gtheta_ad_update, gb_ad_update, gk_ad_update, gl_ad_update, adjust_ad_update, wealth_ad_update,
				V_noad_update, gc_noad_update, gtheta_noad_update, gb_noad_update, gk_noad_update, gl_noad_update, adjust_noad_update, wealth_noad_update,
				V_update, gc_update, gtheta_update, gb_update, gk_update, gl_update, adjust_update, wealth_update,
				gridb, gridtheta, gridz, probz, trans_z, w_guess[i_r][i_w], r_guess[i_r], psi, zeta,
				pdf_B_forward, gridb_pdf, gridtheta_pdf,
				w_implied[i_r][i_w], labor_demand[i_r][i_w], labor_supply[i_r][i_w],
				l);

			double wdf = -1e6;
			for (int t = 0; t < period_trans - 1; t++)
				for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
					wdf = max_2num(wdf, absolute(w_implied[i_r][i_w][t][i_mkt] - w_guess[i_r][i_w][t][i_mkt]));

			if (wdf < wconv || i_w == (temp - 1))
			{
				if (i_r != (maxit_trans_r - 1))
				{
					for (int w_update_period = 0; w_update_period < period_trans; ++w_update_period)
						for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
							w_guess[i_r + 1][0][w_update_period][i_mkt] = w_guess[i_r][i_w][w_update_period][i_mkt];
				}
				// Here, V_update is calculated based on w_implied[i_w]; V_backward is calculated based on w_guess[i_w]
				break;
			}
			else
			{
				for (int w_update_period = 0; w_update_period < period_trans; ++w_update_period)
					for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
						w_guess[i_r][i_w + 1][w_update_period][i_mkt] = update_ratio_w * w_implied[i_r][i_w][w_update_period][i_mkt] + (1 - update_ratio_w) * w_guess[i_r][i_w][w_update_period][i_mkt];
			}
		}

		update_r(Z, 
			theta_candidate,
			V_ad_update, gc_ad_update, gtheta_ad_update, gb_ad_update, gk_ad_update, gl_ad_update, adjust_ad_update, wealth_ad_update,
			V_noad_update, gc_noad_update, gtheta_noad_update, gb_noad_update, gk_noad_update, gl_noad_update, adjust_noad_update, wealth_noad_update,
			V_update, gc_update, gtheta_update, gb_update, gk_update, gl_update, adjust_update, wealth_update,
			gridb, gridtheta, gridz, probz, trans_z, w_guess[i_r][i_w_last], r_guess[i_r], psi, zeta,
			pdf_B_forward, gridb_pdf, gridtheta_pdf,
			r_implied[i_r], capital_demand_sum, capital_supply_sum,
			temp_cap_supply);

		double rdf = -1e6;
		for (int t = 0; t < period_trans - 1; t++)
			rdf = max_2num(rdf, absolute(r_implied[i_r][t] - r_guess[i_r][t]));

		if (rdf < rconv || i_r == (maxit_trans_r - 1))
			break;
		else
		{
			for (int r_update_period = 0; r_update_period < period_trans - 1; ++r_update_period)
				r_guess[i_r + 1][r_update_period] = update_ratio_r * r_implied[i_r][r_update_period] + (1 - update_ratio_r) * r_guess[i_r][r_update_period];
		}
	}


	filename = currentpath + ("sum_pdf_A_trans.txt");
	write2(filename, sum_pdf_A_trans, period_trans, nmkt);
	filename = currentpath + ("error_pdf_A_trans.txt");
	write2(filename, error_pdf_A_trans, period_trans, nmkt);
	filename = currentpath + ("sum_pdf_B_trans.txt");
	write2(filename, sum_pdf_B_trans, period_trans, nmkt);
	filename = currentpath + ("error_pdf_B_trans.txt");
	write2(filename, error_pdf_B_trans, period_trans, nmkt);

	filename = currentpath + ("w_guess.txt");
	write4(filename, w_guess, maxit_trans_r, maxit_trans_w, period_trans, nmkt);
	filename = currentpath + ("w_implied.txt");
	write4(filename, w_implied, maxit_trans_r, maxit_trans_w, period_trans, nmkt);
	filename = currentpath + ("r_guess.txt");
	write2(filename, r_guess, maxit_trans_r, period_trans);
	filename = currentpath + ("r_implied.txt");
	write2(filename, r_implied, maxit_trans_r, period_trans);

	filename = currentpath + ("labor_demand.txt");
	write4(filename, labor_demand, maxit_trans_r, maxit_trans_w, period_trans, nmkt);
	filename = currentpath + ("labor_supply.txt");
	write4(filename, labor_supply, maxit_trans_r, maxit_trans_w, period_trans, nmkt);
	filename = currentpath + ("capital_demand.txt");
	write4(filename, capital_demand, maxit_trans_r, maxit_trans_w, period_trans, nmkt);
	filename = currentpath + ("capital_supply.txt");
	write4(filename, capital_supply, maxit_trans_r, maxit_trans_w, period_trans, nmkt);
	filename = currentpath + ("GDP_output_trans.txt");
	write4(filename, GDP_output_trans, maxit_trans_r, maxit_trans_w, period_trans, nmkt);

	filename = currentpath + ("mig_trans.txt");
	write6(filename, mig_trans, period_trans, nmkt, nz, nb_pdf, ntheta_pdf, nmkt);

	filename = currentpath + ("pdf_A_forward.txt");
	write5(filename, pdf_A_forward, period_trans, nmkt, nz, nb_pdf, ntheta_pdf);
	filename = currentpath + ("pdf_B_forward.txt");
	write5(filename, pdf_B_forward, period_trans, nmkt, nz, nb_pdf, ntheta_pdf);


	filename = currentpath + ("V_ad_backward.txt");
	write5(filename, V_ad_backward, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gc_ad_backward.txt");
	write5(filename, gc_ad_backward, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gtheta_ad_backward.txt");
	write5(filename, gtheta_ad_backward, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gb_ad_backward.txt");
	write5(filename, gb_ad_backward, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("adjust_ad_backward.txt");
	write5(filename, adjust_ad_backward, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("wealth_ad_backward.txt");
	write5(filename, wealth_ad_backward, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gk_ad_backward.txt");
	write5(filename, gk_ad_backward, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gl_ad_backward.txt");
	write5(filename, gl_ad_backward, period_trans, nmkt, nz, nb, ntheta);

	filename = currentpath + ("V_noad_backward.txt");
	write5(filename, V_noad_backward, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gc_noad_backward.txt");
	write5(filename, gc_noad_backward, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gtheta_noad_backward.txt");
	write5(filename, gtheta_noad_backward, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gb_noad_backward.txt");
	write5(filename, gb_noad_backward, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("adjust_noad_backward.txt");
	write5(filename, adjust_noad_backward, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("wealth_noad_backward.txt");
	write5(filename, wealth_noad_backward, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gk_noad_backward.txt");
	write5(filename, gk_noad_backward, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gl_noad_backward.txt");
	write5(filename, gl_noad_backward, period_trans, nmkt, nz, nb, ntheta);

	filename = currentpath + ("V_backward.txt");
	write5(filename, V_backward, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gc_backward.txt");
	write5(filename, gc_backward, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gtheta_backward.txt");
	write5(filename, gtheta_backward, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gb_backward.txt");
	write5(filename, gb_backward, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("adjust_backward.txt");
	write5(filename, adjust_backward, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("wealth_backward.txt");
	write5(filename, wealth_backward, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gk_backward.txt");
	write5(filename, gk_backward, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gl_backward.txt");
	write5(filename, gl_backward, period_trans, nmkt, nz, nb, ntheta);


	filename = currentpath + ("V_ad_update.txt");
	write5(filename, V_ad_update, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gc_ad_update.txt");
	write5(filename, gc_ad_update, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gtheta_ad_update.txt");
	write5(filename, gtheta_ad_update, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gb_ad_update.txt");
	write5(filename, gb_ad_update, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("adjust_ad_update.txt");
	write5(filename, adjust_ad_update, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("wealth_ad_update.txt");
	write5(filename, wealth_ad_update, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gk_ad_update.txt");
	write5(filename, gk_ad_update, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gl_ad_update.txt");
	write5(filename, gl_ad_update, period_trans, nmkt, nz, nb, ntheta);

	filename = currentpath + ("V_noad_update.txt");
	write5(filename, V_noad_update, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gc_noad_update.txt");
	write5(filename, gc_noad_update, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gtheta_noad_update.txt");
	write5(filename, gtheta_noad_update, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gb_noad_update.txt");
	write5(filename, gb_noad_update, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("adjust_noad_update.txt");
	write5(filename, adjust_noad_update, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("wealth_noad_update.txt");
	write5(filename, wealth_noad_update, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gk_noad_update.txt");
	write5(filename, gk_noad_update, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gl_noad_update.txt");
	write5(filename, gl_noad_update, period_trans, nmkt, nz, nb, ntheta);

	filename = currentpath + ("V_update.txt");
	write5(filename, V_update, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gc_update.txt");
	write5(filename, gc_update, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gtheta_update.txt");
	write5(filename, gtheta_update, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gb_update.txt");
	write5(filename, gb_update, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("adjust_update.txt");
	write5(filename, adjust_update, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("wealth_update.txt");
	write5(filename, wealth_update, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gk_update.txt");
	write5(filename, gk_update, period_trans, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gl_update.txt");
	write5(filename, gl_update, period_trans, nmkt, nz, nb, ntheta);

	// generate stats
	double*** income_entre = ini_matrix3(period_trans - 1, nmkt, nincome_pdf);
	double*** income_worker = ini_matrix3(period_trans - 1, nmkt, nincome_pdf);

	double*** noninterest_income_worker = ini_matrix3(period_trans - 1, nmkt, nincome_pdf);
	double*** noninterest_income_entre = ini_matrix3(period_trans - 1, nmkt, nincome_pdf);

	double*** wealth_entre = ini_matrix3(period_trans - 1, nmkt, nb_pdf);
	double*** wealth_worker = ini_matrix3(period_trans - 1, nmkt, nb_pdf);

	double* grid_labor = ini_matrix1(nlabor_pdf);	linspace(grid_labor, min_labor, max_labor, nlabor_pdf);
	double* labor_dist = ini_matrix1(nlabor_pdf);

	double* grid_income = ini_matrix1(nincome_pdf);	linspace(grid_income, min_income, max_income, nincome_pdf);

	double** frac_entre = ini_matrix2(period_trans - 1, nmkt);
	double** frac_credit = ini_matrix2(period_trans - 1, nmkt);

	double** L = ini_matrix2(period_trans - 1, nmkt);
	double** K = ini_matrix2(period_trans - 1, nmkt);
	double** Y = ini_matrix2(period_trans - 1, nmkt);
	double** TFP = ini_matrix2(period_trans - 1, nmkt);
	double** GDP = ini_matrix2(period_trans - 1, nmkt);

	double** capital_demand_stat = ini_matrix2(period_trans - 1, nmkt);
	double** capital_supply_stat = ini_matrix2(period_trans - 1, nmkt);
	double** labor_demand_stat = ini_matrix2(period_trans - 1, nmkt);
	double** labor_supply_stat = ini_matrix2(period_trans - 1, nmkt);

	double** cash = ini_matrix2(period_trans - 1, nmkt);
	double** deposit = ini_matrix2(period_trans - 1, nmkt);
	double** credit = ini_matrix2(period_trans - 1, nmkt);

	double occupation = 0.0;
	double m = 0.0;
	double a = 0.0;
	double k = 0.0;
	double labor = 0.0;
	double y = 0.0;
	double fin = 0.0;
	double b = 0.0;
	double theta = 0.0;
	double borrowing = 0.0;
	double adjust_new_help = 0.0;

	double income = 0.0;
	int income_index = 0;
	double noninterest_income = 0.0;
	int noninterest_income_index = 0;

	int labor_index = 0;
	double width_income = grid_income[1] - grid_income[0];
	double width_labor = grid_labor[1] - grid_labor[0];

	for (int i_t = 0; i_t < period_trans - 1; ++i_t)
		for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
		{
			w[i_mkt] = w_guess[i_r_last][i_w_last][i_t][i_mkt];
			r = r_guess[i_r_last][i_t];

			for (int i_z = 0; i_z < nz; ++i_z)
				for (int i_b = 0; i_b < nb; ++i_b)
					for (int i_theta = 0; i_theta < ntheta; ++i_theta)
					{
						if (adjust_backward[i_t][i_mkt][i_z][i_b][i_theta] <= 4.0)
							temp_fin[i_z][i_b][i_theta] = 0.0;
						else
							temp_fin[i_z][i_b][i_theta] = 1.0;


						if (adjust_backward[i_t][i_mkt][i_z][i_b][i_theta] <= 2.0)
							l[i_z][i_b][i_theta] = 1.0;
						else
							l[i_z][i_b][i_theta] = 0.0;


						if (adjust_backward[i_t][i_mkt][i_z][i_b][i_theta] == 1.0 || adjust_backward[i_t][i_mkt][i_z][i_b][i_theta] == 2.0)
							temp_cap_supply[i_z][i_b][i_theta] = gridb[i_b] * (1 - gridtheta[i_theta]);

						else if (adjust_backward[i_t][i_mkt][i_z][i_b][i_theta] == 3.0 || adjust_backward[i_t][i_mkt][i_z][i_b][i_theta] == 4.0)
							temp_cap_supply[i_z][i_b][i_theta] = gridb[i_b] * (1 - gridtheta[i_theta]) + gk_backward[i_t][i_mkt][i_z][i_b][i_theta];

						else if (adjust_backward[i_t][i_mkt][i_z][i_b][i_theta] == 5.0 || adjust_backward[i_t][i_mkt][i_z][i_b][i_theta] == 6.0)
							temp_cap_supply[i_z][i_b][i_theta] = gridb[i_b] * (1 - gridtheta[i_theta]) + (gridb[i_b] * gridtheta[i_theta] - psi[i_t][i_mkt]);


						if (adjust_backward[i_t][i_mkt][i_z][i_b][i_theta] == 2.0 || adjust_backward[i_t][i_mkt][i_z][i_b][i_theta] == 4.0 || adjust_backward[i_t][i_mkt][i_z][i_b][i_theta] == 6.0)
							temp_adjust[i_z][i_b][i_theta] = 1.0;
						else
							temp_adjust[i_z][i_b][i_theta] = 0.0;
					}

			for (int i_z = 0; i_z < nz; ++i_z)
				for (int i_b = 0; i_b < nb_pdf; ++i_b)
					for (int i_theta = 0; i_theta < ntheta_pdf; ++i_theta)
					{
						if (pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta] != 0.0)
						{
							k = lin_interpo_2(nb, ntheta, gridb, gridtheta, gk_backward[i_t][i_mkt][i_z], gridb_pdf[i_b], gridtheta_pdf[i_theta]);
							labor = lin_interpo_2(nb, ntheta, gridb, gridtheta, gl_backward[i_t][i_mkt][i_z], gridb_pdf[i_b], gridtheta_pdf[i_theta]);

							capital_demand_stat[i_t][i_mkt] += k * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
							capital_supply_stat[i_t][i_mkt] += lin_interpo_2(nb, ntheta, gridb, gridtheta, temp_cap_supply[i_z], gridb_pdf[i_b], gridtheta_pdf[i_theta]) * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
							labor_demand_stat[i_t][i_mkt] += labor * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
							labor_supply_stat[i_t][i_mkt] += lin_interpo_2(nb, ntheta, gridb, gridtheta, l[i_z], gridb_pdf[i_b], gridtheta_pdf[i_theta]) * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];

							L[i_t][i_mkt] += labor * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
							K[i_t][i_mkt] += k * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
							if ((k < 1e-10) || (labor < 1e-10))
							{
								k = 0;
								labor = 0;
								y = 0;
							}
							else if ((k > 0.0) && (labor > 0.0))
							{
								y = Z[i_mkt] * gridz[i_z] * pow(pow(k, alpha) * pow(labor, 1 - alpha), 1 - nu);
								frac_entre[i_t][i_mkt] += pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
							}

							Y[i_t][i_mkt] += y * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];

							fin = lin_interpo_2(nb, ntheta, gridb, gridtheta, temp_fin[i_z], gridb_pdf[i_b], gridtheta_pdf[i_theta]);

							if (fin > 0.5)
								frac_credit[i_t][i_mkt] += pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];

							b = gridb_pdf[i_b];
							theta = gridtheta_pdf[i_theta];
							m = b * theta;
							a = b * (1 - theta);

							if (k > m)
								borrowing = k + psi[i_t][i_mkt] - m;
							else
								borrowing = 0.0;

							if (y > 0.0)
								GDP[i_t][i_mkt] += (y - fin * psi[i_t][i_mkt] - borrowing * chi) * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];

							adjust_new_help = lin_interpo_2(nb, ntheta, gridb, gridtheta, temp_adjust[i_z], gridb_pdf[i_b], gridtheta_pdf[i_theta]);

							if ((k == 0.0) && (labor == 0.0) && adjust_new_help < 0.5)
							{
								cash[i_t][i_mkt] += m * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
								deposit[i_t][i_mkt] += a * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
								credit[i_t][i_mkt] += 0 * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];

								income = (w[i_mkt] + a * r);
								noninterest_income = (w[i_mkt] * 1);
								income_index = sear(nincome_pdf, grid_income, width_income, income);
								noninterest_income_index = sear(nincome_pdf, grid_income, width_income, noninterest_income);

								income_worker[i_t][i_mkt][income_index] += pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
								noninterest_income_worker[i_t][i_mkt][noninterest_income_index] += pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
								wealth_worker[i_t][i_mkt][i_b] += pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
							}
							else if ((k == 0.0) && (labor == 0.0) && adjust_new_help >= 0.5)
							{
								cash[i_t][i_mkt] += m * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
								deposit[i_t][i_mkt] += a * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
								credit[i_t][i_mkt] += 0 * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];

								income = (w[i_mkt] + a * r - zeta[i_t][i_mkt]);
								noninterest_income = (w[i_mkt] * 1);
								income_index = sear(nincome_pdf, grid_income, width_income, income);
								noninterest_income_index = sear(nincome_pdf, grid_income, width_income, noninterest_income);

								income_worker[i_t][i_mkt][income_index] += pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
								noninterest_income_worker[i_t][i_mkt][noninterest_income_index] += pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
								wealth_worker[i_t][i_mkt][i_b] += pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
							}
							else if ((k > 0.0) && (labor > 0.0) && fin < 0.5 && adjust_new_help < 0.5)
							{
								cash[i_t][i_mkt] += m * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
								deposit[i_t][i_mkt] += a * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
								credit[i_t][i_mkt] += 0 * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];

								income = (y - w[i_mkt] * labor + a * r - delta * k);
								noninterest_income = income - a * r;
								income_index = sear(nincome_pdf, grid_income, width_income, income);
								noninterest_income_index = sear(nincome_pdf, grid_income, width_income, noninterest_income);

								income_entre[i_t][i_mkt][income_index] += pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
								noninterest_income_entre[i_t][i_mkt][noninterest_income_index] += pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
								wealth_entre[i_t][i_mkt][i_b] += pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
							}
							else if ((k > 0.0) && (labor > 0.0) && fin < 0.5 && adjust_new_help >= 0.5)
							{
								cash[i_t][i_mkt] += m * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
								deposit[i_t][i_mkt] += a * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
								credit[i_t][i_mkt] += 0 * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];


								income = (y - w[i_mkt] * labor + a * r - delta * k - zeta[i_t][i_mkt]);
								noninterest_income = income - (a * r);
								income_index = sear(nincome_pdf, grid_income, width_income, income);
								noninterest_income_index = sear(nincome_pdf, grid_income, width_income, noninterest_income);

								income_entre[i_t][i_mkt][income_index] += pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
								noninterest_income_entre[i_t][i_mkt][noninterest_income_index] += pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
								wealth_entre[i_t][i_mkt][i_b] += pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
							}
							else if ((k > 0.0) && (labor > 0.0) && fin >= 0.5 && adjust_new_help < 0.5)
							{
								cash[i_t][i_mkt] += m * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
								deposit[i_t][i_mkt] += a * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
								credit[i_t][i_mkt] += (k + psi[i_t][i_mkt] - m) * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];

								income = (y - psi[i_t][i_mkt] - w[i_mkt] * labor - (k + psi[i_t][i_mkt] - m) * (r + chi) - delta * k + a * r);
								noninterest_income = income - a * r;
								income_index = sear(nincome_pdf, grid_income, width_income, income);
								noninterest_income_index = sear(nincome_pdf, grid_income, width_income, noninterest_income);

								income_entre[i_t][i_mkt][income_index] += pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
								noninterest_income_entre[i_t][i_mkt][noninterest_income_index] += pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
								wealth_entre[i_t][i_mkt][i_b] += pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
							}
							else if ((k > 0.0) && (labor > 0.0) && fin >= 0.5 && adjust_new_help >= 0.5)
							{
								cash[i_t][i_mkt] += m * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
								deposit[i_t][i_mkt] += a * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
								credit[i_t][i_mkt] += (k + psi[i_t][i_mkt] - m) * pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];

								income = (y - psi[i_t][i_mkt] - w[i_mkt] * labor - (k + psi[i_t][i_mkt] - m) * (r + chi) - delta * k + a * r - zeta[i_t][i_mkt]);
								noninterest_income = income - a * r;
								income_index = sear(nincome_pdf, grid_income, width_income, income);
								noninterest_income_index = sear(nincome_pdf, grid_income, width_income, noninterest_income);

								income_entre[i_t][i_mkt][income_index] += pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
								noninterest_income_entre[i_t][i_mkt][noninterest_income_index] += pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
								wealth_entre[i_t][i_mkt][i_b] += pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
							}

							if (i_t == 0)
							{
								labor_index = sear(nlabor_pdf, grid_labor, width_labor, labor);
								labor_dist[labor_index] += pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta];
							}

						}
					}


		}
	for (int i_t = 0; i_t < period_trans - 1; ++i_t)
		for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
			TFP[i_t][i_mkt] = Y[i_t][i_mkt] / pow(L[i_t][i_mkt], 1 - alpha) / pow(K[i_t][i_mkt], alpha);

	filename = currentpath + ("income_worker.txt");
	write3(filename, income_worker, period_trans - 1, nmkt, nincome_pdf);
	filename = currentpath + ("income_entre.txt");
	write3(filename, income_entre, period_trans - 1, nmkt, nincome_pdf);

	filename = currentpath + ("noninterest_income_worker.txt");
	write3(filename, noninterest_income_worker, period_trans - 1, nmkt, nincome_pdf);
	filename = currentpath + ("noninterest_income_entre.txt");
	write3(filename, noninterest_income_entre, period_trans - 1, nmkt, nincome_pdf);

	filename = currentpath + ("wealth_entre.txt");
	write3(filename, wealth_entre, period_trans - 1, nmkt, nb_pdf);
	filename = currentpath + ("wealth_worker.txt");
	write3(filename, wealth_worker, period_trans - 1, nmkt, nb_pdf);

	filename = currentpath + ("labor_dist.txt");
	write1(filename, labor_dist, nlabor_pdf);
	filename = currentpath + ("grid_income.txt");
	write1(filename, grid_income, nincome_pdf);
	filename = currentpath + ("grid_labor.txt");
	write1(filename, grid_labor, nlabor_pdf);
	filename = currentpath + ("frac_entre.txt");
	write2(filename, frac_entre, period_trans - 1, nmkt);
	filename = currentpath + ("frac_credit.txt");
	write2(filename, frac_credit, period_trans - 1, nmkt);

	filename = currentpath + ("L.txt");
	write2(filename, L, period_trans - 1, nmkt);
	filename = currentpath + ("K.txt");
	write2(filename, K, period_trans - 1, nmkt);
	filename = currentpath + ("Y.txt");
	write2(filename, Y, period_trans - 1, nmkt);
	filename = currentpath + ("TFP.txt");
	write2(filename, TFP, period_trans - 1, nmkt);
	filename = currentpath + ("GDP.txt");
	write2(filename, GDP, period_trans - 1, nmkt);

	filename = currentpath + ("cash.txt");
	write2(filename, cash, period_trans - 1, nmkt);
	filename = currentpath + ("deposit.txt");
	write2(filename, deposit, period_trans - 1, nmkt);
	filename = currentpath + ("credit.txt");
	write2(filename, credit, period_trans - 1, nmkt);

	filename = currentpath + ("capital_demand_stat.txt");
	write2(filename, capital_demand_stat, period_trans - 1, nmkt);
	filename = currentpath + ("capital_supply_stat.txt");
	write2(filename, capital_supply_stat, period_trans - 1, nmkt);
	filename = currentpath + ("labor_demand_stat.txt");
	write2(filename, labor_demand_stat, period_trans - 1, nmkt);
	filename = currentpath + ("labor_supply_stat.txt");
	write2(filename, labor_supply_stat, period_trans - 1, nmkt);



	return 0;
}
