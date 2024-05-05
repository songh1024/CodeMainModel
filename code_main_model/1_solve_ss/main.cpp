#include "parameters.h"
#include <stdio.h>  /* defines FILENAME_MAX */
#include <direct.h>  //windows
//#include <unistd.h>  //linux
//#include <chrono> 

#define GetCurrentDir _getcwd
//using namespace std::chrono; 
using namespace std;

int main()
{
	char cCurrentPath[FILENAME_MAX];
	GetCurrentDir(cCurrentPath, sizeof(cCurrentPath));
	std::string currentpath = std::string(cCurrentPath);
	std::string inputpath = std::string(cCurrentPath);
	currentpath = currentpath.append("/computation_results/"); //linux
	inputpath = inputpath.append("/input/"); //linux
	std::string filename;

	// Initialize markets
	double* Pi = ini_matrix1(nmkt);	
	filename=inputpath+("Pi.dat");	
	read_dat1(Pi, nmkt, filename);
	double Pi_sum = sum(nmkt, Pi);
	double* Z = ini_matrix1(nmkt);
	filename = inputpath + ("Z.dat");
	read_dat1(Z, nmkt, filename);
	double* d = ini_matrix1(nmkt);
	filename = inputpath + ("d.dat");
	read_dat1(d, nmkt, filename);

	double* zeta=ini_matrix1(nmkt);	
	double* psi=ini_matrix1(nmkt);
	for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
	{
		psi[i_mkt] = conspsi + d[i_mkt] * thetapsi;
		zeta[i_mkt] = conszeta + d[i_mkt] * thetazeta;
		Pi[i_mkt] = Pi[i_mkt] / Pi_sum;
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
	write1(filename, d, nmkt);
	filename=currentpath+("zeta.txt");
	write1(filename, zeta, nmkt);
	filename=currentpath+("psi.txt");
	write1(filename, psi, nmkt);

	// Initialize some commonly used variables in iteration, simulation and GDP accounting
	// iteration
	double* theta_candidate = ini_matrix1(coefficient);
	linspace(theta_candidate, mintheta, maxtheta, coefficient);

	// simulation
	int*** b_AtoB_M_index = ini_int_matrix3(nz, nb_pdf, ntheta_pdf);
	int*** theta_AtoB_M_index=ini_int_matrix3(nz, nb_pdf, ntheta_pdf);
	double***** V_AtoB=ini_matrix5(nmkt, nz, nb_pdf, ntheta_pdf, nmkt);
	double***** step_AtoB_prob=ini_matrix5(nmkt, nz, nb_pdf, ntheta_pdf, nmkt);
	int**** b_BtoA_next_index=ini_int_matrix4(nmkt, nz, nb_pdf, ntheta_pdf);
	int**** theta_BtoA_next_index=ini_int_matrix4(nmkt, nz, nb_pdf, ntheta_pdf);
	double** diff=ini_matrix2(maxit_DP, nmkt);

	// demand_supply
	double*** temp_fin = ini_matrix3(nz, nb, ntheta);
	double*** temp_cap_supply=ini_matrix3(nz, nb, ntheta);
	double*** l=ini_matrix3(nz, nb, ntheta);

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

	double**** V_ad=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gc_ad=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gtheta_ad=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gb_ad=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gk_ad=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gl_ad=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** adjust_ad=ini_matrix4(nmkt, nz, nb, ntheta); // 1 worker, 2 noborrow entrepreneur, 4 borrow entrepreneur
	double**** wealth_ad=ini_matrix4(nmkt, nz, nb, ntheta);

	double**** V_noad=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gc_noad=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gtheta_noad=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gb_noad=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gk_noad=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gl_noad=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** adjust_noad=ini_matrix4(nmkt, nz, nb, ntheta); // 1 worker, 2 noborrow entrepreneur, 4 borrow entrepreneur
	double**** wealth_noad=ini_matrix4(nmkt, nz, nb, ntheta);

	double***** V=ini_matrix5(maxit_DP, nmkt, nz, nb, ntheta);
	double**** gc=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gtheta=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gb=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gk=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** gl=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** adjust=ini_matrix4(nmkt, nz, nb, ntheta);
	double**** wealth=ini_matrix4(nmkt, nz, nb, ntheta);

	// Initialization for Simulation
	double* gridb_pdf=ini_matrix1(nb_pdf);
	linspace(gridb_pdf,minb,maxb,nb_pdf);
	double* gridtheta_pdf=ini_matrix1(ntheta_pdf);
	linspace(gridtheta_pdf,mintheta,maxtheta,ntheta_pdf);

	double***** mig= ini_matrix5(period, nmkt, nz, nb_pdf, ntheta_pdf);
	double***** pdf_A=ini_matrix5(period, nmkt, nz, nb_pdf, ntheta_pdf);
	double***** pdf_B=ini_matrix5(period, nmkt, nz, nb_pdf, ntheta_pdf);
	double** sum_pdf_B=ini_matrix2(period, nmkt);
	double** error_pdf_B=ini_matrix2(period, nmkt);
	
	// Initialization for Bisection Search
	double*** w_iterate=ini_matrix3(maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	double* r_iterate=ini_matrix1(maxit_bisec_ss_r);
	double*** labor_demand_iterate=ini_matrix3(maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	double*** labor_supply_iterate=ini_matrix3(maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	double*** capital_demand_iterate=ini_matrix3(maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	double*** capital_supply_iterate=ini_matrix3(maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	double*** GDP_output=ini_matrix3(maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);

	double* w=ini_matrix1(nmkt);
	double r=0.0;

	GE(Pi, Z, psi, zeta, w, &r,

			 V_ad, gc_ad, gtheta_ad, gb_ad, gk_ad, gl_ad, adjust_ad, wealth_ad, 
			 V_noad, gc_noad, gtheta_noad, gb_noad, gk_noad, gl_noad, adjust_noad, wealth_noad, 
			 V, gc, gtheta, gb, gk, gl, adjust, wealth, 
			 gridb, gridtheta, gridz, probz, trans_z,
			 gridb_pdf, gridtheta_pdf,

			 pdf_A, pdf_B,  
			 w_iterate, r_iterate,
			 labor_demand_iterate, labor_supply_iterate, capital_demand_iterate, capital_supply_iterate, GDP_output,
			 mig,

			 sum_pdf_B, error_pdf_B, 
			 
			 diff, theta_candidate, b_AtoB_M_index, theta_AtoB_M_index, V_AtoB, step_AtoB_prob, b_BtoA_next_index, theta_BtoA_next_index, temp_fin, temp_cap_supply, l, currentpath);
	
	for (int i_mkt=0; i_mkt <nmkt; ++i_mkt)
		w_iterate[maxit_bisec_ss_r-1][maxit_bisec_ss_w-1][i_mkt]=w[i_mkt];
	r_iterate[maxit_bisec_ss_r-1]=r;
	

	filename=currentpath+("gridz.txt");
	write1(filename, gridz, nz);
	filename=currentpath+("probz.txt");
	write1(filename, probz, nz);
	filename=currentpath+("trans_z.txt");
	write2(filename, trans_z, nz, nz);
	filename=currentpath+("gridb.txt");
	write1(filename, gridb, nb);
	filename=currentpath+("gridtheta.txt");
	write1(filename, gridtheta, ntheta);

	filename=currentpath+("gridb_pdf.txt");
	write1(filename, gridb_pdf, nb_pdf);
	filename=currentpath+("gridtheta_pdf.txt");
	write1(filename, gridtheta_pdf, ntheta_pdf);


	filename=currentpath+("sum_pdf_B.txt");
	write2(filename, sum_pdf_B, period, nmkt);
	filename=currentpath+("error_pdf_B.txt");
	write2(filename, error_pdf_B, period, nmkt);

	filename=currentpath+("w_iterate.txt");
	write3(filename, w_iterate, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	filename=currentpath+("r_iterate.txt");
	write1(filename, r_iterate, maxit_bisec_ss_r);

	filename=currentpath+("labor_demand_iterate.txt");
	write3(filename, labor_demand_iterate, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	filename=currentpath+("labor_supply_iterate.txt");
	write3(filename, labor_supply_iterate, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	filename=currentpath+("capital_demand_iterate.txt");
	write3(filename, capital_demand_iterate, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	filename=currentpath+("capital_supply_iterate.txt");
	write3(filename, capital_supply_iterate, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);
	filename=currentpath+("GDP_output.txt");
	write3(filename, GDP_output, maxit_bisec_ss_r, maxit_bisec_ss_w, nmkt);

	filename=currentpath+("pdf_A.txt");
	write4(filename, pdf_A[period-2], nmkt, nz, nb_pdf, ntheta_pdf);
	//write5(filename, pdf_A, period, nmkt, nz, nb_pdf, ntheta_pdf);
	filename=currentpath+("pdf_B.txt");
	write4(filename, pdf_B[period-2], nmkt, nz, nb_pdf, ntheta_pdf);
	//write5(filename, pdf_B, period, nmkt, nz, nb_pdf, ntheta_pdf);
	filename = currentpath + ("mig.txt");
	write5(filename, mig, period, nmkt, nz, nb_pdf, ntheta_pdf);

	filename = currentpath + ("V_ad.txt");
	write4(filename, V_ad, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gc_ad.txt");
	write4(filename, gc_ad, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gtheta_ad.txt");
	write4(filename, gtheta_ad, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gb_ad.txt");
	write4(filename, gb_ad, nmkt, nz, nb, ntheta);
	filename = currentpath + ("adjust_ad.txt");
	write4(filename, adjust_ad, nmkt, nz, nb, ntheta);
	filename = currentpath + ("wealth_ad.txt");
	write4(filename, wealth_ad, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gk_ad.txt");
	write4(filename, gk_ad, nmkt, nz, nb, ntheta);
	filename = currentpath + ("gl_ad.txt");
	write4(filename, gl_ad, nmkt, nz, nb, ntheta);

	filename=currentpath+("V_noad.txt");
	write4(filename, V_noad, nmkt, nz, nb, ntheta);
	filename=currentpath+("gc_noad.txt");
	write4(filename, gc_noad, nmkt, nz, nb, ntheta);
	filename=currentpath+("gtheta_noad.txt");
	write4(filename, gtheta_noad, nmkt, nz, nb, ntheta);
	filename=currentpath+("gb_noad.txt");
	write4(filename, gb_noad, nmkt, nz, nb, ntheta);
	filename=currentpath+("adjust_noad.txt");
	write4(filename, adjust_noad, nmkt, nz, nb, ntheta);
	filename=currentpath+("wealth_noad.txt");
	write4(filename, wealth_noad, nmkt, nz, nb, ntheta);
	filename=currentpath+("gk_noad.txt");
	write4(filename, gk_noad, nmkt, nz, nb, ntheta);
	filename=currentpath+("gl_noad.txt");
	write4(filename, gl_noad, nmkt, nz, nb, ntheta);

	filename=currentpath+("V.txt");
	write5(filename, V, maxit_DP, nmkt, nz, nb, ntheta);
	filename=currentpath+("gc.txt");
	write4(filename, gc, nmkt, nz, nb, ntheta);
	filename=currentpath+("gtheta.txt");
	write4(filename, gtheta, nmkt, nz, nb, ntheta);
	filename=currentpath+("gb.txt");
	write4(filename, gb, nmkt, nz, nb, ntheta);
	filename=currentpath+("adjust.txt");
	write4(filename, adjust, nmkt, nz, nb, ntheta);
	filename=currentpath+("wealth.txt");
	write4(filename, wealth, nmkt, nz, nb, ntheta);
	filename=currentpath+("gk.txt");
	write4(filename, gk, nmkt, nz, nb, ntheta);
	filename=currentpath+("gl.txt");
	write4(filename, gl, nmkt, nz, nb, ntheta);

		
	return 0;
}
