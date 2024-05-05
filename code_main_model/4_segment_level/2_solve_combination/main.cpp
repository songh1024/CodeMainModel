#include "parameters.h"
#include <stdio.h>  /* defines FILENAME_MAX */
#include <direct.h>  //windows
#include <vector>
//#include <unistd.h>  //linux
//#include <chrono> 

#define GetCurrentDir _getcwd
//using namespace std::chrono; 
using namespace std;

double kappa = kappa_new;
double eta = eta_new;

int main()
{
	char cCurrentPath[FILENAME_MAX];
	GetCurrentDir(cCurrentPath, sizeof(cCurrentPath));
	std::string currentpath = std::string(cCurrentPath);
	std::string inputpath = std::string(cCurrentPath);
	std::string inputpath_comb = std::string(cCurrentPath);
	currentpath = currentpath.append("/computation_results/"); //linux
	inputpath = inputpath.append("/input/"); //linux
	inputpath_comb = inputpath_comb.append("/input/combination_input/"); //linux
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
	double*** temp_cap_supply = ini_matrix3(nz, nb, ntheta);

	// prices
	double** w = ini_matrix2(period_trans, nmkt);
	filename = inputpath + ("w.dat");
	read_dat2(w, period_trans, nmkt, filename);
	double* r = ini_matrix1(period_trans);
	filename = inputpath + ("r.dat");
	read_dat1(r, period_trans, filename);

	// Initial distribution and distance
	double**** pdf_ini = ini_matrix4(nmkt, nz, nb_pdf, ntheta_pdf);
	filename = inputpath + ("pdf_ini.dat");
	read_dat4(pdf_ini, nmkt, nz, nb_pdf, ntheta_pdf, filename);

	double** distance = ini_matrix2(nmkt, nmkt);
	filename = inputpath + ("distance.dat");
	read_dat2(distance, nmkt, nmkt, filename);
	double* d_ini = ini_matrix1(nmkt);
	filename = inputpath + ("d_ini.dat");
	read_dat1(d_ini, nmkt, filename);

	// Initialize Markets
	double* Pi = ini_matrix1(nmkt);	
	filename=inputpath+("Pi.dat");	
	read_dat1(Pi, nmkt, filename);

	int* access = ini_int_matrix1(nmkt);
	for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
		access[i_mkt] = (Pi[i_mkt] > h) ? 0 : 1;

	double Pi_sum = sum(nmkt, Pi);
	for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
		Pi[i_mkt] = Pi[i_mkt] / Pi_sum;

	double* Z = ini_matrix1(nmkt);
	filename = inputpath + ("Z.dat");
	read_dat1(Z, nmkt, filename);

	// Initialize Combinations
	int* seg = ini_int_matrix1(nmkt);
	filename = inputpath_comb + ("seg.dat");
	read_int_dat1(seg, nmkt, filename);
	int** seg_open = ini_int_matrix2(nseg, nmkt);
	filename = inputpath_comb + ("seg_open.dat");
	read_int_dat2(seg_open, nseg, nmkt, filename);
	int** seg_nnewbank_bound = ini_int_matrix2(nseg, nmkt);
	filename = inputpath_comb + ("seg_nnewbank_bound.dat");
	read_int_dat2(seg_nnewbank_bound, nseg, 2, filename);


	filename = currentpath + ("parameters.txt");
	write_parameters(filename);
	filename = currentpath + ("parameters_matlab.txt");
	write_parameters_matlab(filename);
	filename = currentpath + ("w.txt");
	write2(filename, w, period_trans, nmkt);
	filename = currentpath + ("r.txt");
	write1(filename, r, period_trans);
	filename = currentpath + ("pdf_ini.txt");
	write4(filename, pdf_ini, nmkt, nz, nb_pdf, ntheta_pdf);
	filename = currentpath + ("distance.txt");
	write2(filename, distance, nmkt, nmkt);
	filename = currentpath + ("d_ini.txt");
	write1(filename, d_ini, nmkt);
	filename = currentpath + ("Pi.txt");
	write1(filename, Pi, nmkt);
	filename = currentpath + ("Z.txt");
	write1(filename, Z, nmkt);
	filename = currentpath + ("access.txt");
	write_int1(filename, access, nmkt);
	filename = currentpath + ("seg.txt");
	write_int1(filename, seg, nmkt);
	filename = currentpath + ("seg_open.txt");
	write_int2(filename, seg_open, nseg, nmkt);
	filename = currentpath + ("seg_nnewbank_bound.txt");
	write_int2(filename, seg_nnewbank_bound, nseg, 2);

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

	// Initialization for New Steady-State
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

	// Initialization for Transitional Dynamics
	double****** mig_trans = ini_matrix6(period_sim, nmkt, nz, nb_pdf, ntheta_pdf, nmkt);
	double** sum_pdf_A_trans = ini_matrix2(period_sim, nmkt);
	double** error_pdf_A_trans = ini_matrix2(period_sim, nmkt);
	double** sum_pdf_B_trans = ini_matrix2(period_sim, nmkt);
	double** error_pdf_B_trans = ini_matrix2(period_sim, nmkt);
	double***** pdf_A_forward = ini_matrix5(period_sim, nmkt, nz, nb_pdf, ntheta_pdf); // 2 is for period 0 and 1
	double***** pdf_B_forward = ini_matrix5(period_sim, nmkt, nz, nb_pdf, ntheta_pdf);
	for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
		for (int i_z = 0; i_z < nz; ++i_z)
			for (int i_b = 0; i_b < nb_pdf; ++i_b)
				for (int i_theta = 0; i_theta < ntheta_pdf; ++i_theta)
					for (int i_t = 1; i_t < period_sim; ++i_t)
					{
						pdf_A_forward[0][i_mkt][i_z][i_b][i_theta] = pdf_ini[i_mkt][i_z][i_b][i_theta];
						pdf_B_forward[0][i_mkt][i_z][i_b][i_theta] = pdf_ini[i_mkt][i_z][i_b][i_theta];
					}

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

	for (int i_seg = 0; i_seg < nseg; ++i_seg)
	{
		/* For a segment, first get the segment level distance matrix, and sure and unsure index */
		vector<int> seg_mktidx;
		vector<int> seg_mktidx_86;
		vector<int> seg_mktidx_sure;
		vector<int> seg_mktidx_unsure;
		for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
		{
			if (seg[i_mkt] == i_seg)
			{
				seg_mktidx.push_back(i_mkt);
				if (seg_open[i_seg][i_mkt] > 0)
					seg_mktidx_sure.push_back(i_mkt);
				if (seg_open[i_seg][i_mkt] == 0)
					seg_mktidx_86.push_back(i_mkt);
				else if (seg_open[i_seg][i_mkt] == -1)
					seg_mktidx_unsure.push_back(i_mkt);
			}
		}
		int seg_nmkt = seg_mktidx.size();
		int seg_nmkt_86 = seg_mktidx_86.size();
		int seg_nmkt_sure = seg_mktidx_sure.size();
		int seg_nmkt_unsure = seg_mktidx_unsure.size();
		filename = currentpath + ("seg_mktidx" + to_string(i_seg) + ".txt");
		write_vecint1(filename, seg_mktidx, seg_nmkt);
		filename = currentpath + ("seg_mktidx_86" + to_string(i_seg) + ".txt");
		write_vecint1(filename, seg_mktidx_86, seg_nmkt_86);
		filename = currentpath + ("seg_mktidx_sure" + to_string(i_seg) + ".txt");
		write_vecint1(filename, seg_mktidx_sure, seg_nmkt_sure);
		filename = currentpath + ("seg_mktidx_unsure" + to_string(i_seg) + ".txt");
		write_vecint1(filename, seg_mktidx_unsure, seg_nmkt_unsure);

		/* Then get combinations */
		int ncase = seg_nnewbank_bound[i_seg][1] - seg_nnewbank_bound[i_seg][0] + 1;
		double* profit = ini_matrix1(ncase);
		int** location = ini_int_matrix2(ncase, nmkt);

		for (int i_case = 0; i_case < ncase; ++i_case)
		{
			int nnewbank = seg_nnewbank_bound[i_seg][0] + i_case;
			int ncomb = num_combination(seg_nmkt_unsure, nnewbank);
			int** mkt_unsure_comb = ini_int_matrix2(ncomb, nnewbank);
			combination(seg_mktidx_unsure, mkt_unsure_comb, ncomb, seg_nmkt_unsure, nnewbank);
			filename = currentpath + ("mkt_unsure_comb.txt");
			write_int2(filename, mkt_unsure_comb, ncomb, nnewbank);

			/* Compute the profit if nnewbank new branch is opened */
			double* profit_comb = ini_matrix1(ncomb);
			for (int i_comb = 0; i_comb < ncomb; ++i_comb)
				profit_comb[i_comb] = compute_profit(currentpath,
					nnewbank, mkt_unsure_comb[i_comb],
					seg_nmkt, seg_mktidx, seg_nmkt_86, seg_mktidx_86, seg_nmkt_sure, seg_mktidx_sure, seg_nmkt_unsure, seg_mktidx_unsure,
					 Pi, Z, distance, w, r, pdf_ini, d_ini,
					 V_ad_new, gc_ad_new, gtheta_ad_new, gb_ad_new, gk_ad_new, gl_ad_new, adjust_ad_new, wealth_ad_new,
					 V_noad_new, gc_noad_new, gtheta_noad_new, gb_noad_new, gk_noad_new, gl_noad_new, adjust_noad_new, wealth_noad_new,
					 V_new, gc_new, gtheta_new, gb_new, gk_new, gl_new, adjust_new, wealth_new,
					 V_ad_backward, gc_ad_backward, gtheta_ad_backward, gb_ad_backward, gk_ad_backward, gl_ad_backward, adjust_ad_backward, wealth_ad_backward,
					 V_noad_backward, gc_noad_backward, gtheta_noad_backward, gb_noad_backward, gk_noad_backward, gl_noad_backward, adjust_noad_backward, wealth_noad_backward,
					 V_backward, gc_backward, gtheta_backward, gb_backward, gk_backward, gl_backward, adjust_backward, wealth_backward,
					 mig_trans, sum_pdf_A_trans, error_pdf_A_trans, sum_pdf_B_trans, error_pdf_B_trans, pdf_A_forward, pdf_B_forward,
					 gridz, probz, trans_z, gridb, gridtheta, gridb_pdf, gridtheta_pdf,
					 theta_candidate, b_AtoB_M_index, theta_AtoB_M_index, V_AtoB, step_AtoB_prob, b_BtoA_next_index, theta_BtoA_next_index, diff, temp_cap_supply,
					access);

			int opt_index = maxindex(ncomb, profit_comb);
			profit[i_case] = profit_comb[opt_index];
			for (int i_mkt = 0; i_mkt < nnewbank; ++i_mkt)
				location[i_case][mkt_unsure_comb[opt_index][i_mkt]] = 1;

			release_int_matrix2(mkt_unsure_comb, ncomb, nnewbank);
		}


		filename = currentpath + ("profit" + to_string(i_seg) + ".txt");
		write1(filename, profit, ncase);
		filename = currentpath + ("location" + to_string(i_seg) + ".txt");
		write_int2(filename, location, ncase, nmkt);
	}


	return 0;
}
