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
	std::string inputpath_knap = std::string(cCurrentPath);
	currentpath = currentpath.append("/computation_results/"); //linux
	inputpath = inputpath.append("/input/"); //linux
	inputpath_knap = inputpath_knap.append("/input/knapsack_input/"); //linux
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
	int* bank86 = ini_int_matrix1(nmkt);
	filename = inputpath_knap + ("bank86.dat");
	read_int_dat1(bank86, nmkt, filename);

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
	int* nnewbank = ini_int_matrix1(nyear);
	filename = inputpath_knap + ("nnewbank.dat");
	read_int_dat1(nnewbank, nyear, filename);
	int* seg = ini_int_matrix1(nmkt);
	filename = inputpath_knap + ("seg.dat");
	read_int_dat1(seg, nmkt, filename);
	int** seg_nnewbank_bound = ini_int_matrix2(nseg, 2);
	filename = inputpath_knap + ("seg_nnewbank_bound.dat");
	read_int_dat2(seg_nnewbank_bound, nseg, 2, filename);
	int** seg_year_newbranch = ini_int_matrix2(nyear, nseg);
	filename = inputpath_knap + ("seg_year_newbranch.dat");
	read_int_dat2(seg_year_newbranch, nyear, nseg, filename);
	int** segment_openorder = ini_int_matrix2(nseg, nmkt);
	filename = inputpath_knap + ("segment_openorder.dat");
	read_int_dat2(segment_openorder, nseg, nmkt, filename);

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
	filename = currentpath + ("bank86.txt");
	write_int1(filename, bank86, nmkt);
	filename = currentpath + ("nnewbank.txt");
	write_int1(filename, nnewbank, nyear);
	filename = currentpath + ("seg.txt");
	write_int1(filename, seg, nmkt);
	filename = currentpath + ("seg_nnewbank_bound.txt");
	write_int2(filename, seg_nnewbank_bound, nseg, 2);
	filename = currentpath + ("seg_year_newbranch.txt");
	write_int2(filename, seg_year_newbranch, nyear, nseg);
	filename = currentpath + ("segment_openorder.txt");
	write_int2(filename, segment_openorder, nseg, nmkt);

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

	int* bank_opt = ini_int_matrix1(nmkt);
	for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
	{
		bank_opt[i_mkt] = 99;
		if (bank86[i_mkt] == 1)
			bank_opt[i_mkt] = 0; // open in year 0
	}
	double** d_opt = ini_matrix2(period_trans, nmkt); // corresponding to optimal branch expansion
	for (int t = 0; t < nyear; ++t)
		for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
			d_opt[t][i_mkt] = d_ini[i_mkt];
	for (int i_t = nyear; i_t < period_trans; ++i_t)
		for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
			d_opt[i_t][i_mkt] = d_opt[nyear - 1][i_mkt];
	int** n_opt = ini_int_matrix2(nyear, nseg);
	int** n_cum_opt = ini_int_matrix2(nyear, nseg);
	
	filename = currentpath + ("bank_opt.txt");
	write_int1(filename, bank_opt, nmkt);
	filename = currentpath + ("d_opt.txt");
	write2(filename, d_opt, period_trans, nmkt);
	filename = currentpath + ("n_opt.txt");
	write_int2(filename, n_opt, nyear, nseg);
	filename = currentpath + ("n_cum_opt.txt");
	write_int2(filename, n_cum_opt, nyear, nseg);

	/* Solve knapsack problem year by year */
	for (int t = 1; t < nyear; ++t)
	{
		int** nnewbank_seg = ini_int_matrix2(nseg, 2 * delta_l + 1);
		double** profit_seg = ini_matrix2(nseg, 2 * delta_l + 1);

		/* First, get profit matrix segment by segment. For each l_{g,t}, open the future branches as soon and many as possible, ranked by the order */
		for (int i_seg = 0; i_seg < nseg; ++i_seg)
		{
			vector<int> seg_mktidx;
			for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
			{
				if (seg[i_mkt] == i_seg)
					seg_mktidx.push_back(i_mkt);
			}
			int seg_nmkt = seg_mktidx.size();

			for (int i = 0; i < 2 * delta_l + 1; ++i)
			{
				int l_gt = seg_year_newbranch[t][i_seg] - delta_l + i;
				nnewbank_seg[i_seg][i] = l_gt;
				if (l_gt < 0 || l_gt + n_cum_opt[t - 1][i_seg] > seg_nnewbank_bound[i_seg][1])
				{
					profit_seg[i_seg][i] = -1e6;
					continue;
				}
					

				/* compute profit_seg[i_seg][i] */
				profit_seg[i_seg][i] = compute_profit(currentpath,
					seg_nmkt, seg_mktidx, i_seg, t, bank_opt, n_cum_opt[t-1][i_seg], l_gt, seg_nnewbank_bound[i_seg][1], segment_openorder[i_seg], seg_year_newbranch,
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
			}

		}
		filename = currentpath + ("profit_seg.txt");
		write2(filename, profit_seg, nseg, 2 * delta_l + 1);
		filename = currentpath + ("nnewbank_seg.txt");
		write_int2(filename, nnewbank_seg, nseg, 2 * delta_l + 1);

		/* Solve the knapsack problem */
		knapsack(n_opt[t], nnewbank[t], nnewbank_seg, profit_seg);

		/* Update optimal solution */

		// update bank_opt, n_cum_opt
		for (int i_seg = 0; i_seg < nseg; ++i_seg)
		{
			for (int i = 0; i < n_opt[t][i_seg]; ++i)
				bank_opt[segment_openorder[i_seg][n_cum_opt[t - 1][i_seg] + i] - 1] = t;
			n_cum_opt[t][i_seg] = n_cum_opt[t - 1][i_seg] + n_opt[t][i_seg];
		}

		// update d_opt
		for (int i_t = 0; i_t < nyear; ++i_t)
			for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
			{
				d_opt[i_t][i_mkt] = 1e10;
				for (int j_mkt = 0; j_mkt < nmkt; ++j_mkt)
				{
					if (bank_opt[j_mkt] <= i_t && (access[j_mkt] == 1 || i_mkt == j_mkt) )
						d_opt[i_t][i_mkt] = min_2num(d_opt[i_t][i_mkt], distance[i_mkt][j_mkt]);
				}
			}
		for (int i_t = nyear; i_t < period_trans; ++i_t)
			for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
				d_opt[i_t][i_mkt] = d_opt[nyear - 1][i_mkt];


		filename = currentpath + ("bank_opt.txt");
		write_int1(filename, bank_opt, nmkt);
		filename = currentpath + ("d_opt.txt");
		write2(filename, d_opt, period_trans, nmkt);
		filename = currentpath + ("n_opt.txt");
		write_int2(filename, n_opt, nyear, nseg);
		filename = currentpath + ("n_cum_opt.txt");
		write_int2(filename, n_cum_opt, nyear, nseg);

		release_int_matrix2(nnewbank_seg, nseg, 2 * delta_l + 1);
		release_matrix2(profit_seg, nseg, 2 * delta_l + 1);
	}


	return 0;
}
