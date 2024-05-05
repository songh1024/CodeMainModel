#include "parameters.h"
#include <stdio.h>  /* defines FILENAME_MAX */
#include <direct.h>  //windows
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

	double** d = ini_matrix2(period_trans, nmkt);
	double** zeta = ini_matrix2(period_trans, nmkt);
	double** psi = ini_matrix2(period_trans, nmkt);
	for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
		d[0][i_mkt] = d_ini[i_mkt];

	double** bank_profit = ini_matrix2(nmkt, nmkt); // When each mkt has a new branch, what is the profit of each of the markets

	bool flag = 0;

	for (int ii_mkt = 0; ii_mkt < nmkt; ++ii_mkt)
	{
		// some markets (ii_mkt) already have branch in 86, we do not calculate highest one-year net profit for them
		if (d_ini[ii_mkt] == 0 && flag == 1)
			continue;
		else if(d_ini[ii_mkt] == 0 && flag == 0)
			flag = 1;

		// Steady State
		for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
			d[1][i_mkt] = d[0][i_mkt];
		d[1][ii_mkt] = 0; // open new branch
		if (access[ii_mkt]==1)
		{
			for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
				d[1][i_mkt] = min_2num(distance[i_mkt][ii_mkt], d[1][i_mkt]);
		}

		for (int t = 2; t < period_trans; ++t)
			for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
				d[t][i_mkt] = d[1][i_mkt];



		for (int t = 0; t < period_trans; ++t)
			for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
			{
				psi[t][i_mkt] = conspsi + d[t][i_mkt] * thetapsi;
				zeta[t][i_mkt] = conszeta + d[t][i_mkt] * thetazeta;
			}


		for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
			optimal_decision(Z[i_mkt],
				adjust_ad_new[i_mkt], wealth_ad_new[i_mkt], gk_ad_new[i_mkt], gl_ad_new[i_mkt],
				adjust_noad_new[i_mkt], wealth_noad_new[i_mkt], gk_noad_new[i_mkt], gl_noad_new[i_mkt],
				gridb, gridtheta, gridz, probz, trans_z, w[period_trans - 1][i_mkt], r[period_trans - 1], psi[period_trans - 1][i_mkt], zeta[period_trans - 1][i_mkt]);

		int it = 0;
		for (it = 1; it < maxit_DP; ++it)
		{
			for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
			{
				diff[it][i_mkt] = itvfun_V(it, i_mkt,
					theta_candidate,
					V_ad_new[i_mkt], gc_ad_new[i_mkt], gtheta_ad_new[i_mkt], gb_ad_new[i_mkt], gk_ad_new[i_mkt], gl_ad_new[i_mkt], adjust_ad_new[i_mkt], wealth_ad_new[i_mkt],
					V_noad_new[i_mkt], gc_noad_new[i_mkt], gtheta_noad_new[i_mkt], gb_noad_new[i_mkt], gk_noad_new[i_mkt], gl_noad_new[i_mkt], adjust_noad_new[i_mkt], wealth_noad_new[i_mkt],
					V_new[it][i_mkt], V_new[it - 1], gc_new[i_mkt], gtheta_new[i_mkt], gb_new[i_mkt], gk_new[i_mkt], gl_new[i_mkt], adjust_new[i_mkt], wealth_new[i_mkt],
					gridb, gridtheta, gridz, probz, trans_z, r[period_trans - 1]);
				//cout<<"iter: "<<it<<std::setw(20)<<"difference: "<<diff[it][i_mkt]<<"     ";
			}
			//cout<<endl;
			if (max(nmkt, diff[it]) < tol_itvfun)
				break;
		}

		it = (it != maxit_DP) ? it : (maxit_DP - 1);
		for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
			for (int i_b = 0; i_b < nb; ++i_b)
				for (int i_z = 0; i_z < nz; ++i_z)
					for (int i_theta = 0; i_theta < ntheta; ++i_theta)
						V_new[0][i_mkt][i_z][i_b][i_theta] = V_new[it][i_mkt][i_z][i_b][i_theta];


		// Transition
		for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
			for (int i_z = 0; i_z < nz; ++i_z)
				for (int i_b = 0; i_b < nb; ++i_b)
					for (int i_theta = 0; i_theta < ntheta; ++i_theta)
						V_backward[period_trans - 1][i_mkt][i_z][i_b][i_theta] = V_new[0][i_mkt][i_z][i_b][i_theta]; // The last value function iteration is stored at the beginning

		for (int i_t = period_trans - 2; i_t >= 0; --i_t)
			for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
				optimal_decision(Z[i_mkt],
					adjust_ad_backward[i_t][i_mkt], wealth_ad_backward[i_t][i_mkt], gk_ad_backward[i_t][i_mkt], gl_ad_backward[i_t][i_mkt],
					adjust_noad_backward[i_t][i_mkt], wealth_noad_backward[i_t][i_mkt], gk_noad_backward[i_t][i_mkt], gl_noad_backward[i_t][i_mkt],
					gridb, gridtheta, gridz, probz, trans_z, w[i_t][i_mkt], r[i_t], psi[i_t][i_mkt], zeta[i_t][i_mkt]);

		double temp_diff = 0;
		for (int i_t = period_trans - 2; i_t >= 0; --i_t) //calculate value functions along transtion path using backward induction
			for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
				temp_diff = backward(i_t, i_mkt,
					theta_candidate,
					V_ad_backward[i_t][i_mkt], gc_ad_backward[i_t][i_mkt], gtheta_ad_backward[i_t][i_mkt], gb_ad_backward[i_t][i_mkt], gk_ad_backward[i_t][i_mkt], gl_ad_backward[i_t][i_mkt], adjust_ad_backward[i_t][i_mkt], wealth_ad_backward[i_t][i_mkt],
					V_noad_backward[i_t][i_mkt], gc_noad_backward[i_t][i_mkt], gtheta_noad_backward[i_t][i_mkt], gb_noad_backward[i_t][i_mkt], gk_noad_backward[i_t][i_mkt], gl_noad_backward[i_t][i_mkt], adjust_noad_backward[i_t][i_mkt], wealth_noad_backward[i_t][i_mkt],
					V_backward[i_t][i_mkt], V_backward[i_t + 1], gc_backward[i_t][i_mkt], gtheta_backward[i_t][i_mkt], gb_backward[i_t][i_mkt], gk_backward[i_t][i_mkt], gl_backward[i_t][i_mkt], adjust_backward[i_t][i_mkt], wealth_backward[i_t][i_mkt],
					gridb, gridtheta, gridz, probz, trans_z, r[i_t]);

		for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
			for (int i_z = 0; i_z < nz; ++i_z)
				for (int i_b = 0; i_b < nb_pdf; ++i_b)
					for (int i_theta = 0; i_theta < ntheta_pdf; ++i_theta)
						for (int i_t = 1; i_t < period_sim; ++i_t)
						{
							pdf_A_forward[i_t][i_mkt][i_z][i_b][i_theta] = 0.0;
							pdf_B_forward[i_t][i_mkt][i_z][i_b][i_theta] = 0.0;
						}

		forward_simulation(mig_trans,
			pdf_A_forward, pdf_B_forward,
			V_backward, gtheta_backward, gb_backward, gridb, gridtheta, gridb_pdf, gridtheta_pdf, gridz, probz, trans_z,
			b_AtoB_M_index, theta_AtoB_M_index, V_AtoB, step_AtoB_prob, b_BtoA_next_index, theta_BtoA_next_index);
		error_simulation_trans(sum_pdf_A_trans, error_pdf_A_trans, pdf_A_forward);
		error_simulation_trans(sum_pdf_B_trans, error_pdf_B_trans, pdf_B_forward);

		for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
			cal_bank_profit(i_mkt, bank_profit[ii_mkt], pdf_B_forward[period_sim - 1][i_mkt], gk_backward[period_sim - 1][i_mkt], adjust_backward[period_sim - 1][i_mkt],
				gridb, gridtheta, gridb_pdf, gridtheta_pdf,
				psi[period_sim - 1][i_mkt],
				temp_cap_supply);

		filename = currentpath + ("d.txt");
		write2(filename, d, period_trans, nmkt);
		filename = currentpath + ("psi.txt");
		write2(filename, psi, period_trans, nmkt);
		filename = currentpath + ("zeta.txt");
		write2(filename, zeta, period_trans, nmkt);
		filename = currentpath + ("bank_profit.txt");
		write2(filename, bank_profit, nmkt, nmkt);

		/*filename = currentpath + ("V_ad_new.txt");
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

		filename = currentpath + ("sum_pdf_A_trans.txt");
		write2(filename, sum_pdf_A_trans, period_sim, nmkt);
		filename = currentpath + ("error_pdf_A_trans.txt");
		write2(filename, error_pdf_A_trans, period_sim, nmkt);
		filename = currentpath + ("sum_pdf_B_trans.txt");
		write2(filename, sum_pdf_B_trans, period_sim, nmkt);
		filename = currentpath + ("error_pdf_B_trans.txt");
		write2(filename, error_pdf_B_trans, period_sim, nmkt);

		filename = currentpath + ("mig_trans.txt");
		write6(filename, mig_trans, period_sim, nmkt, nz, nb_pdf, ntheta_pdf, nmkt);

		filename = currentpath + ("pdf_A_forward.txt");
		write5(filename, pdf_A_forward, period_sim, nmkt, nz, nb_pdf, ntheta_pdf);
		filename = currentpath + ("pdf_B_forward.txt");
		write5(filename, pdf_B_forward, period_sim, nmkt, nz, nb_pdf, ntheta_pdf);


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
		write5(filename, gl_backward, period_trans, nmkt, nz, nb, ntheta);*/
	}

	return 0;
}