#include "parameters.h" 

double compute_profit(std::string currentpath, 
	int nnewbank, int* mkt_unsure_comb,
	int seg_nmkt, std::vector<int> seg_mktidx, int seg_nmkt_86, std::vector<int> seg_mktidx_86, int seg_nmkt_sure, std::vector<int> seg_mktidx_sure, int seg_nmkt_unsure, std::vector<int> seg_mktidx_unsure,
	double* Pi, double* Z, double** distance, double** w, double* r, double**** pdf_ini, double* d_ini,
	double**** V_ad_new, double**** gc_ad_new, double**** gtheta_ad_new, double**** gb_ad_new, double**** gk_ad_new, double**** gl_ad_new, double**** adjust_ad_new, double**** wealth_ad_new,
	double**** V_noad_new, double**** gc_noad_new, double**** gtheta_noad_new, double**** gb_noad_new, double**** gk_noad_new, double**** gl_noad_new, double**** adjust_noad_new, double**** wealth_noad_new,
	double***** V_new, double**** gc_new, double**** gtheta_new, double**** gb_new, double**** gk_new, double**** gl_new, double**** adjust_new, double**** wealth_new,
	double***** V_ad_backward, double***** gc_ad_backward, double***** gtheta_ad_backward, double***** gb_ad_backward, double***** gk_ad_backward, double***** gl_ad_backward, double***** adjust_ad_backward, double***** wealth_ad_backward,
	double***** V_noad_backward, double***** gc_noad_backward, double***** gtheta_noad_backward, double***** gb_noad_backward, double***** gk_noad_backward, double***** gl_noad_backward, double***** adjust_noad_backward, double***** wealth_noad_backward,
	double***** V_backward, double***** gc_backward, double***** gtheta_backward, double***** gb_backward, double***** gk_backward, double***** gl_backward, double***** adjust_backward, double***** wealth_backward,
	double****** mig_trans, double** sum_pdf_A_trans, double** error_pdf_A_trans, double** sum_pdf_B_trans, double** error_pdf_B_trans, double***** pdf_A_forward, double***** pdf_B_forward,
	double* gridz, double* probz, double** trans_z, double* gridb, double* gridtheta, double* gridb_pdf, double* gridtheta_pdf,
	double* theta_candidate, int*** b_AtoB_M_index, int*** theta_AtoB_M_index, double***** V_AtoB, double***** step_AtoB_prob, int**** b_BtoA_next_index, int**** theta_BtoA_next_index, double** diff, double*** temp_cap_supply,
	int* access)
{
	double* profit = ini_matrix1(seg_nmkt);

	std::string filename;

	double** d = ini_matrix2(period_trans, nmkt);
	double** zeta = ini_matrix2(period_trans, nmkt);
	double** psi = ini_matrix2(period_trans, nmkt);

	// 86 bank
	for (int i_mkt = 0; i_mkt < seg_nmkt; ++i_mkt)
		d[0][seg_mktidx[i_mkt]] = d_ini[seg_mktidx[i_mkt]];

	// open new branch
	for (int i_mkt = 0; i_mkt < seg_nmkt; ++i_mkt)
		d[1][seg_mktidx[i_mkt]] = d[0][seg_mktidx[i_mkt]];
	for (int i_mkt = 0; i_mkt < seg_nmkt_sure; ++i_mkt)
		d[1][seg_mktidx_sure[i_mkt]] = 0;
	for (int i_mkt = 0; i_mkt < nnewbank; ++i_mkt)
		d[1][mkt_unsure_comb[i_mkt]] = 0;
	for (int i_mkt = 0; i_mkt < seg_nmkt; ++i_mkt)
	{
		for (int j_mkt = 0; j_mkt < seg_nmkt_sure; ++j_mkt)
		{
			if (access[seg_mktidx_sure[j_mkt]] == 1)
				d[1][seg_mktidx[i_mkt]] = (distance[seg_mktidx[i_mkt]][seg_mktidx_sure[j_mkt]] < d[1][seg_mktidx[i_mkt]]) ? distance[seg_mktidx[i_mkt]][seg_mktidx_sure[j_mkt]] : d[1][seg_mktidx[i_mkt]];
		}
		for (int j_mkt = 0; j_mkt < nnewbank; ++j_mkt)
		{
			if (access[mkt_unsure_comb[j_mkt]] == 1)
				d[1][seg_mktidx[i_mkt]] = (distance[seg_mktidx[i_mkt]][mkt_unsure_comb[j_mkt]] < d[1][seg_mktidx[i_mkt]]) ? distance[seg_mktidx[i_mkt]][mkt_unsure_comb[j_mkt]] : d[1][seg_mktidx[i_mkt]];
		}
	}

	// no new branch starting 88
	for (int t = 2; t < period_trans; ++t)
		for (int i_mkt = 0; i_mkt < seg_nmkt; ++i_mkt)
			d[t][seg_mktidx[i_mkt]] = d[1][seg_mktidx[i_mkt]];

	for (int t = 0; t < period_trans; ++t)
		for (int i_mkt = 0; i_mkt < seg_nmkt; ++i_mkt)
		{
			psi[t][seg_mktidx[i_mkt]] = conspsi + d[t][seg_mktidx[i_mkt]] * thetapsi;
			zeta[t][seg_mktidx[i_mkt]] = conszeta + d[t][seg_mktidx[i_mkt]] * thetazeta;
		}
	
	filename = currentpath + ("d.txt");
	write2(filename, d, period_trans, nmkt);
	filename = currentpath + ("psi.txt");
	write2(filename, psi, period_trans, nmkt);
	filename = currentpath + ("zeta.txt");
	write2(filename, zeta, period_trans, nmkt);



	for (int i_mkt = 0; i_mkt < seg_nmkt; ++i_mkt)
		optimal_decision(Z[seg_mktidx[i_mkt]],
			adjust_ad_new[seg_mktidx[i_mkt]], wealth_ad_new[seg_mktidx[i_mkt]], gk_ad_new[seg_mktidx[i_mkt]], gl_ad_new[seg_mktidx[i_mkt]],
			adjust_noad_new[seg_mktidx[i_mkt]], wealth_noad_new[seg_mktidx[i_mkt]], gk_noad_new[seg_mktidx[i_mkt]], gl_noad_new[seg_mktidx[i_mkt]],
			gridb, gridtheta, gridz, probz, trans_z, w[period_trans - 1][seg_mktidx[i_mkt]], r[period_trans - 1], psi[period_trans - 1][seg_mktidx[i_mkt]], zeta[period_trans - 1][seg_mktidx[i_mkt]]);

	int it = 0;
	for (it = 1; it < maxit_DP; ++it)
	{
		for (int i_mkt = 0; i_mkt < seg_nmkt; ++i_mkt)
		{
			diff[it][seg_mktidx[i_mkt]] = itvfun_V(seg_nmkt, seg_mktidx,
				it, seg_mktidx[i_mkt],
				theta_candidate,
				V_ad_new[seg_mktidx[i_mkt]], gc_ad_new[seg_mktidx[i_mkt]], gtheta_ad_new[seg_mktidx[i_mkt]], gb_ad_new[seg_mktidx[i_mkt]], gk_ad_new[seg_mktidx[i_mkt]], gl_ad_new[seg_mktidx[i_mkt]], adjust_ad_new[seg_mktidx[i_mkt]], wealth_ad_new[seg_mktidx[i_mkt]],
				V_noad_new[seg_mktidx[i_mkt]], gc_noad_new[seg_mktidx[i_mkt]], gtheta_noad_new[seg_mktidx[i_mkt]], gb_noad_new[seg_mktidx[i_mkt]], gk_noad_new[seg_mktidx[i_mkt]], gl_noad_new[seg_mktidx[i_mkt]], adjust_noad_new[seg_mktidx[i_mkt]], wealth_noad_new[seg_mktidx[i_mkt]],
				V_new[it][seg_mktidx[i_mkt]], V_new[it - 1], gc_new[seg_mktidx[i_mkt]], gtheta_new[seg_mktidx[i_mkt]], gb_new[seg_mktidx[i_mkt]], gk_new[seg_mktidx[i_mkt]], gl_new[seg_mktidx[i_mkt]], adjust_new[seg_mktidx[i_mkt]], wealth_new[seg_mktidx[i_mkt]],
				gridb, gridtheta, gridz, probz, trans_z, r[period_trans - 1]);
			//cout<<"iter: "<<it<<std::setw(20)<<"difference: "<<diff[it][seg_mktidx[i_mkt]]<<"     ";
		}
		//cout<<endl;
		if (max(seg_nmkt, diff[it]) < tol_itvfun)
			break;
	}

	it = (it != maxit_DP) ? it : (maxit_DP - 1);
	for (int i_mkt = 0; i_mkt < seg_nmkt; ++i_mkt)
		for (int i_b = 0; i_b < nb; ++i_b)
			for (int i_z = 0; i_z < nz; ++i_z)
				for (int i_theta = 0; i_theta < ntheta; ++i_theta)
					V_new[0][seg_mktidx[i_mkt]][i_z][i_b][i_theta] = V_new[it][seg_mktidx[i_mkt]][i_z][i_b][i_theta];


	// Transition
	for (int i_mkt = 0; i_mkt < seg_nmkt; ++i_mkt)
		for (int i_z = 0; i_z < nz; ++i_z)
			for (int i_b = 0; i_b < nb; ++i_b)
				for (int i_theta = 0; i_theta < ntheta; ++i_theta)
					V_backward[period_trans - 1][seg_mktidx[i_mkt]][i_z][i_b][i_theta] = V_new[0][seg_mktidx[i_mkt]][i_z][i_b][i_theta]; // The last value function iteration is stored at the beginning

	for (int i_t = period_trans - 2; i_t >= 0; --i_t)
		for (int i_mkt = 0; i_mkt < seg_nmkt; ++i_mkt)
			optimal_decision(Z[seg_mktidx[i_mkt]],
				adjust_ad_backward[i_t][seg_mktidx[i_mkt]], wealth_ad_backward[i_t][seg_mktidx[i_mkt]], gk_ad_backward[i_t][seg_mktidx[i_mkt]], gl_ad_backward[i_t][seg_mktidx[i_mkt]],
				adjust_noad_backward[i_t][seg_mktidx[i_mkt]], wealth_noad_backward[i_t][seg_mktidx[i_mkt]], gk_noad_backward[i_t][seg_mktidx[i_mkt]], gl_noad_backward[i_t][seg_mktidx[i_mkt]],
				gridb, gridtheta, gridz, probz, trans_z, w[i_t][seg_mktidx[i_mkt]], r[i_t], psi[i_t][seg_mktidx[i_mkt]], zeta[i_t][seg_mktidx[i_mkt]]);

	double temp_diff = 0;
	for (int i_t = period_trans - 2; i_t >= 0; --i_t) //calculate value functions along transtion path using backward induction
		for (int i_mkt = 0; i_mkt < seg_nmkt; ++i_mkt)
			temp_diff = backward(seg_nmkt, seg_mktidx,
				i_t, seg_mktidx[i_mkt],
				theta_candidate,
				V_ad_backward[i_t][seg_mktidx[i_mkt]], gc_ad_backward[i_t][seg_mktidx[i_mkt]], gtheta_ad_backward[i_t][seg_mktidx[i_mkt]], gb_ad_backward[i_t][seg_mktidx[i_mkt]], gk_ad_backward[i_t][seg_mktidx[i_mkt]], gl_ad_backward[i_t][seg_mktidx[i_mkt]], adjust_ad_backward[i_t][seg_mktidx[i_mkt]], wealth_ad_backward[i_t][seg_mktidx[i_mkt]],
				V_noad_backward[i_t][seg_mktidx[i_mkt]], gc_noad_backward[i_t][seg_mktidx[i_mkt]], gtheta_noad_backward[i_t][seg_mktidx[i_mkt]], gb_noad_backward[i_t][seg_mktidx[i_mkt]], gk_noad_backward[i_t][seg_mktidx[i_mkt]], gl_noad_backward[i_t][seg_mktidx[i_mkt]], adjust_noad_backward[i_t][seg_mktidx[i_mkt]], wealth_noad_backward[i_t][seg_mktidx[i_mkt]],
				V_backward[i_t][seg_mktidx[i_mkt]], V_backward[i_t + 1], gc_backward[i_t][seg_mktidx[i_mkt]], gtheta_backward[i_t][seg_mktidx[i_mkt]], gb_backward[i_t][seg_mktidx[i_mkt]], gk_backward[i_t][seg_mktidx[i_mkt]], gl_backward[i_t][seg_mktidx[i_mkt]], adjust_backward[i_t][seg_mktidx[i_mkt]], wealth_backward[i_t][seg_mktidx[i_mkt]],
				gridb, gridtheta, gridz, probz, trans_z, r[i_t]);

	for (int i_mkt = 0; i_mkt < seg_nmkt; ++i_mkt)
		for (int i_z = 0; i_z < nz; ++i_z)
			for (int i_b = 0; i_b < nb_pdf; ++i_b)
				for (int i_theta = 0; i_theta < ntheta_pdf; ++i_theta)
					for (int i_t = 1; i_t < period_sim; ++i_t)
					{
						pdf_A_forward[i_t][seg_mktidx[i_mkt]][i_z][i_b][i_theta] = 0.0;
						pdf_B_forward[i_t][seg_mktidx[i_mkt]][i_z][i_b][i_theta] = 0.0;
					}

	forward_simulation(seg_nmkt, seg_mktidx,
		mig_trans,
		pdf_A_forward, pdf_B_forward,
		V_backward, gtheta_backward, gb_backward, gridb, gridtheta, gridb_pdf, gridtheta_pdf, gridz, probz, trans_z,
		b_AtoB_M_index, theta_AtoB_M_index, V_AtoB, step_AtoB_prob, b_BtoA_next_index, theta_BtoA_next_index);
	
	for (int i_mkt = 0; i_mkt < seg_nmkt; ++i_mkt)
		cal_bank_profit(i_mkt, profit, pdf_B_forward[period_sim - 1][seg_mktidx[i_mkt]], gk_backward[period_sim - 1][seg_mktidx[i_mkt]], adjust_backward[period_sim - 1][seg_mktidx[i_mkt]],
			gridb, gridtheta, gridb_pdf, gridtheta_pdf,
			psi[period_sim - 1][seg_mktidx[i_mkt]],
			temp_cap_supply);

	double result = sum(seg_nmkt, profit);
	filename = currentpath + ("profit.txt");
	write1(filename, profit, seg_nmkt);

	release_matrix1(profit, seg_nmkt);
	release_matrix2(d, period_trans, nmkt);
	release_matrix2(zeta, period_trans, nmkt);
	release_matrix2(psi, period_trans, nmkt);

	return result;
}