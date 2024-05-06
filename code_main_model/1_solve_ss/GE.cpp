#include "parameters.h" 
using namespace std;

void GE(double* Pi, double* Z, double* psi, double* zeta, double* w_final, double* r_final,

	double**** V_ad, double**** gc_ad, double**** gtheta_ad, double**** gb_ad, double**** gk_ad, double**** gl_ad, double**** adjust_ad, double**** wealth_ad,
	double**** V_noad, double**** gc_noad, double**** gtheta_noad, double**** gb_noad, double**** gk_noad, double**** gl_noad, double**** adjust_noad, double**** wealth_noad,
	double***** V, double**** gc, double**** gtheta, double**** gb, double**** gk, double**** gl, double**** adjust, double**** wealth,
	double* gridb, double* gridtheta, double* gridz, double* probz, double** trans_z,
	double* gridb_pdf, double* gridtheta_pdf,

	double***** pdf_A, double***** pdf_B,
	double*** w_iterate, double* r_iterate,
	double*** labor_demand_iterate, double*** labor_supply_iterate, double*** capital_demand_iterate, double*** capital_supply_iterate, double*** GDP_output,
	double***** mig,

	double** sum_pdf_B, double** error_pdf_B,

	double** diff, double* theta_candidate, int*** b_AtoB_M_index, int*** theta_AtoB_M_index, double***** V_AtoB, double***** step_AtoB_prob, int**** b_BtoA_next_index, int**** theta_BtoA_next_index, double*** temp_fin, double*** temp_cap_supply, double*** l, std::string currentpath)
{
	std::string filename;
	
	double wlow[nmkt];
	double whigh[nmkt];
	double rlow=minr;
	double rhigh=maxr;
	double w[nmkt];
	double r = minr;

	for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
		w_iterate[0][0][i_mkt]= f_w;
	r_iterate[0]= minr;

	double capital_demand_sum;
	double capital_supply_sum;
	double net_labor_demand[nmkt];
	double it_last=0.0;
	int i_r_last=0; int i_w_last=0;
	bool mkt_over;

	// General Eq
	int temp=0;
	for (int i_r = 0; i_r < maxit_bisec_ss_r; ++i_r) {
		i_r_last = i_r;
		r = r_iterate[i_r];
		temp = (i_r < maxit_bisec_ss_r - 1) ? maxit_bisec_ss_w_mid : maxit_bisec_ss_w;
		for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt) {
			wlow[i_mkt]=minw;
			whigh[i_mkt]=maxw;
		}

		for (int i_w = 0; i_w < temp; ++i_w) {	//for (i_w=0; i_w<1; ++i_w)
			filename = currentpath + ("progress.txt");
            write_progress(filename, i_r, i_w);
			i_w_last = i_w;
			for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
				w[i_mkt] = w_iterate[i_r][i_w][i_mkt];

			for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt) {
                optimal_decision(Z[i_mkt],
                                 adjust_ad[i_mkt], wealth_ad[i_mkt], gk_ad[i_mkt], gl_ad[i_mkt],
                                 adjust_noad[i_mkt], wealth_noad[i_mkt], gk_noad[i_mkt], gl_noad[i_mkt],
                                 gridb, gridtheta, gridz, probz, trans_z, w[i_mkt], r, psi[i_mkt], zeta[i_mkt]);
            }

			int it=0;
			for(it = 1; it < maxit_DP; ++it) {
				for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt) {
					diff[it][i_mkt] = itvfun_V(it, i_mkt,
						theta_candidate, 
						V_ad[i_mkt], gc_ad[i_mkt], gtheta_ad[i_mkt], gb_ad[i_mkt], gk_ad[i_mkt], gl_ad[i_mkt], adjust_ad[i_mkt], wealth_ad[i_mkt],
						V_noad[i_mkt], gc_noad[i_mkt], gtheta_noad[i_mkt], gb_noad[i_mkt], gk_noad[i_mkt], gl_noad[i_mkt], adjust_noad[i_mkt], wealth_noad[i_mkt],
						V[it][i_mkt], V[it-1], gc[i_mkt], gtheta[i_mkt], gb[i_mkt], gk[i_mkt], gl[i_mkt], adjust[i_mkt], wealth[i_mkt],
						gridb, gridtheta, gridz, probz, trans_z, r);
					std::cout << "i_r: " << i_r << " ";
                    std::cout << "i_w: " << i_w << " ";
                    std::cout << "iter: " << it << " ";
                    std::cout << "i_mkt: " << i_mkt <<" ";
                    std::cout << "difference: "<< diff[it][i_mkt] << " ";
                    std::cout << std::endl;
				}
					//cout<<endl;
                if (max(nmkt, diff[it])<tol_itvfun)
                    break;
			}

            // 将最后一次迭代的结果保存到V[0]里面
			it = (it != maxit_DP) ? it :(maxit_DP - 1);
			for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
				for (int i_b = 0; i_b < nb; ++i_b)
					for (int i_z = 0; i_z < nz; ++i_z)
						for (int i_theta = 0; i_theta < ntheta; ++i_theta)
							V[0][i_mkt][i_z][i_b][i_theta] = V[it][i_mkt][i_z][i_b][i_theta];

            // 给pdf_A 和 pdf_B 赋初值
			for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
				for (int i_z = 0; i_z < nz; ++i_z)
					for (int i_b = 0; i_b < nb_pdf; ++i_b) {
						for (int t = 0; t < period; ++t)
							for (int i_theta = 0; i_theta < ntheta_pdf; ++i_theta) {
								pdf_A[t][i_mkt][i_z][i_b][i_theta] = 0.0;
								pdf_B[t][i_mkt][i_z][i_b][i_theta] = 0.0;
							}
						pdf_A[0][i_mkt][i_z][i_b][ntheta_pdf - 1] = probz[i_z] / nb_pdf * Pi[i_mkt];
						pdf_B[0][i_mkt][i_z][i_b][ntheta_pdf - 1] = probz[i_z] / nb_pdf * Pi[i_mkt];
					}

            std::cout << "start simulation" << std::endl;
			simulation_permanent(mig,
				pdf_A, pdf_B,
				V[0], gtheta, gb, gridb, gridtheta, gridb_pdf, gridtheta_pdf, gridz, probz, trans_z,
				b_AtoB_M_index, theta_AtoB_M_index, V_AtoB, step_AtoB_prob, b_BtoA_next_index, theta_BtoA_next_index);

			mkt_over=true;
			for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt) {
				labor_demand_supply(pdf_B[period-2][i_mkt], gl[i_mkt], adjust[i_mkt], gridb, gridtheta, gridb_pdf, gridtheta_pdf,
							&labor_demand_iterate[i_r][i_w][i_mkt], &labor_supply_iterate[i_r][i_w][i_mkt], l);
				net_labor_demand[i_mkt] = (labor_demand_iterate[i_r][i_w][i_mkt] - labor_supply_iterate[i_r][i_w][i_mkt])/Pi[i_mkt];
				mkt_over=(mkt_over)&&((absolute(net_labor_demand[i_mkt])<tol_w_ss)||(whigh[i_mkt] - wlow[i_mkt] <tol_w_ss));
				GDP(Z[i_mkt],
					pdf_B[period-2][i_mkt], gk[i_mkt], gl[i_mkt], adjust[i_mkt], gridb, gridtheta, gridb_pdf, gridtheta_pdf, gridz, probz, trans_z,
						&GDP_output[i_r][i_w][i_mkt], psi[i_mkt],
						temp_fin);
				capital_demand_supply(pdf_B[period - 2][i_mkt], gk[i_mkt], adjust[i_mkt], gridb, gridtheta, gridb_pdf, gridtheta_pdf,
					&capital_demand_iterate[i_r][i_w][i_mkt], &capital_supply_iterate[i_r][i_w][i_mkt], psi[i_mkt],
					temp_cap_supply);
            }

			if (mkt_over || i_w == (temp - 1)) {
				if (i_r != (maxit_bisec_ss_r - 1)) {
					for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
						w_iterate[i_r + 1][0][i_mkt] = w_iterate[i_r][i_w][i_mkt];
				}
				break;
			} else {
				for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt) {
					if (net_labor_demand[i_mkt] < -tol_w_ss)
						whigh[i_mkt] = w[i_mkt];
					else if (net_labor_demand[i_mkt] > tol_w_ss)
						wlow[i_mkt] = w[i_mkt];

					w_iterate[i_r][i_w + 1][i_mkt] = (wlow[i_mkt] + whigh[i_mkt]) / 2.0;
				}
			}
		}

		capital_demand_sum = sum(nmkt, capital_demand_iterate[i_r][i_w_last]);
		capital_supply_sum = sum(nmkt, capital_supply_iterate[i_r][i_w_last]);

		double net_capital_demand= capital_demand_sum - capital_supply_sum;
		if (absolute(net_capital_demand) < tol_r_ss* capital_supply_sum || (rhigh - rlow) < tol_r_ss)
			break;
		else {
			if (net_capital_demand>tol_r_ss* capital_supply_sum)
				rlow = r;
			else if (net_capital_demand<-tol_r_ss* capital_supply_sum)
				rhigh = r;
			if (i_r != (maxit_bisec_ss_r-1))
				r_iterate[i_r+1] = (rlow + rhigh) / 2.0;
		}
	}

	*r_final = r_iterate[i_r_last];
	for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
		w_final[i_mkt] = w_iterate[i_r_last][i_w_last][i_mkt];

	error_simulation(sum_pdf_B, error_pdf_B, pdf_B);
}