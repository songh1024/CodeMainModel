#include "parameters.h"
using namespace std;

void update_r(double* Z,
	double* theta_candidate,
	double***** V_ad, double***** gc_ad, double***** gtheta_ad, double***** gb_ad, double***** gk_ad, double***** gl_ad, double***** adjust_ad, double***** wealth_ad,
	double***** V_noad, double***** gc_noad, double***** gtheta_noad, double***** gb_noad, double***** gk_noad, double***** gl_noad, double***** adjust_noad, double***** wealth_noad,
	double***** V, double***** gc, double***** gtheta, double***** gb, double***** gk, double***** gl, double***** adjust, double***** wealth,
	double* gridb, double* gridtheta, double* gridz, double* probz, double** trans_z, double** w_guess, double* r_guess, double** psi, double** zeta,
	double***** pdf_B, double* gridb_pdf, double* gridtheta_pdf,
	double* r_implied, double* capital_demand, double* capital_supply,
	double*** temp_cap_supply)
{
	double* capital_demand_temp = ini_matrix1(nmkt);
	double* capital_supply_temp = ini_matrix1(nmkt);

	double capital_demand_final=0.0;
	double capital_supply_final=0.0;

	double net_capital_demand=0.0;

	double w[nmkt];
	double r=0.0;

	double rlow=0.0;
	double rhigh=0.0;

	double temp_diff=0.0;

	for (int i_t=period_trans-2; i_t>=0; --i_t)
	{
		rlow=minr;
		rhigh=maxr;

		capital_demand_final=capital_demand[i_t];
		capital_supply_final=capital_supply[i_t];

		for (int i_mkt=0; i_mkt <nmkt; ++i_mkt)
			w[i_mkt]=w_guess[i_t][i_mkt];			
		r=r_guess[i_t];
		
		// Initialize V_update
		for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
		{
			optimal_decision(Z[i_mkt],
				adjust_ad[i_t][i_mkt], wealth_ad[i_t][i_mkt], gk_ad[i_t][i_mkt], gl_ad[i_t][i_mkt],
								adjust_noad[i_t][i_mkt], wealth_noad[i_t][i_mkt], gk_noad[i_t][i_mkt], gl_noad[i_t][i_mkt],
								gridb, gridtheta, gridz, probz, trans_z, w[i_mkt], r, psi[i_t][i_mkt], zeta[i_t][i_mkt]);
			temp_diff=backward(i_t, i_mkt,
				theta_candidate, 
				V_ad[i_t][i_mkt], gc_ad[i_t][i_mkt], gtheta_ad[i_t][i_mkt], gb_ad[i_t][i_mkt], gk_ad[i_t][i_mkt], gl_ad[i_t][i_mkt], adjust_ad[i_t][i_mkt], wealth_ad[i_t][i_mkt], 
				V_noad[i_t][i_mkt], gc_noad[i_t][i_mkt], gtheta_noad[i_t][i_mkt], gb_noad[i_t][i_mkt], gk_noad[i_t][i_mkt], gl_noad[i_t][i_mkt], adjust_noad[i_t][i_mkt], wealth_noad[i_t][i_mkt], 
				V[i_t][i_mkt], V[i_t+1], gc[i_t][i_mkt], gtheta[i_t][i_mkt], gb[i_t][i_mkt], gk[i_t][i_mkt], gl[i_t][i_mkt], adjust[i_t][i_mkt], wealth[i_t][i_mkt], 
				gridb, gridtheta, gridz, probz, trans_z, r);
		}

		for (int i_r=0; i_r<maxit_bisec_trans_r; ++i_r)
		{
			net_capital_demand=capital_demand_final-capital_supply_final;
			if (absolute(net_capital_demand)<tol_r_trans* capital_supply_final || (rhigh-rlow)<tol_r_trans || i_r==maxit_bisec_trans_r-1)
			{
				r_implied[i_t]=r;
				break;
			}
			else
			{
				if (net_capital_demand>tol_r_trans* capital_supply_final)
					rlow=r;
				else if (net_capital_demand<-tol_r_trans* capital_supply_final)
					rhigh=r;
				r=(rlow+rhigh)/2;
				// Overall, whenever we have a new guess of r, we update V_update with r. Therefore, V_update is consistent with r_implied

				for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
				{
					optimal_decision(Z[i_mkt],
						adjust_ad[i_t][i_mkt], wealth_ad[i_t][i_mkt], gk_ad[i_t][i_mkt], gl_ad[i_t][i_mkt],
								adjust_noad[i_t][i_mkt], wealth_noad[i_t][i_mkt], gk_noad[i_t][i_mkt], gl_noad[i_t][i_mkt], 
								gridb, gridtheta, gridz, probz, trans_z, w[i_mkt], r, psi[i_t][i_mkt], zeta[i_t][i_mkt]);
					temp_diff=backward(i_t, i_mkt,
						theta_candidate, 
						V_ad[i_t][i_mkt], gc_ad[i_t][i_mkt], gtheta_ad[i_t][i_mkt], gb_ad[i_t][i_mkt], gk_ad[i_t][i_mkt], gl_ad[i_t][i_mkt], adjust_ad[i_t][i_mkt], wealth_ad[i_t][i_mkt], 
						V_noad[i_t][i_mkt], gc_noad[i_t][i_mkt], gtheta_noad[i_t][i_mkt], gb_noad[i_t][i_mkt], gk_noad[i_t][i_mkt], gl_noad[i_t][i_mkt], adjust_noad[i_t][i_mkt], wealth_noad[i_t][i_mkt], 
						V[i_t][i_mkt], V[i_t+1], gc[i_t][i_mkt], gtheta[i_t][i_mkt], gb[i_t][i_mkt], gk[i_t][i_mkt], gl[i_t][i_mkt], adjust[i_t][i_mkt], wealth[i_t][i_mkt], 
						gridb, gridtheta, gridz, probz, trans_z, r);
					
					capital_demand_temp[i_mkt]=0.0;	capital_supply_temp[i_mkt]=0.0;
					capital_demand_supply(pdf_B[i_t][i_mkt], gk[i_t][i_mkt], adjust[i_t][i_mkt], gridb, gridtheta, gridb_pdf, gridtheta_pdf,
						&capital_demand_temp[i_mkt], &capital_supply_temp[i_mkt], psi[i_t][i_mkt],
						temp_cap_supply);
				}
				capital_demand_final=sum(nmkt, capital_demand_temp);
				capital_supply_final=sum(nmkt, capital_supply_temp);
			}
		}
	}

	release_matrix1(capital_demand_temp, nmkt);
	release_matrix1(capital_supply_temp, nmkt);
}

void update_w(double* Z,
	double** sum_B,
	double* theta_candidate,
	double***** V_ad, double***** gc_ad, double***** gtheta_ad, double***** gb_ad, double***** gk_ad, double***** gl_ad, double***** adjust_ad, double***** wealth_ad,
	double***** V_noad, double***** gc_noad, double***** gtheta_noad, double***** gb_noad, double***** gk_noad, double***** gl_noad, double***** adjust_noad, double***** wealth_noad,
	double***** V, double***** gc, double***** gtheta, double***** gb, double***** gk, double***** gl, double***** adjust, double***** wealth,
	double* gridb, double* gridtheta, double* gridz, double* probz, double** trans_z, double** w_guess, double* r_guess, double** psi, double** zeta,
	double***** pdf_B, double* gridb_pdf, double* gridtheta_pdf,
	double** w_implied, double** labor_demand, double** labor_supply,
	double*** l)
{
	double labor_demand_temp[nmkt];
	double labor_supply_temp[nmkt];
	double net_labor_demand[nmkt];

	double w[nmkt];
	double r=0.0;
	double wlow[nmkt];
	double whigh[nmkt];

	double temp_diff=0.0;
	bool mkt_over;

	for (int i_t=period_trans-2; i_t>=0; --i_t)
	{	
		r=r_guess[i_t];

		// Initialize V_update
		for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
		{
			wlow[i_mkt]=minw;
			whigh[i_mkt]=maxw;
			w[i_mkt]=w_guess[i_t][i_mkt];
			labor_demand_temp[i_mkt]=labor_demand[i_t][i_mkt];
			labor_supply_temp[i_mkt]=labor_supply[i_t][i_mkt];

			optimal_decision(Z[i_mkt],
				adjust_ad[i_t][i_mkt], wealth_ad[i_t][i_mkt], gk_ad[i_t][i_mkt], gl_ad[i_t][i_mkt],
								adjust_noad[i_t][i_mkt], wealth_noad[i_t][i_mkt], gk_noad[i_t][i_mkt], gl_noad[i_t][i_mkt], 
								gridb, gridtheta, gridz, probz, trans_z, w[i_mkt], r, psi[i_t][i_mkt], zeta[i_t][i_mkt]);
			temp_diff=backward(i_t, i_mkt,
				theta_candidate, 
				V_ad[i_t][i_mkt], gc_ad[i_t][i_mkt], gtheta_ad[i_t][i_mkt], gb_ad[i_t][i_mkt], gk_ad[i_t][i_mkt], gl_ad[i_t][i_mkt], adjust_ad[i_t][i_mkt], wealth_ad[i_t][i_mkt], 
				V_noad[i_t][i_mkt], gc_noad[i_t][i_mkt], gtheta_noad[i_t][i_mkt], gb_noad[i_t][i_mkt], gk_noad[i_t][i_mkt], gl_noad[i_t][i_mkt], adjust_noad[i_t][i_mkt], wealth_noad[i_t][i_mkt], 
				V[i_t][i_mkt], V[i_t+1], gc[i_t][i_mkt], gtheta[i_t][i_mkt], gb[i_t][i_mkt], gk[i_t][i_mkt], gl[i_t][i_mkt], adjust[i_t][i_mkt], wealth[i_t][i_mkt], 
				gridb, gridtheta, gridz, probz, trans_z, r);
		}

		for (int i_w=0; i_w<maxit_bisec_trans_w; ++i_w)
		{
			mkt_over=true;
			for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
			{
				net_labor_demand[i_mkt] = (labor_demand_temp[i_mkt] - labor_supply_temp[i_mkt]) / sum_B[i_t][i_mkt];
				mkt_over = (mkt_over) && ((absolute(net_labor_demand[i_mkt]) < tol_w_trans) || (whigh[i_mkt] - wlow[i_mkt] < tol_w_trans));
			}
			mkt_over=(mkt_over)||(i_w==maxit_bisec_trans_w-1);
				
			if (mkt_over)
			{
				for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
					w_implied[i_t][i_mkt]=w[i_mkt];
				break;
			}
			else
			{
				for (int i_mkt = 0; i_mkt < nmkt; ++i_mkt)
				{
					if (net_labor_demand[i_mkt]>tol_w_trans)
						wlow[i_mkt]=w[i_mkt];
					else if (net_labor_demand[i_mkt]<-tol_w_trans)
						whigh[i_mkt]=w[i_mkt];

					w[i_mkt]=(wlow[i_mkt]+whigh[i_mkt])/2;
					// Overall, whenever we have a new guess of w, we update V_update with w. Therefore, V_update is consistent with w_implied

					optimal_decision(Z[i_mkt],
						adjust_ad[i_t][i_mkt], wealth_ad[i_t][i_mkt], gk_ad[i_t][i_mkt], gl_ad[i_t][i_mkt],
										adjust_noad[i_t][i_mkt], wealth_noad[i_t][i_mkt], gk_noad[i_t][i_mkt], gl_noad[i_t][i_mkt], 
										gridb, gridtheta, gridz, probz, trans_z, w[i_mkt], r, psi[i_t][i_mkt], zeta[i_t][i_mkt]);
					temp_diff=backward(i_t, i_mkt,
						theta_candidate, 
						V_ad[i_t][i_mkt], gc_ad[i_t][i_mkt], gtheta_ad[i_t][i_mkt], gb_ad[i_t][i_mkt], gk_ad[i_t][i_mkt], gl_ad[i_t][i_mkt], adjust_ad[i_t][i_mkt], wealth_ad[i_t][i_mkt], 
						V_noad[i_t][i_mkt], gc_noad[i_t][i_mkt], gtheta_noad[i_t][i_mkt], gb_noad[i_t][i_mkt], gk_noad[i_t][i_mkt], gl_noad[i_t][i_mkt], adjust_noad[i_t][i_mkt], wealth_noad[i_t][i_mkt], 
						V[i_t][i_mkt], V[i_t+1], gc[i_t][i_mkt], gtheta[i_t][i_mkt], gb[i_t][i_mkt], gk[i_t][i_mkt], gl[i_t][i_mkt], adjust[i_t][i_mkt], wealth[i_t][i_mkt], 
						gridb, gridtheta, gridz, probz, trans_z, r);
				
					labor_demand_temp[i_mkt]=0.0;	labor_supply_temp[i_mkt] =0.0;
					labor_demand_supply(pdf_B[i_t][i_mkt], gl[i_t][i_mkt], adjust[i_t][i_mkt], gridb, gridtheta, gridb_pdf, gridtheta_pdf, 
							&labor_demand_temp[i_mkt], &labor_supply_temp[i_mkt],
							l);
				}		
			}
		}
	}
}
