#include "parameters.h" 
using namespace std;

void optimal_decision(double p, 
	double*** adjust_ad, double*** wealth_ad, double*** gk_ad, double*** gl_ad,
	double*** adjust_noad, double*** wealth_noad, double*** gk_noad, double*** gl_noad,
	double* gridb, double* gridtheta, double* gridz, double* probz, double** trans_z,
	double w, double r, double psi, double zeta) {
	double borrowing_multiplier = 1 / ksi;
	double b = 0.0;
	double theta = 0.0;
	double m = 0.0;
	double a = 0.0;
	double z = 0.0;
	double total_wealth_w_noad = 0.0;
	double total_wealth_w_ad = 0.0;
	double total_wealth_e_nobo_noad = 0.0;
	double total_wealth_e_nobo_ad = 0.0;
	double total_wealth_e_bo_noad = 0.0;
	double total_wealth_e_bo_ad = 0.0;

	double* gk_unconstrained_bo_ad = ini_matrix1(nz);
	double* gk_unconstrained_nobo_ad = ini_matrix1(nz);
	double* gk_unconstrained_bo_noad = ini_matrix1(nz);
	double* gk_unconstrained_nobo_noad = ini_matrix1(nz);
	double* gl_coefficient = ini_matrix1(nz);

	double k_coefficient = pow(((1 - alpha) / w), ((1 - alpha) * (1 - nu) / nu)) * pow((1 - nu), (1 / nu));
	double l_coefficient = pow((w / ((1 - nu) * (1 - alpha))), (1 / (alpha * nu - nu - alpha)));
	for (int i_z = 0; i_z < nz; ++i_z) {
		gk_unconstrained_nobo_ad[i_z] = k_coefficient * pow(((delta) / alpha), ((alpha * nu - nu - alpha) / nu)) * pow(p * gridz[i_z], (1 / nu));
		gk_unconstrained_bo_ad[i_z] = k_coefficient * pow(((delta + r + chi) / alpha), ((alpha * nu - nu - alpha) / nu)) * pow(p * gridz[i_z], (1 / nu));
		gk_unconstrained_nobo_noad[i_z] = gk_unconstrained_nobo_ad[i_z];
		gk_unconstrained_bo_noad[i_z] = gk_unconstrained_bo_ad[i_z];
		gl_coefficient[i_z] = l_coefficient * pow((1 / (p * gridz[i_z])), (1 / (alpha * nu - nu - alpha)));
	}

	double optimal_k_nobo_ad = 0.0;
	double optimal_k_bo_ad = 0.0;
	double optimal_k_nobo_noad = 0.0;
	double optimal_k_bo_noad = 0.0;
	double optimal_l_nobo_ad = 0.0;
	double optimal_l_bo_ad = 0.0;
	double optimal_l_nobo_noad = 0.0;
	double optimal_l_bo_noad = 0.0;

	for (int i_b = 0; i_b < nb; ++i_b)
		for (int i_theta = 0; i_theta < ntheta; ++i_theta) {
			b = gridb[i_b];
			theta = gridtheta[i_theta];
			m = b * theta;
			a = b * (1 - theta);
			total_wealth_w_ad = m + a * (1 + r) + w - zeta;
			for (int i_z = 0; i_z < nz; ++i_z) {
				//if ((i_b == 0 && i_z == 0 && i_theta == 0))
					//cout << endl;
				z = gridz[i_z];

				optimal_k_nobo_ad = (gk_unconstrained_nobo_ad[i_z] > m) ? m : gk_unconstrained_nobo_ad[i_z];
				optimal_l_nobo_ad = gl_coefficient[i_z] * pow((1 / pow(optimal_k_nobo_ad, (alpha * (1 - nu)))), (1 / (alpha * nu - nu - alpha)));
				total_wealth_e_nobo_ad = m - optimal_k_nobo_ad + a * (1 + r) + p * z * pow(pow(optimal_k_nobo_ad, alpha) * pow(optimal_l_nobo_ad, (1 - alpha)), (1 - nu)) + (1 - delta) * optimal_k_nobo_ad - w * optimal_l_nobo_ad - zeta;

				if (m > psi)
				{
					optimal_k_bo_ad = (gk_unconstrained_bo_ad[i_z] > ((m - psi) * borrowing_multiplier)) ? ((m - psi) * borrowing_multiplier) : gk_unconstrained_bo_ad[i_z];
					optimal_l_bo_ad = gl_coefficient[i_z] * pow((1 / pow(optimal_k_bo_ad, (alpha * (1 - nu)))), (1 / (alpha * nu - nu - alpha)));
					if (optimal_k_bo_ad > (m - psi))
						total_wealth_e_bo_ad = (m - optimal_k_bo_ad - psi) * (1 + chi + r) + a * (1 + r) + p * z * pow(pow(optimal_k_bo_ad, alpha) * pow(optimal_l_bo_ad, (1 - alpha)), (1 - nu)) + (1 - delta) * optimal_k_bo_ad - w * optimal_l_bo_ad - zeta;
					else
						total_wealth_e_bo_ad = (m - optimal_k_bo_ad - psi) * 1 + a * (1 + r) + p * z * pow(pow(optimal_k_bo_ad, alpha) * pow(optimal_l_bo_ad, (1 - alpha)), (1 - nu)) + (1 - delta) * optimal_k_bo_ad - w * optimal_l_bo_ad - zeta;
				}
				else
				{
					optimal_k_bo_ad = 0.0;
					optimal_l_bo_ad = 0.0;
					total_wealth_e_bo_ad = -1e7; // if an entrepreneur's liquid wealth < psi, then the entrepreneur cannot enter credit market and borrow
				}

				adjust_ad[i_z][i_b][i_theta] = rank3(total_wealth_w_ad, total_wealth_e_nobo_ad, total_wealth_e_bo_ad);

				if (adjust_ad[i_z][i_b][i_theta] == 1)
				{
					wealth_ad[i_z][i_b][i_theta] = total_wealth_w_ad;
					gk_ad[i_z][i_b][i_theta] = 0.0;
					gl_ad[i_z][i_b][i_theta] = 0.0;
				}
				else if (adjust_ad[i_z][i_b][i_theta] == 2)
				{
					wealth_ad[i_z][i_b][i_theta] = total_wealth_e_nobo_ad;
					gk_ad[i_z][i_b][i_theta] = optimal_k_nobo_ad;
					gl_ad[i_z][i_b][i_theta] = optimal_l_nobo_ad;
				}
				else if (adjust_ad[i_z][i_b][i_theta] == 3)
				{
					wealth_ad[i_z][i_b][i_theta] = total_wealth_e_bo_ad;
					gk_ad[i_z][i_b][i_theta] = optimal_k_bo_ad;
					gl_ad[i_z][i_b][i_theta] = optimal_l_bo_ad;
				}
			}
		}


	for (int i_b = 0; i_b < nb; ++i_b)
		for (int i_theta = 0; i_theta < ntheta; ++i_theta)
		{
			b = gridb[i_b];
			theta = gridtheta[i_theta];
			m = b * theta;
			a = b * (1 - theta);
			total_wealth_w_noad = m + w;
			for (int i_z = 0; i_z < nz; ++i_z)
			{
				//if ((i_b == 0 && i_z == 0 && i_theta == 0))
					//cout << endl;
				z = gridz[i_z];

				optimal_k_nobo_noad = (gk_unconstrained_nobo_noad[i_z] > m) ? m : gk_unconstrained_nobo_noad[i_z];
				optimal_l_nobo_noad = gl_coefficient[i_z] * pow((1 / pow(optimal_k_nobo_noad, (alpha * (1 - nu)))), (1 / (alpha * nu - nu - alpha)));
				total_wealth_e_nobo_noad = m - optimal_k_nobo_noad + p * z * pow(pow(optimal_k_nobo_noad, alpha) * pow(optimal_l_nobo_noad, (1 - alpha)), (1 - nu)) + (1 - delta) * optimal_k_nobo_noad - w * optimal_l_nobo_noad;

				if (m > psi)
				{
					optimal_k_bo_noad = (gk_unconstrained_bo_noad[i_z] > ((m - psi) * borrowing_multiplier)) ? ((m - psi) * borrowing_multiplier) : gk_unconstrained_bo_noad[i_z];
					optimal_l_bo_noad = gl_coefficient[i_z] * pow((1 / pow(optimal_k_bo_noad, (alpha * (1 - nu)))), (1 / (alpha * nu - nu - alpha)));
					if (optimal_k_bo_noad > (m - psi))
						total_wealth_e_bo_noad = (m - optimal_k_bo_noad - psi) * (1 + chi + r) + p * z * pow(pow(optimal_k_bo_noad, alpha) * pow(optimal_l_bo_noad, (1 - alpha)), (1 - nu)) + (1 - delta) * optimal_k_bo_noad - w * optimal_l_bo_noad;
					else
						total_wealth_e_bo_noad = (m - optimal_k_bo_noad - psi) * 1 + p * z * pow(pow(optimal_k_bo_noad, alpha) * pow(optimal_l_bo_noad, (1 - alpha)), (1 - nu)) + (1 - delta) * optimal_k_bo_noad - w * optimal_l_bo_noad;
				}
				else
				{
					optimal_k_bo_noad = 0.0;
					optimal_l_bo_noad = 0.0;
					total_wealth_e_bo_noad = -1e7; // if an entrepreneur's liquid wealth < psi, then the entrepreneur cannot enter credit market and borrow
				}

				adjust_noad[i_z][i_b][i_theta] = rank3(total_wealth_w_noad, total_wealth_e_nobo_noad, total_wealth_e_bo_noad);

				if (adjust_noad[i_z][i_b][i_theta] == 1)
				{
					wealth_noad[i_z][i_b][i_theta] = total_wealth_w_noad;
					gk_noad[i_z][i_b][i_theta] = 0.0;
					gl_noad[i_z][i_b][i_theta] = 0.0;
				}
				else if (adjust_noad[i_z][i_b][i_theta] == 2)
				{
					wealth_noad[i_z][i_b][i_theta] = total_wealth_e_nobo_noad;
					gk_noad[i_z][i_b][i_theta] = optimal_k_nobo_noad;
					gl_noad[i_z][i_b][i_theta] = optimal_l_nobo_noad;
				}
				else if (adjust_noad[i_z][i_b][i_theta] == 3)
				{
					wealth_noad[i_z][i_b][i_theta] = total_wealth_e_bo_noad;
					gk_noad[i_z][i_b][i_theta] = optimal_k_bo_noad;
					gl_noad[i_z][i_b][i_theta] = optimal_l_bo_noad;
				}
			}
		}

}


int rank3(double t_w, double t_e_nobo, double t_e_bo)
{
	double result=t_w;
	int occupation=1;
	if (result<t_e_nobo)
	{
		result=t_e_nobo;
		occupation=2;
	}
	if (result<t_e_bo)
	{
		result=t_e_bo;
		occupation=3;
	}
	return occupation;

}