#include "parameters.h"

void forward_simulation(int seg_nmkt, std::vector<int> seg_mktidx,
	double****** mig,
	double***** pdf_A, double***** pdf_B,
	double***** V, double***** gtheta, double***** gb, double* gridb, double* gridtheta, double* gridb_pdf, double* gridtheta_pdf, double* gridz_pdf, double* probz_pdf, double** trans_z_pdf,
	int*** b_AtoB_M_index, int*** theta_AtoB_M_index, double***** V_AtoB, double***** step_AtoB_prob, int**** b_BtoA_next_index, int**** theta_BtoA_next_index)
{
	double width_b = gridb_pdf[1] - gridb_pdf[0];
	double width_theta = gridtheta_pdf[1] - gridtheta_pdf[0];

	double m;	double m_M;
	double b_M;
	double theta_M;
	double V_sum;
	double b_next;
	double theta_next;

	// get b_AtoB_M_index, theta_AtoB_M_index
	for (int i_z = 0; i_z < nz; ++i_z)
		for (int i_b = 0; i_b < nb_pdf; ++i_b)
			for (int i_theta = 0; i_theta < ntheta_pdf; ++i_theta)
			{
				m = gridb_pdf[i_b] * gridtheta_pdf[i_theta];
				m_M = m - kappa;
				if (m >= kappa)
				{
					b_M = gridb_pdf[i_b] - kappa;
					theta_M = (b_M < 1e-10) ? 1.0 : (m_M / b_M);
					b_AtoB_M_index[i_z][i_b][i_theta] = min2((int)floor(b_M / width_b), nb_pdf - 1);
					theta_AtoB_M_index[i_z][i_b][i_theta] = min2((int)floor(theta_M / width_theta), ntheta_pdf - 1);
				}
			}

	int t;
	for (t = 0; t < period_sim - 1; ++t)
	{
		// get b_BtoA_next_index, theta_BtoA_next_index
		for (int i_mkt = 0; i_mkt < seg_nmkt; ++i_mkt)
			for (int i_z = 0; i_z < nz; ++i_z)
				for (int i_b = 0; i_b < nb_pdf; ++i_b)
					for (int i_theta = 0; i_theta < ntheta_pdf; ++i_theta)
					{
						b_next = lin_interpo_2(nb, ntheta, gridb, gridtheta, gb[t][seg_mktidx[i_mkt]][i_z], gridb_pdf[i_b], gridtheta_pdf[i_theta]);
						theta_next = lin_interpo_2(nb, ntheta, gridb, gridtheta, gtheta[t][seg_mktidx[i_mkt]][i_z], gridb_pdf[i_b], gridtheta_pdf[i_theta]);
						b_BtoA_next_index[seg_mktidx[i_mkt]][i_z][i_b][i_theta] = min2((int)floor(b_next / width_b), nb_pdf - 1);
						theta_BtoA_next_index[seg_mktidx[i_mkt]][i_z][i_b][i_theta] = min2((int)floor(theta_next / width_theta), ntheta_pdf - 1);
					}

		// from pdf_B to next period pdf_A
		for (int i_mkt = 0; i_mkt < seg_nmkt; ++i_mkt)
			for (int i_z = 0; i_z < nz; ++i_z)
				for (int i_b = 0; i_b < nb_pdf; ++i_b)
					for (int i_theta = 0; i_theta < ntheta_pdf; ++i_theta)
					{
						if (pdf_B[t][seg_mktidx[i_mkt]][i_z][i_b][i_theta] != 0.0)
						{
							for (int i_z_next = 0; i_z_next < nz; ++i_z_next)
								pdf_A[t + 1][seg_mktidx[i_mkt]][i_z_next][b_BtoA_next_index[seg_mktidx[i_mkt]][i_z][i_b][i_theta]][theta_BtoA_next_index[seg_mktidx[i_mkt]][i_z][i_b][i_theta]] +=
								trans_z_pdf[i_z][i_z_next] * pdf_B[t][seg_mktidx[i_mkt]][i_z][i_b][i_theta];
						}
					}

		// get V_AtoB, step_AtoB_prob
		for (int i_mkt = 0; i_mkt < seg_nmkt; ++i_mkt)
			for (int i_z = 0; i_z < nz; ++i_z)
				for (int i_b = 0; i_b < nb_pdf; ++i_b)
					for (int i_theta = 0; i_theta < ntheta_pdf; ++i_theta)
					{
						m = gridb_pdf[i_b] * gridtheta_pdf[i_theta];
						m_M = m - kappa;
						if (m >= kappa)
						{
							b_M = gridb_pdf[i_b] - kappa;
							theta_M = (b_M < 1e-10) ? 1.0 : (m_M / b_M);
							for (int j_mkt = 0; j_mkt < seg_nmkt; ++j_mkt)
							{
								if (j_mkt != i_mkt)
									V_AtoB[seg_mktidx[i_mkt]][i_z][i_b][i_theta][seg_mktidx[j_mkt]] = exp(1 / eta * lin_interpo_2(nb, ntheta, gridb, gridtheta, V[t + 1][seg_mktidx[j_mkt]][i_z], b_M, theta_M));
								else
									V_AtoB[seg_mktidx[i_mkt]][i_z][i_b][i_theta][seg_mktidx[j_mkt]] = exp(1 / eta * lin_interpo_2(nb, ntheta, gridb, gridtheta, V[t + 1][seg_mktidx[j_mkt]][i_z], gridb_pdf[i_b], gridtheta_pdf[i_theta]));
							}
							
							V_sum = sum_select(seg_nmkt, seg_mktidx, V_AtoB[seg_mktidx[i_mkt]][i_z][i_b][i_theta]);
							
							for (int j_mkt = 0; j_mkt < seg_nmkt; ++j_mkt)
								step_AtoB_prob[seg_mktidx[i_mkt]][i_z][i_b][i_theta][seg_mktidx[j_mkt]] = V_AtoB[seg_mktidx[i_mkt]][i_z][i_b][i_theta][seg_mktidx[j_mkt]] / V_sum;
						}
						else
						{
							for (int j_mkt = 0; j_mkt < seg_nmkt; ++j_mkt)
							{
								if (j_mkt != i_mkt)
									step_AtoB_prob[seg_mktidx[i_mkt]][i_z][i_b][i_theta][seg_mktidx[j_mkt]] = 0.0;
								else
									step_AtoB_prob[seg_mktidx[i_mkt]][i_z][i_b][i_theta][seg_mktidx[j_mkt]] = 1.0;
							}
						}
					}

		// from pdf_A to pdf_B, migration
		for (int i_b = 0; i_b < nb_pdf; ++i_b)
			for (int i_z = 0; i_z < nz; ++i_z)
				for (int i_theta = 0; i_theta < ntheta_pdf; ++i_theta)
					for (int i_mkt = 0; i_mkt < seg_nmkt; ++i_mkt)
					{
						for (int j_mkt = 0; j_mkt < seg_nmkt; ++j_mkt)
							mig[t + 1][seg_mktidx[i_mkt]][i_z][i_b][i_theta][seg_mktidx[j_mkt]] = 0.0;
						if (pdf_A[t + 1][seg_mktidx[i_mkt]][i_z][i_b][i_theta] != 0.0)
						{
							m = gridb_pdf[i_b] * gridtheta_pdf[i_theta];
							if (m >= kappa)
							{
								for (int j_mkt = 0; j_mkt < seg_nmkt; ++j_mkt)
								{
									if (j_mkt != i_mkt)
									{
										mig[t + 1][seg_mktidx[i_mkt]][i_z][i_b][i_theta][seg_mktidx[j_mkt]] = pdf_A[t + 1][seg_mktidx[i_mkt]][i_z][i_b][i_theta] * step_AtoB_prob[seg_mktidx[i_mkt]][i_z][i_b][i_theta][seg_mktidx[j_mkt]];
										pdf_B[t + 1][seg_mktidx[j_mkt]][i_z][b_AtoB_M_index[i_z][i_b][i_theta]][theta_AtoB_M_index[i_z][i_b][i_theta]] +=
											mig[t + 1][seg_mktidx[i_mkt]][i_z][i_b][i_theta][seg_mktidx[j_mkt]];
									}
									else
										pdf_B[t + 1][seg_mktidx[j_mkt]][i_z][i_b][i_theta] += pdf_A[t + 1][seg_mktidx[i_mkt]][i_z][i_b][i_theta] * step_AtoB_prob[seg_mktidx[i_mkt]][i_z][i_b][i_theta][seg_mktidx[j_mkt]];
								}
							}
							else
								pdf_B[t + 1][seg_mktidx[i_mkt]][i_z][i_b][i_theta] += pdf_A[t + 1][seg_mktidx[i_mkt]][i_z][i_b][i_theta];

						}
					}


	}


}