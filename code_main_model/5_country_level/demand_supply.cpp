#include "parameters.h" 

void cal_bank_profit(int t, int i_mkt, double** bank_profit, double*** pdf_ss, double*** gk, double*** adjust, 
	double* gridb, double* gridtheta, double* gridb_pdf, double* gridtheta_pdf,
	double psi,
	double*** loan_demand)
{
	for (int i_z = 0; i_z < nz; ++i_z)
		for (int i_b = 0; i_b < nb; ++i_b)
			for (int i_theta = 0; i_theta < ntheta; ++i_theta)
			{
				loan_demand[i_z][i_b][i_theta] = 0;
				if (adjust[i_z][i_b][i_theta] == 5.0 || adjust[i_z][i_b][i_theta] == 6.0)
					loan_demand[i_z][i_b][i_theta] = gk[i_z][i_b][i_theta] + psi - gridb[i_b] * gridtheta[i_theta];
			}

	for (int i_z = 0; i_z < nz; ++i_z)
		for (int i_b = 0; i_b < nb_pdf; ++i_b)
			for (int i_theta = 0; i_theta < ntheta_pdf; ++i_theta)
			{
				if (pdf_ss[i_z][i_b][i_theta] != 0.0)
					bank_profit[t][i_mkt] += lin_interpo_2(nb, ntheta, gridb, gridtheta, loan_demand[i_z], gridb_pdf[i_b], gridtheta_pdf[i_theta]) * pdf_ss[i_z][i_b][i_theta] * chi * pow(beta, t);

			}
}

void capital_demand_supply(double*** pdf_ss, double*** gk, double*** adjust, double* gridb, double* gridtheta,
				   double* gridb_pdf, double* gridtheta_pdf,
					double* capital_demand, double* capital_supply, double psi,
					double*** temp_cap_supply)
{
	
	for (int i_z=0; i_z<nz; ++i_z)
		for (int i_b=0; i_b<nb; ++i_b)
			for (int i_theta=0; i_theta<ntheta; ++i_theta)
			{
				if (adjust[i_z][i_b][i_theta] == 1.0 || adjust[i_z][i_b][i_theta] == 2.0)
					temp_cap_supply[i_z][i_b][i_theta] = gridb[i_b] * (1 - gridtheta[i_theta]);

				else if (adjust[i_z][i_b][i_theta] == 3.0 || adjust[i_z][i_b][i_theta] == 4.0)
					temp_cap_supply[i_z][i_b][i_theta] = gridb[i_b] * (1 - gridtheta[i_theta]) + gk[i_z][i_b][i_theta];

				else if (adjust[i_z][i_b][i_theta] == 5.0 || adjust[i_z][i_b][i_theta] == 6.0)
					temp_cap_supply[i_z][i_b][i_theta] = gridb[i_b] * (1 - gridtheta[i_theta]) + (gridb[i_b] * gridtheta[i_theta] - psi);
			}

	for (int i_z=0; i_z<nz; ++i_z)
		for (int i_b=0; i_b<nb_pdf; ++i_b)
			for (int i_theta=0; i_theta<ntheta_pdf; ++i_theta)
			{
				if (pdf_ss[i_z][i_b][i_theta]!=0.0)
				{
					*capital_demand += lin_interpo_2(nb, ntheta, gridb, gridtheta, gk[i_z], gridb_pdf[i_b], gridtheta_pdf[i_theta]) * pdf_ss[i_z][i_b][i_theta];
					*capital_supply += lin_interpo_2(nb, ntheta, gridb, gridtheta, temp_cap_supply[i_z], gridb_pdf[i_b], gridtheta_pdf[i_theta]) * pdf_ss[i_z][i_b][i_theta];
				}

			}
}

void labor_demand_supply(double*** pdf_ss, double*** gl, double*** adjust, double* gridb, double* gridtheta, 
				   double* gridb_pdf, double* gridtheta_pdf,
					double* labor_demand, double* labor_supply, 
					double*** l)
{
	for (int i_z=0; i_z<nz; ++i_z)
		for (int i_b=0; i_b<nb; ++i_b)
			for (int i_theta=0; i_theta<ntheta; ++i_theta)
			{
				if (adjust[i_z][i_b][i_theta]<=2.0+1e-5)
					l[i_z][i_b][i_theta]=1.0;

				else if (adjust[i_z][i_b][i_theta]>=2.0+1e-5)
					l[i_z][i_b][i_theta]=0.0;
			}

	for (int i_z=0; i_z<nz; ++i_z)
		for (int i_b=0; i_b<nb_pdf; ++i_b)
			for (int i_theta=0; i_theta<ntheta_pdf; ++i_theta)
			{
				if (pdf_ss[i_z][i_b][i_theta]!=0.0)
				{
					*labor_demand+=lin_interpo_2(nb, ntheta, gridb, gridtheta, gl[i_z], gridb_pdf[i_b], gridtheta_pdf[i_theta])*pdf_ss[i_z][i_b][i_theta];
					*labor_supply+=lin_interpo_2(nb, ntheta, gridb, gridtheta, l[i_z], gridb_pdf[i_b], gridtheta_pdf[i_theta])*pdf_ss[i_z][i_b][i_theta];
				}

			}
	
}

void GDP(double p,
	double*** pdf_ss, double*** gk, double*** gl, double*** adjust, double* gridb, double* gridtheta,
	double* gridb_pdf, double* gridtheta_pdf, double* gridz_pdf, double* probz_pdf, double** trans_z_pdf,
	double* GDP_output, double psi,
	double*** temp_fin)
{
	double b=0.0;
	double z=0.0;
	double m = 0.0;
	double k=0.0;
	double l=0.0;
	double y=0.0;
	double fin=0.0;
	double borrow=0.0;
	for (int i_z=0; i_z<nz; ++i_z)
		for (int i_b=0; i_b<nb; ++i_b)
			for (int i_theta=0; i_theta<ntheta; ++i_theta)
			{
				if (adjust[i_z][i_b][i_theta]<=4.0+1e-5)
					temp_fin[i_z][i_b][i_theta]=0.0;
				else
					temp_fin[i_z][i_b][i_theta]=1.0;
			}

	for (int i_z=0; i_z<nz; ++i_z)
		for (int i_b=0; i_b<nb_pdf; ++i_b)
			for (int i_theta=0; i_theta<ntheta_pdf; ++i_theta)
			{
				if (pdf_ss[i_z][i_b][i_theta]!=0.0)
				{
					b = gridb_pdf[i_b];
					z = gridz_pdf[i_z];
					m = b * gridtheta_pdf[i_theta];
					k = lin_interpo_2(nb, ntheta, gridb, gridtheta, gk[i_z], gridb_pdf[i_b], gridtheta_pdf[i_theta]);
					l = lin_interpo_2(nb, ntheta, gridb, gridtheta, gl[i_z], gridb_pdf[i_b], gridtheta_pdf[i_theta]);
					k = (k < 0) ? 0 : k;
					l = (l < 0) ? 0 : l;
				
					fin = lin_interpo_2(nb, ntheta, gridb, gridtheta, temp_fin[i_z], gridb_pdf[i_b], gridtheta_pdf[i_theta]) * psi;
					y = p * z * pow((pow(k, alpha) * pow(l, (1 - alpha))), (1 - nu));
					borrow=(k>m)?(k+psi-m):0.0;
					*GDP_output += (y - fin - borrow * chi) * pdf_ss[i_z][i_b][i_theta];
				}
		}
}