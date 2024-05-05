#include "parameters.h" 

using namespace std;

void write3_bool(std::string filename, bool*** v, int dim1, int dim2, int dim3)
{
	std::ofstream f;
	//f.open (filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
	f.open(filename.c_str(), std::ios::out | std::ios::trunc);

	for (int i = 0; i < dim1; i++)
		for (int j = 0; j < dim2; j++)
			for (int k = 0; k < dim3; k++)
			{
				if (i == dim1 - 1 && j == dim2 - 1 && k == dim3 - 1)
					f << v[i][j][k];
				else
					f << v[i][j][k] << ",";
			}
	f.close();

}

void write_int3(std::string filename, int*** v, int dim1, int dim2, int dim3)
{
	std::ofstream f;
	//f.open (filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
	f.open(filename.c_str(), std::ios::out | std::ios::trunc);

	for (int i = 0; i < dim1; i++)
		for (int j = 0; j < dim2; j++)
			for (int k = 0; k < dim3; k++)
			{
				if (i == dim1 - 1 && j == dim2 - 1 && k == dim3 - 1)
					f << setprecision(writing_precision) << v[i][j][k];
				else
					f << setprecision(writing_precision) << v[i][j][k] << ",";
			}
	f.close();

}

void write8(std::string filename, double******** v, int dim1, int dim2, int dim3, int dim4, int dim5, int dim6, int dim7, int dim8)
{
	std::ofstream f;
	//f.open (filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
	f.open(filename.c_str(), std::ios::out | std::ios::trunc);
	for (int g = 0; g < dim1; ++g)
		for (int h = 0; h < dim2; ++h)
			for (int i = 0; i < dim3; ++i)
				for (int j = 0; j < dim4; ++j)
					for (int k = 0; k < dim5; ++k)
						for (int l = 0; l < dim6; ++l)
							for (int m = 0; m < dim7; ++m)
								for (int n = 0; n < dim8; ++n)

								{
									if (g == dim1 - 1 && h == dim2 - 1 && i == dim3 - 1 && j == dim4 - 1 && k == dim5 - 1 && l == dim6 - 1 && m == dim7 - 1 && n == dim8 - 1)
										f << setprecision(writing_precision) << v[g][h][i][j][k][l][m][n];
									else
										f << setprecision(writing_precision) << v[g][h][i][j][k][l][m][n] << ",";
								}

	f.close();

}

void write7(std::string filename, double******* v, int dim1, int dim2, int dim3, int dim4, int dim5, int dim6, int dim7)
{
	std::ofstream f;
	//f.open (filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
	f.open(filename.c_str(), std::ios::out | std::ios::trunc);
	for (int g = 0; g < dim1; ++g)
		for (int h = 0; h < dim2; ++h)
			for (int i = 0; i < dim3; ++i)
				for (int j = 0; j < dim4; ++j)
					for (int k = 0; k < dim5; ++k)
						for (int l = 0; l < dim6; ++l)
							for (int m = 0; m < dim7; ++m)

							{
								if (g == dim1 - 1 && h == dim2 - 1 && i == dim3 - 1 && j == dim4 - 1 && k == dim5 - 1 && l == dim6 - 1 && m == dim7 - 1)
									f << setprecision(writing_precision) << v[g][h][i][j][k][l][m];
								else
									f << setprecision(writing_precision) << v[g][h][i][j][k][l][m] << ",";
							}

	f.close();

}

void write6(std::string filename, double****** v, int dim1, int dim2, int dim3, int dim4, int dim5, int dim6)
{
	std::ofstream f;
	//f.open (filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
	f.open(filename.c_str(), std::ios::out | std::ios::trunc);
	for (int g = 0; g < dim1; ++g)
		for (int h = 0; h < dim2; ++h)
			for (int i = 0; i < dim3; ++i)
				for (int j = 0; j < dim4; ++j)
					for (int k = 0; k < dim5; ++k)
						for (int l = 0; l < dim6; ++l)

						{
							if (g == dim1 - 1 && h == dim2 - 1 && i == dim3 - 1 && j == dim4 - 1 && k == dim5 - 1 && l == dim6 - 1)
								f << setprecision(writing_precision) << v[g][h][i][j][k][l];
							else
								f << setprecision(writing_precision) << v[g][h][i][j][k][l] << ",";
						}

	f.close();

}

void write5(std::string filename, double***** v, int dim1, int dim2, int dim3, int dim4, int dim5)
{
	std::ofstream f;
	//f.open (filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
	f.open(filename.c_str(), std::ios::out | std::ios::trunc);
	for (int g = 0; g < dim1; ++g)
		for (int h = 0; h < dim2; ++h)
			for (int i = 0; i < dim3; ++i)
				for (int j = 0; j < dim4; ++j)
					for (int k = 0; k < dim5; ++k)


					{
						if (g == dim1 - 1 && h == dim2 - 1 && i == dim3 - 1 && j == dim4 - 1 && k == dim5 - 1)
							f << setprecision(writing_precision) << v[g][h][i][j][k];
						else
							f << setprecision(writing_precision) << v[g][h][i][j][k] << ",";
					}

	f.close();

}


void write4(std::string filename, double**** v, int dim1, int dim2, int dim3, int dim4)
{
	std::ofstream f;
	//f.open (filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
	f.open(filename.c_str(), std::ios::out | std::ios::trunc);
	for (int g = 0; g < dim1; ++g)
		for (int h = 0; h < dim2; ++h)
			for (int i = 0; i < dim3; ++i)
				for (int j = 0; j < dim4; ++j)

				{
					if (g == dim1 - 1 && h == dim2 - 1 && i == dim3 - 1 && j == dim4 - 1)
						f << setprecision(writing_precision) << v[g][h][i][j];
					else
						f << setprecision(writing_precision) << v[g][h][i][j] << ",";
				}

	f.close();

}

void write3(std::string filename, double*** v, int dim1, int dim2, int dim3)
{
	std::ofstream f;
	//f.open (filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
	f.open(filename.c_str(), std::ios::out | std::ios::trunc);

	for (int i = 0; i < dim1; i++)
		for (int j = 0; j < dim2; j++)
			for (int k = 0; k < dim3; k++)
			{
				if (i == dim1 - 1 && j == dim2 - 1 && k == dim3 - 1)
					f << setprecision(writing_precision) << v[i][j][k];
				else
					f << setprecision(writing_precision) << v[i][j][k] << ",";
			}
	f.close();

}


void write2(std::string filename, double** v, int dim1, int dim2) //write grids to file
{
	std::ofstream f;
	//f.open (filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
	f.open(filename.c_str(), std::ios::out | std::ios::trunc);
	for (int i = 0; i < dim1; ++i)
		for (int j = 0; j < dim2; ++j)
		{
			if (i == dim1 - 1 && j == dim2 - 1)
				f << setprecision(writing_precision) << v[i][j];
			else
				f << setprecision(writing_precision) << v[i][j] << ",";
		}
	f.close();
}

void write1(std::string filename, double* v, int dim1) //write grids to file
{
	std::ofstream f;
	//f.open (filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
	f.open(filename.c_str(), std::ios::out | std::ios::trunc);
	for (int i = 0; i < dim1; ++i)
	{
		if (i == dim1 - 1)
			f << setprecision(writing_precision) << v[i];
		else
			f << setprecision(writing_precision) << v[i] << ",";
	}
	f.close();
}

void write_int2(std::string filename, int** v, int dim1, int dim2) //write grids to file
{
	std::ofstream f;
	//f.open (filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
	f.open(filename.c_str(), std::ios::out | std::ios::trunc);
	for (int i = 0; i < dim1; ++i)
		for (int j = 0; j < dim2; ++j)
		{
			if (i == dim1 - 1 && j == dim2 - 1)
				f << setprecision(writing_precision) << v[i][j];
			else
				f << setprecision(writing_precision) << v[i][j] << ",";
		}
	f.close();
}

void write_int1(std::string filename, int* v, int dim1) //write grids to file
{
	std::ofstream f;
	//f.open (filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
	f.open(filename.c_str(), std::ios::out | std::ios::trunc);
	for (int i = 0; i < dim1; ++i)
	{
		if (i == dim1 - 1)
			f << setprecision(writing_precision) << v[i];
		else
			f << setprecision(writing_precision) << v[i] << ",";
	}
	f.close();
}

void write_vecint1(std::string filename, std::vector<int>& v, int dim1)
{
	std::ofstream f;
	//f.open (filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
	f.open(filename.c_str(), std::ios::out | std::ios::trunc);
	for (int i = 0; i < dim1; ++i)
	{
		if (i == dim1 - 1)
			f << setprecision(writing_precision) << v[i];
		else
			f << setprecision(writing_precision) << v[i] << ",";
	}
	f.close();
}



void write_parameters(std::string filename)
{
	std::ofstream f;
	//f.open (filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
	f.open(filename.c_str(), std::ios::out | std::ios::trunc);



	f << "beta\t" << beta << std::endl;
	f << "sigma\t" << sigma << std::endl;

	f << "shape\t" << shape << std::endl;
	f << "scale\t" << scale << std::endl;
	f << "cutoff\t" << cutoff << std::endl;
	f << "grammar\t" << grammar << std::endl;

	f << "time_preference_rate\t" << time_preference_rate << std::endl;

	f << "delta\t" << delta << std::endl;
	f << "alpha\t" << alpha << std::endl;
	f << "nu\t" << nu << std::endl;

	f << "ksi\t" << ksi << std::endl;
	f << "chi\t" << chi << std::endl;
	f << "conspsi\t" << conspsi << std::endl;
	f << "thetapsi\t" << thetapsi << std::endl;
	f << "conszeta\t" << conszeta << std::endl;
	f << "thetazeta\t" << thetazeta << std::endl;

	f << "kappa_ini\t" << kappa_ini << std::endl;
	f << "eta_ini\t" << eta_ini << std::endl;
	f << "kappa_new\t" << kappa_new << std::endl;
	f << "eta_new\t" << eta_new << std::endl;

	f << "nb\t" << nb << std::endl;
	f << "ntheta\t" << ntheta << std::endl;
	f << "nz\t" << nz << std::endl;

	f << "minb\t" << minb << std::endl;
	f << "maxb\t" << maxb << std::endl;
	f << "mintheta\t" << mintheta << std::endl;
	f << "maxtheta\t" << maxtheta << std::endl;

	f << "maxit_DP\t" << maxit_DP << std::endl;
	f << "maxit_gss\t" << maxit_gss << std::endl;
	f << "tol_itvfun\t" << tol_itvfun << std::endl;
	f << "tol_gss\t" << tol_gss << std::endl;
	f << "golden_ratio\t" << golden_ratio << std::endl;

	f << "coefficient\t" << coefficient << std::endl;

	f << "writing_precision\t" << writing_precision << std::endl;

	f << "nb_pdf\t" << nb_pdf << std::endl;
	f << "ntheta_pdf\t" << ntheta_pdf << std::endl;

	f << "initial_saving\t" << initial_saving << std::endl;
	f << "period\t" << period << std::endl;

	f << "nmkt\t" << nmkt << std::endl;
	f << "f_w\t" << f_w << std::endl;

	f << "h\t" << h << std::endl;

	f << "period_trans\t" << period_trans << std::endl;
	f << "period_sim\t" << period_sim << std::endl;

	f << "nseg\t" << nseg << std::endl;

	f.close();
}

void write_parameters_matlab(std::string filename)
{
	std::ofstream f;
	//f.open (filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
	f.open(filename.c_str(), std::ios::out | std::ios::trunc);

	f << "beta," << beta << ",";
	f << "sigma," << sigma << ",";

	f << "shape," << shape << ",";
	f << "scale," << scale << ",";
	f << "cutoff," << cutoff << ",";
	f << "grammar," << grammar << ",";

	f << "time_preference_rate," << time_preference_rate << ",";

	f << "delta," << delta << ",";
	f << "alpha," << alpha << ",";
	f << "nu," << nu << ",";

	f << "ksi," << ksi << ",";
	f << "chi," << chi << ",";
	f << "conspsi," << conspsi << ",";
	f << "thetapsi," << thetapsi << ",";
	f << "conszeta," << conszeta << ",";
	f << "thetazeta," << thetazeta << ",";

	f << "kappa_ini," << kappa_ini << ",";
	f << "eta_ini," << eta_ini << ",";
	f << "kappa_new," << kappa_new << ",";
	f << "eta_new," << eta_new << ",";

	f << "nb," << nb << ",";
	f << "ntheta," << ntheta << ",";
	f << "nz," << nz << ",";

	f << "minb," << minb << ",";
	f << "maxb," << maxb << ",";
	f << "mintheta," << mintheta << ",";
	f << "maxtheta," << maxtheta << ",";

	f << "maxit_DP," << maxit_DP << ",";
	f << "maxit_gss," << maxit_gss << ",";
	f << "tol_itvfun," << tol_itvfun << ",";
	f << "tol_gss," << tol_gss << ",";
	f << "golden_ratio," << golden_ratio << ",";

	f << "coefficient," << coefficient << ",";

	f << "writing_precision," << writing_precision << ",";

	f << "nb_pdf," << nb_pdf << ",";
	f << "ntheta_pdf," << ntheta_pdf << ",";

	f << "initial_saving," << initial_saving << ",";
	f << "period," << period << ",";

	f << "nmkt," << nmkt << ",";
	f << "f_w," << f_w << ",";

	f << "h," << h << ",";
	
	f << "period_trans," << period_trans << ",";
	f << "period_sim," << period_sim << ",";

	f << "nseg," << nseg;

	f.close();
}

void write_progress(std::string filename, int i_r, int i_w) //write grids to file
{
	std::ofstream f;
	//f.open (filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
	f.open(filename.c_str(), std::ios::out | std::ios::trunc);
	f << "i_r : " << i_r << std::endl << "i_w : " << i_w << endl;

	f.close();
}

void read_int_dat1(int* a, int length, std::string filename)
{
	std::ifstream f;
	//f.open (filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
	f.open(filename.c_str(), std::ios::in);
	for (int i = 0; i < length; ++i)
	{

		f >> a[i];
	}
	f.close();
}

void read_int_dat2(int** a, int nx, int ny, std::string filename)
{
	std::ifstream f;
	//f.open (filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
	f.open(filename.c_str(), std::ios::in);
	for (int i = 0; i < nx; ++i)
	{
		for (int j = 0; j < ny; ++j)
		{

			f >> a[i][j];
		}
	}
	f.close();
}

void read_dat1(double* a, int length, std::string filename)
{
	std::ifstream f;
	//f.open (filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
	f.open(filename.c_str(), std::ios::in);
	for (int i = 0; i < length; ++i)
	{

		f >> a[i];
	}
	f.close();
}

void read_dat2(double** a, int nx, int ny, std::string filename)
{
	std::ifstream f;
	//f.open (filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
	f.open(filename.c_str(), std::ios::in);
	for (int i = 0; i < nx; ++i)
	{
		for (int j = 0; j < ny; ++j)
		{

			f >> a[i][j];
		}
	}
	f.close();
}

void read_dat3(double*** a, int nx, int ny, int nz, std::string filename)
{
	std::ifstream f;
	//f.open (filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
	f.open(filename.c_str(), std::ios::in);
	for (int i = 0; i < nx; ++i)
	{
		for (int j = 0; j < ny; ++j)
		{
			for (int k = 0; k < nz; ++k)
			{
				f >> a[i][j][k];
			}
		}
	}
	f.close();
}

void read_dat4(double**** a, int nx, int ny, int nz, int np, std::string filename)
{
	std::ifstream f;
	//f.open (filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
	f.open(filename.c_str(), std::ios::in);
	for (int i = 0; i < nx; ++i)
	{
		for (int j = 0; j < ny; ++j)
		{
			for (int k = 0; k < nz; ++k)
			{
				for (int m = 0; m < np; ++m)
				{
					f >> a[i][j][k][m];
				}
			}
		}
	}
	f.close();
}

void read_dat5(double***** a, int nx, int ny, int nz, int np, int nq, std::string filename)
{
	std::ifstream f;
	//f.open (filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
	f.open(filename.c_str(), std::ios::in);
	for (int i = 0; i < nx; ++i)
	{
		for (int j = 0; j < ny; ++j)
		{
			for (int k = 0; k < nz; ++k)
			{
				for (int m = 0; m < np; ++m)
				{
					for (int n = 0; n < nq; ++n)
					{
						f >> a[i][j][k][m][n];
					}
				}
			}
		}
	}
	f.close();
}
