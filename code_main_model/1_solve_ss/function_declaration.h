#include <string>

#ifndef FUNCTION_DECLARATION_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define FUNCTION_DECLARATION_H

double absolute(double a);
double* absolute_vec(int length, double* a);
double* addition(int length, double* a, double* b);
double* subtraction(int length, double* a, double* b);
double max(int length, double* a);
double max2(int dim1, int dim2, double** a);
double min(int length, double* a);
double sum(int dim1, double* a);
double sum2(int dim1, int dim2, double** a);
double max_2num(double a, double b);
int min2(int a, int b);

bool*** ini_bool_matrix3(int dim1, int dim2, int dim3);
int* ini_int_matrix1(int dim1);
int** ini_int_matrix2(int dim1, int dim2);
int*** ini_int_matrix3(int dim1, int dim2, int dim3);
int**** ini_int_matrix4(int dim1, int dim2, int dim3, int dim4);
int***** ini_int_matrix5(int dim1, int dim2, int dim3, int dim4, int dim5);
double* ini_matrix1(int dim1);
double** ini_matrix2(int dim1, int dim2);
double*** ini_matrix3(int dim1, int dim2, int dim3);
double**** ini_matrix4(int dim1, int dim2, int dim3, int dim4);
double***** ini_matrix5(int dim1, int dim2, int dim3, int dim4, int dim5);
double****** ini_matrix6(int dim1, int dim2, int dim3, int dim4, int dim5, int dim6);
double******* ini_matrix7(int dim1, int dim2, int dim3, int dim4, int dim5, int dim6, int dim7);
double******** ini_matrix8(int dim1, int dim2, int dim3, int dim4, int dim5, int dim6, int dim7, int dim8);

double dot(int length, double* a, double* b);
int sear(int length, double* a, double width, double value);
int sear2(int length, double* a, double value);
double** transpose(int rownum, int colnum, double** m);
void linspace(double* grid, double min, double max, int length);
double lin_interpo_1(int nx, double* gridx, double* v, double x);
double lin_interpo_2_most_robust(int nx, int ny, double* gridx, double* gridy, double** vfun, double x, double y);
double lin_interpo_2(int nx, int ny, double* gridx, double* gridy, double** vfun, double x, double y);

void z_generator(double* gridz, double* probz, double** trans_z);

void read_int_dat1(int* a, int length, std::string filename);
void read_dat1(double* a, int length, std::string filename);
void read_dat2(double** a, int nx, int ny, std::string filename);
void read_dat3(double*** a, int nx, int ny, int nz, std::string filename);
void read_dat4(double**** a, int nx, int ny, int nz, int np, std::string filename);
void read_dat5(double***** a, int nx, int ny, int nz, int np, int nq, std::string filename);
void write3_bool(std::string filename, bool ***v, int dim1, int dim2, int dim3); 
void write8(std::string filename, double******** v, int dim1, int dim2, int dim3, int dim4, int dim5, int dim6, int dim7, int dim8);
void write7(std::string filename, double******* v, int dim1, int dim2, int dim3, int dim4, int dim5, int dim6, int dim7);
void write6(std::string filename, double ******v, int dim1, int dim2, int dim3, int dim4, int dim5, int dim6);
void write5(std::string filename, double *****v, int dim1, int dim2, int dim3, int dim4, int dim5);
void write4(std::string filename, double ****v, int dim1, int dim2, int dim3, int dim4);
void write3(std::string filename, double ***v, int dim1, int dim2, int dim3);
void write2(std::string filename, double** v, int dim1, int dim2);
void write1(std::string filename, double* v, int dim1);
void write_int3(std::string filename, int*** v, int dim1, int dim2, int dim3);
void write_parameters(std::string filename);
void write_parameters_matlab(std::string filename);
void write_progress(std::string filename, int i_r, int i_w);

void clear_matrix5(double***** matrix, int dim1, int dim2, int dim3, int dim4, int dim5);
void release_bool_matrix3(bool*** matrix, int dim1, int dim2, int dim3);
void release_int_matrix3(int*** matrix, int dim1, int dim2, int dim3);
void release_matrix5(double***** matrix, int dim1, int dim2, int dim3, int dim4, int dim5);
void release_matrix4(double**** matrix, int dim1, int dim2, int dim3, int dim4);
void release_matrix3(double*** matrix, int dim1, int dim2, int dim3);
void release_matrix2(double** matrix, int dim1, int dim2);
void release_matrix1(double* matrix, int dim1);

void optimal_decision(double p,
	double*** adjust_ad, double*** wealth_ad, double*** gk_ad, double*** gl_ad, 
					  double*** adjust_noad, double*** wealth_noad, double*** gk_noad, double*** gl_noad,
					  double* gridb, double* gridtheta, double* gridz, double* probz, double** trans_z, 
					  double w, double r, double psi, double zeta);
int rank3(double t_w, double t_e_nobo, double t_e_bo);


void golden_section_search_ad(int it, int i_mkt, int state_z, int state_b, int state_theta, double* theta_candidate, double*** V, double*** gc, double*** gtheta, double*** gb, double**** V_prv, double*** wealth,
						   double* gridb, double* gridtheta, double* gridz, double* probz, double** trans_z);
void golden_section_search_noad(int it, int i_mkt, int state_z, int state_b, int state_theta, double*** V, double*** gc, double*** gtheta, double*** gb, double**** V_prv, double*** wealth,
						   double* gridb, double* gridtheta, double* gridz, double* probz, double** trans_z, double r);

double value_ad(int i_mkt, int state_b, int state_z, double wealth, double consumption, double theta, double**** V_prv, double* gridb, double* gridtheta, double** trans_z);
double value_noad(int i_mkt, int state_b, int state_z, double wealth, double deposit, double consumption, double**** V_prv, double* gridb, double* gridtheta, double** trans_z);

void simulation_permanent(double***** mig,
	double***** pdf_A, double***** pdf_B,
	double**** V, double**** gtheta, double**** gb, double* gridb, double* gridtheta,
	double* gridb_pdf, double* gridtheta_pdf, double* gridz_pdf, double* probz_pdf, double** trans_z_pdf,
	int*** b_AtoB_M_index, int*** theta_AtoB_M_index, double***** V_AtoB, double***** step_AtoB_prob, int**** b_BtoA_next_index, int**** theta_BtoA_next_index);

void error_simulation(double** sum, double** error, double***** pdf);

void capital_demand_supply(double*** pdf_ss, double*** gk, double*** adjust, double* gridb, double* gridtheta,
				   double* gridb_pdf, double* gridtheta_pdf,
					double* capital_demand, double* capital_supply, double psi,
					double*** temp_cap_supply);

void labor_demand_supply(double*** pdf_ss, double*** gl, double*** adjust, double* gridb, double* gridtheta, 
				   double* gridb_pdf, double* gridtheta_pdf,
					double* labor_demand, double* labor_supply, 
					double*** l);

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

	double** diff, double* theta_candidate, int*** b_AtoB_M_index, int*** theta_AtoB_M_index, double***** V_AtoB, double***** step_AtoB_prob, int**** b_BtoA_next_index, int**** theta_BtoA_next_index, double*** temp_fin, double*** temp_cap_supply, double*** l, std::string currentpath);
void GDP(double p, 
	double*** pdf_ss, double*** gk, double*** gl, double*** adjust, double* gridb, double* gridtheta,
				   double* gridb_pdf, double* gridtheta_pdf, double* gridz_pdf, double* probz_pdf, double** trans_z_pdf,
					double* GDP_output, double psi,
					double*** temp_fin);

double itvfun_V(int it, int i_mkt, 
				double* theta_candidate, 
				double*** V_ad, double*** gc_ad, double*** gtheta_ad, double*** gb_ad, double*** gk_ad, double*** gl_ad, double*** adjust_ad, double*** wealth_ad,
				double*** V_noad, double*** gc_noad, double*** gtheta_noad, double*** gb_noad, double*** gk_noad, double*** gl_noad, double*** adjust_noad, double*** wealth_noad,
				double*** V_next, double**** V_prv, double*** gc, double*** gtheta, double*** gb, double*** gk, double*** gl, double*** adjust, double*** wealth,
				double* gridb, double* gridtheta, double* gridz, double* probz, double** trans_z, double r);

#endif