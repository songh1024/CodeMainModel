#include <string>

#ifndef FUNCTION_DECLARATION_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define FUNCTION_DECLARATION_H

double absolute(double a);
double* absolute_vec(int length, double* a);
double* addition(int length, double* a, double* b);
double* subtraction(int length, double* a, double* b);
double max(int length, double* a);
int maxindex(int length, double* a);
double max2(int dim1, int dim2, double** a);
double min(int length, double* a);
double sum(int dim1, double* a);
double sum2(int dim1, int dim2, double** a);
double sum_select(int dim1, std::vector<int> v, double* a);
double max_2num(double a, double b);
double min_2num(double a, double b);
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
void read_int_dat2(int** a, int nx, int ny, std::string filename);
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
void write_int1(std::string filename, int* v, int dim1);
void write_int2(std::string filename, int** v, int dim1, int dim2);
void write_parameters(std::string filename);
void write_parameters_matlab(std::string filename);
void write_progress(std::string filename, int i_r, int i_w);
void write_vecint1(std::string filename, std::vector<int>& v, int dim1);

void clear_matrix5(double***** matrix, int dim1, int dim2, int dim3, int dim4, int dim5);
void release_bool_matrix3(bool*** matrix, int dim1, int dim2, int dim3);
void release_int_matrix3(int*** matrix, int dim1, int dim2, int dim3);
void release_int_matrix2(int** matrix, int dim1, int dim2);
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


void golden_section_search_ad(int seg_nmkt, std::vector<int> seg_mktidx,
							int it, int i_mkt, int state_z, int state_b, int state_theta, double* theta_candidate, double*** V, double*** gc, double*** gtheta, double*** gb, double**** V_prv, double*** wealth,
						   double* gridb, double* gridtheta, double* gridz, double* probz, double** trans_z);
void golden_section_search_noad(int seg_nmkt, std::vector<int> seg_mktidx, 
							int it, int i_mkt, int state_z, int state_b, int state_theta, double*** V, double*** gc, double*** gtheta, double*** gb, double**** V_prv, double*** wealth,
						   double* gridb, double* gridtheta, double* gridz, double* probz, double** trans_z, double r);

double value_ad(int seg_nmkt, std::vector<int> seg_mktidx, 
							int i_mkt, int state_b, int state_z, double wealth, double consumption, double theta, double**** V_prv, double* gridb, double* gridtheta, double** trans_z);
double value_noad(int seg_nmkt, std::vector<int> seg_mktidx,
							int i_mkt, int state_b, int state_z, double wealth, double deposit, double consumption, double**** V_prv, double* gridb, double* gridtheta, double** trans_z);

void capital_demand_supply(double*** pdf_ss, double*** gk, double*** adjust, double* gridb, double* gridtheta,
				   double* gridb_pdf, double* gridtheta_pdf,
					double* capital_demand, double* capital_supply, double psi,
					double*** temp_cap_supply);

void labor_demand_supply(double*** pdf_ss, double*** gl, double*** adjust, double* gridb, double* gridtheta, 
				   double* gridb_pdf, double* gridtheta_pdf,
					double* labor_demand, double* labor_supply, 
					double*** l);

void GDP(double p, 
	double*** pdf_ss, double*** gk, double*** gl, double*** adjust, double* gridb, double* gridtheta,
				   double* gridb_pdf, double* gridtheta_pdf, double* gridz_pdf, double* probz_pdf, double** trans_z_pdf,
					double* GDP_output, double psi,
					double*** temp_fin);

double itvfun_V(int seg_nmkt, std::vector<int> seg_mktidx, 
				int it, int i_mkt,
				double* theta_candidate, 
				double*** V_ad, double*** gc_ad, double*** gtheta_ad, double*** gb_ad, double*** gk_ad, double*** gl_ad, double*** adjust_ad, double*** wealth_ad,
				double*** V_noad, double*** gc_noad, double*** gtheta_noad, double*** gb_noad, double*** gk_noad, double*** gl_noad, double*** adjust_noad, double*** wealth_noad,
				double*** V_next, double**** V_prv, double*** gc, double*** gtheta, double*** gb, double*** gk, double*** gl, double*** adjust, double*** wealth,
				double* gridb, double* gridtheta, double* gridz, double* probz, double** trans_z, double r);

double backward(int seg_nmkt, std::vector<int> seg_mktidx, 
	int it, int i_mkt,
	double* theta_candidate,
	double*** V_ad, double*** gc_ad, double*** gtheta_ad, double*** gb_ad, double*** gk_ad, double*** gl_ad, double*** adjust_ad, double*** wealth_ad,
	double*** V_noad, double*** gc_noad, double*** gtheta_noad, double*** gb_noad, double*** gk_noad, double*** gl_noad, double*** adjust_noad, double*** wealth_noad,
	double*** V_next, double**** V_prv, double*** gc, double*** gtheta, double*** gb, double*** gk, double*** gl, double*** adjust, double*** wealth,
	double* gridb, double* gridtheta, double* gridz, double* probz, double** trans_z, double r);

void forward_simulation(int seg_nmkt, std::vector<int> seg_mktidx,
	double****** mig,
	double***** pdf_A, double***** pdf_B,
	double***** V, double***** gtheta, double***** gb, double* gridb, double* gridtheta, double* gridb_pdf, double* gridtheta_pdf, double* gridz_pdf, double* probz_pdf, double** trans_z_pdf,
	int*** b_AtoB_M_index, int*** theta_AtoB_M_index, double***** V_AtoB, double***** step_AtoB_prob, int**** b_BtoA_next_index, int**** theta_BtoA_next_index);

void cal_bank_profit(int i_mkt, double* bank_profit, double*** pdf_ss, double*** gk, double*** adjust,
	double* gridb, double* gridtheta, double* gridb_pdf, double* gridtheta_pdf,
	double psi,
	double*** loan_demand);

int num_combination(int n, int k);
void combination(std::vector<int> v, int** comb, int ncomb, int N, int K);

double compute_profit(std::string currentpath,
	int nnewbank, int* mkt_unsure_comb,
	int seg_nmkt, std::vector<int> seg_mktidx, int seg_nmkt_86, std::vector<int> seg_mktidx_86, int seg_nmkt_sure, std::vector<int> seg_mktidx_sure, int seg_nmkt_unsure, std::vector<int> seg_mktidx_unsure,
	double* Pi, double* Z, double** distance, double** w, double* r, double**** pdf_ini, double* d_ini,
	double**** V_ad_new, double**** gc_ad_new, double**** gtheta_ad_new, double**** gb_ad_new, double**** gk_ad_new, double**** gl_ad_new, double**** adjust_ad_new, double**** wealth_ad_new,
	double**** V_noad_new, double**** gc_noad_new, double**** gtheta_noad_new, double**** gb_noad_new, double**** gk_noad_new, double**** gl_noad_new, double**** noadjust_noad_new, double**** wealth_noad_new,
	double***** V_new, double**** gc_new, double**** gtheta_new, double**** gb_new, double**** gk_new, double**** gl_new, double**** adjust_new, double**** wealth_new,
	double***** V_ad_backward, double***** gc_ad_backward, double***** gtheta_ad_backward, double***** gb_ad_backward, double***** gk_ad_backward, double***** gl_ad_backward, double***** adjust_ad_backward, double***** wealth_ad_backward,
	double***** V_noad_backward, double***** gc_noad_backward, double***** gtheta_noad_backward, double***** gb_noad_backward, double***** gk_noad_backward, double***** gl_noad_backward, double***** noadjust_noad_backward, double***** wealth_noad_backward,
	double***** V_backward, double***** gc_backward, double***** gtheta_backward, double***** gb_backward, double***** gk_backward, double***** gl_backward, double***** adjust_backward, double***** wealth_backward,
	double****** mig_trans, double** sum_pdf_A_trans, double** error_pdf_A_trans, double** sum_pdf_B_trans, double** error_pdf_B_trans, double***** pdf_A_forward, double***** pdf_B_forward,
	double* gridz, double* probz, double** trans_z, double* gridb, double* gridtheta, double* gridb_pdf, double* gridtheta_pdf, 
	double* theta_candidate, int*** b_AtoB_M_index, int*** theta_AtoB_M_index, double***** V_AtoB, double***** step_AtoB_prob, int**** b_BtoA_next_index, int**** theta_BtoA_next_index, double** diff, double*** temp_cap_supply,
	int* access);

#endif