#include "parameters.h" 
using namespace std;
double itvfun_V(int it, int i_mkt,
	double* theta_candidate,
	double*** V_ad, double*** gc_ad, double*** gtheta_ad, double*** gb_ad, double*** gk_ad, double*** gl_ad, double*** adjust_ad, double*** wealth_ad,
	double*** V_noad, double*** gc_noad, double*** gtheta_noad, double*** gb_noad, double*** gk_noad, double*** gl_noad, double*** adjust_noad, double*** wealth_noad,
	double*** V_next, double**** V_prv, double*** gc, double*** gtheta, double*** gb, double*** gk, double*** gl, double*** adjust, double*** wealth,
	double* gridb, double* gridtheta, double* gridz, double* probz, double** trans_z, double r)
{
	#pragma omp parallel for schedule(dynamic, 1)
	for (int i_b=0; i_b<nb; ++i_b)
		for (int i_z=0; i_z<nz; ++i_z)
		{
			for (int i_theta=0; i_theta<ntheta; ++i_theta)
			{
				//if ((it == 1 && i_b == 0 && i_z == 0 && i_theta == 0))
					//cout << endl;
				golden_section_search_ad(it, i_mkt, i_z, i_b, i_theta, theta_candidate, V_ad, gc_ad, gtheta_ad, gb_ad, V_prv, wealth_ad, gridb, gridtheta, gridz, probz, trans_z);
				golden_section_search_noad(it, i_mkt, i_z, i_b, i_theta, V_noad, gc_noad, gtheta_noad, gb_noad, V_prv, wealth_noad, gridb, gridtheta, gridz, probz, trans_z, r);
				
				if (V_ad[i_z][i_b][i_theta] >= V_noad[i_z][i_b][i_theta])
				{
					V_next[i_z][i_b][i_theta] = V_ad[i_z][i_b][i_theta];
					gc[i_z][i_b][i_theta] = gc_ad[i_z][i_b][i_theta];
					gtheta[i_z][i_b][i_theta] = gtheta_ad[i_z][i_b][i_theta];
					gb[i_z][i_b][i_theta] = gb_ad[i_z][i_b][i_theta];
					gk[i_z][i_b][i_theta] = gk_ad[i_z][i_b][i_theta];
					gl[i_z][i_b][i_theta] = gl_ad[i_z][i_b][i_theta];
					wealth[i_z][i_b][i_theta] = wealth_ad[i_z][i_b][i_theta];
					adjust[i_z][i_b][i_theta] = (2*adjust_ad[i_z][i_b][i_theta]);
				}
				else
				{
					V_next[i_z][i_b][i_theta]=V_noad[i_z][i_b][i_theta];
					gc[i_z][i_b][i_theta]=gc_noad[i_z][i_b][i_theta];
					gtheta[i_z][i_b][i_theta]=gtheta_noad[i_z][i_b][i_theta];
					gb[i_z][i_b][i_theta]=gb_noad[i_z][i_b][i_theta];
					gk[i_z][i_b][i_theta]=gk_noad[i_z][i_b][i_theta];
					gl[i_z][i_b][i_theta]=gl_noad[i_z][i_b][i_theta];
					wealth[i_z][i_b][i_theta]=wealth_noad[i_z][i_b][i_theta];
					adjust[i_z][i_b][i_theta]=(2*adjust_noad[i_z][i_b][i_theta]-1);
				}
			}

		}

	//evaluating convergence
	double diff=-1e6; 
	for (int i_b=0; i_b<nb; ++i_b)
		for (int i_z=0; i_z<nz; ++i_z)
			for (int i_theta=0; i_theta<ntheta; ++i_theta)
				diff=max_2num(diff, absolute(V_next[i_z][i_b][i_theta]-V_prv[i_mkt][i_z][i_b][i_theta])); 

	return diff;
}


void golden_section_search_ad(int it, int i_mkt, int state_z, int state_b, int state_theta, double* theta_candidate, double*** V, double*** gc, double*** gtheta, double*** gb, double**** V_prv, double*** wealth,
	double* gridb, double* gridtheta, double* gridz, double* probz, double** trans_z)
{
	double total_wealth=wealth[state_z][state_b][state_theta];

	double consumption=0.0;
	double clow=0.0;
	double chigh=total_wealth;
	double diffc=chigh-clow;

	double temp_clow=clow+(1-golden_ratio)*diffc;
	double temp_chigh=chigh-(1-golden_ratio)*diffc;

	double v_clow= -1e7;
	double v_chigh= -1e7;

	int itc=0;
	int flag_c=0;

	int thetaindex=0;
	double theta=0.0;

	double v_opt=-1e7;
	double c_opt=0.0;
	double theta_opt=0.0;

	for (int index=0; index<coefficient; ++index) 
	{
		theta=theta_candidate[index];
		clow=0.0;
		chigh=total_wealth;
		itc=0;
		diffc=chigh-clow;
		flag_c=0;
		temp_clow=clow+(1-golden_ratio)*diffc;
	    temp_chigh=chigh-(1-golden_ratio)*diffc;
		while (itc<maxit_gss && diffc>=tol_gss)
		{
			++itc;
			if (flag_c!=2)
			{
				consumption=temp_clow;
				v_clow=value_ad(i_mkt, state_b, state_z, total_wealth, consumption, theta, V_prv, gridb, gridtheta, trans_z);
			}

			if (flag_c!=1)
			{
				consumption=temp_chigh;
				v_chigh=value_ad(i_mkt, state_b, state_z, total_wealth, consumption, theta, V_prv, gridb, gridtheta, trans_z);
			}

			if (v_clow>v_chigh)
			{
				flag_c=1;
				chigh=temp_chigh;
				temp_chigh=temp_clow;
				v_chigh=v_clow;
				temp_clow=clow+(1-golden_ratio)*(chigh-clow);
				diffc=temp_chigh-temp_clow;
			}
			else
			{
				flag_c=2;
				clow=temp_clow;
				temp_clow=temp_chigh;
				v_clow=v_chigh;
				temp_chigh=chigh-(1-golden_ratio)*(chigh-clow);
				diffc=temp_chigh-temp_clow;
			}
			

		}
		if (v_clow>v_opt)
		{
			v_opt=v_clow; // The final v_clow is not necessary the value evaluated at clow, could be chigh
			c_opt=clow;
			theta_opt=theta;
		}

	}

	//store optimal policies and values
	V[state_z][state_b][state_theta] = v_opt;
	gc[state_z][state_b][state_theta] = c_opt;
	gtheta[state_z][state_b][state_theta] = theta_opt;
	gb[state_z][state_b][state_theta] = total_wealth - c_opt;
	

}

void golden_section_search_noad(int it, int i_mkt, int state_z, int state_b, int state_theta, double*** V, double*** gc, double*** gtheta, double*** gb, double**** V_prv, double*** wealth,
	double* gridb, double* gridtheta, double* gridz, double* probz, double** trans_z, double r)
{
	double total_wealth=wealth[state_z][state_b][state_theta];

	double consumption=0.0;
	double clow=0.0;
	double chigh=total_wealth;
	double diffc=chigh-clow;

	double temp_clow=clow+(1-golden_ratio)*diffc;
	double temp_chigh=chigh-(1-golden_ratio)*diffc;

	double v_clow= -1e7;
	double v_chigh= -1e7;

	int itc=0;
	int flag_c=0;

	double deposit=gridb[state_b]*(1-gridtheta[state_theta])*(1+r);

	while (itc<maxit_gss && diffc>=tol_gss)
	{
		++itc;
		if (flag_c!=2)
		{
			consumption=temp_clow;
			v_clow=value_noad(i_mkt, state_b, state_z, total_wealth, deposit, consumption, V_prv, gridb, gridtheta, trans_z);
		}

		if (flag_c!=1)
		{
			consumption=temp_chigh;
			v_chigh=value_noad(i_mkt, state_b, state_z, total_wealth, deposit, consumption, V_prv, gridb, gridtheta, trans_z);
		}

		if (v_clow>v_chigh)
		{
			flag_c=1;
			chigh=temp_chigh;
			temp_chigh=temp_clow;
			v_chigh=v_clow;
			temp_clow=clow+(1-golden_ratio)*(chigh-clow);
			diffc=temp_chigh-temp_clow;
		}
		else
		{
			flag_c=2;
			clow=temp_clow;
			temp_clow=temp_chigh;
			v_clow=v_chigh;
			temp_chigh=chigh-(1-golden_ratio)*(chigh-clow);
			diffc=temp_chigh-temp_clow;
		}

	}

	//store optimal policies and values
	consumption=clow;
	double m_next_NM=total_wealth-consumption;
	double b_next_NM=m_next_NM+deposit;
	V[state_z][state_b][state_theta]=v_clow;
	gc[state_z][state_b][state_theta]=consumption;
	gb[state_z][state_b][state_theta]=b_next_NM;
	gtheta[state_z][state_b][state_theta]=(b_next_NM<1e-10) ? 1.0 : (m_next_NM/b_next_NM);

	consumption=total_wealth;
	v_clow=value_noad(i_mkt, state_b, state_z, total_wealth, deposit, consumption, V_prv, gridb, gridtheta, trans_z);
	if (v_clow>V[state_z][state_b][state_theta])
	{
		V[state_z][state_b][state_theta]=v_clow;
		gc[state_z][state_b][state_theta]=consumption;
		gb[state_z][state_b][state_theta]=deposit;
		gtheta[state_z][state_b][state_theta]=(deposit<1e-10) ? 1.0 : 0.0;
	}
}

double value_ad(int i_mkt, int state_b, int state_z, double wealth, double consumption, double theta, double**** V_prv, double* gridb, double* gridtheta, double** trans_z)
{
	double value_next=0.0;
	double temp=0.0;
	double b_next_NM=wealth-consumption;
	double theta_next_NM=theta;
	double m_next_NM=b_next_NM*theta_next_NM;

	double b_next_M=wealth-consumption-kappa;
	double m_next_M=m_next_NM-kappa;
	double theta_next_M=(absolute(b_next_M)<1e-10) ? 1.0 : (m_next_M/b_next_M);

	if (m_next_NM<kappa)
	{
		for (int i_z_next=0; i_z_next<nz; ++i_z_next)
			value_next+=trans_z[state_z][i_z_next]*lin_interpo_2(nb, ntheta, gridb, gridtheta, V_prv[i_mkt][i_z_next], b_next_NM, theta_next_NM);
	}
	else
	{
		for (int i_z_next = 0; i_z_next < nz; ++i_z_next)
		{
			for (int j_mkt = 0; j_mkt < nmkt; ++j_mkt)
			{
				if (j_mkt == i_mkt)
					temp += exp(1 / eta * lin_interpo_2(nb, ntheta, gridb, gridtheta, V_prv[j_mkt][i_z_next], b_next_NM, theta_next_NM));
				else
					temp += exp(1 / eta * lin_interpo_2(nb, ntheta, gridb, gridtheta, V_prv[j_mkt][i_z_next], b_next_M, theta_next_M));
			}
			value_next += trans_z[state_z][i_z_next] * eta * log(temp);
			temp = 0.0;
		}
	}

	double value=pow(consumption, 1-sigma)/(1-sigma)+beta*value_next;

	return value;
}

double value_noad(int i_mkt, int state_b, int state_z, double wealth, double deposit, double consumption, double**** V_prv, double* gridb, double* gridtheta, double** trans_z)
{
	double value_next=0.0;
	double temp=0.0;
	double m_next_NM=wealth-consumption;
	double b_next_NM=m_next_NM+deposit;
	double theta_next_NM=(absolute(b_next_NM)<1e-10) ? 1.0 : (m_next_NM/b_next_NM);
	double m_next_M=wealth-consumption-kappa;
	double b_next_M=m_next_M+deposit;
	double theta_next_M=(absolute(b_next_M)<1e-10) ? 1.0 : (m_next_M/b_next_M);

	if (m_next_NM<kappa)
	{
		for (int i_z_next=0; i_z_next<nz; ++i_z_next)
			value_next+=trans_z[state_z][i_z_next]*lin_interpo_2(nb, ntheta, gridb, gridtheta, V_prv[i_mkt][i_z_next], b_next_NM, theta_next_NM);
	}
	else
	{
		for (int i_z_next = 0; i_z_next < nz; ++i_z_next)
		{
			for (int j_mkt = 0; j_mkt < nmkt; ++j_mkt)
			{
				if (j_mkt == i_mkt)
					temp += exp(1 / eta * lin_interpo_2(nb, ntheta, gridb, gridtheta, V_prv[j_mkt][i_z_next], b_next_NM, theta_next_NM));
				else
					temp += exp(1 / eta * lin_interpo_2(nb, ntheta, gridb, gridtheta, V_prv[j_mkt][i_z_next], b_next_M, theta_next_M));
			}
			value_next += trans_z[state_z][i_z_next] * eta * log(temp);
			temp = 0.0;
		}
	}

	double value=pow(consumption, 1-sigma)/(1-sigma)+beta*value_next;

	return value;
}
