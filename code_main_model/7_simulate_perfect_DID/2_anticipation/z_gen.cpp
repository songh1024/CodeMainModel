#include "parameters.h" 

void z_generator(double* gridz, double* probz, double** trans_z)
{
	double minz=scale;
	double maxz=scale/pow((1-cutoff),(1.0/shape));
	double* gridz_temp=ini_matrix1(nz+1);
	linspace(gridz_temp, minz, maxz, (nz+1));

	/*double sum=0.0;
	for (int i=0; i<nz; ++i)
	{
		gridz[i]=(gridz_temp[i]+gridz_temp[i+1])/2.0;
		probz[i]=(1-pow((scale/gridz_temp[i+1]),shape))-(1-pow((scale/gridz_temp[i]),shape));
		sum+=probz[i];
	}
	for (int i=0; i<nz; ++i)
	{
		probz[i]/=sum;
	}

	release_matrix1(gridz_temp, nz+1);
	for (int i=0; i<nz; ++i)
		for (int j=0; j<nz; ++j)
			trans_z[i][j]=((i==j)*(1-grammar)+probz[j]*grammar);*/

	gridz[nz-1]=maxz;
	for (int i=0; i<nz; i++)
		gridz[i]=1+(i+1.0)/nz*(gridz[nz-1]-1);
	double CP_z[nz];
	for (int i=0; i<nz; i++)
		CP_z[i]=1-pow(gridz[i], -shape);
	probz[0]=CP_z[0];
	for (int i=1; i<nz; i++)
		probz[i]=CP_z[i]-CP_z[i-1];
	double total_prob_z=0;
	for (int i=0; i<nz; i++)
		total_prob_z+=probz[i];
	for (int i=0; i<nz; i++)
		probz[i]=probz[i]/total_prob_z;
	for (int i=0; i<nz; ++i)
		for (int j=0; j<nz; ++j)
			trans_z[i][j]=((i==j)*(1-grammar)+probz[j]*grammar);

	release_matrix1(gridz_temp, nz+1);
	
	
}
	
