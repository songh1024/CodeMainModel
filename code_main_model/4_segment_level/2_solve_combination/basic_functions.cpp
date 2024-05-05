#include "parameters.h" 

double absolute(double a)
{
	double b=0.0;
	b=(a>=0)?a:-a;
	return b;
}

double* absolute_vec(int length, double* a)
{
	double* b=new double[length];
	for(int i=0; i<length; ++i)
		b[i]=(a[i]>=0)?a[i]:-a[i];
	return b;
}


double* addition(int length, double* a, double* b)
{
	double* c=new double[length];
	for(int i=0; i<length; ++i)
		c[i]=a[i]+b[i];
	return c;
}

double dot(int length, double* a, double* b) // dot product
{
	double sum=0;
	for(int i=0; i<length; ++i)
	{
		sum+=a[i]*b[i];
	}

	return sum;

}
double max_2num(double a, double b)
{
	
	return (a>b)?a:b;
}

double min_2num(double a, double b)
{

	return (a > b) ? b : a;
}

int min2(int a, int b)
{
	
	return (a>b)?b:a;
}

double max(int length, double* a) // find the maximum value in a vector
{
	double opt;
	opt=a[0];
	for(int i=1; i<length; ++i)
	{
		opt=(a[i]>opt)?a[i]:opt;
	}
	return opt;
}

int maxindex(int length, double* a)
{
	int opt = 0;
	for (int i = 1; i < length; ++i)
	{
		opt = (a[i] > a[opt]) ? i : opt;
	}
	return opt;
}

double max2(int dim1, int dim2, double** a) // find the maximum value in a vector
{
	double opt;
	opt = a[0][0];
	for (int i = 0; i < dim1; ++i)
		for (int j=0; j < dim2; ++j)
		{
			opt = (a[i][j] > opt) ? a[i][j] : opt;
		}
	return opt;
}

double min(int length, double* a) // find the minimum value in a vector
{
	double opt;
	opt=a[0];
	for(int i=1; i<length; ++i)
	{
		opt=(a[i]<opt)?a[i]:opt;
	}

	return opt;
}

int sear(int length, double* a, double width, double value)
{

	if (value<=a[0])
		return 0;

	else if (value>=a[length-1])
		return (length-1);

	/*else
		return ((int)floor(value/width));*/

	/*for (int i=0; i<length; ++i)
	{
		if ((a[i]-value)<=width/2 && (a[i]-value)>=-width/2) // this is because C++ and MATLAB have different accuracies
			return i;
	}*/
	int index=(int)floor(value/width);
	if (absolute(a[index]-value)>absolute(a[index+1]-value))
		return index+1;
	else
		return index;

	//return 99999;

}

int sear2(int length, double* a, double value)
{

	if (value == a[length - 1]) // Note that this takes care of the case where value==a[length - 1]
		return (length - 2);

	else
	{
		for (int i = 0; i < length-1; ++i)
			if (value >= a[i] && value < a[i+1])
				return i;
	}

	//return 99999;

}

double* subtraction(int length, double* a, double* b)
{
	double* c=new double[length];
	for(int i=0; i<length; ++i)
		c[i]=a[i]-b[i];
	return c;
}

double sum2(int dim1, int dim2, double** a)
{
	double opt=0;
	for (int i = 0; i < dim1; ++i)
		for (int j = 0; j < dim2; ++j)
			opt += a[i][j];
		
	return opt;
}

double sum(int dim1, double* a)
{
	double opt=0;
	for(int i=0; i<dim1; ++i)
		opt+=a[i];
	
	return opt;
}

double sum_select(int dim1, std::vector<int> v, double* a)
{
	double opt = 0;
	for (int i = 0; i < dim1; ++i)
		opt += a[v[i]];

	return opt;
}


double** transpose(int rownum, int colnum, double** m)
{
	double** a=ini_matrix2(colnum, rownum);

	for(int i=0; i<colnum; ++i)
		for(int j=0; j<rownum; ++j)
			a[i][j]=m[j][i];

	return a;


}