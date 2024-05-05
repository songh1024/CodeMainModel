#include "parameters.h" 

bool*** ini_bool_matrix3(int dim1, int dim2, int dim3)
{
	bool*** result=new bool** [dim1];
	for (int i=0; i<dim1; ++i)
		result[i]=new bool*[dim2];
	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			result[i][j]=new bool[dim3];

	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			for (int k=0; k<dim3; ++k)
				result[i][j][k]=0;
	return result;
}

int* ini_int_matrix1(int dim1)
{
	int* result = new int [dim1];
	
	return result;
}

int** ini_int_matrix2(int dim1, int dim2)
{
	int** result=new int* [dim1];
	for (int i=0; i<dim1; ++i)
		result[i]=new int[dim2];
	
	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			result[i][j]=0;
	return result;
}

int*** ini_int_matrix3(int dim1, int dim2, int dim3)
{
	int*** result=new int** [dim1];
	for (int i=0; i<dim1; ++i)
		result[i]=new int*[dim2];
	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			result[i][j]=new int[dim3];

	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			for (int k=0; k<dim3; ++k)
				result[i][j][k]=0;
	return result;
}

int**** ini_int_matrix4(int dim1, int dim2, int dim3, int dim4)
{
	int**** result = new int*** [dim1];

	for (int i = 0; i < dim1; ++i)
		result[i] = new int** [dim2];

	for (int i = 0; i < dim1; ++i)
		for (int j = 0; j < dim2; ++j)
			result[i][j] = new int* [dim3];

	for (int i = 0; i < dim1; ++i)
		for (int j = 0; j < dim2; ++j)
			for (int k = 0; k < dim3; ++k)
				result[i][j][k] = new int[dim4];

	for (int i = 0; i < dim1; ++i)
		for (int j = 0; j < dim2; ++j)
			for (int k = 0; k < dim3; ++k)
				for (int l = 0; l < dim4; ++l)
					result[i][j][k][l] = 0;
	return result;
}

int***** ini_int_matrix5(int dim1, int dim2, int dim3, int dim4, int dim5)
{
	int***** result = new int**** [dim1];

	for (int i = 0; i < dim1; ++i)
		result[i] = new int*** [dim2];

	for (int i = 0; i < dim1; ++i)
		for (int j = 0; j < dim2; ++j)
			result[i][j] = new int** [dim3];

	for (int i = 0; i < dim1; ++i)
		for (int j = 0; j < dim2; ++j)
			for (int k = 0; k < dim3; ++k)
				result[i][j][k] = new int* [dim4];

	for (int i = 0; i < dim1; ++i)
		for (int j = 0; j < dim2; ++j)
			for (int k = 0; k < dim3; ++k)
				for (int l = 0; l < dim4; ++l)
					result[i][j][k][l] = new int[dim5];

	for (int i = 0; i < dim1; ++i)
		for (int j = 0; j < dim2; ++j)
			for (int k = 0; k < dim3; ++k)
				for (int l = 0; l < dim4; ++l)
					for (int m = 0; m < dim5; ++m)
						result[i][j][k][l][m] = 0;

	return result;
}

double* ini_matrix1(int dim1)
{
	double* result=new double [dim1];
	for (int i=0; i<dim1; ++i)
		result[i]=0.0;

	return result;
}

double** ini_matrix2(int dim1, int dim2)
{
	double** result=new double* [dim1];
	for (int i=0; i<dim1; ++i)
		result[i]=new double[dim2];

	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			result[i][j]=0.0;
	return result;
}

				
double*** ini_matrix3(int dim1, int dim2, int dim3)
{
	double*** result=new double** [dim1];
	for (int i=0; i<dim1; ++i)
		result[i]=new double*[dim2];
	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			result[i][j]=new double[dim3];

	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			for (int k=0; k<dim3; ++k)
				result[i][j][k]=0.0;
	return result;
}

double**** ini_matrix4(int dim1, int dim2, int dim3, int dim4)
{
	double**** result=new double*** [dim1];

	for (int i=0; i<dim1; ++i)
		result[i]=new double**[dim2];

	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			result[i][j]=new double*[dim3];

	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			for (int k=0; k<dim3; ++k)
				result[i][j][k]=new double[dim4];

	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			for (int k=0; k<dim3; ++k)
				for (int l=0; l<dim4; ++l)
					result[i][j][k][l]=0.0;
	return result;
}

double***** ini_matrix5(int dim1, int dim2, int dim3, int dim4, int dim5)
{
	double***** result=new double**** [dim1];

	for (int i=0; i<dim1; ++i)
		result[i]=new double***[dim2];

	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			result[i][j]=new double**[dim3];

	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			for (int k=0; k<dim3; ++k)
				result[i][j][k]=new double*[dim4];

	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			for (int k=0; k<dim3; ++k)
				for (int l=0; l<dim4; ++l)
					result[i][j][k][l]=new double[dim5];

	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			for (int k=0; k<dim3; ++k)
				for (int l=0; l<dim4; ++l)
					for (int m=0; m<dim5; ++m)
						result[i][j][k][l][m]=0.0;

	return result;
}

double****** ini_matrix6(int dim1, int dim2, int dim3, int dim4, int dim5, int dim6)
{
	double****** result=new double***** [dim1];

	for (int i=0; i<dim1; ++i)
		result[i]=new double****[dim2];

	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			result[i][j]=new double***[dim3];

	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			for (int k=0; k<dim3; ++k)
				result[i][j][k]=new double**[dim4];

	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			for (int k=0; k<dim3; ++k)
				for (int l=0; l<dim4; ++l)
					result[i][j][k][l]=new double*[dim5];

	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			for (int k=0; k<dim3; ++k)
				for (int l=0; l<dim4; ++l)
					for (int m=0; m<dim5; ++m)
						result[i][j][k][l][m]=new double [dim6];

	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			for (int k=0; k<dim3; ++k)
				for (int l=0; l<dim4; ++l)
					for (int m=0; m<dim5; ++m)
						for (int n=0; n<dim6; ++n)
							result[i][j][k][l][m][n]=0.0;

	return result;
}

double******* ini_matrix7(int dim1, int dim2, int dim3, int dim4, int dim5, int dim6, int dim7)
{
	double******* result = new double****** [dim1];

	for (int i = 0; i < dim1; ++i)
		result[i] = new double***** [dim2];

	for (int i = 0; i < dim1; ++i)
		for (int j = 0; j < dim2; ++j)
			result[i][j] = new double**** [dim3];

	for (int i = 0; i < dim1; ++i)
		for (int j = 0; j < dim2; ++j)
			for (int k = 0; k < dim3; ++k)
				result[i][j][k] = new double*** [dim4];

	for (int i = 0; i < dim1; ++i)
		for (int j = 0; j < dim2; ++j)
			for (int k = 0; k < dim3; ++k)
				for (int l = 0; l < dim4; ++l)
					result[i][j][k][l] = new double** [dim5];

	for (int i = 0; i < dim1; ++i)
		for (int j = 0; j < dim2; ++j)
			for (int k = 0; k < dim3; ++k)
				for (int l = 0; l < dim4; ++l)
					for (int m = 0; m < dim5; ++m)
						result[i][j][k][l][m] = new double* [dim6];

	for (int i = 0; i < dim1; ++i)
		for (int j = 0; j < dim2; ++j)
			for (int k = 0; k < dim3; ++k)
				for (int l = 0; l < dim4; ++l)
					for (int m = 0; m < dim5; ++m)
						for (int n = 0; n < dim6; ++n)
							result[i][j][k][l][m][n] = new double [dim7];

	for (int i = 0; i < dim1; ++i)
		for (int j = 0; j < dim2; ++j)
			for (int k = 0; k < dim3; ++k)
				for (int l = 0; l < dim4; ++l)
					for (int m = 0; m < dim5; ++m)
						for (int n = 0; n < dim6; ++n)
							for (int o = 0; o < dim7; ++o)
								result[i][j][k][l][m][n][o] = 0.0;

	return result;
}

double******** ini_matrix8(int dim1, int dim2, int dim3, int dim4, int dim5, int dim6, int dim7, int dim8)
{
	double******** result = new double******* [dim1];

	for (int i = 0; i < dim1; ++i)
		result[i] = new double****** [dim2];

	for (int i = 0; i < dim1; ++i)
		for (int j = 0; j < dim2; ++j)
			result[i][j] = new double***** [dim3];

	for (int i = 0; i < dim1; ++i)
		for (int j = 0; j < dim2; ++j)
			for (int k = 0; k < dim3; ++k)
				result[i][j][k] = new double**** [dim4];

	for (int i = 0; i < dim1; ++i)
		for (int j = 0; j < dim2; ++j)
			for (int k = 0; k < dim3; ++k)
				for (int l = 0; l < dim4; ++l)
					result[i][j][k][l] = new double*** [dim5];

	for (int i = 0; i < dim1; ++i)
		for (int j = 0; j < dim2; ++j)
			for (int k = 0; k < dim3; ++k)
				for (int l = 0; l < dim4; ++l)
					for (int m = 0; m < dim5; ++m)
						result[i][j][k][l][m] = new double** [dim6];

	for (int i = 0; i < dim1; ++i)
		for (int j = 0; j < dim2; ++j)
			for (int k = 0; k < dim3; ++k)
				for (int l = 0; l < dim4; ++l)
					for (int m = 0; m < dim5; ++m)
						for (int n = 0; n < dim6; ++n)
							result[i][j][k][l][m][n] = new double* [dim7];

	for (int i = 0; i < dim1; ++i)
		for (int j = 0; j < dim2; ++j)
			for (int k = 0; k < dim3; ++k)
				for (int l = 0; l < dim4; ++l)
					for (int m = 0; m < dim5; ++m)
						for (int n = 0; n < dim6; ++n)
							for (int o = 0; o < dim7; ++o)
								result[i][j][k][l][m][n][o] = new double [dim8];

	for (int i = 0; i < dim1; ++i)
		for (int j = 0; j < dim2; ++j)
			for (int k = 0; k < dim3; ++k)
				for (int l = 0; l < dim4; ++l)
					for (int m = 0; m < dim5; ++m)
						for (int n = 0; n < dim6; ++n)
							for (int o = 0; o < dim7; ++o)
								for (int p = 0; p < dim8; ++p)
									result[i][j][k][l][m][n][o][p] = 0.0;

	return result;
}


void release_bool_matrix3(bool*** matrix, int dim1, int dim2, int dim3)
{
	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			delete[] matrix[i][j];

	for (int i=0; i<dim1; ++i)
		delete[] matrix[i];

	delete[] matrix;
}

void release_int_matrix3(int*** matrix, int dim1, int dim2, int dim3)
{
	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			delete[] matrix[i][j];

	for (int i=0; i<dim1; ++i)
		delete[] matrix[i];

	delete[] matrix;
}

void release_int_matrix2(int** matrix, int dim1, int dim2)
{

	for (int i = 0; i < dim1; ++i)
		delete[] matrix[i];

	delete[] matrix;
}

void release_matrix5(double***** matrix, int dim1, int dim2, int dim3, int dim4, int dim5)
{
	
	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			for (int k=0; k<dim3; ++k)
				for (int l=0; l<dim4; ++l)
					delete[] matrix[i][j][k][l];

	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			for (int k=0; k<dim3; ++k)
				delete[] matrix[i][j][k];

	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			delete[] matrix[i][j];

	for (int i=0; i<dim1; ++i)
		delete[] matrix[i];

	delete[] matrix;

}

void release_matrix4(double**** matrix, int dim1, int dim2, int dim3, int dim4)
{
	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			for (int k=0; k<dim3; ++k)
				delete[] matrix[i][j][k];

	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			delete[] matrix[i][j];

	for (int i=0; i<dim1; ++i)
		delete[] matrix[i];

	delete[] matrix;

}

void release_matrix3(double*** matrix, int dim1, int dim2, int dim3)
{
	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			delete[] matrix[i][j];

	for (int i=0; i<dim1; ++i)
		delete[] matrix[i];

	delete[] matrix;
}

void release_matrix2(double** matrix, int dim1, int dim2)
{

	for (int i=0; i<dim1; ++i)
		delete[] matrix[i];

	delete[] matrix;
}

void release_matrix1(double* matrix, int dim1)
{

	delete[] matrix;
}

void clear_matrix5(double***** matrix, int dim1, int dim2, int dim3, int dim4, int dim5)
{
	for (int i=0; i<dim1; ++i)
		for (int j=0; j<dim2; ++j)
			for (int k=0; k<dim3; ++k)
				for (int m=0; m<dim4; ++m)
					for (int n=0; n<dim5; ++n)
						matrix[i][j][k][m][n]=0.0;;


}		
				