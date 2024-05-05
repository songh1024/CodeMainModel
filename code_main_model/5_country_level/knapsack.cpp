#include "parameters.h" 
using namespace std;
void knapsack(int* n_opt, int nnewbank, int** nnewbank_seg, double** profit_seg)
{
	double** RecurF = ini_matrix2(nseg, nnewbank + 1);

	// Initialize RecurF[0]
	for (int i = 0; i < nnewbank + 1; ++i)
		RecurF[0][i] = -1e6;
	for (int i = 0; i < 2 * delta_l + 1; ++i)
	{
		if (nnewbank_seg[0][i] >= 0)
			RecurF[0][nnewbank_seg[0][i]] = profit_seg[0][i];
	}

	//cout << 0 << ": " << endl;
	//for (int i = 0; i < nnewbank + 1; ++i)
	//	cout << RecurF[0][i] << "  ";
	//cout << endl;


	/* Phi^k (x) = max_{b = 0, ..., x} {Theta^k (b) + Phi^{k-1} (x-b)} */
	int from_k, from_prev;
	double* temp = ini_matrix1(2 * delta_l + 1);
	for (int i_k = 1; i_k < nseg; ++i_k)
	{
		for (int i_x = 0; i_x < nnewbank + 1; ++i_x)
		{
			for (int i_b = 0; i_b < 2 * delta_l + 1; ++i_b) // Note that this i_b does not directly translates to the "b" in the equation, nnewbank_seg[i_k][i_b] does
			{
				from_k = nnewbank_seg[i_k][i_b];
				from_prev = i_x - from_k;

				if (from_k < 0 || from_prev < 0)
					temp[i_b] = -1e6;
				else
					temp[i_b] = profit_seg[i_k][i_b] + RecurF[i_k - 1][from_prev];
			}
			RecurF[i_k][i_x] = max(2 * delta_l + 1, temp);
		}
		//cout << i_k << ": " << endl;
		//for (int i = 0; i < nnewbank + 1; ++i)
		//	cout << RecurF[i_k][i] << "  ";
		//cout << endl;
	}

	int* n_cum_opt = ini_int_matrix1(nseg);
	n_cum_opt[nseg - 1] = nnewbank;
	double* profit_cum_opt = ini_matrix1(nseg);
	profit_cum_opt[nseg - 1] = RecurF[nseg - 1][nnewbank];

	for (int i_k = nseg - 1; i_k > 0; --i_k)
		for (int i_b = 0; i_b < 2 * delta_l + 1; ++i_b)
		{
			from_k = nnewbank_seg[i_k][i_b];
			from_prev = n_cum_opt[i_k] - from_k;
			if (from_k >= 0 && from_prev >= 0 && profit_cum_opt[i_k] == profit_seg[i_k][i_b] + RecurF[i_k - 1][from_prev])
			{
				n_opt[i_k] = from_k;
				n_cum_opt[i_k - 1] = n_cum_opt[i_k] - from_k;
				profit_cum_opt[i_k - 1] = RecurF[i_k - 1][from_prev];
				break;
			}
		}
	n_opt[0] = n_cum_opt[0];
	//for (int i_seg = 0; i_seg < nseg; ++i_seg)
	//	cout << n_opt[i_seg] << "  ";


	release_int_matrix1(n_cum_opt, nseg);
	release_matrix1(profit_cum_opt, nseg);
	release_matrix1(temp, 2 * delta_l + 1);
	release_matrix2(RecurF, nseg, 2 * delta_l + 1);
}

