#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
 
int num_combination(int n, int k)
{
	int result=1;
	for (int i=n; i>n-k; --i)
		result*=i;
	for (int i=k; i>0; --i)
		result/=i;
	return result;
}

void combination(std::vector<int> v, int** comb, int ncomb, int N, int K)
{
    std::string bitmask(K, 1); // K leading 1's
    bitmask.resize(N, 0); // N-K trailing 0's
 
    // print integers and permute bitmask

	int i_comb=0;
    int i_cum = 0;
    do {
        i_cum = 0;
        for (int i = 0; i < N; ++i) // [0..N-1] integers
        {
            if (bitmask[i])
            {
                comb[i_comb][i_cum] = v[i];
                ++i_cum;
            }
        }
		++i_comb;
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
}
