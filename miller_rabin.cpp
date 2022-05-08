#include <iostream>
#include <fstream>
#include <math.h>
#include <ctime>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <random>
#include <algorithm>
#include <numeric>
#include <vector>
#include <chrono>
#include <omp.h>

std::mt19937 rng(time(0));
std::uniform_real_distribution<double> dist;

// ------------------------------------------------------------------------------------
unsigned long long fastExp(unsigned long long b, unsigned long long e, unsigned long long m)
{
	unsigned long long result = 1;
	if (1 & e)
		result = b;
	while (1) {
		if (!e) break;
		e >>= 1;
		b = (b * b) % m;
		if (e & 1)
			result = (result * b) % m;
	}
	return result;
}

int miller_rabin(int n, int confidence) {

    unsigned long long q = 1;
    unsigned long long  k = 0;
    unsigned long long f = n - 1;
    while (1) {
        if (f % 2 == 0) {
            k += 1;
            f = f >> 1;
        }
        else {
            q = (n - 1)/(1 << k);
            break;
        }
    }

    for (int i = 0; i < confidence; ++i) {

        a = static_cast<int>((n - 2)*dist(rng) + 1);
        unsigned long long fastexp_aq = fastExp(a, q, n);
        if (fastexp_aq == 1) {
            continue;
        }
        flag_outer_loop = 0;
        for (int j = 0; j < k; ++j) {
            if (fastExp(fastexp_aq, (1 << j), n) == n - 1) {
                flag_outer_loop = 1;
                break;
            }
        }
        if (flag_outer_loop == 1) {
            continue;
        }

        return 0;
    }

    return 1;
}


int how_many_primes(int n_min_, int n_max, int confidence) {

    int n_min = n_min_;
    if (n_min_ < 7) {
        n_min = 7;
    }

    int size_sieve = 22;
    int initial_sieve[size_sieve] = {11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 67, 71, 73, 79, 83, 89, 97};
    int count_array[omp_get_max_threads()] = {0};
    int count = 2 + size_sieve;
    int s_min = static_cast<int>(static_cast<float>(n_min)/6);
    int s_max = static_cast<int>(static_cast<float>(n_max)/6);

    if (miller_rabin(6*s_min + 1, confidence) == 1) {
        count += 1;
    }

    if (miller_rabin(6*s_max - 1, confidence) == 1) {
        count += 1;
    }

    if (6*s_min - 1 >= n_min) {
        if (miller_rabin(6*s_min - 1, confidence) == 1) {
        count += 1;
        }
    }

    if (6*s_max + 1 <= n_max) {
        if (miller_rabin(6*s_max + 1, confidence) == 1) {
        count += 1;
        }
    }

    #pragma omp parallel num_threads(omp_get_max_threads())
    #pragma omp for
    for (int i = s_min + 1; i < s_max; ++i){

        int flag_loop = 0;
        for (int j = 0; j < size_sieve; ++j) {
            if ((6*i - 1)%initial_sieve[j] == 0) {
                flag_loop = 1;
                break;
            }
        }
        if (flag_loop == 1) {
            continue;
        }

        if (miller_rabin(6*i - 1, confidence) == 1) {
            count_array[omp_get_thread_num()] += 1;
        }
    }
    
    #pragma omp parallel num_threads(omp_get_max_threads())
    #pragma omp for
    for (int i = s_min + 1; i < s_max; ++i){
        
        int flag_loop = 0;
        for (int j = 0; j < size_sieve; ++j) {
            if ((6*i + 1)%initial_sieve[j] == 0) {
                flag_loop = 1;
                break;
            }
        }
        if (flag_loop == 1) {
            continue;
        }


        if (miller_rabin(6*i + 1, confidence) == 1) {
            count_array[omp_get_thread_num()] += 1;
        }
    }

    #pragma critical
    for (int j = 0; j < omp_get_max_threads(); j++) {
        count += count_array[j];
    }

    return count;
}

int main(int argc, char** argv) {

    int n_min = atoi(argv[1]);
    int n_max = atoi(argv[2]);
    int confidence = atoi(argv[3]);

    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "Number of primes found in the interval " << "[" << n_min << ", " << n_max << "[" << " - " << how_many_primes(n_min, n_max, confidence) << "\n";
    auto stop = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << "Took " << elapsed_time.count() << " milliseconds" << "\n"; 
    std::cout << "Average of " << n_max/static_cast<float>(elapsed_time.count()*0.001)/1000000 << " Million Numbers/s";
    return 1;

}
