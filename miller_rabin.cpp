#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <chrono>
#include <omp.h>
#include </home/darkness/Desktop/ClusterLibraries/eigen-3.4.0/Eigen/Dense>

std::mt19937 rng(time(0));
std::uniform_real_distribution<double> dist;

// ------------------------------------------------------------------------------------
unsigned int size_initial_sieve = 10;
unsigned int size_sufficient_witnesses = 3;
unsigned long long initial_sieve[10] = {0};
unsigned long long sufficient_witnesses[3] = {2, 7, 61}; // for n < 4,759,123,141
//unsigned long long count = 3 + size_initial_sieve;

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

unsigned long long miller_rabin(int n, int confidence) {

    unsigned long long q = 1;
    unsigned long long k = 0;
    unsigned long long f = n - 1;

    int k_pow = f & (~(f - 1));
    q = f/k_pow;

    while (k_pow) {
        k_pow = k_pow >> 1;
        k += 1;
    }

    k = k - 1;
    //std::cout << q << "\n";
    for (int i = 0; i < size_sufficient_witnesses; ++i) {

        //unsigned long long a = static_cast<int>((n - 2)*dist(rng) + 1);
        unsigned long long fastexp_aq = fastExp(sufficient_witnesses[i], q, n);

        if (fastexp_aq == 1 or fastexp_aq == n - 1)
            continue;

        unsigned long long fastexp_aq_j = fastexp_aq;
        for (int j = 1; j < k; ++j) {
            unsigned long long fastexp_aq_j_s = (fastexp_aq_j % n); 
            fastexp_aq_j = (fastexp_aq_j_s*fastexp_aq_j_s) % n;
            if (fastexp_aq_j == n - 1)
                goto jump1;
        }
        return 0;
        jump1: 1;
    }
    return 1;
}

void set_eratostenes_sieve() {
    int k_prime = 10;
    for (int j = 0; j < size_initial_sieve; ++j) {
        for (int k = k_prime + 1; k < 1000; ++k) {
            if (miller_rabin(k, 5) == 1) {
                k_prime = k;
                initial_sieve[j] = k_prime;
                break;
            }
        }
    }
}

int how_many_primes(int n_min_, int n_max, int confidence) {

    int n_min = n_min_;
    if (n_min_ < 7)
        n_min = 7;

    //unsigned long long count_array[omp_get_max_threads()] = {0};
    unsigned long long count = 0;
    unsigned long long s_min = static_cast<int>(static_cast<float>(n_min)/6);
    unsigned long long s_max = static_cast<int>(static_cast<float>(n_max)/6);

    if (miller_rabin(6*s_min + 1, confidence) == 1)
        count += 1;

    if (miller_rabin(6*s_max - 1, confidence) == 1)
        count += 1;

    if (6*s_min - 1 >= n_min) {
        if (miller_rabin(6*s_min - 1, confidence) == 1)
        count += 1;
    }

    if (6*s_max + 1 <= n_max) {
        if (miller_rabin(6*s_max + 1, confidence) == 1)
        count += 1;
    }

    #pragma omp parallel num_threads(omp_get_max_threads())
    #pragma omp for reduction(+:count)
    for (int i = s_min + 1; i < s_max; ++i){
        unsigned long long number = 6*i - 1;
        
        for (int j = 0; j < size_initial_sieve; ++j) {
            if (number%initial_sieve[j] == 0)
                goto jump2;    
        }

        if (miller_rabin(number, confidence) == 1)
            count += 1;
        jump2: 1;
    }
    
    #pragma omp parallel num_threads(omp_get_max_threads())
    #pragma omp for reduction(+:count)
    for (int i = s_min + 1; i < s_max; ++i){
        unsigned long long number = 6*i + 1;
    
        for (int j = 0; j < size_initial_sieve; ++j) {
            if (number%initial_sieve[j] == 0)
                goto jump3;
        }

        if (miller_rabin(number, confidence) == 1)
            count += 1;
        jump3: 1;
    }
    
    //for (int j = 0; j < omp_get_max_threads(); j++)
    //    count += count_array[j];

    return count + 5 + size_initial_sieve;
}

int main(int argc, char** argv) {

    std::ostringstream fn;
    std::ofstream myfile;
	fn << "results.txt";
	myfile.open(fn.str());

    int n_min = atoi(argv[1]);
    int n_max = atoi(argv[2]);
    int confidence = atoi(argv[3]);
    
    set_eratostenes_sieve();
    
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    int count = how_many_primes(n_min, n_max, confidence);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    long elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();

    std::cout << count << "\n";
    std::cout << static_cast<float>(n_max)/(static_cast<float>(elapsed_time)*0.000000001)/1000000.0 << "\n";

    //myfile << count << "\n";
    //myfile << static_cast<float>(n_max)/(static_cast<float>(elapsed_time)*0.000000001)/1000000.0 << "\n";

    myfile.close();
}
