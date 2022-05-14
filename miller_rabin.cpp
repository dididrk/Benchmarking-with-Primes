#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <chrono>
#include <omp.h>

std::mt19937 rng(time(0));
std::uniform_real_distribution<double> dist;

// ------------------------------------------------------------------------------------
unsigned int num_threads = 1;
unsigned int size_initial_sieve = 10;
unsigned int size_sufficient_witnesses = 7;
unsigned long long initial_sieve[10] = {0};
unsigned long long sufficient_witnesses[7] = {2, 3, 5, 7, 11, 13, 17}; // for n < 341,550,071,728,321
// 7,999,252,175,582,851

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

unsigned long long miller_rabin(unsigned long long n, unsigned long long confidence) {

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
  int k_starter = 10;
  for (int j = 0; j < size_initial_sieve; ++j) {
      for (int k = k_starter; k < 1000; ++k) {
          if (miller_rabin(k, 5) == 1) {
              k_starter = k + 1;
              initial_sieve[j] = k;
              break;
          }
      }
  }
}

unsigned int how_many_primes(unsigned long long n_min_, unsigned long long n_max, unsigned long long confidence) {

  unsigned long long n_min = n_min_;
  if (n_min_ < 7)
      n_min = 7;

  unsigned long long count = 0;
  unsigned long long s_min = n_min/6;
  unsigned long long s_max = n_max/6;

  if (6*s_min - 1 >= n_min) 
  {
    if (miller_rabin(6*s_min - 1, confidence) == 1)
      count += 1;
  }
  if (6*s_min + 1 >= n_min)
  {
    if (miller_rabin(6*s_min + 1, confidence) == 1)
      count += 1;
  }

  if (6*s_max + 1 < n_max) 
  {
    if (miller_rabin(6*s_max + 1, confidence) == 1)
      count += 1;
  }
  if (6*s_max - 1 < n_max)
  {
    if (miller_rabin(6*s_max - 1, confidence) == 1)
      count += 1;
  }

  #pragma omp parallel num_threads(num_threads)
  #pragma omp for reduction(+:count)
  for (unsigned long long i = s_min + 1; i < s_max; ++i){
      unsigned long long number = 6*i - 1;
      
      for (int j = 0; j < size_initial_sieve; ++j) {
          if (number%initial_sieve[j] == 0)
              goto jump2;    
      }

      if (miller_rabin(number, confidence) == 1)
      {
        count += 1;
      }
      jump2: 1;
  }

  #pragma omp parallel num_threads(num_threads)
  #pragma omp for reduction(+:count)
  for (unsigned long long i = s_min + 1; i < s_max; ++i){
      unsigned long long number = 6*i + 1;
  
      for (int j = 0; j < size_initial_sieve; ++j) {
          if (number%initial_sieve[j] == 0)
              goto jump3;
      }

      if (miller_rabin(number, confidence) == 1)
      {
        count += 1;
      }
      jump3: 1;
  }

  if (n_min <= 7)
    return count + size_sufficient_witnesses + size_initial_sieve;

  else 
    return count;
}

void continuous_benchmark()
{
  unsigned long long n_low = 1;
  unsigned long long n_high = 1000000;
  unsigned long long count_total = 0;
  unsigned long long number_of_rounds = 1;
  unsigned int n_averages = 1;
  unsigned int number_of_averages = 100;
  double moving_average[number_of_averages] = {0};
  double moving_variance[number_of_averages] = {0};
  double average_speed = 0;
  double std_speed = 0;

  while(1)
  {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    int count = how_many_primes(n_low, n_high, 5);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    long elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();

    count_total += count;
    float speed = static_cast<float>(count)/(static_cast<float>(elapsed_time)*0.000000001)/1000000.0;

    if (number_of_rounds > number_of_averages)
    {
      average_speed -= moving_average[(number_of_rounds - 1)%number_of_averages];
      std_speed -= moving_variance[(number_of_rounds - 1)%number_of_averages];
    }

    average_speed += speed;
    std_speed += speed*speed;
    moving_average[(number_of_rounds - 1)%number_of_averages] = speed;
    moving_variance[(number_of_rounds - 1)%number_of_averages] = speed*speed;
    if (n_averages > number_of_averages)
    {
      n_averages = number_of_averages;
    }

    std::cout << "\033[2J\033[1;1H";
    std::cout << "Primes between 1 and " << static_cast<unsigned long long>(n_high/1000000) << " Million ---- " << count_total << "\n";
    //std::cout << "Instantaneous speed" << std::setprecision(3) << std::fixed << speed << " Million Primes/s" << "\n";
    std::cout << "Speed " << std::setprecision(3) << std::fixed << average_speed/n_averages << " +- " << sqrt((std_speed - pow(average_speed, 2)/n_averages)/n_averages) <<" Million Primes/s" << "\n";
    
    if (n_high >= 341550071728321)
      break;

    n_low = n_high;
    n_high += 1000000;
    number_of_rounds += 1;
    n_averages += 1;
  }
}

int main(int argc, char** argv) {

  std::string nthreads;
  std::cout << "Multithreaded/Single Thread/Specific Threads Test (m/s/n): ";
  std::cin >> nthreads;

  if (nthreads == "m")
    num_threads = omp_get_max_threads();

  if (nthreads == "s")
    num_threads = 1;

  if (nthreads != "m" and nthreads != "s") 
  {
    int nthreadss;
    std::cout << "\nSpecific Number of Threads: ";
    std::cin >> nthreadss;
    num_threads = nthreadss;
  }

  set_eratostenes_sieve();
  continuous_benchmark();
}
