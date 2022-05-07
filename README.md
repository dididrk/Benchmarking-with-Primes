# Benchmarking-with-Primes
This is a simple cpu benchmarking program which computes primes up to n using the Miller Rabin algorithm and sieve of Eratosthenes

It only needs 3 arguments in the terminal ./name.exe (or ./name.out in linux) arg1 arg2 arg3. 
arg1 is the initial number which the program starts counting.
arg2 the final number.
arg3 the confidence number loop, since Miller Rabin is a probabilistic algorithm the change of a false prime is 1 in 4^confidence. I suggest confidence = 4 if up to 100 million.
