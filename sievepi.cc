// sievepi.cc
//
// This program counts (and possibly prints) the prime numbers up to the
// integer command line argument using a blocked sieve of Eratosthenes

// It works up to 100,000,000,000 - I couldn't be bothered to wait beyond 
// that because it takes a few minutes. checked against powers of 10 in 
// https://mathworld.wolfram.com/PrimeCountingFunction.html

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

// tune this value based on your machine. 512K seems pretty fast

#ifndef BLOCK_SIZE
#define BLOCK_SIZE	(1<<20)
#endif

// "ENOUGH" array entries to memoize prime numbers. will eventually not
// be enough, not tested beyond 10^11 but should be good enough for pi(n)
// where n is within 64 bits

#define ENOUGH		65536 

// uncomment to actually print the prime numbers; takes a lot longer than
// generating them. why not make this a command line argument? because we
// want the compiler to optimize away the code in an inner loop that decides 
// whether to print

//#define PRINT_PRIMES

// print a uint64_t
void print (uint64_t x) {
#ifdef PRINT_PRIMES
	printf ("%lu", (unsigned long int) x);
#endif
}

// print a uint64_t and then a newline
void println (uint64_t x) {
#ifdef PRINT_PRIMES
	print (x);
	printf ("\n");
#endif
}

// A is a "block" of bits containing the sieve. it's divided into two halves: one for numbers that are 6 * i + 1 and the other for 6 * 1 + 5

bool A[2][BLOCK_SIZE];

// P2 is an array filled with the first ENOUGH prime numbers. they are
// 32-bit ints because we will never need more than that in the range of
// possible 64-bit primes. we keep two lists of primes, one for multiples of
// 6 * i + 1 and the other for multiples of 6 * i + 5. Z[i] counts how many
// entries are currently in the P[i] list of primes.

uint32_t P[2][ENOUGH], Z[2] = { 0, 0 };

// the pi function. counts primes up to v

uint64_t pi (uint64_t v) {
	// we will only sieve numbers that are multiples of 6 * i + 1 or 5

	const int C[2] = { 1, 5 };

	// useful variable names (I am old)

	uint64_t a, b, c, f, i, j, k, l, m, n, o, p, r, w, y;

	// start the count of primes at 0

	n = 0;

	// w is the number of entries we actually need to sieve, representing
	// numbers that are 6 * i + 1 or 5

	w = 1 + v / 6;

	// count and remember first few primes as special cases

	if (v >= 2) { println (2); n++; }
	if (v >= 3) { println (3); n++; }
	if (v >= 5) { n++; P[1][Z[1]++] = 5; println (5); }
	if (v >= 7) { n++; P[0][Z[0]++] = 7; println (7); }
	if (v >= 11) { n++; P[1][Z[1]++] = 11; println (11); }

	// sieve the first block of numbers

	// r is the end of the first block; either BLOCK_SIZE or w, whichever is smaller

	r = BLOCK_SIZE;
	if (r > w) r = w;

	// we don't need to check for factors beyond sqrt(r/6)

	l = 1 + (unsigned int) sqrt (r / 6);

	// reset all elements in the block

	memset (A, 1, sizeof (A));

	// sieve all of the multiples of 6 * i + 1 or 5 in the first block

	for (m=0; m<2; m++) {
		c = C[m];
		for (p=0; p<2; p++) {
			a = C[p];

			// mark off all multiples of 6 * i + c

			for (i=0; i<l; i++) {

				f = i * 6 + c;
				if (f == 1) continue;

				// find k, the place to start sieving, as the first multiple of i * 6 + a that isn't i

				for (j=f; ; j+=f) {
					k = (j - a) / 6;
					if (k * 6 + a == j)
						if (k != i) break;
				}

				// mark off the multiples of 6 * i + c

				for (j=k; j<r; j+=f) A[p][j] = 0;
			}
		}
	}

	// collect the results for the first block

	for (p=0; p<2; p++) {
		c = C[p];

		// figure out how high to go and adjust if necessary

		k = (v - c) / 6 + 1;
		if (r > k) r = k;

		// count all the 1s in A

		for (i=2; i<r; i++) {
			if (A[p][i]) {

				// if we haven't collected ENOUGH primes of the form 6 * i + c, insert this prime into the P array

				if (Z[p] < ENOUGH) P[p][Z[p]++] = i * 6 + c;

				// one more prime number; print it and count it

				println (i * 6 + c);
				n++;
			}
		}
	}

	// sieve all the remaining numbers block by block

	// b starts out at the end of the last section of numbers (the last "block") and increments by BLOCK_SIZE

	for (b=BLOCK_SIZE; b<=w; b+=BLOCK_SIZE) {

		// r is the end of this block, or up to the last number we're supposed to count to

		r = b + BLOCK_SIZE;
		if (r > w) r = w;

		// we don't need to consider factors more than sqrt(r)

		l = 1 + (unsigned int) (sqrt (r / 6));

		// reset the buffer for this block

		memset (A, 1, sizeof (A));

		// sieve through the numbers in this range using known small(ish) prime numbers

		for (m=0; m<2; m++) for (p=0; p<2; p++) {
			a = C[p];

			// for each prime number that's a multiple of 6 * i + a, use it to sieve

			for (i=0; i<Z[p]; i++) {
				f = P[m][i];
				y = f / 6;

				// stop if we have investigated enough primes

				if (y >= l) break;

				// find where to start sieving

				o = b * 6;
				o += f - (o % f);
				for (j=o; ;j+=f) {
					k = (j - a) / 6;
					if (k * 6 + a == j)
						if (k != y) break;
				}

				// mark off all the multiples of f

				for (j=k; j<r; j+=f) A[p][j-b] = 0;
			}
		}

		// collect the results for the current block

		for (p=0; p<2; p++) {

			// figure out the end of this block

			c = C[p];
			k = (v - c) / 6 + 1;
			if (r > k) r = k;

			// start at the beginning of the block represented in A, which is the end of this block range minus the beginning of the range

			l = r - b;
			for (i=0; i<l; i++) {
				if (A[p][i]) println ((i + b) * 6 + c);
				n += A[p][i];
			}
		}
	}
	return n;
}

int main (int argc, char *argv[]) {
	uint64_t n;
	if (argc != 2) {
		printf ("Usage: %s <n> to compute pi(n)\n", argv[0]);
		return 1;
	}
	sscanf (argv[1], "%lu", &n);
	printf ("pi(%lu) = %lu\n", (unsigned long int) n, (unsigned long int) pi (n));
	return 0;
}

