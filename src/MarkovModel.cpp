#include "MarkovModel.h"
#include "MemoryAllocators.h"
#include "Random.h"

#include <stdio.h>
#include <string.h>
#include <alloca.h>

#undef VERSION
#include <cblas.h>

#define FREE_ARRAY(ptr)    { if ((ptr) != NULL) delete [] ptr; ptr = NULL; }

#define ALIGNED(x) (typeof(x))__builtin_assume_aligned(x,16)

MarkovModel::MarkovModel()
   {
   states = 0;

   backgroundError = 1e-5;

   matrix = NULL;
   }

MarkovModel::~MarkovModel()
   {
   FreeMemory();
   }

void ffill(float * dest, float value, size_t n)
   {
   size_t size = 16;
   float * __restrict__ d = ALIGNED(dest);
   for (int i = 0; i < 16; i++)
      d[i] = value;
   for (; 2*size < n; size *= 2)
      cblas_scopy(size, dest, 1, dest+size, 1);
   cblas_scopy(n-size, dest, 1, dest+size, 1);
   }

void MarkovModel::Transpose(float * from, float * to, double r, bool quick)
   {
   if (r == 0)
      {
      for (int i = 0; i < states; i++)
         to[i] = from[i];
      return;
      }

#ifndef FAST
   if (quick)
#endif
      {
      float flipRate = r * empiricalFlipRate;

      // This is sum(abs(x)) but probabilities are all positive
      float sum = cblas_sasum(states, from, 1);

      sum *= r * (1.0 - empiricalFlipRate) / states;
      float complement = 1. - r;
      // avoid underflows
      if (sum < 1e-10)
         {
         sum *= 1e15;
         flipRate *= 1e15;
         complement *= 1e15;
         }

      ffill(to, sum, states);
      cblas_saxpy(states, complement, from, 1, to, 1);
      if (quick)
	 {
	 cblas_saxpy(states/2, flipRate, from+states/2, 1, to, 1);
	 cblas_saxpy(states/2, flipRate, from, 1, to+states/2, 1);
	 }
      else
	 {
	 cblas_saxpy(states/2, flipRate, from+1, 2, to, 2);
	 cblas_saxpy(states/2, flipRate, from, 2, to+1, 2);
	 }

      return;
      }

#ifndef FAST
   double flipRate = r * empiricalFlipRate;
   double sum = 0.0;
   for (int i = 0; i < states; i++)
      sum += from[i];
   sum *= r * (1.0 - empiricalFlipRate) / states;
   double complement = 1. - r;

   // avoid underflows
   if (sum < 1e-10)
      {
      sum *= 1e15;
      flipRate *= 1e15;
      complement *= 1e15;
      }

   for (int i = 0; i < states; i++)
      to[i] = from[i] * complement + from[i^1] * flipRate + sum;
#endif
   }

#define WIDTH 16 /* size of cached slice of haplotype matrix */

void MarkovModel::Condition(float * vector, char ** haplotypes, int position,
                            char observed, double e, double freq, bool flipped)
   {
   if (observed == 0) return;

   double prandom = e * freq + backgroundError;
   double pmatch = (1. - e) + prandom;

   // Here, a lookup table is faster than a conditional expression
   // If the haplotype data is in SNP-major orientation (during imputation),
   // we take the fast path using floats; otherwise we use doubles.
   if (flipped)
      {
      float P[8] = { prandom, prandom, prandom, prandom,
		     prandom, prandom, prandom, prandom };
      P[observed] = pmatch;
      for (int i = 0; i < states; i++)
	 vector[i] *= P[haplotypes[position][i]];
      }
   else
      {
      static char *cache = NULL;
      static int cache_pos = -1;
#pragma omp threadprivate(cache)
#pragma omp threadprivate(cache_pos)
      if (!cache)
	 cache = (char *)malloc((states+1)*WIDTH);
      if (cache_pos != position/WIDTH)
	 {
	 cache_pos = position/WIDTH;
	 for (int i = 0; i < states; i++)
	    memcpy(cache+i*WIDTH, &haplotypes[i][cache_pos*WIDTH], WIDTH);
	 }
#ifdef FAST
      float P[8] = { prandom, prandom, prandom, prandom,
		     prandom, prandom, prandom, prandom };
#else
      double P[8] = { prandom, prandom, prandom, prandom,
		      prandom, prandom, prandom, prandom };
#endif
      P[observed] = pmatch;
      for (int i = 0; i < states; i++)
	 vector[i] *= P[cache[i*WIDTH + position%WIDTH]];
      }
   }

void MarkovModel::FreeMemory()
   {
   if (matrix != NULL)
       FreeFloatMatrix(matrix, markers);
   }

void MarkovModel::Allocate(int MARKERS, int STATES)
   {
   if (markers != MARKERS || states != STATES)
      {
      FreeMemory();

      MarkovParameters::Allocate(MARKERS);

      states = STATES;
      matrix = AllocateFloatMatrix(markers, states & 1 ? states + 1 : states);

      if (states & 1)
         for (int i = 0; i < markers; i++)
            matrix[i][states] = 0.0;

      imputedHap.Dimension(markers);
      imputedDose.Dimension(markers);
      leaveOneOut.Dimension(markers);

      imputedAlleles.Dimension(markers);
      }
   }

void MarkovModel::WalkLeft(char * observed, char ** haplotypes, float ** freqs, bool flipped)
   {
   // Initialize likelihoods at first position
   for (int i = 0; i < states; i++)
      matrix[0][i] = 1.;

   // Scan along chromosome
   for (int i = 0; i < markers - 1; i++)
      {
      if (observed[i])
         Condition(matrix[i], haplotypes, i, observed[i], E[i], freqs[observed[i]][i], flipped);
      Transpose(matrix[i], matrix[i+1], R[i], flipped);
      }

   if (observed[markers - 1])
      Condition(matrix[markers - 1], haplotypes, markers - 1, observed[markers - 1], E[markers - 1], freqs[observed[markers - 1]][markers - 1], flipped);
   }

void MarkovModel::Impute(char * major, char * observed, char ** haplotypes, float ** freqs)
   {
   float * swap;
   float * __restrict__ vector = new float [states];
   float * __restrict__ extra = new float [states];
   vector = ALIGNED(vector);
   extra = ALIGNED(extra);

   // Initialize likelihoods at first position
   for (int i = 0; i < states; i++)
      vector[i] = 1.;

   // Scan along chromosome
   for (int i = markers - 1; i > 0; i--)
      {
      float * __restrict__ m = ALIGNED(matrix[i]);
      for (int j = 0; j < states; j++)
         extra[j] = vector[j] * m[j];

      Impute(major, observed, extra, haplotypes, freqs, i);

      if (observed[i])
         Condition(vector, haplotypes, i, observed[i], E[i], freqs[observed[i]][i], true);
      Transpose(vector, extra, R[i - 1], true);

      swap = vector; vector = extra; extra = swap;
      }

   if (observed[0])
      Condition(vector, haplotypes, 0, observed[0], E[0], freqs[observed[0]][0], true);
   Impute(major, observed, vector, haplotypes, freqs, 0);

   delete [] vector;
   delete [] extra;
   }


void MarkovModel::Impute(char * major, char * observed, float * probs,
                         char ** haplotypes, float ** freqs, int position)
   {
   float P[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
#if 0
   // trouble: repeatedly updating a computed element is very slow!
   for (int i = 0; i < states; i++)
      P[haplotypes[i][position]] += probs[i];
#endif
   float * __restrict__ mask = (float *)alloca(states*sizeof(float));
   char * __restrict__ hp = ALIGNED(haplotypes[position]);
   float ptotal = cblas_sasum(states, probs, 1), pmajor = ptotal;
   for (int j = 1; j < 8; j++)
      if ((j != major[position]) && (freqs[j][position] > 0))
	 {
	 for (int i = 0; i < states; i++)
	    {
	    // make gcc's autovectorization happy
	    mask[i] = hp[i];
	    mask[i] = (mask[i] == j);
	    }
	 P[j] = cblas_sdot(states,mask,1,probs,1);
	 if (P[j] > pmajor)
	    P[j] = pmajor;
	 pmajor -= P[j];
	 }
   P[major[position]] = pmajor;

   int mle = 1;
   for (int i = 2; i <= 7; i++)
      if (P[i] >= P[mle])
         mle = i;

   char labels[] = {0, 'a', 'c', 'g', 't', 'd', 'i', 'r'};

   imputedDose[position] += imputedHap[position] = (pmajor / ptotal);
   imputedAlleles[position] = labels[mle];

   double fmatch = 1.0 / (1. - E[position] + E[position] * freqs[major[position]][position] + backgroundError);
   double fmismatch = 1.0 / (E[position] * freqs[major[position]][position] + backgroundError);

   for (int i = 1; i <= 7; i++)
      if (observed[position] == i)
         P[i] *= fmatch;
      else
         P[i] *= fmismatch;

   ptotal = P[1] + P[2] + P[3] + P[4] + P[5] + P[6] + P[7];
   pmajor = P[major[position]];

   leaveOneOut[position] = pmajor / ptotal;
   }

void MarkovModel::ClearImputedDose()
   {
   imputedDose.Zero();
   }

void MarkovModel::ProfileModel(char * observed, char ** haplotypes, float ** freqs)
   {
   if (markers == 0) return;

   // Cumulative probability
   double sum = 0.0;

   // Sample state at the first position
   for (int i = 0; i < states; i++)
      sum += matrix[markers - 1][i];

   double r = globalRandom.Next() * sum;
   int    state = 0;

   for ( sum = 0.0 ; state < states - 1 && sum < r; state++)
      sum = sum + matrix[markers - 1][state];

   if (observed[markers - 1])
      empE[markers - 1] += CountErrors(haplotypes[state][markers - 1], observed[markers - 1], E[markers - 1], freqs[observed[markers - 1]][markers - 1]);
   else
      empE[markers - 1] += E[markers - 1];

   for (int m = markers - 2; m >= 0; m--)
      {
      double sum = 0.0;

      for (int i = 0; i < states; i++)
         sum += matrix[m][i];

      double norec = matrix[m][state] * (1.0 - R[m]);
      double flip = matrix[m][state ^ 1] * R[m] * empiricalFlipRate;
      double rec = sum * R[m] * (1.0 - empiricalFlipRate) / states;

      sum = norec + flip + rec;

      double r = globalRandom.Next() * sum;

      if (r > norec)
         if (r > norec + flip)
            {
            empR[m]++;

            r -= norec - flip;
            r *= states / (R[m] * (1.0 - empiricalFlipRate));

            state = 0;
            for ( sum = 0.0 ; state < states - 1; state++)
               if ( (sum += matrix[m][state]) > r)
                  break;
            }
         else
            {
            empR[m]++;
            empiricalFlips++;

            state ^= 1;
            }

      if (observed[m])
         empE[m] += CountErrors(haplotypes[state][m], observed[m], E[m], freqs[observed[m]][m]);
      else
         empE[m] += E[m];
      }

   empiricalCount++;
   }

double MarkovModel::CountErrors(char copied, char observed, double e, double freq)
   {
   if (observed == 0)
      return e;

   if (observed == copied)
      return e * freq / (1.0 - e + e * freq + backgroundError);

   return e * freq / (e * freq + backgroundError);
   }

double MarkovModel::CountErrors(float * vector, char ** haplotypes, int position, char observed, double e, double freq)
   {
   if (observed == 0)
      return e;

   double match = 0;
   double mismatch = 0;
   double background = 0;

   for (int i = 0; i < states; i++)
      if (haplotypes[i][position] == observed)
         match += vector[i];
      else
         mismatch += vector[i];

   background = (match + mismatch) * backgroundError;
   mismatch = (match + mismatch) * e * freq;
   match *= 1.0 - e;

   return mismatch / (mismatch + match + background);
   }

double MarkovModel::CountRecombinants(float * from, float * to, double r)
   {
   if (r == 0)
      return 0.0;

   double sum = 0.0;
   for (int i = 0; i < states; i++)
      sum += from[i];

   double rsum = 0.0;
   double fsum = 0.0;
   double nrsum = 0.0;

   for (int i = 0; i < states; i++)
      {
      rsum += to[i];
      fsum += from[i] * to[i^1];
      nrsum += from[i] * to[i];
      }

   fsum *= r * empiricalFlipRate;
   rsum *= sum * r * (1.0 - empiricalFlipRate) / states;
   nrsum *= 1.0 - r;

   double total = fsum + rsum + nrsum;

   empiricalFlips += fsum / total;

   return (rsum + fsum) / total;
   }

void MarkovModel::CountExpected(char * observed, char ** haplotypes, float ** freq)
   {
   float * swap;
   float * vector = new float [states];
   float * extra = new float [states];

   // Initialize likelihoods at first position
   for (int i = 0; i < states; i++)
      vector[i] = 1.;

   // Scan along chromosome
   for (int i = markers - 1; i > 0; i--)
      {
      for (int j = 0; j < states; j++)
         extra[j] = vector[j] * matrix[i][j];

      if (observed[i])
         {
         empE[i] += CountErrors(extra, haplotypes, i, observed[i], E[i], freq[observed[i]][i]);
         Condition(vector, haplotypes, i, observed[i], E[i], freq[observed[i]][i], false);
         }
      else
         empE[i] += E[i];

      Transpose(vector, extra, R[i - 1], false);

      empR[i-1] += CountRecombinants(vector, matrix[i - 1], R[i-1]);

      swap = vector; vector = extra; extra = swap;
      }

   if (observed[0])
      {
      Condition(vector, haplotypes, 0, observed[0], E[0], freq[observed[0]][0], false);
      empE[0] += CountErrors(vector, haplotypes, 0, observed[0], E[0], freq[observed[0]][0]);
      }
   else
      empE[0] += E[0];

   empiricalCount++;

   delete [] vector;
   delete [] extra;
   }

/* Local variables: */
/* c-file-style: "minimac" */
/* End: */

