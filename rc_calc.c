/* --------------------------------------------------------------- */
/*                   r/c calculator scripting                      */
/* ----------------------------------------------------------------*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/macros.h"

#define _MAIN_ 1
#include "eseries.h"

#define MAXCOMP 10 /* Maximum number of components */


void resCombination(REAL *res, int num_res, int num_comb, int index, REAL *comb, REAL target);



REAL closest_val; /* closest value so far */
REAL closest_diff = 100000000.00; /* diff of val and target */
REAL closest[MAXCOMP]; /* array detailing values of components */
BOOL ser_par_config[MAXCOMP]; /* array detailing serial/parallel */
int  nComp = 0;

/*
Retrieves and prints the best resistor configuration
* target - target resistance value
* numComp - total number of resistors allowed to be used to achieve target val
* compVals - array of resistor values
* num_res - number of resistor values
*/
void resistor(REAL target, int numComp, REAL *compVals, int num_res)
{
   int i;
   /* run through all possible number of components */
   for (i=1; i<=numComp; i++)
   {
      REAL data[MAXCOMP];
      resCombination(compVals, num_res, i, 0, data, target);
   }

   /* print results */
   fprintf(stdout, "Closest Value: %.3f\n", closest_val);
   fprintf(stdout, "Difference:    %.3f\n", closest_diff);
   fprintf(stdout, "Resistor Configuration: ");
   for (i=0; i<numComp; i++)
   {
      if (i<numComp)
      {
         fprintf(stdout, "R%d %f ", i, closest[i]);
         if (i+1<numComp)
         {
            if (ser_par_config[i+1])
               fprintf(stdout, "|| ");
            else
               fprintf(stdout, "+ ");
         }
         else
         {
            break;
         }
      }
   }
}

/*
Calculates the best combination of resistors to achieve a target value.
	* res[] - input array of resistor values
	* num_res	- size of input array of resistor values
	* num_comb	- number of resistors allowed
	* index - index of comb[]
	* comb[] - array of current combination
	* target - the target value
	* No return value - passes current best combination to global values
*/
void resCombination(REAL *res, int num_res, int num_comb, int index, REAL *comb, REAL target)
{
   int i;
   
   /* current combination is complete */
   if (index == num_comb)
   {
      int j, k;
      
      REAL ser_par_size = pow(2,num_comb); /* 2^(number of components) */
      BOOL *ser_par = NULL;
      REAL calc; /* calculated equivalent resistance value */

      if((ser_par = (BOOL *)malloc(ser_par_size * sizeof(BOOL)))==NULL)
         return;
      
      /* step through every possible series/parallel config of current combination */
      for (j=0; j<ser_par_size; j++)
      {
         calc = 0.0;
         /* creates a boolean array of 0s & 1s for all possible combinations */
         for (k=0; k<num_comb; k++)
            ser_par[k] = (j >> k) & 1;

         /* do the calculations for the combination based on series/parallel combo */
         for (k=0; k<num_comb; k++)
         {
            /* first number, just add */
            if (k==0)
            {
               calc = comb[k];
            }
            /* zero means series, add resistance values */
            else if (!ser_par[k])
            {
               calc += comb[k];
            }
            /* one means parallel, inverse of the sum of reciprocals */
            else if (ser_par[k])
            {
               calc = (calc*comb[k])/(calc+comb[k]);
            }

            /* check to see if difference is less than previous best */
            if (fabs(calc - target) < closest_diff)
            {
               /* it is less, so update global values */
               closest_val = calc;
               closest_diff = fabs(calc - target);
               /* clear to zero */
               for (k=0; k<num_comb; k++)
                  closest[k] = 0;

               /* update closest value & series/parallel arrays */
               for (k=0; k<num_comb; k++)
               {
                  closest[k] = comb[k];
                  ser_par_config[k] = ser_par[k];
               }
            }
         }
         FREE(ser_par);
         return;
      }
      FREE(ser_par);
   }

   /* recursively call and replace the index with all possible values */
   for (i=0; i<=num_res && num_res-i+1 >= num_comb-index; i++)
   {
      comb[index] = res[i];
      resCombination(res, num_res, num_comb, index+1, comb, target);
   }
}

int main(int argc, char **argv)
{
   REAL *compVals = NULL;
   REAL target = 5.0;
   int nComp = 4;
   int num_res;

   if((compVals = PopulateESeries(e6Base, 6, &num_res, 1, 6))==NULL)
   {
      return(1);
   }
   resistor(target, nComp, compVals, num_res);
   FREE(compVals);
   return(0);
}
