#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/macros.h"

REAL e3Base[]  = {1.0, 2.2, 4.7};
REAL e6Base[]  = {1.0, 1.5, 2.2, 3.3, 4.7, 6.8};
REAL e12Base[] = {1.0, 1.2, 1.5, 1.8, 2.2, 2.7,
                  3.3, 3.9, 4.7, 5.6, 6.8, 8.2};
REAL e24Base[] = {1.0, 1.1, 1.2, 1.3, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4,
                  2.7, 3.0, 3.3, 3.6, 3.9, 4.3, 4.7, 5.1, 5.6, 6.2,
                  6.8, 7.5, 8.2, 9.1};
REAL e48Base[] = {1.00, 1.05, 1.10, 1.15, 1.21, 1.27, 1.33, 1.40, 1.47,
                  1.54, 1.62, 1.69, 1.78, 1.87, 1.96, 2.05, 2.15, 2.26,
                  2.37, 2.49, 2.61, 2.74, 2.87, 3.01, 3.16, 3.32, 3.48,
                  3.65, 3.83, 4.02, 4.22, 4.42, 4.64, 4.87, 5.11, 5.36,
                  5.62, 5.90, 6.19, 6.49, 6.81, 7.15, 7.50, 7.87, 8.25,
                  8.66, 9.09, 9.53};
REAL e96Base[] = {1.00, 1.02, 1.05, 1.07, 1.10, 1.13, 1.15, 1.18, 1.21,
                  1.24, 1.27, 1.30, 1.33, 1.37, 1.40, 1.43, 1.47, 1.50,
                  1.54, 1.58, 1.62, 1.65, 1.69, 1.74, 1.78, 1.82, 1.87,
                  1.91, 1.96, 2.00, 2.05, 2.10, 2.15, 2.21, 2.26, 2.32,
                  2.37, 2.43, 2.49, 2.55, 2.61, 2.67, 2.74, 2.80, 2.87,
                  2.94, 3.01, 3.09, 3.16, 3.24, 3.32, 3.40, 3.48, 3.57,
                  3.65, 3.74, 3.83, 3.92, 4.02, 4.12, 4.22, 4.32, 4.42,
                  4.53, 4.64, 4.75, 4.87, 4.99, 5.11, 5.23, 5.36, 5.49,
                  5.62, 5.76, 5.90, 6.04, 6.19, 6.34, 6.49, 6.65, 6.81,
                  6.98, 7.15, 7.32, 7.50, 7.68, 7.87, 8.06, 8.25, 8.45,
                  8.66, 8.87, 9.09, 9.31, 9.53, 9.76};
REAL e192Base[] = {1.00, 1.01, 1.02, 1.04, 1.05, 1.06, 1.07, 1.09, 1.10,
                   1.11, 1.13, 1.14, 1.15, 1.17, 1.18, 1.20, 1.21, 1.23,
                   1.24, 1.26, 1.27, 1.29, 1.30, 1.32, 1.33, 1.35, 1.37,
                   1.38, 1.40, 1.42, 1.43, 1.45, 1.47, 1.49, 1.50, 1.52,
                   1.54, 1.56, 1.58, 1.60, 1.62, 1.64, 1.65, 1.67, 1.69,
                   1.72, 1.74, 1.76, 1.78, 1.80, 1.82, 1.84, 1.87, 1.89,
                   1.91, 1.93, 1.96, 1.98, 2.00, 2.03, 2.05, 2.08, 2.10,
                   2.13, 2.15, 2.18, 2.21, 2.23, 2.26, 2.29, 2.32, 2.34,
                   2.37, 2.40, 2.43, 2.46, 2.49, 2.52, 2.55, 2.58, 2.61,
                   2.64, 2.67, 2.71, 2.74, 2.77, 2.80, 2.84, 2.87, 2.91,
                   2.94, 2.98, 3.01, 3.05, 3.09, 3.12, 3.16, 3.20, 3.24,
                   3.28, 3.32, 3.36, 3.40, 3.44, 3.48, 3.52, 3.57, 3.61,
                   3.65, 3.70, 3.74, 3.79, 3.83, 3.88, 3.92, 3.97, 4.02,
                   4.07, 4.12, 4.17, 4.22, 4.27, 4.32, 4.37, 4.42, 4.48,
                   4.53, 4.59, 4.64, 4.70, 4.75, 4.81, 4.87, 4.93, 4.99,
                   5.05, 5.11, 5.17, 5.23, 5.30, 5.36, 5.42, 5.49, 5.56,
                   5.62, 5.69, 5.76, 5.83, 5.90, 5.97, 6.04, 6.12, 6.19,
                   6.26, 6.34, 6.42, 6.49, 6.57, 6.65, 6.73, 6.81, 6.90,
                   6.98, 7.06, 7.15, 7.23, 7.32, 7.41, 7.50, 7.59, 7.68,
                   7.77, 7.87, 7.96, 8.06, 8.16, 8.25, 8.35, 8.45, 8.56,
                   8.66, 8.76, 8.87, 8.98, 9.09, 9.20, 9.31, 9.42, 9.53,
                   9.65, 9.76, 9.88};


typedef struct _gene
{
   REAL *values;
   int  *operators; /* Note, the first operator is ignored */
   int  NComp;
} GENE;

typedef struct _eval
{
   REAL value,
        error,
        percentageError,
        compDifference,
        score;
} EVAL;
      
#define OP_PARALLEL 0
#define OP_SERIES   1
#define TYPE_CAP    0
#define TYPE_RES    1

/***********************************************************************/
REAL PickRandomEValue(REAL *values, int NValues)
{
   int offset = rand() % NValues;
   return(values[offset]);
}

int PickRandomOperator(void)
{
   int offset = rand() % 2;
   return(offset?OP_SERIES:OP_PARALLEL);
}

int PickRandomComponentNumber(int minNum, int maxNum)
{
   return(rand() % (maxNum + 1 - minNum) + minNum);
}

/***********************************************************************/
GENE *InitializePopulation(int NGenes, int minComponents,
                           int maxComponents,
                           REAL *values, int NValues)
{
   GENE *genes = NULL;
   BOOL error = FALSE;
   int  i;
   
   /* Allocate memory for the array of genes */
   if((genes = (GENE *)malloc(NGenes * sizeof(GENE)))==NULL)
      return(NULL);

   /* Allocate memory for each gene */
   for(i=0; i<NGenes; i++)
   {
      if((genes[i].values =
          (REAL *)malloc(maxComponents * sizeof(REAL)))==NULL)
      {
         error=TRUE;
         break;
      }

      if((genes[i].operators =
          (int *)malloc(maxComponents * sizeof(int)))==NULL)
      {
         error=TRUE;
         break;
      }
   }

   if(error)
   {
      /* TODO: free up memory of items in the genes on error */
      FREE(genes);
      return(NULL);
   }

   /* Populate the genes with random values */
   for(i=0; i<NGenes; i++)
   {
      int j;

      genes[i].NComp = PickRandomComponentNumber(minComponents,
                                                 maxComponents);

      for(j=0; j<genes[i].NComp; j++)
      {
         genes[i].values[j]    = PickRandomEValue(values, NValues);
         genes[i].operators[j] = PickRandomOperator();
      }
   }

   return(genes);
}
   

/***********************************************************************/
REAL *PopulateESeries(REAL *eBase, int numberInSeries, int *NValues,
                      REAL lowPower, REAL highPower)
{
   int i, j, power;
   int NRanges = (highPower-lowPower)+1;
   REAL *eValues = NULL;
   
   /* Calculate values in the whole series between 1R and 10M */
   *NValues = 1 + (NRanges * numberInSeries);
   /* Allocate space */
   if((eValues = (REAL *)malloc(*NValues * sizeof(REAL)))==NULL)
      return(NULL);
   /* Populate */
   i=0;
   for(power=lowPower; power<=highPower; power++)
   {
      for(j=0; j<numberInSeries; j++)
      {
         eValues[i++] = eBase[j] * pow(10,power);
      }
   }
   /* Add the first value above the last range
      (e.g. 10M for resistors)
   */
   eValues[i++] = pow(10, (highPower+1));

   return(eValues);
}

EVAL EvaluateGene(GENE *gene, int type, REAL target)
{
   int i;
   REAL value = gene->values[0];
   EVAL eval;
   REAL compDiff = 0.0;
   int  NPairs = 0;

   for(i=1; i<gene->NComp; i++)
   {
      /* Resistors in series or caps in parallel                       */
      if(((type == TYPE_RES) && (gene->operators[i] == OP_SERIES)) ||
         ((type == TYPE_CAP) && (gene->operators[i] == OP_PARALLEL)))
      {
         value += gene->values[i];
      }
      else /* Resistors in parallel or caps in series                  */
      {
         value = 1/((1/value) + (1/gene->values[i]));
      }
   }

   for(i=0; i<gene->NComp; i++)
   {
      int j;
      for(j=i+1; j<gene->NComp; j++)
      {
         compDiff += fabs(gene->values[i] - gene->values[j]) /
            MAX(gene->values[i], gene->values[j]);
         NPairs++;
      }
   }
   
   eval.value = value;
   eval.error = fabs(value - target);
   eval.percentageError = 100*(eval.error/target);
   if(gene->NComp == 1)
      eval.compDifference = 0;
   else
      eval.compDifference = compDiff / NPairs;

   /* To optimize.... */
   eval.score = 10 * eval.percentageError + eval.compDifference;

   return(eval);
}


/* holds the address of the array of which the sorted index
 * order needs to be found
 */
int * base_arr = NULL;
/* Note how the compare function compares the values of the
 * array to be sorted. The passed value to this function
 * by `qsort' are actually the `idx' array elements.
 */

static REAL *scores = NULL;

int compareScores(const void *a, const void *b)
{
    int aa = *((int *) a),
        bb = *((int *) b);

    if (scores[aa] < scores[bb])
    {
       return -1;
    }
    else if (scores[aa] == scores[bb])
    {
       return 0;
    }

    return 1;
}

int *RankPopulation(GENE *genes, int NGenes, int type, REAL target)
{
   int *idx = NULL,
      i;

   /* Allocate memory for the index                                    */
   if((idx = (int *)malloc(NGenes * sizeof(int)))==NULL)
      return(NULL);

   /* Allocate memory for the scores                                   */
   if((scores = (REAL *)malloc(NGenes * sizeof(REAL)))==NULL)
   {
      FREE(idx);
      return(NULL);
   }

   /* Populate the scores and initialize the index                     */
   for(i=0; i<NGenes; i++)
   {
      EVAL eval;
      eval = EvaluateGene(&(genes[i]), type, target);
      scores[i] = eval.score;
      idx[i] = i;
   }

   /* Perform the sort on the index                                    */
   qsort(idx, NGenes, sizeof(int), compareScores);
   FREE(scores);

   return idx;
}

/***********************************************************************/
int main(int argc, char **argv)
{
   int NValues, i, j, k;
   int NGenes = 100;
   int minComponents = 1;
   int maxComponents = 5;
   int type = TYPE_RES;
   REAL *values = PopulateESeries(e3Base, 3, &NValues, 0, 6);
   GENE *genes = NULL;
   EVAL eval;
   REAL target = 100.0;
   int  *rank = NULL;

   srand(time(NULL));

/*
   for(i=0; i<NValues; i++)
      printf("%.2f\n", values[i]);
*/

   genes = InitializePopulation(NGenes, minComponents, maxComponents,
                                values, NValues);

   rank = RankPopulation(genes, NGenes, type, target);


   for(i=0; i<NGenes; i++)
   {
      j=rank[i];
      printf("%3d: ", j);
      for(k=0; k<genes[j].NComp; k++)
      {
         printf("[%.1f %s] ",
                genes[j].values[k],
                (genes[j].operators[k]==OP_SERIES?"PAR":"SER"));
      }
      printf("\n");
      
      eval = EvaluateGene(&(genes[j]), type, target);
      
      printf("   EVAL: V=%.1f S=%.1f E=%.1f P=%.1f%% D=%.1f\n",
             eval.value, eval.score, eval.error, eval.percentageError,
             eval.compDifference);
   }

   FREE(rank);
   
   return(0);
}
