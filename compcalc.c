#define _GNU_SOURCE 1
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/macros.h"

#include "eseries.h"

typedef struct _gene
{
   REAL *values;
   int  *operators; /* Note, the first operator is ignored             */
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

/***********************************************************************/
int PickRandomOperator(void)
{
   int offset = rand() % 2;
   return(offset?OP_SERIES:OP_PARALLEL);
}

/***********************************************************************/
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
   
   /* Allocate memory for the array of genes                           */
   if((genes = (GENE *)malloc(NGenes * sizeof(GENE)))==NULL)
      return(NULL);

   /* Allocate memory for each gene                                    */
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
      /* TODO: free up memory of items in the genes on error           */
      FREE(genes);
      return(NULL);
   }

   /* Populate the genes with random values                            */
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
   int  i, j,
        power,
        NRanges  = (highPower-lowPower)+1;
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


/***********************************************************************/
EVAL EvaluateGene(GENE *gene, int type, REAL target)
{
   int  i,
        NPairs   = 0;
   REAL value    = gene->values[0],
        compDiff;
   EVAL eval;

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

   /* Calculate diversity in the component values                      */
   compDiff = 0.0;
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

   /* To optimize....                                                  */
   eval.score = 10 * eval.percentageError + eval.compDifference;

   return(eval);
}


/***********************************************************************/
#ifdef _GNU_SOURCE
int compareScores(const void *a, const void *b, void *vScores)
#else
static REAL *scores = NULL;
int compareScores(const void *a, const void *b)
#endif
{
    int aa = *((int *) a),
        bb = *((int *) b);
#ifdef _GNU_SOURCE
    REAL *scores = (REAL *)vScores;
#endif
    
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


/***********************************************************************/
int *RankPopulation(GENE *genes, int NGenes, int type, REAL target)
{
   int *idx = NULL,
        i;
#ifdef _GNU_SOURCE
   REAL *scores = NULL;
#endif

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

      eval      = EvaluateGene(&(genes[i]), type, target);
      scores[i] = eval.score;
      idx[i]    = i;
   }

   /* Perform the sort on the index                                    */
#ifdef _GNU_SOURCE
   qsort_r(idx, NGenes, sizeof(int), compareScores, scores);
#else
   qsort(idx, NGenes, sizeof(int), compareScores);
#endif
   
   FREE(scores);

   return idx;
}

/***********************************************************************/
int main(int argc, char **argv)
{
   int  NValues,
        i, j, k,
        *rank = NULL,
        NGenes = 100,
        minComponents = 1,
        maxComponents = 5,
        type = TYPE_RES;
   REAL *values = NULL,
        target = 100.0;
   GENE *genes = NULL;
   EVAL eval;

   srand(time(NULL));

   if((values = PopulateESeries(e3Base, 3, &NValues, 0, 6))==NULL)
   {
      return(1);
   }
   
#ifdef DEBUG
   for(i=0; i<NValues; i++)
      fprintf(stderr, "%.2f ", values[i]);
   fprintf(stderr,"\n");
#endif

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
