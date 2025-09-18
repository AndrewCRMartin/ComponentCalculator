#include "compcalc.h"


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
GENE *InitializePopulation(ULONG NGenes, int minComponents,
                           int maxComponents,
                           REAL *values, int NValues,
                           BOOL setVals)
{
   GENE *genes = NULL;
   BOOL error = FALSE;
   ULONG  i;
   
   /* Allocate memory for the array of genes                           */
   if((genes = (GENE *)malloc(NGenes * sizeof(GENE)))==NULL)
      return(NULL);

   /* Allocate memory for each gene                                    */
   for(i=0; i<NGenes; i++)
   {
      genes[i].minComp = minComponents;
      genes[i].maxComp = maxComponents;
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
         if(setVals)
            genes[i].values[j] = PickRandomEValue(values, NValues);

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
EVAL EvaluateGene(GENE *gene, int type, REAL target, int compTarget)
{
   int  i,
        NPairs   = 0;
   REAL value,
        compDiff;
   EVAL eval;

   value = gene->values[0];
#ifdef DEBUG
   fprintf(stderr, "Type = %s\n", (type==TYPE_RES)?"Res":"Cap");
   fprintf(stderr, "Value[0] = %f\n", value);
#endif
   for(i=1; i<gene->NComp; i++)
   {
#ifdef DEBUG
      fprintf(stderr, "Value[%d] = %f ", i, gene->values[i]);
      fprintf(stderr, "in %s\n",
              (gene->operators[i]==OP_SERIES)?"Series":"Parallel");
#endif
      
      /* Resistors in series or caps in parallel                       */
      if(((type == TYPE_RES) && (gene->operators[i] == OP_SERIES)) ||
         ((type == TYPE_CAP) && (gene->operators[i] == OP_PARALLEL)))
      {
#ifdef DEBUG
         fprintf(stderr, "Adding series value %f to %f\n",
                 gene->values[i], value);
#endif
         if(gene->values[i] > 0.0)
            value += gene->values[i];
      }
      else /* Resistors in parallel or caps in series                  */
      {
         if(gene->values[i] > 0.0)
            value = 1/((1/value) + (1/gene->values[i]));
      }
   }

   /* Calculate diversity in the component values                      */
   compDiff = 0.0;
   for(i=0; i<gene->NComp; i++)
   {
      int j;

      if(gene->values[i] > 0.0)
      {
         for(j=i+1; j<gene->NComp; j++)
         {
            if(gene->values[j] > 0.0)
            {
               compDiff += fabs(gene->values[i] - gene->values[j]) /
                  MAX(gene->values[i], gene->values[j]);
               NPairs++;
            }
         }
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
   switch(compTarget)
   {
   case CT_LOW:
      eval.score += gene->NComp - gene->minComp;
      break;
   case CT_HIGH:
      eval.score += gene->maxComp - gene->NComp;
      break;
   default:
      break;
   }

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
ULONG *RankPopulation(GENE *genes, ULONG NGenes, int type, REAL target,
                      int compTarget)
{
   ULONG *idx = NULL,
         i;
#ifdef _GNU_SOURCE
   REAL *scores = NULL;
#endif

   /* Allocate memory for the index                                    */
   if((idx = (ULONG *)malloc(NGenes * sizeof(int)))==NULL)
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

      eval      = EvaluateGene(&(genes[i]), type, target, compTarget);
      scores[i] = eval.score;
      idx[i]    = i;
   }

   /* Perform the sort on the index                                    */
#ifdef _GNU_SOURCE
   qsort_r(idx, NGenes, sizeof(ULONG), compareScores, scores);
#else
   qsort(idx, NGenes, sizeof(ULONG), compareScores);
#endif
   
   FREE(scores);

   return idx;
}

void MutateValue(GENE gene, REAL *values, int NValues)
{
   int component = rand() % gene.NComp;
   gene.values[component] = PickRandomEValue(values, NValues);
}

void MutateOperator(GENE gene)
{
   /* 10% chance                                                       */
   if((rand() % 100) > 89)
   {
      /* Swap the operator for a random component                      */
      int component = rand() % gene.NComp;
      gene.operators[component] = gene.operators[component]==0?1:0;
   }
}

void  MutateNumberOfComponents(GENE gene, int minComp, int maxComp)
{
   /* 5% chance                                                        */
   if((rand() % 100) > 94)
   {
      if((gene.NComp == minComp) && (gene.NComp != maxComp))
      {
         gene.NComp++;
      }
      else if((gene.NComp == maxComp) && (gene.NComp != minComp))
      {
         gene.NComp--;
      }
      else
      {
         if(rand() % 2)
            gene.NComp++;
         else
            gene.NComp--;
      }
   }
}

/***********************************************************************/
void MutatePopulation(GENE *genes, ULONG NGenes, int *rank,
                      int minComp, int maxComp,
                      REAL *values, int NValues)
{
   ULONG i, j,
      replaceFrom = NGenes/2 - 1;

   /* Copy the best half to the worst half                             */
   for(i=0; i<replaceFrom; i++)
      genes[rank[i+replaceFrom]] = genes[rank[i]];

   /* Mutate the second half                                           */
   for(i=replaceFrom; i<NGenes; i++)
   {
      j=rank[i];
      MutateValue(genes[j], values, NValues);
      MutateOperator(genes[j]);
      MutateNumberOfComponents(genes[j], minComp, maxComp);
   }
}

