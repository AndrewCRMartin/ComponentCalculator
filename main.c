#define SCAN_ALL 1
#include "compcalc.h"
#include "eseries.h"

void PopulateGenes(GENE *genes, ULONG NGenes, int minComponents, int maxComponents, REAL *values, int NValues)
{
   int i, j;
   ULONG geneNum = 0;

   for(i=0; i<maxComponents; i++)
   {
      geneNum = 0;
      for(j=0; j<NValues-1; j++)
         genes[geneNum++].values[i] = values[j];

      genes[geneNum++].values[i] = 0.0;

      if(geneNum > NGenes)
      {
         printf("Logic Error!!!\n");
         exit(1);
      }
   }
}


GENE *CreateAll(ULONG *NGenes, int minComponents, int maxComponents, REAL *values, int NValues)
{
   GENE *genes = NULL;
   ULONG i;

   *NGenes = (ULONG)pow((NValues+1), maxComponents);
   printf("Calculated we need %ld genes\n", *NGenes);
   if((genes = InitializePopulation(*NGenes, minComponents, maxComponents, values, NValues, FALSE))==NULL)
   {
      fprintf(stderr, "No memory for genes\n");
      return(NULL);
   }

   fprintf(stderr, "Created genes\n");

   /* We use all components and a value of 0 to indicate an unused component   */
   for(i=0; i<*NGenes; i++)
      genes[i].NComp = maxComponents;

   PopulateGenes(genes, *NGenes, minComponents, maxComponents, values, NValues);
   fprintf(stderr, "Populated genes\n");

   return(genes);
}


void ScanAllCombinations(int minComponents, int maxComponents, REAL *values, int NValues, REAL target, int compTarget)
{
   int compNum, *valNums;
   GENE gene;

   valNums = (int *)malloc(maxComponents * sizeof(int));  /* TODO */
   for(compNum=0; compNum<maxComponents; compNum++)
      valNums[compNum] = 0;
   
   for(compNum=0; compNum<maxComponents; compNum++)
   {
      gene.values[compNum] = values[valNums[compNum]];
      valsNums[compNum]++;

      if(geneNum > NGenes)
      {
         printf("Logic Error!!!\n");
         exit(1);
      }
   }

   free(valNums);
   
}


/***********************************************************************/
int main(int argc, char **argv)
{
   int  NValues,
        j, k, iter,
        minComponents = 1,
        maxComponents = 3,
        type = TYPE_RES,
        compTarget = CT_LOW;
   REAL *values = NULL,
        target = 105.0;
   GENE *genes = NULL;
   EVAL eval;
   ULONG NGenes = (ULONG)POPSIZE,
         *rank = NULL;

   
   srand(time(NULL));

   if((values = PopulateESeries(e6Base, 6, &NValues, 0, 6))==NULL)
   {
      return(1);
   }
   
#ifdef DEBUG
   for(i=0; i<NValues; i++)
      fprintf(stderr, "%.2f ", values[i]);
   fprintf(stderr,"\n");
#endif

#ifdef SCAN_ALL
   ScanAllCombinations(minComponents, maxComponents, values, NValues, target, compTarget);
#endif
   
#ifdef ALL_VALUES
   genes = CreateAll(&NGenes, minComponents, maxComponents, values, NValues);
   if((rank = RankPopulation(genes, NGenes, type, target, compTarget))==NULL)
   {
      return(1);
   }
   else
   {
      eval = EvaluateGene(&(genes[rank[0]]), type, target,
                          compTarget);
      if((eval.value == target) &&
         (eval.compDifference < 0.5))
      {
         printf("No good solution found\n");
      }
      printf("Best solution is:\n");
      j=rank[0];
      printf("%3d: ", j);
      for(k=0; k<genes[j].NComp; k++)
      {
         if(k==0)
         {
            printf("[%.1f] ",
                   genes[j].values[k]);
         }
         else
         {
            printf("[%.1f %s] ",
                   genes[j].values[k],
                   (genes[j].operators[k]==OP_SERIES?"SER":"PAR"));
         }
      }
      printf("\n");
   }
#endif
   
   
#ifdef GA
   genes = InitializePopulation(NGenes, minComponents, maxComponents,
                                values, NValues, TRUE);

   for(iter=0; iter<MAXITER; iter++)
   {
      if((rank = RankPopulation(genes, NGenes, type, target,
                                compTarget))!=NULL)
      {
         /* Evaluate the best one                                      */
         eval = EvaluateGene(&(genes[rank[0]]), type, target,
                             compTarget);
         if((eval.value == target) &&
            (eval.compDifference < 0.5))
            break;

         if(iter > MAXITER/2)
         {
            if((eval.percentageError < 1) &&
               (eval.compDifference < 0.5))
               break;
         }

         if(iter > 2*MAXITER/3)
         {
            if(eval.percentageError < 1)
               break;
         }
         MutatePopulation(genes, NGenes, rank,
                          minComponents, maxComponents,
                          values, NValues);
         FREE(rank);
      }
      else
      {
         fprintf(stderr, "Error\n");
         exit(1);
      }
   }
#endif
   
   if(rank == NULL)
   {
      if((rank = RankPopulation(genes, NGenes, type, target,
                                compTarget))==NULL)
      {
         fprintf(stderr, "Error\n");
         exit(1);
      }
   }

   j=rank[0];
   printf("%3d: ", j);
   for(k=0; k<genes[j].NComp; k++)
   {
      if(k==0)
      {
         printf("[%.1f] ",
                genes[j].values[k]);
      }
      else
      {
         printf("[%.1f %s] ",
                genes[j].values[k],
                (genes[j].operators[k]==OP_SERIES?"SER":"PAR"));
      }
   }
   printf("\n");
      
   eval = EvaluateGene(&(genes[j]), type, target, compTarget);
      
   printf("   EVAL: V=%.1f S=%.1f E=%.1f P=%.1f%% D=%.1f\n",
          eval.value, eval.score, eval.error, eval.percentageError,
          eval.compDifference);
   
   FREE(rank);
   
   return(0);
}


void PrintAllGenes(GENE *genes, int NGenes, int *rank,
                   int type, REAL target, int compTarget)
{
   int i, j, k;
   EVAL eval;
   
   for(i=0; i<NGenes; i++)
   {
      j=rank[i];
      printf("%3d: ", j);
      for(k=0; k<genes[j].NComp; k++)
      {
         if(k==0)
         {
            printf("[%.1f] ",
                   genes[j].values[k]);
         }
         else
         {
            printf("[%.1f %s] ",
                   genes[j].values[k],
                   (genes[j].operators[k]==OP_SERIES?"SER":"PAR"));
         }
      }
      printf("\n");
      
      eval = EvaluateGene(&(genes[j]), type, target, compTarget);
      
      printf("   EVAL: V=%.1f S=%.1f E=%.1f P=%.1f%% D=%.1f\n",
             eval.value, eval.score, eval.error, eval.percentageError,
             eval.compDifference);
   }
}
