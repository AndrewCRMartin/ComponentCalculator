#include "compcalc.h"
#include "eseries.h"

/***********************************************************************/
int main(int argc, char **argv)
{
   int  NValues,
        j, k, iter,
        *rank = NULL,
        NGenes = POPSIZE,
        minComponents = 1,
        maxComponents = 5,
        type = TYPE_RES,
        compTarget = CT_LOW;
   REAL *values = NULL,
        target = 105.0;
   GENE *genes = NULL;
   EVAL eval;

   srand(time(NULL));

   if((values = PopulateESeries(e24Base, 24, &NValues, 0, 6))==NULL)
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
