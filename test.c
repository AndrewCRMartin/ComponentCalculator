#include "compcalc.h"
#include "eseries.h"

/***********************************************************************/
void TestEvaluate(void)
{
   EVAL eval;
   GENE gene;
   REAL *values = NULL,
        target = 105.0;
   int  i, NValues,
        minComponents = 2,
        maxComponents = 3,
        type = TYPE_RES;

   values = PopulateESeries(e24Base, 24, &NValues, 0, 6);
   gene.values    = (REAL *)malloc(maxComponents * sizeof(REAL));
   gene.operators = (int  *)malloc(maxComponents * sizeof(int));

   gene.NComp = PickRandomComponentNumber(minComponents,
                                          maxComponents);

   for(i=0; i<gene.NComp; i++)
   {
      gene.values[i]    = PickRandomEValue(values, NValues);
      gene.operators[i] = PickRandomOperator();
   }

   for(i=0; i<gene.NComp; i++)
   {
      if(i==0)
      {
         printf("[%.1f] ",
                gene.values[i]);
      }
      else
      {
         printf("[%.1f %s] ",
                gene.values[i],
                (gene.operators[i]==OP_SERIES?"SER":"PAR"));
      }
   }
   printf("\n");
   
   eval = EvaluateGene(&gene, type, target, CT_NONE);
   printf("   EVAL: V=%.1f S=%.1f E=%.1f P=%.1f%% D=%.1f\n",
          eval.value, eval.score, eval.error, eval.percentageError,
          eval.compDifference);
}

/***********************************************************************/
int main(int argc, char **argv)
{
   srand(time(NULL));
   
   TestEvaluate();
   return(0);
}
