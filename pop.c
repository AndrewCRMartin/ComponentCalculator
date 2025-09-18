#include "compcalc.h"
#include "eseries.h"


void ScanAllCombinations(int minComponents, int maxComponents, REAL *values, int NValues,
                         int type, REAL target, int compTarget)
{
   int compNum, *valNums;
   int bumpComp = 0,
      i;
   GENE gene;
   REAL *components;

   valNums    = (int  *)malloc(maxComponents * sizeof(int));  /* TODO */
   components = (REAL *)malloc(maxComponents * sizeof(REAL));  /* TODO */

   for(compNum=0; compNum<maxComponents; compNum++)
      valNums[compNum] = 0;

   for(bumpComp = 0; bumpComp<maxComponents; bumpComp++)
   {
      while(valNums[bumpComp] < NValues)
      {
         for(compNum=0; compNum<maxComponents; compNum++)
            components[compNum] = values[valNums[compNum]];

         for(i=0; i<maxComponents; i++)
            printf("%10.1f ", components[i]);
         printf("\n");
         
         
         valNums[bumpComp]++;
      }
   }
   
   free(valNums);
   
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



#ifdef TEST
int main(int argc, char **argv)
{
   int  NValues,
        minComponents = 1,
        maxComponents = 3,
        type = TYPE_RES,
        compTarget = CT_LOW;
   REAL *values = NULL,
        target = 105.0;

   if((values = PopulateESeries(e3Base, 3, &NValues, 0, 6))==NULL)
   {
      return(1);
   }

   
   ScanAllCombinations(minComponents, maxComponents, values, NValues, type, target, compTarget);
   return(0);
}
#endif
