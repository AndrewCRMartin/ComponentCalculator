#define _GNU_SOURCE 1
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/macros.h"

typedef struct _gene
{
   REAL *values;
   int  *operators; /* Note, the first operator is ignored             */
   int  NComp;
   int  minComp;  /* These should be in a master simulation struct     */
   int  maxComp;
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
#define CT_NONE     0
#define CT_LOW      1
#define CT_HIGH     2

#define MAXITER     100000
#define POPSIZE     1000

REAL PickRandomEValue(REAL *values, int NValues);
int PickRandomOperator(void);
int PickRandomComponentNumber(int minNum, int maxNum);
GENE *InitializePopulation(int NGenes, int minComponents,
                           int maxComponents,
                           REAL *values, int NValues);
REAL *PopulateESeries(REAL *eBase, int numberInSeries, int *NValues,
                      REAL lowPower, REAL highPower);
EVAL EvaluateGene(GENE *gene, int type, REAL target, int compTarget);
#ifdef _GNU_SOURCE
int compareScores(const void *a, const void *b, void *vScores);
#else
int compareScores(const void *a, const void *b);
#endif
int *RankPopulation(GENE *genes, int NGenes, int type, REAL target,
                    int compTarget);
void MutateValue(GENE gene, REAL *values, int NValues);
void MutateOperator(GENE gene);
void  MutateNumberOfComponents(GENE gene, int minComp, int maxComp);
void MutatePopulation(GENE *genes, int NGenes, int *rank,
                      int minComp, int maxComp,
                      REAL *values, int NValues);
