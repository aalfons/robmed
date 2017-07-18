/*------------------------------------------------------*/
/* Support functions for lmRob()                        */
/* Author: Jeffrey Wang                                 */
/* Date  : 12/20/1999                                   */
/*------------------------------------------------------*/

#include <R.h>
#include <time.h>

#if defined(S_newio_IS_STDIO_FILE)
#include "newredef.h"
#endif

void F77_NAME(mminitclk)(Sint *ct)
{
  *ct = (Sint) clock();
}

void F77_NAME(mmprint)(Sint *nrep, Sint *itmp, Sint *ct, 
                  Sint *elap, Sint *ninc)
{
  Sint lapt, nleft, ihr;
  Sfloat tmp;

  nleft=(*nrep-(*ninc)*(*itmp))/(*ninc)+1;
  *elap += clock()-*ct;
  *ct = (Sint) clock();

  tmp = (Sfloat) nleft /CLOCKS_PER_SEC;
  lapt = (Sint) ((*elap)/(*itmp)*tmp); 
  
  if (lapt < 60)
    Rprintf("00:00:%02d left\n", lapt);
//    printf("00:00:%02ld left\n", lapt);
  else if (lapt < 360) {
    nleft = lapt/60;
    lapt = lapt % 60;
    Rprintf("00:%02ld:%02d left\n", nleft, lapt);
//    printf("00:%02ld:%02ld left\n", nleft, lapt);
  }
  else {
    ihr = lapt/360;
    lapt = lapt % 360;
    nleft = lapt/60;
    lapt = lapt % 60;
    Rprintf("%ld:%02ld:%02d left\n", ihr, nleft, lapt);
//    printf("%ld:%02ld:%02ld left\n", ihr, nleft, lapt);
  }
}
