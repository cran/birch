#include "ll.h"

// Bring in the global variables
extern int DIM;
extern double RADIUS;

vector<int> which_min(double distances[], int nrx, bool min){
  mdebug1("In which_min\n");

  vector<int> keypair(2);
  keypair[0] = 0;   keypair[1] = 0;
  double mDistance = distances[0];
  for (int i=1; i < (nrx-1); i++)
    for (int j=(i+1); j < nrx; j++){
      if ((min && distances[i + nrx*j] < mDistance) || (!min && distances[i + nrx*j] > mDistance)) {
	keypair[0] = i;
	keypair[1] = j;
	mDistance = distances[i+nrx*j];
      }
    }
  return keypair;
}

