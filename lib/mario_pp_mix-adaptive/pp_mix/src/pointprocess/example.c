/***
Code for generating perfect samples of  a Strauss point point process.
Is compiled using: make -f example.makefile
Produces executable: example
Run using: ./example 1234   (where 1234 is a user specified random seed)
***/

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits>

/*GSL library*/
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "point2pattern.h"
#include "pointprocess.h"
#include "sampler.h"
#include "strauss.h"

int main(int argc, char *argv[]){
  FILE *out;
  
  // Initialise random number generator
  gsl_rng *r;
  r = gsl_rng_alloc(gsl_rng_taus); 
  gsl_rng_set(r, 12312312);

  // Initialise point pattern
  Point2Pattern ExamplePattern(0,200,0,200, 9,9);
  // initialise sampler
  Sampler ExampleSimulator;

  StraussProcess ExampleProcess(0,10,  // range of x
				0,5,   // range of y
				5,     // beta
				0.5,   // gamma
				0.25); // R

  // Generate perfect sample of Strauss process
  ExampleSimulator.Sim(&ExamplePattern,&ExampleProcess,r);

  std::cout << ExamplePattern.to_eigen_mat() << std::endl;
  // Write point pattern to disk
  ExamplePattern.DumpToFile("example.dat");

  // Print summary of parameters to screen
  ExampleProcess.Print();

}
