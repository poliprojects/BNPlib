#ifndef GUARD_sampler_h
#define GUARD_sampler_h

#include "point2pattern.h"
#include "pointprocess.h"
#include <math.h>

class Sampler{
 public:
  //PointProcess *PP;
  //Point2Pattern *P2P;
  long int GeneratedPoints, LivingPoints, NoP;
  struct Point{ long int No; double X; double Y; double *Mark; float R; 
    struct Point *next; }; 
  struct Point3{ char Case; char XCell; char YCell; struct Point3 *next; }; 
  //long int UpperLiving[2];
  //Sampler(PointProcess *p){ PP = p;}
  //Sampler();
  ~Sampler(){}
  void Sim(Point2Pattern *p2p, PointProcess *pp, gsl_rng *r);
  void Sim(Point2Pattern *p2p, PointProcess *pp, 
		      long int **summary, gsl_rng *r);
  void Complement(Point2Pattern *p2p, Point2Pattern *inputp2p, long int*,
		  PointProcess *pp, gsl_rng*);
  void Complement(Point2Pattern *p2p, 
		  Point2Pattern *inputp2p, 
		  Point2Pattern *returnp2p,
		  PointProcess *pp, gsl_rng*);
  void GenerateDominatingPoisson(Point *headPoint,
				 long int *GeneratedPoints,
				 long int *LivingPoints,
				 long int *NoP,
				 PointProcess *pp,
				 gsl_rng *r);
  long int BirthDeath(long int TimeStep,
		      struct Point *headLiving,
		      struct Point *headDeleted,
		      struct Point3 *headTransition,
		      Point2Pattern *p2p,
		      PointProcess *pp,
		      gsl_rng *r);
  void Forward(long int TS, long int TT, int TX, int TY,
	       struct Point *Proposal, long int *DDD,
	       Point2Pattern *p2p, PointProcess *pp);
  long int ForwardBD(double Time, Point2Pattern *p2p, PointProcess *pp,
		 gsl_rng *r);
  void ForwardMH(long int Steps, Point2Pattern *p2p, PointProcess *pp,
		 gsl_rng *r);
};

#endif
