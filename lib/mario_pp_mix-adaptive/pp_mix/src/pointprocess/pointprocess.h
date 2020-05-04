#ifndef GUARD_pointprocess_h
#define GUARD_pointprocess_h

#include "point2pattern.h"
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


class PointProcess {
 public:
  double Xmin, Xmax, Ymin, Ymax, TotalBirthRate, InteractionRange;
  int MarkLength;
  long int InteractionRangeUpdateFrequency;
  PointProcess(double xmin, double xmax, double ymin, double ymax){
    Xmin = xmin; Xmax = xmax;
    Ymin = ymin; Ymax = ymax;
    MarkLength = 0;
    InteractionRangeUpdateFrequency = -1;
  }
  ~PointProcess(){}
  virtual void NewEvent(double *x, double *y, double **mark,
			char *InWindow, gsl_rng *r)=0;
  virtual double logNewEventDensity(double x, double y, double *mark){
    printf("NewEventDenisty is not defined\n");
    exit(0);
  };
  /*
  virtual void GeneratePoisson(Point *headPoint, 
			       long int *GeneratedPoints,
			       long int *LivingPoints,
			       long int *NoP,
			       gsl_rng *r)=0;
  */
  //virtual double Interaction(double r)=0;
  virtual double Interaction(struct Point2 *Point1, struct Point2 *Point2)=0;
  virtual void CalcBeta(long int xsidepomm, long int ysidepomm, 
		   double *betapomm){ 
    printf("Define CalcBeta...\n");
  }
  virtual void CheckBeta(long int xsidepomm, long int ysidepomm, 
		   double *betapomm){ 
    printf("Define CheckBeta...\n");
  }
  virtual void CalcBeta(Point2Pattern *p2p){
    printf("Define CalcBeta...\n");
  }
  virtual void CalcTempBeta(Point2Pattern *p2p){
    printf("Define CalcTempBeta...\n");
  }
  virtual double lnCondInt(struct Point2 *TempCell, Point2Pattern *p2p)
  { printf("lnCondInt not defined for this process\n"); exit(0);};
  virtual double lnCondIntRatio(struct Point2 *TempCell, Point2Pattern *p2p)
  { printf("lnCondIntRatio not defined for this process\n"); exit(0);};
  virtual double lnDens(Point2Pattern *p2p);
  virtual void Beta(struct Point2 *TempCell){
    TempCell->Beta = 0;
    printf("Define Beta...\n");
    exit(0);
  };
  virtual void Print(){
    printf("Define process...\n");};
  void GenerateDominatingPoisson(Point2Pattern *p2p, gsl_rng *r);
  virtual void ResetInteractionRange(void);
  virtual void UpdateInteractionRange(Point2Pattern *p2p);
  virtual void UpdateInteractionRange(struct Point2 *Point);
  virtual double InteractionRangeFct(struct Point2 *Point);
};

#endif 
