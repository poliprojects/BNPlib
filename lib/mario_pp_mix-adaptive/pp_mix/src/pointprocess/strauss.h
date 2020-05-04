#ifndef GUARD_strauss_h
#define GUARD_strauss_h

#include "point2pattern.h"
#include "pointprocess.h"

class StraussProcess : public PointProcess {
 public:
  StraussProcess(double xmin, double xmax, double ymin, double ymax, 
		double b, double g, double Ri);
  ~StraussProcess(){}
  void NewEvent(double *x, double *y, double **mark, 
		char *InWindow, gsl_rng *r);
  //double Interaction(double r);
  double Interaction(struct Point2 *Point1, struct Point2 *Point2);
  void CalcBeta(long int xsidepomm, long int ysidepomm, 
	   double *betapomm);
  void CheckBeta(long int xsidepomm, long int ysidepomm, 
		 double *betapomm);
  double lnCondInt(struct Point2 *TempCell, Point2Pattern *p2p);
  double lnCondIntRatio(struct Point2 *TempCell, Point2Pattern *p2p);
  void Beta(struct Point2 *TempCell);
  long int Pairs(Point2Pattern *p2p);
  void CalcBeta(Point2Pattern *p2p);
  double Getbeta(){return(beta);};
  double Getgamma(){return(gamma);};
  double GetR(){return(R);};
  void SetParameters(double b, double g, double Ri);
  void Print();
  protected:
  double beta, gamma, R;
};

class StraussProcessInhom : public PointProcess {
 public:
  StraussProcessInhom(double xmin, double xmax, double ymin, double ymax, 
		 double b, double g, double Ri,
		 double ntrcpt, double slpx, double slpx2,
		 double slpy, double slpy2, double slpxy);
  ~StraussProcessInhom(){}
  void NewEvent(double *x, double *y, double **mark, 
		char *InWindow, gsl_rng *r);
  //double Interaction(double r);
  double Interaction(struct Point2 *Point1, struct Point2 *Point2);
  void CalcBeta(long int xsidepomm, long int ysidepomm, 
	   double *betapomm);
  void CheckBeta(long int xsidepomm, long int ysidepomm, 
		 double *betapomm);
  double lnCondInt(struct Point2 *TempCell, Point2Pattern *p2p);
  double lnCondIntRatio(struct Point2 *TempCell, Point2Pattern *p2p);
  void Beta(struct Point2 *TempCell);
  double Beta(double x, double y);
  void CalcBeta(Point2Pattern *p2p);
  void SetParameters(double b, double g, double Ri);
  void SetParameters(double b, double g, double Ri,
		     double ntrcpt, double slpx, double slpx2,
		     double slpy, double slpy2, double slpxy);
  void Print();
  void SummaryToFile(const char *FileName);
  protected:
  double beta, gamma, R;
  double intercept, slopex, slopex2, slopey, slopey2, slopexy;
  double betaMax, betaMin;
};

#endif
