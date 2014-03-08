#include <R.h>
#include <Rmath.h>


extern "C" {
  
  void teststatistic(double *xValues,double *yValues,int *xLength,int*yLength,double *lambda,double *result)
  {
	  double m=*xLength;
	  double n=*yLength;
	  double l=*lambda;

	  double T_m=0;
	  for(int i=0; i<m; i++)
	  {
		  for(int j=0; j<m; j++)
		  {
			  
			  T_m=bessel_i((2*sqrt(xValues[j]*xValues[i]))/l,0,1)*exp( (-xValues[j]-xValues[i] )/l)+T_m;
		  }
	  } 
	  T_m=(1/(m*m))*T_m;
	  
	  double T_n=0;
	  for(int i=0; i<n; i++)
	  {
		  for(int j=0; j<n; j++)
		  {
			  T_n=bessel_i((2*sqrt(yValues[j]*yValues[i]))/l,0,1)*exp( (-yValues[j]-yValues[i] )/l)+T_n;
		  }
	  } 
	  T_n=(1/(n*n))*T_n;
	  
	  double T_mixed=0;
	  for(int i=0; i<m; i++)
	  {
		  for(int j=0; j<n; j++)
		  {
			  T_mixed=bessel_i((2*sqrt(yValues[j]*xValues[i]))/l,0,1)*exp( (-yValues[j]-xValues[i] )/l)+T_mixed;
		  }
	  } 
	  T_mixed=(2/(m*n))*T_mixed;
	  
	  *result=( ((m*n)/(m+n))*(T_n+T_m-T_mixed) );

  }
  
    
  void teststatisticStandardized(double *xValues,double *yValues,int *xLength,int*yLength,double *lambda,double *meanXY,double *result)
  {
	  double m=*xLength;
	  double n=*yLength;	
	  double l=*lambda;
	  double XY_mean=*meanXY;
	  

	  double T_m=0;
	  for(int i=0; i<m; i++)
	  {
		  for(int j=0; j<m; j++)
		  {
			  
			  T_m=bessel_i(((2/(XY_mean*l))*sqrt(xValues[j]*xValues[i])),0,1)*exp((-xValues[j]-xValues[i])/(XY_mean*l))+T_m;
		  }
	  } 
	  T_m=(1/(m*m))*T_m;
	  
	  double T_n=0;
	  for(int i=0; i<n; i++)
	  {
		  for(int j=0; j<n; j++)
		  {
			  T_n=bessel_i(((2/(XY_mean*l))*sqrt(yValues[j]*yValues[i])),0,1)*exp((-yValues[j]-yValues[i])/(XY_mean*l))+T_n;
		  }
	  } 
	  T_n=(1/(n*n))*T_n;
	  
	  double T_mixed=0;
	  for(int i=0; i<m; i++)
	  {
		  for(int j=0; j<n; j++)
		  {
			  T_mixed=bessel_i(((2/(XY_mean*l))*sqrt((yValues[j]*xValues[i]))),0,1)*exp(-(yValues[j]/(XY_mean*l))-(xValues[i]/(XY_mean*l)))+T_mixed;
		  }
	  } 
	  T_mixed=(2/(m*n))*T_mixed;
	  
	  *result=( ((m*n)/(m+n))*(T_n+T_m-T_mixed) );

  }



  
}
