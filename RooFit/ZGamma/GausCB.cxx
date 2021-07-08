/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "GausCB.h" 
#include "Math/DistFunc.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 

ClassImp(GausCB); 

 GausCB::GausCB(const char *name, const char *title, 
                        RooAbsReal& _m,
                        RooAbsReal& _mean,
                        RooAbsReal& _gaus_std,
                        RooAbsReal& _cb_std,
                        RooAbsReal& _cb_alpha,
                        RooAbsReal& _cb_n,
                        RooAbsReal& _cb_frac) :
   RooAbsPdf(name,title), 
   m("m","m",this,_m),
   mean("mean","mean",this,_mean),
   gaus_std("gaus_std","gaus_std",this,_gaus_std),
   cb_std("cb_std","cb_std",this,_cb_std),
   cb_alpha("cb_alpha","cb_alpha",this,_cb_alpha),
   cb_n("cb_n","cb_n",this,_cb_n),
   cb_frac("cb_frac","cb_frac",this,_cb_frac)
 { 
 } 


 GausCB::GausCB(const GausCB& other, const char* name) :  
   RooAbsPdf(other,name), 
   m("m",this,other.m),
   mean("mean",this,other.mean),
   gaus_std("gaus_std",this,other.gaus_std),
   cb_std("cb_std",this,other.cb_std),
   cb_alpha("cb_alpha",this,other.cb_alpha),
   cb_n("cb_n",this,other.cb_n),
   cb_frac("cb_frac",this,other.cb_frac)
 { 
 } 



 Double_t GausCB::evaluate() const 
 { 
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
   double g = ROOT::Math::gaussian_pdf(m,gaus_std,mean);
   double cb = ROOT::Math::crystalball_pdf(m,cb_alpha, cb_n, cb_std, mean);
   return (1-cb_frac)*g+cb_frac*cb;
 } 


