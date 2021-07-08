/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "ModGaus01.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 

ClassImp(ModGaus01); 

 ModGaus01::ModGaus01(const char *name, const char *title, 
                        RooAbsReal& _m,
                        RooAbsReal& _m0,
                        RooAbsReal& _v0,
                        RooAbsReal& _s0,
                        RooAbsReal& _s1) :
   RooAbsPdf(name,title), 
   m("m","m",this,_m),
   m0("m0","m0",this,_m0),
   v0("v0","v0",this,_v0),
   s0("s0","s0",this,_s0),
   s1("s1","s1",this,_s1)
 { 
 } 


 ModGaus01::ModGaus01(const ModGaus01& other, const char* name) :  
   RooAbsPdf(other,name), 
   m("m",this,other.m),
   m0("m0",this,other.m0),
   v0("v0",this,other.v0),
   s0("s0",this,other.s0),
   s1("s1",this,other.s1)
 { 
 } 

 Double_t ModGaus01::evaluate() const 
 { 
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
  Double_t width, power;
  power = v0; 
  width = s0 + s1*(m - 100)/(180-100);
  return exp(-1*pow(fabs((m-m0)/width), power));
 } 
