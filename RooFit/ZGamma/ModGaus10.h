/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef MODGAUS10
#define MODGAUS10

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class ModGaus10 : public RooAbsPdf {
public:
  ModGaus10() {} ; 
  ModGaus10(const char *name, const char *title,
	      RooAbsReal& _m,
	      RooAbsReal& _m0,
	      RooAbsReal& _v0,
	      RooAbsReal& _v1,
	      RooAbsReal& _s0);
  ModGaus10(const ModGaus10& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new ModGaus10(*this,newname); }
  inline virtual ~ModGaus10() { }

protected:

  RooRealProxy m ;
  RooRealProxy m0 ;
  RooRealProxy v0 ;
  RooRealProxy v1 ;
  RooRealProxy s0 ;
  
  Double_t evaluate() const ;

private:

  ClassDef(ModGaus10,1) // Your description goes here...
};
 
#endif
