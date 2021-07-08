#include "FitCategory.c"
#include "TString.h"

void runCatFits() {
  gROOT->SetBatch();
  vector<TString> cats = {"total" , "dijet1", "dijet2", "dijet3",
                          "untag1", "untag2", "untag3", "untag4"};
  vector<TString> types = {"ModGaus","ModGaus11"};
  for(size_t it(0); it < types.size(); it++) 
    for(size_t ic(0); ic < cats.size(); ic++) 
      FitCategory(cats.at(ic),types.at(it));
  return;
}
