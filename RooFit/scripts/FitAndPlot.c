#include "WeightedSBFit.c"
#include "PlotFits.c"

void FitAndPlot() {
  vector<int> mus = {0};
  vector<int> lumis = {42};
  for(size_t im = 0; im < mus.size(); im++) {
    for(size_t il = 0; il < lumis.size(); il++) {
      WeightedFit(mus.at(im),lumis.at(il));
      PlotFit(mus.at(im), lumis.at(il));
    }
  }
  return;
}

