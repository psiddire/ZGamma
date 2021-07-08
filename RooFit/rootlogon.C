void rootlogon() {
  string processLine = string(".include") + string(getenv("ANALYSIS_DIR")) + "/inc";
  gROOT->ProcessLine(processLine.c_str());
}
