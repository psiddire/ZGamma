double gauss(double *x, double *par) {
  return ROOT::Math::gaussian_pdf(x[0],par[0],par[1]);
}
double crys(double *x, double *par) {
  return ROOT::Math::crystalball_pdf(x[0],2,20,3,10);
}
double b_pdf(double *x, double *par) {
  double vl(2.2), vr(0.7), s0(11), sl(8.2), sh(50), m0(112.6);
  Double_t width, power;
  power = vl + vr*(x[0] - 100)/(180-100);
  if(x[0] <= m0)
    width = sl + (s0 - sl)*(x[0] - 100)/(m0-100);
  else
    width = s0 + (sh - s0)*(x[0] - m0)/(180 - m0);
  return exp(-1*pow(fabs((x[0]-m0)/width), power));
}

void IntegrationTest() {
  TF1 *g = new TF1("g",gauss,-10,10,2); g->SetParameters(2,2);
  TF1 *c = new TF1("c",crys,-100,100,0);
  TF1 *p = new TF1("p",b_pdf,100,180,0);
  cout << "Gaussian Int: " << g->Integral(-10,10) << endl;
  cout << "Crystal Ball Int: " << c->Integral(-100,100) << endl;
  cout << "Background Int: " << p->Integral(100,180) << endl;
  return;
}
