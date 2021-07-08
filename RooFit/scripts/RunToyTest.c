#include "ToyTest.c"

void RunToyTest(){
  int ntoys(100);
  int lumi(42);
// Using best-fit function to generate toys
//   ToyTest( 0,lumi,ntoys); 
//   ToyTest( 1,lumi,ntoys); 
//   ToyTest( 2,lumi,ntoys); 
//   ToyTest( 3,lumi,ntoys);
//   ToyTest( 4,lumi,ntoys);
//   ToyTest( 5,lumi,ntoys);
//   ToyTest( 6,lumi,ntoys);
//   ToyTest( 7,lumi,ntoys);
//   ToyTest( 8,lumi,ntoys);
//   ToyTest( 9,lumi,ntoys);
//   ToyTest(10,lumi,ntoys);
// Using MC to generate toys
  ToyTest( 0,lumi,ntoys,true); 
//   ToyTest( 1,lumi,ntoys,true); 
//   ToyTest( 2,lumi,ntoys,true); 
//   ToyTest( 3,lumi,ntoys,true);
//   ToyTest( 4,lumi,ntoys,true);
//   ToyTest( 5,lumi,ntoys,true);
//   ToyTest( 6,lumi,ntoys,true); 
//   ToyTest( 7,lumi,ntoys,true); 
//   ToyTest( 8,lumi,ntoys,true); 
//   ToyTest( 9,lumi,ntoys,true);
//   ToyTest(10,lumi,ntoys,true);
  return;
}
