## Modified code from HZZ->4l analysis with Z1 edition
 
 Kept the FSR candidate 

# Method to implement in code
  
  **Instalation**
  - put in CMSSW/src, do scram b -j 8 **Notice that some of CMSSW version could not compile it**
  - gSystem->LoadMacro("xx/xx/KinZfitter.C")&Helperfunction.C
  - compile it when you run it 

  ``` root[0] .L xx/xx/KinZfitter.C ``` --> compile in root first
  ``` root[1] .L xxx.C ``` --> run your code for the fitter


  **run it**

  0. (Add the package into your BuildFile.xml)
  
  1. include the head file  

  ``` #include "KinZfitter/KinZfitter/interface/KinZfitter.h"```
  
  2. Declare and then initialize the KinZfitter class when initializing your analyzer 

  ``` KinZfitter *kinZfitter; kinZfitter = new KinZfitter(isData); //(In data, (isData=true). In mc (isData=false))```

  3. Prepare inputs after Higgs candidate is formed

  leptons: Suppose Lep_Z1_1,Lep_Z1_2 in TLorentzVector

  ``` TLorentzVector lep(pt,eta,phi,m)```

  4. Setup, refit and get the refitted results:

  In your analyzer, do:

  ```cpp
  kinZfitter->Setup(selectedLeptons, selectedFsrMap);
  kinZfitter->KinRefitZ();
  
  // refit mass4l
  double mass4lREFIT = kinZfitter->GetRefitM4l();
  // four 4-vectors after refitting order by Z1_1,Z1_2,Z2_1,Z2_2
  vector < TLorentzVector > p4 = kinZfitter->GetRefitP4s(); 
  ```

  5.Support functions
  
  **refitZmass**

  ``` double massZ1REFIT = kinZfitter->GetRefitMZ1();```

  **refit lepton four vector**

  ``` vector<TLorentzVector> refitlep = kinZfitter->GetRefitP4s()```