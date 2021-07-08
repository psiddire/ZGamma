H->Zgamma spectrum fitting
========

Code for fitting the mllg spectrum with RooFit, using a variety of custom functions.

## Setup

Custom libraries were created for the pdfs used to fit the signal and background shapes. In order to save these objects in a RooWorkspace, the class code must be included in the workspace:

~~~~bash
source set_env.sh
./compile_pdf.sh
~~~~

## Fitting methods

### Background fits

The code in `FitCategory.c` performs a background-only fit in the specified analysis category (currently the categories available are those currently being used in the CMS H->Zgamma analysis, consisting of 3 dijet categories and 4 untagged categories).

### Signal


### Signal + Background 

