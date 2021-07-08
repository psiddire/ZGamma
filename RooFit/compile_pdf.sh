# Make rootDict file(src/customPdfDict.cxx), rootmap lib file(lib/libcustomPdf.rootmap), and pcm file
rootcling -f src/customPdfDict.cxx -rml lib/libcustomPdf.dylib -rmf lib/libcustomPdf.rootmap -Iinc GausCB.h ModGaus.h ModGaus11.h ModGaus01.h ModThreeGaus.h Linkdef.h
# Move pcm to location it can by found by LD_LIBRARY_PATH
mv src/customPdfDict_rdict.pcm lib
# Make lib/libcustomPdf.so
g++ -fPIC src/GausCB.cxx src/ModGaus.cxx src/ModGaus11.cxx src/ModGaus01.cxx src/ModThreeGaus.cxx src/customPdfDict.cxx -shared -o lib/libcustomPdf.dylib -g -Iinc `root-config --cflags` `root-config --glibs` `root-config --ldflags` -lRooFit -lMathCore -lRooFitCore
