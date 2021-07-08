# Make rootDict file(src/customPdfDict.cxx), rootmap lib file(lib/libcustomPdf.rootmap), and pcm file
rootcling -f src/customPdfDict.cxx -rml lib/libcustomPdf.so -rmf lib/libcustomPdf.rootmap -Iinc GausCB.h ModGaus.h Linkdef.h
# Move pcm to location it can by found by LD_LIBRARY_PATH
mv src/customPdfDict_rdict.pcm lib
# Make lib/libcustomPdf.so
g++ -fPIC src/GausCB.cxx src/ModGaus.cxx src/customPdfDict.cxx -shared -o lib/libcustomPdf.so -g -Iinc `root-config --cflags` `root-config --glibs` `root-config --ldflags` -lRooFit -lMathCore
