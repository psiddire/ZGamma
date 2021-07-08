/* Reader of ggNtuplizer's and other TTrees.
 *
 * Works both with old as well as with new ntuples. In particular, for tree
 * branches containing vector<...> objects of elementary data types as well as
 * two-dimensional vector<vector<float> > or vector<vector<int> > arrays,
 * underlying arrays are returned (vector::front()) by TreeReader::Get*() call
 * family.
 *
 * NOTE: vector<bool> cannot be handled in this way since booleans are packed
 * as bits into bytes.
 *
 * NOTE: for branches containing variable length arrays, corresponding counter
 * branches must be read off first.
 *
 * Usage of int, long, ... elementary data types for retrieving branch contents
 * is not portable and should be avoided. This is due to the fact that
 * sizeof(int) and sizeof(long) may differ on 32bit and 64bit architectures.
 * Therefore, it is advised to use the following conventions for types of branch
 * contents:
 *   1. Int_t instead of int;
 *   2. UInt_t instead of unsigned int;
 *   3. Long64_t instead of long;
 *   4. ULong64_t instead of unsigned long, etc.
 *
 * For instance, it is suggested to use
 *   Int_t nEle = reader.GetInt("nEle");
 * instead of
 *   int nEle = reader.GetInt("nEle");
 *
 * It is generally safe, however, to use float and double in place of Float_t
 * and Double_t.
 *
 * Some notes:
 * 1. In-place modifications of arrays (within array boundaries) are allowed and
 *    encouraged.
 * 2. The code loads only requested tree branches. Therefore, in general, an
 *    analysis will run much faster.
 * 3. Branch contents are accessed directly by accessing leafs. Thus, it is
 *    assumed that all leafs have unique names (equal to corresponding names of
 *    branches) which is indeed the case when each branch contains only one
 *    leaf.
 * 4. In the code, no protection against invalid code usage is provided.
 * 5. For integers, it is assumed that char = signed char = Char_t.
 * 6. Tree friends are not supported.
 */

#ifndef UNTUPLIZER_H
#define UNTUPLIZER_H

#include <map>
#include <string>
#include <vector>
#include <fnmatch.h>

#include <TFile.h>
#include <TList.h>
#include <TChain.h>
#include <TLeafB.h>
#include <TLeafC.h>
#include <TLeafD.h>
#include <TLeafF.h>
#include <TLeafI.h>
#include <TLeafL.h>
#include <TLeafO.h>
#include <TLeafS.h>
#include <TString.h>
#include <TSystem.h>
#include <TSystemFile.h>
#include <TLeafObject.h>
#include <TLeafElement.h>
#include <TSystemDirectory.h>

// prints a message and exits gracefully
#ifndef FATAL
#define FATAL(msg) do { fprintf(stderr, "FATAL: %s\n", msg); gSystem->Exit(1); } while (0)
#endif

class TreeReader {

 public:

   // types of branch/leaf contents
   enum ETypes {
      kBool,             // 1-byte boolean (TLeafO)
      kChar,             // 1-byte integer (TLeafB)
      kShort,            // 2-byte integer (TLeafS)
      kInt,              // 4-byte integer (TLeafI)
      kLong64,           // 8-byte integer (TLeafL)
      kFloat,            // 4-byte float (TLeafF)
      kDouble,           // 8-byte float (TLeafD)
      kString,           // variable length string (TLeafC)
      kArrBool,          // array of 1-byte booleans (TLeafO only)
      // NOTE: for integers, it is assumed that char = signed char = Char_t
      kArrChar,          // array of 1-byte integers (TLeafB or vector<char>)
      kArrCharTLeaf,     // array of 1-byte integers (TLeafB)
      kArrCharVector,    // array of 1-byte integers (vector<char>)
      kArrUCharVector,   // array of 1-byte unsigned integers (vector<unsigned char>)
      kArrShort,         // array of 2-byte integers (TLeafS or vector<short>)
      kArrShortTLeaf,    // array of 2-byte integers (TLeafS)
      kArrShortVector,   // array of 2-byte integers (vector<short>)
      kArrUShortVector,  // array of 2-byte unsigned integers (vector<unsigned short>)
      kArrInt,           // array of 4-byte integers (TLeafI or vector<int>)
      kArrIntTLeaf,      // array of 4-byte integers (TLeafI)
      kArrIntVector,     // array of 4-byte integers (vector<int>)
      kArrUIntVector,    // array of 4-byte unsigned integers (vector<unsigned int>)
      kArrLong64,        // array of 8-byte integers (TLeafL or vector<long>)
      kArrLong64TLeaf,   // array of 8-byte integers (TLeafL)
      kArrLong64Vector,  // array of 8-byte integers (vector<long>)
      kArrULong64Vector, // array of 8-byte unsigned integers (vector<unsigned long>)
      kArrFloat,         // array of 4-byte floats (TLeafF or vector<float>)
      kArrFloatTLeaf,    // array of 4-byte floats (TLeafF)
      kArrFloatVector,   // array of 4-byte floats (vector<float>)
      kArrDouble,        // array of 8-byte floats (TLeafD or vector<double>)
      kArrDoubleTLeaf,   // array of 8-byte floats (TLeafD)
      kArrDoubleVector,  // array of 8-byte floats (vector<double>)
      kArrVectInt,       // array of vector<int> (vector<vector<int> >)
      kArrVectFloat,     // array of vector<float> (vector<vector<float> >)
      kTObject,          // general object inherited from TObject
      kVoidPtr           // all other data types
   };

   // TTree/TChain initializers
   TreeReader(TTree* tree);
   TreeReader(const char* patt,               const char* treename = "ggNtuplizer/EventTree");
   TreeReader(const char** paths, int npaths, const char* treename = "ggNtuplizer/EventTree");
   TreeReader(std::vector<std::string> paths, const char* treename = "ggNtuplizer/EventTree");

   virtual ~TreeReader();

   // simple getters
   TTree*   GetTree()        { return fTree; }
   Long64_t GetEntriesFast() { return fTree->GetEntriesFast(); }
   Bool_t   HasMC()          { return fkMC; }

   // useful to determine which type of variable to use for which branch
   void Print();

   // sets event number to retrieve next time TreeReader::Get*() called
   void GetEntry(Long64_t entry);

   // returns either a pointer to first element of one- or multi-dimensional
   // array or a pointer to an unprocessed object.
   void* GetPtr(const char* branch_name, ETypes cktype = kVoidPtr, Int_t* nsize = NULL);

   // specializations of the call above
   // NOTE: for unsigned data types, use same methods with type casting
   Char_t*   GetPtrChar   (const char* bname) { return (Char_t*)    GetPtr(bname, kArrChar);   }
   Short_t*  GetPtrShort  (const char* bname) { return (Short_t*)   GetPtr(bname, kArrShort);  }
   Int_t*    GetPtrInt    (const char* bname) { return (Int_t*)     GetPtr(bname, kArrInt);    }
   Long64_t* GetPtrLong64 (const char* bname) { return (Long64_t*)  GetPtr(bname, kArrLong64); }
   Float_t*  GetPtrFloat  (const char* bname) { return (Float_t*)   GetPtr(bname, kArrFloat);  }
   Double_t* GetPtrDouble (const char* bname) { return (Double_t*)  GetPtr(bname, kArrDouble); }
   char*     GetPtrString (const char* bname) { return (char*)      GetPtr(bname, kString);    }
   TObject*  GetPtrTObject(const char* bname) { return (TObject*)   GetPtr(bname, kTObject);   }

   // NOTE: this works only for TLeafO (and not for vector<bool>) branches
   Bool_t*   GetPtrBool   (const char* bname) { return (Bool_t*) GetPtr(bname, kArrBool); }

   // return branch values for elementary types
   Bool_t   GetBool  (const char* bname) { return ((Bool_t*)   GetPtr(bname, kBool))  [0]; }
   Char_t   GetChar  (const char* bname) { return ((Char_t*)   GetPtr(bname, kChar))  [0]; }
   Short_t  GetShort (const char* bname) { return ((Short_t*)  GetPtr(bname, kShort)) [0]; }
   Int_t    GetInt   (const char* bname) { return ((Int_t*)    GetPtr(bname, kInt))   [0]; }
   Long64_t GetLong64(const char* bname) { return ((Long64_t*) GetPtr(bname, kLong64))[0]; }
   Float_t  GetFloat (const char* bname) { return ((Float_t*)  GetPtr(bname, kFloat)) [0]; }
   Double_t GetDouble(const char* bname) { return ((Double_t*) GetPtr(bname, kDouble))[0]; }

   // vector<vector<float> > and vector<vector<int> > tree branches
   std::vector<Float_t>* GetPtrVectorFloat(const char* bname, Int_t &nsize) {
      return (std::vector<Float_t>*) GetPtr(bname, kArrVectFloat, &nsize);
   }
   std::vector<Float_t>* GetPtrVectorFloat(const char* bname) {
      return (std::vector<Float_t>*) GetPtr(bname, kArrVectFloat);
   }
   std::vector<Int_t>* GetPtrVectorInt(const char* bname, Int_t &nsize) {
      return (std::vector<Int_t>*) GetPtr(bname, kArrVectInt, &nsize);
   }
   std::vector<Int_t>* GetPtrVectorInt(const char* bname) {
      return (std::vector<Int_t>*) GetPtr(bname, kArrVectInt);
   }

 protected:

   std::vector<std::string> FindFiles(const char* patt);

   void  InitTreeChain(std::vector<std::string> paths, const char* treename);
   void  FindLeaf(const char* bname);

   TFile*   fFile;     // file handle associated with fTree
   TTree*   fTree;     // pointer to TTree or TChain
   Int_t    fTreeNum;  // (for TChain) current tree number in a list of TTrees
   Bool_t   fkMC;      // if kTRUE, MC truth info is available
   Long64_t fEntry;    // fTree entry number to read with Get*() methods

   // caching
   std::map<std::string,size_t> fLeafIdx;    // leaf name => index in arrays below
   std::vector<TLeaf*>          fLeafAddr;   // cached leaf address
   std::vector<ETypes>          fLeafType;   // cached type of leaf contents
   std::vector<void*>           fLeafValue;  // cached payload address
};

//______________________________________________________________________________
TreeReader::TreeReader(TTree* tree) :
   fFile(0),
   fTree(0),
   fTreeNum(-1),
   fkMC(kFALSE),
   fEntry(-1)
{
   /* Associates an external TTree or TChain with this class.
    *
    * tree = reference to TTree or TChain. The caller is the owner of the
    * object.
    */

   // verify sizes of elementary data types
   if (sizeof(short) != 2 || sizeof(int) != 4 || sizeof(long) != 8 ||
       sizeof(float) != 4 || sizeof(double) != 8)
      FATAL("short/int/long/float/double of unsupported size");

   fTree = tree;

   // find out availability of MC truth info (check existence of "nMC" branch)
   fkMC = fTree->GetBranch("nMC") ? kTRUE : kFALSE;
}

//______________________________________________________________________________
TreeReader::TreeReader(const char* patt, const char* treename) :
   fFile(0),
   fTree(0),
   fTreeNum(-1),
   fkMC(kFALSE),
   fEntry(-1)
{
   /* Makes TTree or TChain based on the path pattern given.
    *
    * patt may be:
    *   - path to an existing root file (e.g. "/my/dir/file.root");
    *   - path to an existing directory (e.g "/my/dir" or "/my/dir/");
    *   - path wildcard to directories with root files (e.g. "/my/d*r");
    *   - path wildcard to root files (e.g. "/my/d*r/2*11/hgg*.root").
    *
    * Note: in the last case, patt must always end with ".root".
    *
    * Paths may be absolute or relative. Any root-supported file protocols
    * should work, e.g. "rfio:/home/2*1/123*3/abc*.root" is a valid path.
    *
    * All wildcard patterns from "man 7 glob" are supported. In short:
    *   - "*" matches any string, including the empty string;
    *   - "?" matches any single character;
    *   - "[...]" matches any character or ranges of characters inside the
    *     brackets, e.g. "[abc]" or "[a-Z]";
    *   - "[!...]" negates the matching above.
    *
    * Note: "*" and "?" can never match the "/" character, e.g. "foo*bar" never
    * matches to "foo/bar".
    *
    * patt is inspected in the following order:
    *
    * 1. If patt ends with ".root" and the corresponding path exists, it is
    *    taken as a TTree input.
    *
    * 2. Otherwise all files in all subdirectories matching the pattern are
    *    searched for and are taken as TTree/TChain input.
    *
    * Note: if patt is a directory which ends with ".root", add "/" at the end
    * to distinguish it from the case 1, e.g. "/path/to/dir.root/".
    *
    * Note: in patt, the very first directory component must not be a pattern.
    * E.g. the wildcard "12*34" is invalid, while "./12*34" is valid. Another
    * example: "rfio:reldir/file*.root" and "rfio:/abs*e/file*.root" are valid
    * wildcards, while "rfio:rel*e/file*.root" is not.
    */

   InitTreeChain(FindFiles(patt), treename);
}

//______________________________________________________________________________
TreeReader::TreeReader(const char** paths, int npaths, const char* treename) :
   fFile(0),
   fTree(0),
   fTreeNum(-1),
   fkMC(kFALSE),
   fEntry(-1)
{
   /* Makes TChain of trees from several root files.
    *
    * If input contains only one path, this path is considered as a wildcard,
    * see documentation of the TreeReader(const char*, const char*) constructor.
    *
    * paths = array of any root-supported paths to root files with TTrees;
    * npaths = number of elements in the array above.
    */

   std::vector<std::string> paths_;

   if (npaths == 1)
      paths_ = FindFiles(paths[0]);
   else
      for (int i = 0; i < npaths; i++)
         paths_.push_back(paths[i]);

   InitTreeChain(paths_, treename);
}

//______________________________________________________________________________
TreeReader::TreeReader(std::vector<std::string> paths, const char* treename) :
   fFile(0),
   fTree(0),
   fTreeNum(-1),
   fkMC(kFALSE),
   fEntry(-1)
{
   /* Makes TChain of trees from several root files.
    *
    * If input contains only one path, this path is considered as a wildcard,
    * see documentation of the TreeReader(const char*, const char*) constructor.
    *
    * paths = array of any root-supported paths to root files with TTrees.
    */

   if (paths.size() == 1)
      paths = FindFiles(paths[0].c_str());

   InitTreeChain(paths, treename);
}

//______________________________________________________________________________
TreeReader::~TreeReader()
{
   /* Frees memory.
    */

   if (fTree) delete fTree;
   if (fFile) delete fFile;
}

//______________________________________________________________________________
void TreeReader::Print()
{
   /* Prints to stdout branch names together with descriptions of their content
    * as to be used in an analysis code.
    *
    * Useful to determine which type of variable to use for which branch.
    *
    * NOTE: for branches containing vector<...> objects of elementary data types
    * or vector<...> objects of vector<int> or vector<float> data types, types
    * of underlying arrays are shown (as to be returned by TreeReader::Get*()
    * call family).
    */

   Printf("Branch names with content descriptions:");

   // access leafs directly
   TObjArray* leafs = fTree->GetListOfLeaves();

   // leaf loop
   for (int i = 0; i < leafs->GetEntriesFast(); i++) {
      TLeaf* leaf = dynamic_cast<TLeaf*>(leafs->At(i));
      if (!leaf)
         FATAL("leaf is NULL");

      // vector<...> and other objects not inherited from TObject
      if (leaf->IsA() == TLeafElement::Class()) {
         // content descriptor
         std::string descr(leaf->GetBranch()->GetClassName());

         // known one- and two-dimensional vector<...> arrays
         const char* types[] = {
            "float", "Float_t",
            "char", "Char_t", "signed char", "unsigned char", "UChar_t",
            "short", "short int", "signed short", "signed short int", "Short_t",
            "unsigned short", "unsigned short int", "UShort_t",
            "int", "signed", "signed int", "Int_t",
            "unsigned int", "unsigned", "UInt_t",
            "long", "long int", "signed long", "signed long int", "Long64_t", "Long_t",
            "unsigned long", "unsigned long int", "ULong64_t", "ULong_t",
            "vector<float> ", "vector<int> "};
         const char* type_descr[] = {
            "float", "float",
            "Char_t", "Char_t", "Char_t", "UChar_t", "UChar_t",
            "Short_t", "Short_t", "Short_t", "Short_t", "Short_t",
            "UShort_t", "UShort_t", "UShort_t",
            "Int_t", "Int_t", "Int_t", "Int_t",
            "UInt_t", "UInt_t", "UInt_t",
            "Long64_t", "Long64_t", "Long64_t", "Long64_t", "Long64_t", "Long64_t",
            "ULong64_t", "ULong64_t", "ULong64_t", "Long64_t",
            "vector<float>", "vector<Int_t>"};

         for (int c = 0; c < 34; c++)
            if (descr.compare(Form("vector<%s>", types[c])) == 0) {
               descr = Form("%s*", type_descr[c]);
               break;
            }

         Printf("  %-36s: %s", leaf->GetName(), descr.data());
      }

      // objects inherited from TObject
      else if (leaf->IsA() == TLeafObject::Class()) {
         Printf("  %-36s: %s", leaf->GetName(),
                               leaf->GetBranch()->GetClassName());
      }

      // elementary data types; fixed/variable length arrays of elementary types
      else {
         std::string descr;

         if (leaf->IsA() == TLeafB::Class())
            descr = "Char_t";
         else if (leaf->IsA() == TLeafD::Class())
            descr = "double";
         else if (leaf->IsA() == TLeafF::Class())
            descr = "float";
         else if (leaf->IsA() == TLeafI::Class())
            descr = "Int_t";
         else if (leaf->IsA() == TLeafL::Class())
            descr = "Long64_t";
         else if (leaf->IsA() == TLeafO::Class())
            descr = "Bool_t";
         else if (leaf->IsA() == TLeafS::Class())
            descr = "Short_t";
         else
            FATAL(Form("unsupported leaf of class %s", leaf->ClassName()));

         // variable length string
         if (leaf->IsA() == TLeafC::Class())
            Printf("  %-36s: char* (variable length string)", leaf->GetName());

         // single variable of elementary data type
         else if (!leaf->GetLeafCount() && leaf->GetLenStatic() == 1)
            Printf("  %-36s: %s", leaf->GetName(), descr.data());

         // fixed/variable length array
         else {
            // remove leaf name from the description
            std::string title(leaf->GetTitle());
            size_t ind = title.find('[');
            if (ind >= title.npos)
               FATAL("string::find('[') failed");

            Printf("  %-36s: %s%s", leaf->GetName(), descr.data(), &title[ind]);
         }
      }

   } // leaf loop
}

//______________________________________________________________________________
void TreeReader::GetEntry(Long64_t entry)
{
   /* Sets event number to retrieve next time TreeReader::Get*() called.
    */

   if (fTree->IsA() != TChain::Class())
      fEntry = entry;

   // TChain requires special treatment
   else {
      fEntry = ((TChain*)fTree)->LoadTree(entry);

      // reset caches on switching to new TTree
      if (fTreeNum != ((TChain*)fTree)->GetTreeNumber()) {
         fLeafIdx.clear();
         fLeafType.clear();
         fLeafAddr.clear();
         fLeafValue.clear();

         fTreeNum = ((TChain*)fTree)->GetTreeNumber();
      }
   }

   // reset cache of addresses to leaf payloads
   for (size_t i = 0; i < fLeafValue.size(); i++)
      fLeafValue[i] = NULL;
}

//______________________________________________________________________________
void* TreeReader::GetPtr(const char* branch_name, ETypes cktype, Int_t* nsize)
{
   /* Returns either a pointer to first element of one- or multi-dimensional
    * array or a pointer to an unprocessed object.
    *
    * In particular, for vector<...> arrays of elementary data types as well as
    * for vector<vector<int> > and vector<vector<float> > arrays,
    * vector::front() is returned. In this way, the code will work both with
    * TLeaf* tree branches as well as with vector<...> tree branches (provided
    * that underlying data types are the same in both cases).
    *
    * branch_name = name of branch to search for;
    * cktype = if not kVoidPtr, an additional type verification is performed;
    * nsize = if not NULL, filled with vector::size() for vector<...> branches.
    */

   // entry number in fLeafAddr, fLeafType and fLeafValue
   size_t i;

   // find the entry number
   std::map<std::string,size_t>::const_iterator got = fLeafIdx.find(branch_name);

   if (got != fLeafIdx.end())
      i = got->second;
   else {
      // this code is executed once per branch and per root file
      FindLeaf(branch_name);
      i = fLeafType.size() - 1;
   }

   // verify leaf type, if requested
   if (cktype != kVoidPtr) {
      if (cktype == kArrFloat) {
         if (fLeafType[i] != kArrFloatTLeaf &&
             fLeafType[i] != kArrFloatVector)
            FATAL(Form("branch is not of type float*: %s", branch_name));
      } else if (cktype == kArrInt) {
         if (fLeafType[i] != kArrIntTLeaf &&
             fLeafType[i] != kArrIntVector &&
             fLeafType[i] != kArrUIntVector)
            FATAL(Form("branch is not of type (U)Int_t*: %s", branch_name));
      } else if (cktype == kArrChar) {
         if (fLeafType[i] != kArrCharTLeaf &&
             fLeafType[i] != kArrCharVector &&
             fLeafType[i] != kArrUCharVector)
            FATAL(Form("branch is not of type (U)Char_t*: %s", branch_name));
      } else if (cktype == kArrShort) {
         if (fLeafType[i] != kArrShortTLeaf &&
             fLeafType[i] != kArrShortVector &&
             fLeafType[i] != kArrUShortVector)
            FATAL(Form("branch is not of type (U)Short_t*: %s", branch_name));
      } else if (cktype == kArrLong64) {
         if (fLeafType[i] != kArrLong64TLeaf &&
             fLeafType[i] != kArrLong64Vector &&
             fLeafType[i] != kArrULong64Vector)
            FATAL(Form("branch is not of type (U)Long64_t*: %s", branch_name));
      } else if (cktype == kArrDouble) {
         if (fLeafType[i] != kArrDoubleTLeaf &&
             fLeafType[i] != kArrDoubleVector)
            FATAL(Form("branch is not of type double*: %s", branch_name));
      } else if (cktype == kTObject) {
         if (fLeafType[i] != kTObject)
            FATAL(Form("branch content is not inherited from TObject: %s", branch_name));
      } else if (cktype == kString) {
         if (fLeafType[i] != kString)
            FATAL(Form("branch is not a variable length string: %s", branch_name));
      } else
         if (cktype != fLeafType[i])
            FATAL(Form("invalid branch type requested: %s", branch_name));
   }

   // load branch contents into memory, if necessary
   if (!fLeafValue[i]) {
      // NOTE: it is assumed here that the corresponding branch is not disabled
      TBranch* br = fLeafAddr[i]->GetBranch();
      if (!br)
         FATAL(Form("TLeaf::GetBranch() failed: %s", branch_name));
      if (br->GetEntry(fEntry) < 0)
         FATAL(Form("TBranch::GetEntry() failed: %s", branch_name));

      // pointer to actual leaf payload
      void* ptr = fLeafAddr[i]->GetValuePointer();
      if (fLeafType[i] == kTObject)
         ptr = *((void**)ptr);

      // cache address to payload
      if (fLeafType[i] == kArrFloatVector)
         fLeafValue[i] = &((std::vector<float>*)ptr)->front();
      else if (fLeafType[i] == kArrIntVector)
         fLeafValue[i] = &((std::vector<int>*)ptr)->front();
      else if (fLeafType[i] == kArrUIntVector)
         fLeafValue[i] = &((std::vector<unsigned int>*)ptr)->front();
      else if (fLeafType[i] == kArrCharVector)
         fLeafValue[i] = &((std::vector<char>*)ptr)->front();
      else if (fLeafType[i] == kArrUCharVector)
         fLeafValue[i] = &((std::vector<unsigned char>*)ptr)->front();
      else if (fLeafType[i] == kArrShortVector)
         fLeafValue[i] = &((std::vector<short>*)ptr)->front();
      else if (fLeafType[i] == kArrUShortVector)
         fLeafValue[i] = &((std::vector<unsigned short>*)ptr)->front();
      else if (fLeafType[i] == kArrLong64Vector)
         fLeafValue[i] = &((std::vector<long>*)ptr)->front();
      else if (fLeafType[i] == kArrULong64Vector)
         fLeafValue[i] = &((std::vector<unsigned long>*)ptr)->front();
      else if (fLeafType[i] == kArrDoubleVector)
         fLeafValue[i] = &((std::vector<double>*)ptr)->front();
      else
         fLeafValue[i] = ptr;
   }

   // special case of vector<vector<float> > and vector<vector<int> >
   if (fLeafType[i] == kArrVectFloat) {
      if (nsize)
         *nsize = (Int_t) ((std::vector<std::vector<float> >*)fLeafValue[i])->size();
      return &((std::vector<std::vector<float> >*)fLeafValue[i])->front();
   }
   else if (fLeafType[i] == kArrVectInt) {
      if (nsize)
         *nsize = (Int_t) ((std::vector<std::vector<int> >*)fLeafValue[i])->size();
      return &((std::vector<std::vector<int> >*)fLeafValue[i])->front();
   }

   return fLeafValue[i];
}

//______________________________________________________________________________
std::vector<std::string> TreeReader::FindFiles(const char* patt)
{
   // Returns list of paths (to root files) matching the path pattern given.

   std::vector<std::string> paths;

   TString p(patt);

   if (p.EndsWith(".root")) {
      // check whether target path already exists
      if (!gSystem->AccessPathName(patt)) {
         paths.push_back(patt);
         return paths;
      }
   } else
      // append the file search pattern
      p += (p.EndsWith("/") ? "*.root" : "/*.root");

   // find existing base directory
   TString base = p;
   TList* files;
   do {
      Ssiz_t i = TString(base(0, base.Length() - 1)).Last('/');
      if (i == kNPOS)
         FATAL("no input files found");

      base = TString(base(0, i + 1));  // keep all "/" at the end
      TSystemDirectory dir(base, base);
      files = dir.GetListOfFiles();
   } while (base.Length() > 0 && !files);

   if (base.Length() < 1)
      FATAL("no input files found");

   // get rid of all "//" in p beyond the base directory
   TString trail = p(base.Length(), p.Length());
   Ssiz_t siz = trail.Length();
   while (trail.ReplaceAll("//", "/").Length() != siz)
      siz = trail.Length();

   // NOTE: trail never starts with "/"
   p = base + trail;

   std::vector<std::string> dirs;

   dirs.push_back(base.Data());

   // walk though matching subdirectories, search for matching files
   while (dirs.size() > 0) {
      TString dirpath = dirs.back();
      dirs.pop_back();

      TSystemDirectory dir(dirpath, dirpath);
      TList* files = dir.GetListOfFiles();
      if (!files) continue;

      TIter next(files);
      while (TSystemFile* file = (TSystemFile*) next()) {
         TString fname = file->GetName();
         TString path = TString(dirpath) + fname;

         if (fname == "." || fname == "..") continue;

         if (file->IsDirectory()) {
            // search for a match
            Ssiz_t i = p.Last('/');
            while (i != kNPOS) {
               if (fnmatch(TString(p(0, i)), path, FNM_PATHNAME) == 0) {
                  dirs.push_back((path + "/").Data());
                  break;
               }
               i = TString(p(0, i)).Last('/');
            }

         } else if (fname.EndsWith(".root")) {
            if (fnmatch(p, path, FNM_PATHNAME) == 0) {
//                Printf("TreeReader: adding %s", path.Data());
               paths.push_back(path.Data());
            }
         }
      }

      delete files;
   }

   if (paths.size() < 1)
      FATAL("no input files found");

   Printf("TreeReader: %lu file%s found", paths.size(),
                                          paths.size() == 1 ? "" : "s");
   return paths;
}

//______________________________________________________________________________
void TreeReader::InitTreeChain(std::vector<std::string> paths, const char* treename)
{
   /* Makes TTree or TChain based on the input given.
    */

   // verify sizes of elementary data types
   if (sizeof(short) != 2 || sizeof(int) != 4 || sizeof(long) != 8 ||
       sizeof(float) != 4 || sizeof(double) != 8)
      FATAL("short/int/long/float/double of unsupported size");

   if (paths.size() == 1) {
      // current working directory
      TDirectory* wd = gDirectory;

      // open file with TTree
      fFile = TFile::Open(paths[0].c_str());
      if (!fFile || fFile->IsZombie())
         FATAL("TFile::Open() failed");

      // cd back into previous current working directory
      if (wd) wd->cd();
      else gDirectory = 0;

      // get tree
      fTree = dynamic_cast<TTree*>(fFile->Get(treename));
      if (!fTree)
         FATAL(Form("TTree \"%s\" not found:", treename));

      // be 100% sure: check explicitly object's class
      if (((TObject*)fTree)->IsA() != TTree::Class())
         FATAL(Form("\"%s\" is not a TTree", treename));

   } else {
      fTree = new TChain(treename);

      // add root files with TTrees, reading the number of entries in each file
      for (size_t i = 0; i < paths.size(); i++)
         if (((TChain*)fTree)->AddFile(paths[i].c_str(), 0) != 1)
            FATAL("TChain::AddFile() failed");
   }

   // find out availability of MC truth info (check existence of "nMC" branch)
   fkMC = fTree->GetBranch("nMC") ? kTRUE : kFALSE;
}

//______________________________________________________________________________
void TreeReader::FindLeaf(const char* bname)
{
   /* Finds leaf, determines type of its contents and fills the cache of leafs
    * accordingly.
    */

   // find leaf
   TLeaf* leaf = fTree->FindLeaf(bname);
   if (!leaf)
      FATAL(Form("leaf not found: %s", bname));

   // actual data type of the leaf's payload
   ETypes type;

   // vector<...> arrays and other objects not inherited from TObject
   if (leaf->IsA() == TLeafElement::Class()) {
      std::string descr(leaf->GetBranch()->GetClassName());

      if (descr.compare("vector<float>") == 0 ||
          descr.compare("vector<Float_t>") == 0)
         type = kArrFloatVector;
      else if (descr.compare("vector<double>") == 0 ||
          descr.compare("vector<Double_t>") == 0)
         type = kArrDoubleVector;
      else if (descr.compare("vector<int>") == 0 ||
               descr.compare("vector<signed>") == 0 ||
               descr.compare("vector<signed int>") == 0 ||
               descr.compare("vector<Int_t>") == 0)
         type = kArrIntVector;
      else if (descr.compare("vector<unsigned>") == 0 ||
               descr.compare("vector<unsigned int>") == 0 ||
               descr.compare("vector<UInt_t>") == 0)
         type = kArrUIntVector;
      else if (descr.compare("vector<char>") == 0 ||
               descr.compare("vector<Char_t>") == 0)
         type = kArrCharVector;
      else if (descr.compare("vector<unsigned char>") == 0 ||
               descr.compare("vector<UChar_t>") == 0)
         type = kArrUCharVector;
      else if (descr.compare("vector<short>") == 0 ||
               descr.compare("vector<short int>") == 0 ||
               descr.compare("vector<signed short>") == 0 ||
               descr.compare("vector<signed short int>") == 0 ||
               descr.compare("vector<Short_t>") == 0)
         type = kArrShortVector;
      else if (descr.compare("vector<unsigned short>") == 0 ||
               descr.compare("vector<unsigned short int>") == 0 ||
               descr.compare("vector<UShort_t>") == 0)
         type = kArrUShortVector;
      else if (descr.compare("vector<long>") == 0 ||
               descr.compare("vector<long int>") == 0 ||
               descr.compare("vector<signed long>") == 0 ||
               descr.compare("vector<signed long int>") == 0 ||
               descr.compare("vector<Long64_t>") == 0 ||
               descr.compare("vector<Long>") == 0)
         type = kArrLong64Vector;
      else if (descr.compare("vector<unsigned long>") == 0 ||
               descr.compare("vector<unsigned long int>") == 0 ||
               descr.compare("vector<ULong64_t>") == 0 ||
               descr.compare("vector<ULong_t>") == 0)
         type = kArrULong64Vector;
      else if (descr.compare("vector<vector<float> >") == 0 ||
               descr.compare("vector<vector<Float_t> >") == 0)
         type = kArrVectFloat;
      else if (descr.compare("vector<vector<int> >") == 0 ||
               descr.compare("vector<vector<signed> >") == 0 ||
               descr.compare("vector<vector<signed int> >") == 0 ||
               descr.compare("vector<vector<Int_t> >") == 0)
         type = kArrVectInt;
      else
         type = kVoidPtr;
   } // TLeafElement

   // objects inherited from TObject
   else if (leaf->IsA() == TLeafObject::Class())
      type = kTObject;

   // variable length strings
   else if (leaf->IsA() == TLeafC::Class())
      type = kString;

   // fixed/variable length arrays of elementary data types
   else if (leaf->GetLeafCount() || leaf->GetLenStatic() > 1) {
      if (leaf->IsA() == TLeafF::Class())
         type = kArrFloatTLeaf;
      else if (leaf->IsA() == TLeafI::Class())
         type = kArrIntTLeaf;
      else if (leaf->IsA() == TLeafB::Class())
         type = kArrCharTLeaf;
      else if (leaf->IsA() == TLeafS::Class())
         type = kArrShortTLeaf;
      else if (leaf->IsA() == TLeafL::Class())
         type = kArrLong64TLeaf;
      else if (leaf->IsA() == TLeafO::Class())
         type = kArrBool;
      else if (leaf->IsA() == TLeafD::Class())
         type = kArrDoubleTLeaf;
      else
         FATAL(Form("branch contains an unknown data type: %s", bname));
   }

   // single variables of elementary data types
   else if (!leaf->GetLeafCount() && leaf->GetLenStatic() == 1) {
      if (leaf->IsA() == TLeafF::Class())
         type = kFloat;
      else if (leaf->IsA() == TLeafI::Class())
         type = kInt;
      else if (leaf->IsA() == TLeafB::Class())
         type = kChar;
      else if (leaf->IsA() == TLeafS::Class())
         type = kShort;
      else if (leaf->IsA() == TLeafD::Class())
         type = kDouble;
      else if (leaf->IsA() == TLeafL::Class())
         type = kLong64;
      else if (leaf->IsA() == TLeafO::Class())
         type = kBool;
      else
         FATAL(Form("branch contains an unknown data type: %s", bname));
   }

   // unknown TLeaf variant
   else
      FATAL(Form("branch contains an unknown data type: %s", bname));

   // update the cache
   fLeafIdx[bname] = fLeafAddr.size();
   fLeafAddr.push_back(leaf);
   fLeafType.push_back(type);
   fLeafValue.push_back(NULL);
}

#endif
