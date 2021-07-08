// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME srcdIcustomPdfDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "GausCB.h"
#include "ModGaus.h"
#include "ModGaus11.h"
#include "ModGaus01.h"
#include "ModThreeGaus.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_GausCB(void *p = 0);
   static void *newArray_GausCB(Long_t size, void *p);
   static void delete_GausCB(void *p);
   static void deleteArray_GausCB(void *p);
   static void destruct_GausCB(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::GausCB*)
   {
      ::GausCB *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::GausCB >(0);
      static ::ROOT::TGenericClassInfo 
         instance("GausCB", ::GausCB::Class_Version(), "GausCB.h", 16,
                  typeid(::GausCB), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::GausCB::Dictionary, isa_proxy, 4,
                  sizeof(::GausCB) );
      instance.SetNew(&new_GausCB);
      instance.SetNewArray(&newArray_GausCB);
      instance.SetDelete(&delete_GausCB);
      instance.SetDeleteArray(&deleteArray_GausCB);
      instance.SetDestructor(&destruct_GausCB);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::GausCB*)
   {
      return GenerateInitInstanceLocal((::GausCB*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::GausCB*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_ModGaus(void *p = 0);
   static void *newArray_ModGaus(Long_t size, void *p);
   static void delete_ModGaus(void *p);
   static void deleteArray_ModGaus(void *p);
   static void destruct_ModGaus(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ModGaus*)
   {
      ::ModGaus *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ModGaus >(0);
      static ::ROOT::TGenericClassInfo 
         instance("ModGaus", ::ModGaus::Class_Version(), "ModGaus.h", 16,
                  typeid(::ModGaus), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ModGaus::Dictionary, isa_proxy, 4,
                  sizeof(::ModGaus) );
      instance.SetNew(&new_ModGaus);
      instance.SetNewArray(&newArray_ModGaus);
      instance.SetDelete(&delete_ModGaus);
      instance.SetDeleteArray(&deleteArray_ModGaus);
      instance.SetDestructor(&destruct_ModGaus);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ModGaus*)
   {
      return GenerateInitInstanceLocal((::ModGaus*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::ModGaus*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_ModGaus11(void *p = 0);
   static void *newArray_ModGaus11(Long_t size, void *p);
   static void delete_ModGaus11(void *p);
   static void deleteArray_ModGaus11(void *p);
   static void destruct_ModGaus11(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ModGaus11*)
   {
      ::ModGaus11 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ModGaus11 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("ModGaus11", ::ModGaus11::Class_Version(), "ModGaus11.h", 16,
                  typeid(::ModGaus11), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ModGaus11::Dictionary, isa_proxy, 4,
                  sizeof(::ModGaus11) );
      instance.SetNew(&new_ModGaus11);
      instance.SetNewArray(&newArray_ModGaus11);
      instance.SetDelete(&delete_ModGaus11);
      instance.SetDeleteArray(&deleteArray_ModGaus11);
      instance.SetDestructor(&destruct_ModGaus11);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ModGaus11*)
   {
      return GenerateInitInstanceLocal((::ModGaus11*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::ModGaus11*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_ModGaus01(void *p = 0);
   static void *newArray_ModGaus01(Long_t size, void *p);
   static void delete_ModGaus01(void *p);
   static void deleteArray_ModGaus01(void *p);
   static void destruct_ModGaus01(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ModGaus01*)
   {
      ::ModGaus01 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ModGaus01 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("ModGaus01", ::ModGaus01::Class_Version(), "ModGaus01.h", 16,
                  typeid(::ModGaus01), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ModGaus01::Dictionary, isa_proxy, 4,
                  sizeof(::ModGaus01) );
      instance.SetNew(&new_ModGaus01);
      instance.SetNewArray(&newArray_ModGaus01);
      instance.SetDelete(&delete_ModGaus01);
      instance.SetDeleteArray(&deleteArray_ModGaus01);
      instance.SetDestructor(&destruct_ModGaus01);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ModGaus01*)
   {
      return GenerateInitInstanceLocal((::ModGaus01*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::ModGaus01*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_ModThreeGaus(void *p = 0);
   static void *newArray_ModThreeGaus(Long_t size, void *p);
   static void delete_ModThreeGaus(void *p);
   static void deleteArray_ModThreeGaus(void *p);
   static void destruct_ModThreeGaus(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ModThreeGaus*)
   {
      ::ModThreeGaus *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ModThreeGaus >(0);
      static ::ROOT::TGenericClassInfo 
         instance("ModThreeGaus", ::ModThreeGaus::Class_Version(), "ModThreeGaus.h", 16,
                  typeid(::ModThreeGaus), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ModThreeGaus::Dictionary, isa_proxy, 4,
                  sizeof(::ModThreeGaus) );
      instance.SetNew(&new_ModThreeGaus);
      instance.SetNewArray(&newArray_ModThreeGaus);
      instance.SetDelete(&delete_ModThreeGaus);
      instance.SetDeleteArray(&deleteArray_ModThreeGaus);
      instance.SetDestructor(&destruct_ModThreeGaus);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ModThreeGaus*)
   {
      return GenerateInitInstanceLocal((::ModThreeGaus*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::ModThreeGaus*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr GausCB::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *GausCB::Class_Name()
{
   return "GausCB";
}

//______________________________________________________________________________
const char *GausCB::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::GausCB*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int GausCB::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::GausCB*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *GausCB::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::GausCB*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *GausCB::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::GausCB*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ModGaus::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ModGaus::Class_Name()
{
   return "ModGaus";
}

//______________________________________________________________________________
const char *ModGaus::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ModGaus*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ModGaus::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ModGaus*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ModGaus::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ModGaus*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ModGaus::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ModGaus*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ModGaus11::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ModGaus11::Class_Name()
{
   return "ModGaus11";
}

//______________________________________________________________________________
const char *ModGaus11::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ModGaus11*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ModGaus11::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ModGaus11*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ModGaus11::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ModGaus11*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ModGaus11::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ModGaus11*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ModGaus01::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ModGaus01::Class_Name()
{
   return "ModGaus01";
}

//______________________________________________________________________________
const char *ModGaus01::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ModGaus01*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ModGaus01::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ModGaus01*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ModGaus01::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ModGaus01*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ModGaus01::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ModGaus01*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ModThreeGaus::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ModThreeGaus::Class_Name()
{
   return "ModThreeGaus";
}

//______________________________________________________________________________
const char *ModThreeGaus::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ModThreeGaus*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ModThreeGaus::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ModThreeGaus*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ModThreeGaus::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ModThreeGaus*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ModThreeGaus::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ModThreeGaus*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void GausCB::Streamer(TBuffer &R__b)
{
   // Stream an object of class GausCB.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(GausCB::Class(),this);
   } else {
      R__b.WriteClassBuffer(GausCB::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_GausCB(void *p) {
      return  p ? new(p) ::GausCB : new ::GausCB;
   }
   static void *newArray_GausCB(Long_t nElements, void *p) {
      return p ? new(p) ::GausCB[nElements] : new ::GausCB[nElements];
   }
   // Wrapper around operator delete
   static void delete_GausCB(void *p) {
      delete ((::GausCB*)p);
   }
   static void deleteArray_GausCB(void *p) {
      delete [] ((::GausCB*)p);
   }
   static void destruct_GausCB(void *p) {
      typedef ::GausCB current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::GausCB

//______________________________________________________________________________
void ModGaus::Streamer(TBuffer &R__b)
{
   // Stream an object of class ModGaus.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ModGaus::Class(),this);
   } else {
      R__b.WriteClassBuffer(ModGaus::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ModGaus(void *p) {
      return  p ? new(p) ::ModGaus : new ::ModGaus;
   }
   static void *newArray_ModGaus(Long_t nElements, void *p) {
      return p ? new(p) ::ModGaus[nElements] : new ::ModGaus[nElements];
   }
   // Wrapper around operator delete
   static void delete_ModGaus(void *p) {
      delete ((::ModGaus*)p);
   }
   static void deleteArray_ModGaus(void *p) {
      delete [] ((::ModGaus*)p);
   }
   static void destruct_ModGaus(void *p) {
      typedef ::ModGaus current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ModGaus

//______________________________________________________________________________
void ModGaus11::Streamer(TBuffer &R__b)
{
   // Stream an object of class ModGaus11.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ModGaus11::Class(),this);
   } else {
      R__b.WriteClassBuffer(ModGaus11::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ModGaus11(void *p) {
      return  p ? new(p) ::ModGaus11 : new ::ModGaus11;
   }
   static void *newArray_ModGaus11(Long_t nElements, void *p) {
      return p ? new(p) ::ModGaus11[nElements] : new ::ModGaus11[nElements];
   }
   // Wrapper around operator delete
   static void delete_ModGaus11(void *p) {
      delete ((::ModGaus11*)p);
   }
   static void deleteArray_ModGaus11(void *p) {
      delete [] ((::ModGaus11*)p);
   }
   static void destruct_ModGaus11(void *p) {
      typedef ::ModGaus11 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ModGaus11

//______________________________________________________________________________
void ModGaus01::Streamer(TBuffer &R__b)
{
   // Stream an object of class ModGaus01.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ModGaus01::Class(),this);
   } else {
      R__b.WriteClassBuffer(ModGaus01::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ModGaus01(void *p) {
      return  p ? new(p) ::ModGaus01 : new ::ModGaus01;
   }
   static void *newArray_ModGaus01(Long_t nElements, void *p) {
      return p ? new(p) ::ModGaus01[nElements] : new ::ModGaus01[nElements];
   }
   // Wrapper around operator delete
   static void delete_ModGaus01(void *p) {
      delete ((::ModGaus01*)p);
   }
   static void deleteArray_ModGaus01(void *p) {
      delete [] ((::ModGaus01*)p);
   }
   static void destruct_ModGaus01(void *p) {
      typedef ::ModGaus01 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ModGaus01

//______________________________________________________________________________
void ModThreeGaus::Streamer(TBuffer &R__b)
{
   // Stream an object of class ModThreeGaus.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ModThreeGaus::Class(),this);
   } else {
      R__b.WriteClassBuffer(ModThreeGaus::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ModThreeGaus(void *p) {
      return  p ? new(p) ::ModThreeGaus : new ::ModThreeGaus;
   }
   static void *newArray_ModThreeGaus(Long_t nElements, void *p) {
      return p ? new(p) ::ModThreeGaus[nElements] : new ::ModThreeGaus[nElements];
   }
   // Wrapper around operator delete
   static void delete_ModThreeGaus(void *p) {
      delete ((::ModThreeGaus*)p);
   }
   static void deleteArray_ModThreeGaus(void *p) {
      delete [] ((::ModThreeGaus*)p);
   }
   static void destruct_ModThreeGaus(void *p) {
      typedef ::ModThreeGaus current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ModThreeGaus

namespace {
  void TriggerDictionaryInitialization_customPdfDict_Impl() {
    static const char* headers[] = {
"GausCB.h",
"ModGaus.h",
"ModGaus11.h",
"ModGaus01.h",
"ModThreeGaus.h",
0
    };
    static const char* includePaths[] = {
"inc",
"/Applications/root_v6.12.04/include",
"/Users/adorsett/Desktop/CERNbox/Zgamma/draw_pico/RooFit/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "customPdfDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(Your description goes here...)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$GausCB.h")))  GausCB;
class __attribute__((annotate(R"ATTRDUMP(Your description goes here...)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$ModGaus.h")))  ModGaus;
class __attribute__((annotate(R"ATTRDUMP(Your description goes here...)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$ModGaus11.h")))  ModGaus11;
class __attribute__((annotate(R"ATTRDUMP(Your description goes here...)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$ModGaus01.h")))  ModGaus01;
class __attribute__((annotate(R"ATTRDUMP(Your description goes here...)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$ModThreeGaus.h")))  ModThreeGaus;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "customPdfDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "GausCB.h"
#include "ModGaus.h"
#include "ModGaus11.h"
#include "ModGaus01.h"
#include "ModThreeGaus.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"GausCB", payloadCode, "@",
"ModGaus", payloadCode, "@",
"ModGaus01", payloadCode, "@",
"ModGaus11", payloadCode, "@",
"ModThreeGaus", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("customPdfDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_customPdfDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_customPdfDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_customPdfDict() {
  TriggerDictionaryInitialization_customPdfDict_Impl();
}
