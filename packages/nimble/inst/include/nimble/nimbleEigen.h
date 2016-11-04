#ifndef __NIMBLE_EIGEN
#define __NIMBLE_EIGEN

#include <iostream>
#include <cstdlib>
#include<Rmath.h>

// should the template arguments (types) be Map<MatrixXd> instead of MatrixXd?

// For more general: note each arg will need same scalar type which should be cast upward from logical->integer->double
// Could make a template class for each of 2-5 (say) arguments with a single type for each
// need to wrap constants

// concatenation, c(A1, A2)
template<bool useLinearAccess, typename result_type, typename eigenType, typename Index>
struct nimble_eigen_coeff_impl;

template<typename result_type, typename eigenType, typename Index>
struct nimble_eigen_coeff_impl<true, result_type, eigenType, Index> {
  static result_type getCoeff(const eigenType &Arg, Index i) {return Arg.coeff(i);}
};

template<typename result_type, typename eigenType, typename Index>
struct nimble_eigen_coeff_impl<false, result_type, eigenType, Index> {
  static result_type getCoeff(const eigenType &Arg, Index i) {
    std::div_t divRes = div(i, Arg.rows());
    return Arg.coeff(divRes.rem, floor(divRes.quot));
  }
};

template<typename Derived1, typename Derived2>
class concatenateClass {
 public:
  const Derived1 &Arg1;
  const Derived2 &Arg2;
  int size1, size2;
  typedef double result_type;
 concatenateClass(const Derived1 &A1, const Derived2 &A2) : Arg1(A1), Arg2(A2) {
    size1 = Arg1.size();
    size2 = Arg2.size();
  };
  // Index types should come in as template argument based on DerivedReturn
  typedef typename Eigen::internal::traits<Derived1>::Index Index1;
  typedef typename Eigen::internal::traits<Derived2>::Index Index2;
  
  result_type operator()(Index1 i) const //Eigen::DenseIndex //Assume Index1 type and Index2 type will always be the same, or cast-able.
  {
    std::cout<<"IN 1\n";
    if(i < size1)
      return nimble_eigen_coeff_impl< Eigen::internal::traits<Derived1>::Flags & LinearAccessBit, result_type, Derived1, Index1 >::getCoeff(Arg1, i); //generalization of Arg1(i) or Arg1.coeff(i) 
    else
      return nimble_eigen_coeff_impl< Eigen::internal::traits<Derived2>::Flags & LinearAccessBit, result_type, Derived2, Index1 >::getCoeff(Arg2, i - size1); //Arg2(i - size1);
  }

  result_type operator()(Index1 i, Index2 j) const // I don't think this should normally be called, but if it does, act like a vector
  {
    std::cout<<"IN 2\n";
    return operator()(i);
  }
};


namespace Eigen{
  namespace internal{
    template<typename Derived1, typename Derived2>
    // 3 options for linear_access:
      //1. turn it off: struct functor_has_linear_access<gConcatenateClass<Derived1, Derived2> > { enum { ret = 0 }; }; 
      //2. turn it on: (and let it be resolved by nimble_eigen_coeff_impl>:
    struct functor_has_linear_access<concatenateClass<Derived1, Derived2> > { enum { ret = 1}; }; 
      //3. Set it once according to arguments: (problem here is that if it is off for this expression, that may make expressions using this one forbidden from using linear access, which is a bit harsh: struct functor_has_linear_access<gConcatenateClass<Derived1, Derived2> > { enum { ret = traits<Derived1>::Flags & traits<Derived2>::Flags & LinearAccessBit }; };
      template<typename Derived1, typename Derived2>
      struct functor_traits<concatenateClass<Derived1, Derived2> >
      {
	enum
	{
	  Cost = 10, // there are templated costs available to pick up for this
	  PacketAccess = false, // happy to keep this false for now
	  IsRepeatable = true // default was false. 
	};
      };    
  }
}

template<typename returnDerived>
struct concatenate_impl {
  template<typename Derived1, typename Derived2>
    static CwiseNullaryOp<concatenateClass<const Derived1, const Derived2>, returnDerived > concatenate(const Derived1 &A1, const Derived2 &A2) {
    concatenateClass<const Derived1, const Derived2> c(A1.derived(), A2.derived());
    return(CwiseNullaryOp<concatenateClass<const Derived1, const Derived2>, returnDerived >(A1.size() + A2.size(), 1, c));
  }
};

#define nimCd concatenate_impl<MatrixXd>::concatenate
#define nimCi concatenate_impl<MatrixXi>::concatenate
#define nimCb concatenate_imple<MatrixXb>::concatenate

// rep, rep(x, times, each)
template<typename Derived1>
class repClass {
public:
  const Derived1 &Arg1;
  int times, each;
  typedef enum{useEach, useTimes, useBoth} eachTimesCaseType;
  eachTimesCaseType eachTimesCase;
  int sizeArg1;
  typedef double result_type;
 repClass(const Derived1 &A1, int timesIn, int eachIn) :
  Arg1(A1),
    times(timesIn),
    each(eachIn) {
    // add more complete checking of row vs. col vector here
    // assume col vector for now
    sizeArg1 = Arg1.size();
    if(each > 1)
      if(times > 1)
	eachTimesCase = useBoth;
      else
	eachTimesCase = useEach;
    else
      eachTimesCase = useTimes;
    // may need to handle times = 0 or each = 0
  };
  typedef typename Eigen::internal::traits<Derived1>::Index Index1;

  result_type operator()(Index1 i) const //Eigen::DenseIndex
  {
    Index1 iUse;
    std::cout<<"IN 1 rep\n";
    switch(eachTimesCase) {
    case useEach:
      iUse = floor(i/each);
      break;
    case useTimes:
      iUse = i % sizeArg1;
      break;
    case useBoth:
      iUse = static_cast<int>(floor(i / each)) % sizeArg1;
      break;
    default:
      iUse = 0; //error
    }
    return nimble_eigen_coeff_impl< Eigen::internal::traits<Derived1>::Flags & LinearAccessBit, result_type, Derived1, Index1 >::getCoeff(Arg1, iUse);
  }

  result_type operator()(Index1 i, Index1 j) const
  {
    std::cout<<"IN 2 rep\n";
    return operator()(i);
  }
};

namespace Eigen{
  namespace internal{
    template<typename Derived1>
    struct functor_has_linear_access<repClass<Derived1> > { enum { ret = 1}; }; 
      template<typename Derived1>
      struct functor_traits<repClass<Derived1> >
      {
	enum
	{
	  Cost = 10, // there are templated costs available to pick up for this
	  PacketAccess = false, // happy to keep this false for now
	  IsRepeatable = true // default was false. 
	};
      };    
  }
}

template<typename returnDerived>
struct rep_impl {
  template<typename Derived1>
  static CwiseNullaryOp<repClass<const Derived1>, returnDerived > rep(const Derived1 &A1, int reps, int each) {
    repClass<const Derived1> repObj(A1, reps, each);
    return(CwiseNullaryOp<repClass<const Derived1>, returnDerived >(A1.rows() * reps * each, 1, repObj));
  }
};

#define nimRepd rep_impl<MatrixXd>::rep
#define nimRepi rep_impl<MatrixXi>::rep
#define nimRepb rep_impl<MatrixXb>::rep

// sequences, seq(from, to, by, length.out)  not implementing along.with for now
// not implementing this for `:` because we already have special treatment `:` in for-loop context, so want to tread carefully there

typedef enum{useBy, useLength} byOrLength;

// to-do: make native integer handling where appropriate
// right now from and to are doubles
template<typename DerivedOut, typename scalarFrom, typename scalarTo, typename scalarBy, byOrLength>
  class seqClass;


template<typename DerivedOut, typename scalarFrom, typename scalarTo, typename scalarBy>
  class seqClass<DerivedOut, scalarFrom, scalarTo, scalarBy, useBy> {
public:
  scalarFrom from;
  scalarBy by;
  unsigned int length_out;
  //  typedef enum{useBy, useLength} byOrLength;
  
  typedef double result_type; // need to pull this from DerivedOut
 seqClass(scalarFrom fromIn, scalarTo toIn, scalarBy byIn) :
  from(fromIn),
    by(byIn) {
    printf("Add some checking to seqClass constructor and deal with inconsistent scalar types\n");
      length_out = 1 + static_cast<int>(floor(static_cast<double>(toIn) - static_cast<double>(from)) / static_cast<double>(byIn));
    };
  
  typedef typename Eigen::internal::traits<DerivedOut>::Index Index;
  result_type operator()(Index i) const //Eigen::DenseIndex
  {
    std::cout<<"IN 1 seq\n";
    return( from + static_cast<int>(i) * by );
  }
  result_type operator()(Index i, Index j) const
  {
    std::cout<<"IN 2 seq\n";
    if(j != 0) printf("Problem calling seq in C++ with two indices\n");
    return( from + static_cast<int>(i) * by);
  }
};

template<typename DerivedOut, typename scalarFrom, typename scalarTo>
  class seqClass<DerivedOut, scalarFrom, scalarTo, int, useLength> {
public:
  scalarFrom from;
  double by;
  unsigned int length_out;
  //  typedef enum{useBy, useLength} byOrLength;
  
  typedef double result_type; // need to pull this from DerivedOut
 seqClass(scalarFrom fromIn, scalarTo toIn, unsigned int length_outIn) :
  from(fromIn),
    length_out(length_outIn)
    {
      printf("Add some checking to seqClass constructor and deal with inconsistent scalar types\n");
      if(length_out == 1) {
	std::cout<<"setting by to 0\n";
	by = 0;
      } else {
	std::cout<<static_cast<double>(toIn)<<" "<<static_cast<double>(from)<<" "<<static_cast<double>(length_out)<<"\n";
	by = (static_cast<double>(toIn) - static_cast<double>(from)) / (static_cast<double>(length_out) - 1.);
      }
    };

  
  typedef typename Eigen::internal::traits<DerivedOut>::Index Index;
  result_type operator()(Index i) const //Eigen::DenseIndex
  {
    std::cout<<"IN 1 seq\n";
    return( from + static_cast<int>(i) * by );
  }
  result_type operator()(Index i, Index j) const
  {
    std::cout<<"IN 2 seq\n";
    if(j != 0) printf("Problem calling seq in C++ with two indices\n");
    return( from + static_cast<int>(i) * by);
  }
};


namespace Eigen{
  namespace internal{
    template<typename DerivedOut, typename scalarFrom, typename scalarTo, typename scalarBy>
      struct functor_has_linear_access<seqClass<DerivedOut, scalarFrom, scalarTo, scalarBy, useBy> > { enum { ret = 1}; }; 
    template<typename DerivedOut, typename scalarFrom, typename scalarTo, typename scalarBy>
      struct functor_traits<seqClass<DerivedOut, scalarFrom, scalarTo, scalarBy, useBy> >
      {
	enum
	{
	  Cost = 10, // there are templated costs available to pick up for this
	  PacketAccess = false, // happy to keep this false for now
	  IsRepeatable = true // default was false. 
	};
      };    
  
  template<typename DerivedOut, typename scalarFrom, typename scalarTo, typename scalarBy>
    struct functor_has_linear_access<seqClass<DerivedOut, scalarFrom, scalarTo, scalarBy, useLength> > { enum { ret = 1}; }; 
  template<typename DerivedOut, typename scalarFrom, typename scalarTo, typename scalarBy>
    struct functor_traits<seqClass<DerivedOut, scalarFrom, scalarTo, scalarBy, useLength> >
      {
	enum
	{
	  Cost = 10, // there are templated costs available to pick up for this
	  PacketAccess = false, // happy to keep this false for now
	  IsRepeatable = true // default was false. 
	};
      };    
  }
}


template<typename DerivedOut>
struct seq_impl {
  template<typename scalarFrom, typename scalarTo, typename scalarBy>
    static CwiseNullaryOp<seqClass<DerivedOut, scalarFrom, scalarTo, scalarBy, useBy>, DerivedOut > seqBy(scalarFrom from, scalarTo to, scalarBy by, unsigned int len) {
    seqClass<DerivedOut, scalarFrom, scalarTo, scalarBy, useBy> seqObj(from, to, by);
    return(CwiseNullaryOp<seqClass<DerivedOut, scalarFrom, scalarTo, scalarBy, useBy> , DerivedOut >(seqObj.length_out, 1, seqObj));
  }
  template<typename scalarFrom, typename scalarTo, typename scalarBy>
    static CwiseNullaryOp<seqClass<DerivedOut, scalarFrom, scalarTo, scalarBy, useLength>, DerivedOut > seqLen(scalarFrom from, scalarTo to, scalarBy by, unsigned int len) {
    seqClass<DerivedOut, scalarFrom, scalarTo, scalarBy, useLength> seqObj(from, to, len);
    return(CwiseNullaryOp<seqClass<DerivedOut, scalarFrom, scalarTo, scalarBy, useLength> , DerivedOut >(seqObj.length_out, 1, seqObj));
  }
};

#define nimSeqByD seq_impl<MatrixXd>::seqBy
#define nimSeqLenD seq_impl<MatrixXd>::seqLen
#define nimSeqByI seq_impl<MatrixXi>::seqBy
#define nimSeqLenI seq_impl<MatrixXi>::seqLen

// nonseqIndexed

template<typename DerivedObj, typename DerivedI1, typename DerivedI2>
class nonseqIndexedClass {
 public:
  const DerivedObj &obj;
  const DerivedI1 &index1;
  const DerivedI2 &index2;
  int dim1, dim2;
  typedef double result_type;
 nonseqIndexedClass(const DerivedObj &s, const DerivedI1 &i1, const DerivedI2 &i2) :
  obj(s),
    index1(i1),
    index2(i2) {
      dim1 = i1.size();
      dim2 = i2.size();
    }
  typedef typename Eigen::internal::traits<DerivedObj>::Index IndexObj;

  result_type operator()(IndexObj i) const //Eigen::DenseIndex
  {
    std::cout<<"IN 1\n";
    std::div_t divRes = div(i, dim1);
    return obj.coeff(nimble_eigen_coeff_impl< Eigen::internal::traits<DerivedI1>::Flags & LinearAccessBit, result_type, DerivedI1, typename Eigen::internal::traits<DerivedI2>::Index >::getCoeff(index1, divRes.rem) - 1,
		     nimble_eigen_coeff_impl< Eigen::internal::traits<DerivedI2>::Flags & LinearAccessBit, result_type, DerivedI2, typename Eigen::internal::traits<DerivedI2>::Index >::getCoeff(index2, floor(divRes.quot)) - 1); // This type of the index argument is confusing.  What is being passed is a type from std::div_t, which ought to be castable to any Eigen Index type I hope.
    //index1(divRes.rem)-1, index2(floor(divRes.quot))-1);
  }
  result_type operator()(IndexObj i, IndexObj j) const
  {
    std::cout<<"IN 2\n";
    return obj.coeff(nimble_eigen_coeff_impl< Eigen::internal::traits<DerivedI1>::Flags & LinearAccessBit, result_type, DerivedI1, typename Eigen::internal::traits<DerivedI2>::Index >::getCoeff(index1, i) - 1,
		     nimble_eigen_coeff_impl< Eigen::internal::traits<DerivedI2>::Flags & LinearAccessBit, result_type, DerivedI2, typename Eigen::internal::traits<DerivedI2>::Index >::getCoeff(index2, j) - 1);

    //return obj.coeff(index1(i)-1,
    //		     index2(j)-1);
  }
};

namespace Eigen{
  namespace internal{
    template<typename DerivedObj, typename Derived1, typename Derived2>
    // 3 options for linear_access:
      //1. turn it off: struct functor_has_linear_access<gConcatenateClass<Derived1, Derived2> > { enum { ret = 0 }; }; 
      //2. turn it on: (and let it be resolved by nimble_eigen_coeff_impl>:
      struct functor_has_linear_access<nonseqIndexedClass<DerivedObj, Derived1, Derived2> > { enum { ret = 1}; }; 
    //3. Set it once according to arguments: (problem here is that if it is off for this expression, that may make expressions using this one forbidden from using linear access, which is a bit harsh: struct functor_has_linear_access<gConcatenateClass<Derived1, Derived2> > { enum { ret = traits<Derived1>::Flags & traits<Derived2>::Flags & LinearAccessBit }; };
    template<typename DerivedObj, typename Derived1, typename Derived2>
      struct functor_traits<nonseqIndexedClass<DerivedObj, Derived1, Derived2> >
      {
	enum
	{
	  Cost = 10, // there are templated costs available to pick up for this
	  PacketAccess = false, // happy to keep this false for now
	  IsRepeatable = true // default was false. 
	};
      };    
  }
}

template<typename returnDerived>
struct nonseqIndexed_impl {
  template<typename DerivedObj, typename DerivedI1, typename DerivedI2>
    static CwiseNullaryOp<nonseqIndexedClass<DerivedObj, DerivedI1, DerivedI2 >, returnDerived > nonseqIndexed(const DerivedObj &s, const DerivedI1 &i1, const DerivedI2 &i2) {
    nonseqIndexedClass<DerivedObj, DerivedI1, DerivedI2 > nonseqIndexedObj(s.derived(), i1.derived(), i2.derived());
    return(CwiseNullaryOp<nonseqIndexedClass<DerivedObj, DerivedI1, DerivedI2 >, returnDerived >(nonseqIndexedObj.dim1, nonseqIndexedObj.dim2, nonseqIndexedObj));
  }
};

#define nimNonseqIndexedd nonseqIndexed_impl<MatrixXd>::nonseqIndexed
#define nimNonseqIndexedi nonseqIndexed_impl<MatrixXi>::nonseqIndexed
#define nimNonseqIndexedb nonseqIndexed_imple<MatrixXb>::nonseqIndexed

// vectorization of any scalar function: R's so-called "Recycling Rule"

// wrap access to Eigen's traits::..::LinearAccessBit so we can proxy it with true for a scalar type (double, int, bool)
template<typename T>
struct nimble_eigen_traits {
  enum {LinearAccessBit = int(Eigen::internal::traits<T>::Flags & LinearAccessBit)};
};

template<>
struct nimble_eigen_traits<double> {
  enum {LinearAccessBit = int(1)};
};

template<>
struct nimble_eigen_traits<int> {
  enum {LinearAccessBit = int(1)};
};

template<>
struct nimble_eigen_traits<bool> {
  enum {LinearAccessBit = int(1)};
};

// put the call to arg.size() in a struct so we can proxy it with "1" for a scalar type
template<typename T>
struct nimble_size_impl {
  static unsigned int getSize(const T &arg) {return arg.size();}
};

template<>
struct nimble_size_impl<double> {
  static unsigned int getSize(const double &arg) {return 1;}
};

template<>
struct nimble_size_impl<int> {
  static unsigned int getSize(const int &arg) {return 1;}
};

template<>
struct nimble_size_impl<bool> {
  static unsigned int getSize(const bool &arg) {return 1;}
};

// put the call to arg.coeff(...) in a struct so we can proxy it when LinearAccessBit is false and we can proxy it for a scalar 
template<bool useLinearAccess, typename result_type, typename eigenType, typename Index>
struct nimble_eigen_coeff_mod_impl;

template<typename result_type, typename eigenType, typename Index>
struct nimble_eigen_coeff_mod_impl<true, result_type, eigenType, Index> {
  static result_type getCoeff(const eigenType &Arg, Index i, unsigned int size) {return Arg.coeff(i % size);}
};

template<typename result_type, typename eigenType, typename Index>
struct nimble_eigen_coeff_mod_impl<false, result_type, eigenType, Index> {
  static result_type getCoeff(const eigenType &Arg, Index i, unsigned int size) {
    std::div_t divRes = div(i % size, Arg.rows());
    return Arg.coeff(divRes.rem, floor(divRes.quot));
  }  
};

template<typename result_type, typename Index>
struct nimble_eigen_coeff_mod_impl<true, result_type, double, Index> {
  static result_type getCoeff(const double Arg, Index i, unsigned int size) {return Arg;}
};

template<typename result_type, typename Index>
struct nimble_eigen_coeff_mod_impl<true, result_type, int, Index> {
  static result_type getCoeff(const int Arg, Index i, unsigned int size) {return Arg;}
};

template<typename result_type, typename Index>
struct nimble_eigen_coeff_mod_impl<true, result_type, bool, Index> {
  static result_type getCoeff(const bool Arg, Index i, unsigned int size) {return Arg;}
};

// Here is a test function
double RRtest_add(double a1, double a2) {std::cout<<a1<<" "<<a2<<" "<<a1+a2<<"\n"; return a1 + a2;}

// Here is the large macro for creating a functor class to be used in a NullaryExpr
#define MAKE_RECYCLING_RULE_CLASS2(FUNNAME, RETURNSCALARTYPE) \
template<typename Index, typename DerivedA1, typename DerivedA2> \
class FUNNAME ## RecyclingRuleClass { \
public: \
  const DerivedA1 &Arg1;\
  const DerivedA2 &Arg2;\
  unsigned int size1, size2, outputSize;\
  FUNNAME ## RecyclingRuleClass(const DerivedA1 &A1, const DerivedA2 &A2 ) :\
  Arg1(A1), Arg2(A2)\
{\
  outputSize = size1 = nimble_size_impl<DerivedA1>::getSize(Arg1); \
  size2 = nimble_size_impl<DerivedA2>::getSize(Arg2); \
  if(size2 > outputSize) outputSize = size2; \
  } \
  RETURNSCALARTYPE operator()(Index i) const { \
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::LinearAccessBit), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::LinearAccessBit), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2)); \
  }\
  RETURNSCALARTYPE operator()(Index i, Index j) const {\
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::LinearAccessBit), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::LinearAccessBit), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2)); \
  }\
}; \
\
 namespace Eigen{\
  namespace internal{\
    template<typename Index, typename Derived1, typename Derived2>	\
    struct functor_has_linear_access<FUNNAME ## RecyclingRuleClass<Index, Derived1, Derived2> > { enum { ret = 1}; }; \
    template<typename Index, typename Derived1, typename Derived2>\
      struct functor_traits<FUNNAME ## RecyclingRuleClass<Index, Derived1, Derived2> > { enum {Cost = 10, PacketAccess = false, IsRepeatable = true}; }; \
  }\
}\
\
template<typename DerivedReturn>\
struct FUNNAME ## _RR_impl {\
  typedef typename Eigen::internal::traits<DerivedReturn>::Index IndexReturn;\
  template<typename Derived1, typename Derived2>\
  static CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2>, DerivedReturn >\
  FUNNAME ## _RecyclingRule(const Derived1 &A1, const Derived2 &A2) {\
    FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2> obj(A1, A2);\
    return(CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2>, DerivedReturn >(obj.outputSize, 1, obj));\
  }\
};

// Here is the large macro for creating a functor class to be used in a NullaryExpr
#define MAKE_RECYCLING_RULE_CLASS3(FUNNAME, RETURNSCALARTYPE) \
  template<typename Index, typename DerivedA1, typename DerivedA2, typename DerivedA3>	\
class FUNNAME ## RecyclingRuleClass { \
public: \
  const DerivedA1 &Arg1;\
  const DerivedA2 &Arg2;\
  const DerivedA3 &Arg3;\
  unsigned int size1, size2, size3, outputSize;				\
  FUNNAME ## RecyclingRuleClass(const DerivedA1 &A1, const DerivedA2 &A2, const DerivedA3 &A3 ) : \
  Arg1(A1), Arg2(A2), Arg3(A3) \
{\
  outputSize = size1 = nimble_size_impl<DerivedA1>::getSize(Arg1); \
  size2 = nimble_size_impl<DerivedA2>::getSize(Arg2); \
  size3 = nimble_size_impl<DerivedA3>::getSize(Arg3); \
  if(size2 > outputSize) outputSize = size2; \
  if(size3 > outputSize) outputSize = size3; \
  } \
  RETURNSCALARTYPE operator()(Index i) const { \
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::LinearAccessBit), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::LinearAccessBit), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA3>::LinearAccessBit), RETURNSCALARTYPE, DerivedA3, Index >::getCoeff(Arg3, i, size3)); \
  }\
  RETURNSCALARTYPE operator()(Index i, Index j) const {\
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::LinearAccessBit), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::LinearAccessBit), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA3>::LinearAccessBit), RETURNSCALARTYPE, DerivedA3, Index >::getCoeff(Arg3, i, size3)); \
  }\
}; \
\
 namespace Eigen{\
  namespace internal{\
    template<typename Index, typename Derived1, typename Derived2, typename Derived3>	\
      struct functor_has_linear_access<FUNNAME ## RecyclingRuleClass<Index, Derived1, Derived2, Derived3> > { enum { ret = 1}; }; \
    template<typename Index, typename Derived1, typename Derived2, typename Derived3>	\
      struct functor_traits<FUNNAME ## RecyclingRuleClass<Index, Derived1, Derived2, Derived3> > { enum {Cost = 10, PacketAccess = false, IsRepeatable = true}; }; \
  }\
}\
\
template<typename DerivedReturn>\
struct FUNNAME ## _RR_impl {\
  typedef typename Eigen::internal::traits<DerivedReturn>::Index IndexReturn;\
  template<typename Derived1, typename Derived2, typename Derived3>			\
    static CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3>, DerivedReturn > \
    FUNNAME ## _RecyclingRule(const Derived1 &A1, const Derived2 &A2, const Derived3 &A3) { \
    FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3> obj(A1, A2, A3); \
    return(CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3>, DerivedReturn >(obj.outputSize, 1, obj)); \
  }\
};

// Here is the large macro for creating a functor class to be used in a NullaryExpr
#define MAKE_RECYCLING_RULE_CLASS4(FUNNAME, RETURNSCALARTYPE) \
  template<typename Index, typename DerivedA1, typename DerivedA2, typename DerivedA3, typename DerivedA4>	\
class FUNNAME ## RecyclingRuleClass { \
public: \
  const DerivedA1 &Arg1;\
  const DerivedA2 &Arg2;\
  const DerivedA3 &Arg3;\
  const DerivedA4 &Arg4;\
  unsigned int size1, size2, size3, size4, outputSize;			\
  FUNNAME ## RecyclingRuleClass(const DerivedA1 &A1, const DerivedA2 &A2, const DerivedA3 &A3, const DerivedA4 &A4 ) : \
  Arg1(A1), Arg2(A2), Arg3(A3), Arg4(A4) \
{\
  outputSize = size1 = nimble_size_impl<DerivedA1>::getSize(Arg1); \
  size2 = nimble_size_impl<DerivedA2>::getSize(Arg2); \
  size3 = nimble_size_impl<DerivedA3>::getSize(Arg3); \
  size4 = nimble_size_impl<DerivedA4>::getSize(Arg4); \
  if(size2 > outputSize) outputSize = size2; \
  if(size3 > outputSize) outputSize = size3; \
  if(size4 > outputSize) outputSize = size4; \
  } \
  RETURNSCALARTYPE operator()(Index i) const { \
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::LinearAccessBit), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::LinearAccessBit), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA3>::LinearAccessBit), RETURNSCALARTYPE, DerivedA3, Index >::getCoeff(Arg3, i, size3), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA4>::LinearAccessBit), RETURNSCALARTYPE, DerivedA4, Index >::getCoeff(Arg4, i, size4)); \
  }\
  RETURNSCALARTYPE operator()(Index i, Index j) const {\
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::LinearAccessBit), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::LinearAccessBit), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA3>::LinearAccessBit), RETURNSCALARTYPE, DerivedA3, Index >::getCoeff(Arg3, i, size3), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA4>::LinearAccessBit), RETURNSCALARTYPE, DerivedA4, Index >::getCoeff(Arg4, i, size4)); \
  }\
}; \
\
 namespace Eigen{\
  namespace internal{\
    template<typename Index, typename Derived1, typename Derived2, typename Derived3, typename Derived4> \
      struct functor_has_linear_access<FUNNAME ## RecyclingRuleClass<Index, Derived1, Derived2, Derived3, Derived4> > { enum { ret = 1}; }; \
    template<typename Index, typename Derived1, typename Derived2, typename Derived3, typename Derived4> \
      struct functor_traits<FUNNAME ## RecyclingRuleClass<Index, Derived1, Derived2, Derived3, Derived4> > { enum {Cost = 10, PacketAccess = false, IsRepeatable = true}; }; \
  }\
}\
\
template<typename DerivedReturn>\
struct FUNNAME ## _RR_impl {\
  typedef typename Eigen::internal::traits<DerivedReturn>::Index IndexReturn;\
  template<typename Derived1, typename Derived2, typename Derived3, typename Derived4>	\
    static CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3, Derived4>, DerivedReturn > \
    FUNNAME ## _RecyclingRule(const Derived1 &A1, const Derived2 &A2, const Derived3 &A3, const Derived4 &A4) { \
    FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3, Derived4> obj(A1, A2, A3, A4); \
    return(CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3, Derived4>, DerivedReturn >(obj.outputSize, 1, obj)); \
  }\
};


MAKE_RECYCLING_RULE_CLASS2(RRtest_add, double) // only return type is needed here, and correct number of arguments as last digit of macro name

// number of arguments is always 2 more than number of distribution parameters (1 for x and 1 for log)
MAKE_RECYCLING_RULE_CLASS4(dbinom, double)
MAKE_RECYCLING_RULE_CLASS3(dexp_nimble, double)
MAKE_RECYCLING_RULE_CLASS4(dnbinom, double)
MAKE_RECYCLING_RULE_CLASS3(dpois, double)
MAKE_RECYCLING_RULE_CLASS3(dchisq, double)
MAKE_RECYCLING_RULE_CLASS4(dbeta, double)
MAKE_RECYCLING_RULE_CLASS4(dnorm, double)
MAKE_RECYCLING_RULE_CLASS4(dgamma, double)
MAKE_RECYCLING_RULE_CLASS4(dlnorm, double)
MAKE_RECYCLING_RULE_CLASS4(dunif, double)
MAKE_RECYCLING_RULE_CLASS4(dweibull, double)
MAKE_RECYCLING_RULE_CLASS4(dt_nonstandard, double)


#endif
