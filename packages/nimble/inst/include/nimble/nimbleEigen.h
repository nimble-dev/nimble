#ifndef __NIMBLE_EIGEN
#define __NIMBLE_EIGEN

#include <iostream>
#include <cstdlib>

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

// to-do: make native integer handling where appropriate
// right now from and to are doubles
class seqClass {
public:
  double from, by;
  int length_out;
  //  typedef enum{useEach, useTimes, useBoth} eachTimesCaseType;
  //eachTimesCaseType eachTimesCase;
  
  typedef double result_type;
  seqClass(double fromIn, double toIn, int length_outIn) :
    from(fromIn),
      length_out(length_outIn) {
	if(length_out == 1) by = 0;
	else by = (toIn - from) / (length_out - 1);
	// may need to handle invalid arguments or return of length 0
  };
    // type of the 3rd argument distinguises by vs. length_out
 seqClass(double fromIn, double toIn, double byIn) :
    from(fromIn),
      by(byIn) {
	length_out = 1 + static_cast<int>(floor(toIn - from) / byIn);
    // may need to handle invalid arguments or return of length 0
  };

    
  typedef Eigen::internal::traits<MatrixXd>::Index Index;
  result_type operator()() const
  {
    std::cout<<"IN 0 seq\n";
    return 0;
  }


  result_type operator()(Index i) const //Eigen::DenseIndex
  {
    std::cout<<"IN 1 seq\n";
    // add check of whether i >= length_out
    return( from + static_cast<int>(i) * by );
  }

  // added this
  result_type operator()(Index i, Index j) const
  {
    std::cout<<"IN 2 seq\n";
    return 0;
  }
};

template<typename Derived1>
CwiseNullaryOp<seqClass, MatrixXd > seqBy(double from, double to, double by, int NOTUSED) {
  seqClass seqObj(from, to, by);
  return(CwiseNullaryOp<seqClass, MatrixXd >(seqObj.length_out, 1, seqObj));
}

template<typename Derived1>
CwiseNullaryOp<seqClass, MatrixXd > seqLen(double from, double to, double NOTUSED, int length_out) {
  seqClass seqObj(from, to, length_out);
  return(CwiseNullaryOp<seqClass, MatrixXd >(seqObj.length_out, 1, seqObj));
}

#define nimSeqBy seqBy<Eigen::MatrixXd>
#define nimSeqLen seqLen<Eigen::MatrixXd>

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

// vectorization of any scalar function:

/* // need to make a different template for each number of arguments with matching number of template arguments (plus a F type and return type). */
/* // need to consider vector vs. matrix return (maybe handled all from nimble types) */
/* // see here for getting template types inferred: http://stackoverflow.com/questions/797594/when-a-compiler-can-infer-a-template-parameter */
/* // 2 arguments: */
/* template<typename DerivedA1, typename DerivedA2, typename F> */
/* class vectorizedFun { */
/*  public: */
/*   const DerivedA1 &arg1; */
/*   const DerivedA2 &arg2; */
/*   int size1, size2, outputSize; */
/*   F fun; */
/*   typedef Eigen::internal::traits<MatrixXd>::Index Index; */
  
/*  vectorizedFun(const DerivedA1 &A1, const DerivedA2 &A2, int oSize, F f) : */
/*   arg1(A1), */
/*     arg2(A2), */
/*     outputSize(oSize), */
/*     fun(f) { */
/*       size1 = A1.nrow(); */
/*       size2 = A2.nrow(); */
/*     } */
/*   result_type operator()(Index i) const //Eigen::DenseIndex */
/*   { */
/*     std::cout<<"IN 1\n"; */
/*     return fun(A1((i % size1) - 1), */
/* 	       A2((i % size2) - 1));  */
/*   } */
/* }; */

/* template<typename DerivedA1, typename DerivedA2, typename F> */
/*   vectorizedFun<DerivedA1, DerivedA2, F> newVectorizedFun(const DerivedA1 &A1, const DerivedA2 &A2, F f) {return vectorizedFun(A1, A2, f);} */

#endif
