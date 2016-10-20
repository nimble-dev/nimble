#ifndef __NIMBLE_EIGEN
#define __NIMBLE_EIGEN

#include <iostream>
#include <cstdlib>

// should the template arguments (types) be Map<MatrixXd> instead of MatrixXd?

// concatenation, c(A1, A2)
template<typename Derived1, typename Derived2>
class concatenateClass {
public:
  const Derived1 &Arg1;
  const Derived2 &Arg2;
  int size1, size2;
  typedef double result_type;
  concatenateClass(const Derived1 &A1, const Derived2 &A2) : Arg1(A1), Arg2(A2) {
    // add more complete checking of row vs. col vector here
    // assume col vector for now
    size1 = Arg1.rows();
    size2 = Arg2.rows();
  };
  typedef Eigen::internal::traits<MatrixXd>::Index Index;
  result_type operator()() const
  {
    std::cout<<"IN 0\n";
    return Arg1.coeff(0);
  }


  result_type operator()(Index i) const //Eigen::DenseIndex
  {
    std::cout<<"IN 1\n";
    if(i < size1)
      return Arg1.coeff(i);
    else
      return Arg2.coeff(i - size1);
  }

  // added this
  result_type operator()(Index i, Index j) const
  {
    std::cout<<"IN 2\n";
    return Arg1.coeff(0, 0);
  }
};

template<typename Derived1, typename Derived2>
CwiseNullaryOp<concatenateClass<Derived1, Derived2>, MatrixXd > concatenate(const Derived1 &A1, const Derived1 &A2) {
  concatenateClass<Derived1, Derived2> c(A1, A2);
  return(CwiseNullaryOp<concatenateClass<Derived1, Derived2>, MatrixXd >(A1.rows() + A2.rows(), 1, c));
}

#define nimC concatenate<Eigen::MatrixXd, Eigen::MatrixXd>


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
    sizeArg1 = Arg1.rows();
    if(each > 1)
      if(times > 1)
	eachTimesCase = useBoth;
      else
	eachTimesCase = useEach;
    else
      eachTimesCase = useTimes;
    // may need to handle times = 0 or each = 0
  };
  typedef Eigen::internal::traits<MatrixXd>::Index Index;
  result_type operator()() const
  {
    std::cout<<"IN 0 rep\n";
    return Arg1.coeff(0);
  }


  result_type operator()(Index i) const //Eigen::DenseIndex
  {
    std::cout<<"IN 1 rep\n";
    switch(eachTimesCase) {
    case useEach:
      return(Arg1.coeff(floor(i/each)));
      break;
    case useTimes:
      return(Arg1.coeff(i % sizeArg1));
      break;
    case useBoth:
      return(Arg1.coeff(static_cast<int>(floor(i / each)) % sizeArg1));
    default:
      return(0); //error
    }
    
  }

  // added this
  result_type operator()(Index i, Index j) const
  {
    std::cout<<"IN 2 rep\n";
    return Arg1.coeff(0, 0);
  }
};

template<typename Derived1>
CwiseNullaryOp<repClass<Derived1>, MatrixXd > rep(const Derived1 &A1, int reps, int each) {
repClass<Derived1> repObj(A1, reps, each);
return(CwiseNullaryOp<repClass<Derived1>, MatrixXd >(A1.rows() * reps * each, 1, repObj));
}

#define nimRep rep<Eigen::MatrixXd>


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
  const DerivedObj &source;
  const DerivedI1 &index1;
  const DerivedI2 &index2;
  int dim1, dim2;
  typedef double result_type;
 nonseqIndexedClass(const DerivedObj &s, const DerivedI1 &i1, const DerivedI2 &i2) :
  source(s),
    index1(i1),
    index2(i2) {
      dim1 = i1.rows();
      dim2 = i2.rows();
    }
  typedef Eigen::internal::traits<MatrixXd>::Index Index;

  result_type operator()() const
  {
    std::cout<<"IN 0\n";
    return source.coeff(0);
  }


  result_type operator()(Index i) const //Eigen::DenseIndex
  {
    std::cout<<"IN 1\n";
    std::div_t divRes = div(i, dim1);
    return source.coeff(index1(divRes.rem)-1, index2(floor(divRes.quot))-1);
  }

  
  result_type operator()(Index i, Index j) const
  {
    std::cout<<"IN 2\n";
    return source.coeff(index1(i)-1, index2(j)-1);
  }

  
};

template<typename DerivedObj, typename DerivedI1, typename DerivedI2>
  CwiseNullaryOp<nonseqIndexedClass<DerivedObj, DerivedI1, DerivedI2 >, DerivedObj > nonseqIndexed(const DerivedObj &s, const DerivedI1 &i1, const DerivedI2 &i2) {
  nonseqIndexedClass<DerivedObj, DerivedI1, DerivedI2 > nonseqIndexedObj(s, i1, i2);
  return(CwiseNullaryOp<nonseqIndexedClass<DerivedObj, DerivedI1, DerivedI2 >, DerivedObj >(nonseqIndexedObj.dim1, nonseqIndexedObj.dim2, nonseqIndexedObj));
}

#define nimNonseqIndexed nonseqIndexed<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXi>

#endif
