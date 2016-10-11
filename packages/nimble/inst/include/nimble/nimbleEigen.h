#ifndef __NIMBLE_EIGEN
#define __NIMBLE_EIGEN

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

#define concatenateM concatenate<Eigen::MatrixXd, Eigen::MatrixXd>


// rep
template<typename Derived1, typename Derived2>
class repClass {
public:
  const Derived1 &Arg1;
  bool each;
  int sizeArg1, numReps;
  typedef double result_type;
 concatenateClass(const Derived1 &A1, int numRepsIn, bool eachIn) :
  Arg1(A1),
    numReps(numRepsIn),
    each(eachIn) {
    // add more complete checking of row vs. col vector here
    // assume col vector for now
    sizeArg1 = Arg1.rows();
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
    if(each) {
      return(Arg1.coeff(floor(i/numReps)));
    } else {
      return(Arg1.coeff(i %% sizeArg1));
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
CwiseNullaryOp<concatenateClass<Derived1, Derived2>, MatrixXd > rep(const Derived1 &A1, int reps, bool each) {
repClass<Derived1, Derived2> rep(A1, reps, each);
return(CwiseNullaryOp<repClass<Derived1>, MatrixXd >(A1.rows() * reps, 1, rep));
}

#define repM rep<Eigen::MatrixXd>
  
#endif
