/*
 * NIMBLE: an R package for programming with BUGS models.
 * Copyright (C) 2014-2017 Perry de Valpine, Christopher Paciorek,
 * Daniel Turek, Clifford Anderson-Bergman, Nick Michaud, Fritz Obermeyer,
 * Duncan Temple Lang.
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, a copy is available at
 * https://www.R-project.org/Licenses/
 */

#ifndef __NIMBLE_EIGEN
#define __NIMBLE_EIGEN

#include <iostream>
#include <stdlib.h>
#include<Rmath.h>
#include<vector>
#include<cstdlib>
#include "dists.h"


// a utility function used by nimSeq and generated size expressions to determine the length of a sequence
template<typename fromT, typename toT, typename byT>
  int calcSeqLength(fromT from, toT to, byT by) { // we need this function because of imprecision issues
  double doubleLength = (static_cast<double>(to) - static_cast<double>(from))/static_cast<double>(by);
  return(1 + floorOrEquivalent(doubleLength));
}

// a utility function used to determine missing nrow or ncol for a matrix
template<typename totLenT, typename knownDimT>
  int calcMissingMatrixSize(totLenT totLen, knownDimT knownDim) {
  double doubleLength = (static_cast<double>(totLen) - 1.) / static_cast<double>(knownDim);
  return(1 + floorOrEquivalent(doubleLength));
}

// put the call to arg.size() in a struct so we can proxy it with "1" for a scalar type
// wrap access to Eigen's traits::..::LinearAccessBit so we can proxy it with true for a scalar type (double, int, bool)
template<typename T>
struct nimble_eigen_traits {
  enum {nimbleUseLinearAccess = int(Eigen::internal::traits<T>::Flags & LinearAccessBit)};
  typedef typename Eigen::internal::traits<T>::Scalar Scalar;
};

template<>
struct nimble_eigen_traits<double> {
  enum {nimbleUseLinearAccess = int(1)};
  typedef double Scalar;
};

template<>
struct nimble_eigen_traits<int> {
  enum {nimbleUseLinearAccess = int(1)};
  typedef int Scalar;
};

template<>
struct nimble_eigen_traits<bool> {
  enum {nimbleUseLinearAccess = int(1)};
  typedef bool Scalar;
};

template<typename V>
struct nimble_eigen_traits<std::vector<V> > {
  enum {nimbleUseLinearAccess = int(1)};
  typedef V Scalar;
};

template<typename T>
struct nimble_size_impl {
  static unsigned int getSize(const T &arg) {return arg.size();}
  static unsigned int getRows(const T &arg) {return arg.rows();}
  static unsigned int getCols(const T &arg) {return arg.cols();}
  static unsigned int getDiagRows(const T &arg) {return arg.size();}
  static unsigned int getDiagCols(const T &arg) {return arg.size();}
};

template<>
struct nimble_size_impl<double> {
  static unsigned int getSize(const double &arg) {return 1;}
  static unsigned int getRows(const double &arg) {return 1;}
  static unsigned int getCols(const double &arg) {return 1;}
  static unsigned int getDiagRows(const double &arg) {return floor(arg);}
  static unsigned int getDiagCols(const double &arg) {return floor(arg);}
};

template<>
struct nimble_size_impl<int> {
  static unsigned int getSize(const int &arg) {return 1;}
  static unsigned int getRows(const int &arg) {return 1;}
  static unsigned int getCols(const int &arg) {return 1;}
  static unsigned int getDiagRows(const int &arg) {return arg;}
  static unsigned int getDiagCols(const int &arg) {return arg;}
};

template<>
struct nimble_size_impl<bool> {
  static unsigned int getSize(const bool &arg) {return 1;}
  static unsigned int getRows(const bool &arg) {return 1;}
  static unsigned int getCols(const bool &arg) {return 1;}
  static unsigned int getDiagRows(const bool &arg) {return arg;}
  static unsigned int getDiagCols(const bool &arg) {return arg;}
};


template<bool useLinearAccess, typename result_type, typename eigenType, typename Index>
struct nimble_eigen_coeff_impl;

template<typename result_type, typename eigenType, typename Index>
struct nimble_eigen_coeff_impl<true, result_type, eigenType, Index> {
  static result_type getCoeff(const eigenType &Arg, Index i) {return Arg.coeff(i);}
  static result_type getDiagCoeff(const eigenType &Arg, Index i) {return Arg.coeff(i);}
};

template<typename result_type, typename eigenType, typename Index>
struct nimble_eigen_coeff_impl<false, result_type, eigenType, Index> {
  static result_type getCoeff(const eigenType &Arg, Index i) {
    std::div_t divRes = div(static_cast<int>(i), static_cast<int>(Arg.rows())); // some compilers don't like seeing div(long, int)
    return Arg.coeff(divRes.rem, divRes.quot);
  }
  static result_type getDiagCoeff(const eigenType &Arg, Index i) {
    std::div_t divRes = div(static_cast<int>(i), static_cast<int>(Arg.rows()));
    return Arg.coeff(divRes.rem, divRes.quot);
  }
};

template<typename result_type, typename Index>
struct nimble_eigen_coeff_impl<true, result_type, double, Index> {
  static result_type getCoeff(const double Arg, Index i) {return Arg;}
  static result_type getDiagCoeff(const double Arg, Index i) {return 1.0;}
};

template<typename result_type, typename Index>
struct nimble_eigen_coeff_impl<true, result_type, int, Index> {
  static result_type getCoeff(const int Arg, Index i) {return Arg;}
  static result_type getDiagCoeff(const int Arg, Index i) {return 1;}
};

template<typename result_type, typename Index>
struct nimble_eigen_coeff_impl<true, result_type, bool, Index> {
  static result_type getCoeff(const bool Arg, Index i) {return Arg;}
  static result_type getDiagCoeff(const bool Arg, Index i) {return true;}
};

template<typename result_type, typename Index, typename V>
  struct nimble_eigen_coeff_impl<true, result_type, std::vector<V>, Index> {
  static result_type getCoeff(const std::vector<V> &Arg, Index i) {return Arg[i];}
  static result_type getDiagCoeff(const std::vector<V> &Arg, Index i) {return Arg[i];} // should never be called
};


#include "nimbleEigenNimArr.h"
// diagonal, cases diag(scalar) and diag(vector) to create something behaving like a matrix

template<typename DerivedIndex, typename DerivedSource, typename result_type = double>
  class diagonalClass {
 public:
  const DerivedSource &src;
  int dim1, dim2;
  //  typedef double result_type;
 diagonalClass(const DerivedSource &s) : src(s) {
    dim1 = nimble_size_impl<DerivedSource>::getDiagRows(src);
    dim2 = nimble_size_impl<DerivedSource>::getDiagCols(src);
  }
    result_type operator()(DerivedIndex i) const //Eigen::DenseIndex
  {
    //    std::cout<<"IN 1\n";
    std::div_t divRes = div(static_cast<int>(i), dim1);
    // iRow = divRes.rem
    // iCol = floor(divRes.quot)
    if(divRes.rem == divRes.quot) { // on diagonal
      return nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedSource>::nimbleUseLinearAccess), result_type, DerivedSource, DerivedIndex >::getDiagCoeff(src, divRes.rem);
    }
    return 0; // off diagonal
  }
  result_type operator()(DerivedIndex i, DerivedIndex j) const
  {
    //std::cout<<"IN 2\n";
    if(i == j) { // on diagonal
      return nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedSource>::nimbleUseLinearAccess), result_type, DerivedSource, DerivedIndex >::getDiagCoeff(src, i);
    }
    return 0; // off diagonal
  }
};

namespace Eigen{
  namespace internal{
    template<typename DerivedIndex, typename DerivedSource, typename result_type>
      struct functor_has_linear_access<diagonalClass<DerivedIndex, DerivedSource, result_type> > { enum { ret = 1}; }; 
    template<typename DerivedIndex, typename DerivedSource, typename result_type>
      struct functor_traits<diagonalClass<DerivedIndex, DerivedSource, result_type> > { enum {Cost = 10, PacketAccess = false, IsRepeatable = true }; }; 
  }
}

template<typename returnDerived, typename return_scalar = double>
struct diagonal_impl {
  typedef Eigen::Index IndexReturn;
  template<typename DerivedSource>
  static CwiseNullaryOp<diagonalClass<IndexReturn, DerivedSource, return_scalar >, returnDerived > diagonal(const DerivedSource &s) {
    diagonalClass<IndexReturn, DerivedSource, return_scalar > obj(s);
    return(CwiseNullaryOp<diagonalClass<IndexReturn, DerivedSource, return_scalar >, returnDerived >(obj.dim1, obj.dim2, obj));
  }
};

#define nimDiagonalD diagonal_impl<MatrixXd>::diagonal
#define nimDiagonalI diagonal_impl<MatrixXi>::diagonal  // We will always return a MatrixXd from diag(), so really these could be one case
#define nimDiagonalB diagonal_impl<MatrixXb>::diagonal


//concatenation c()
template<typename Index, typename Derived1, typename result_type = double>
class concatenate1Class {
  /* We need this for something like c(2) or c(2,2).  It also converts from scalar to vector nicely. */
 public:
  const Derived1 &Arg1;
  int size1, totalLength;
  //  typedef double result_type;
 concatenate1Class(const Derived1 &A1) : Arg1(A1) {
    size1 = nimble_size_impl<Derived1>::getSize(Arg1);
    totalLength = size1;
  };
  
  result_type operator()(Index i) const //Eigen::DenseIndex //Assume Index1 type and Index2 type will always be the same, or cast-able.
  {
    //    std::cout<<"IN 1\n";
    return nimble_eigen_coeff_impl< nimble_eigen_traits<Derived1>::nimbleUseLinearAccess, result_type, Derived1, Index >::getCoeff(Arg1, i);
  }

  result_type operator()(Index i, Index j) const // I don't think this should normally be called, but if it does, act like a vector
  {
    //std::cout<<"IN 2\n";
    return operator()(i);
  }
};


template<typename Index, typename Derived1, typename Derived2, typename result_type = double>
class concatenateClass {
 public:
  const Derived1 &Arg1;
  const Derived2 &Arg2;
  int size1, size2, totalLength;
 // double is default for simplicity.  nimDerivs version will use CppAD::AD<double>
  //  typedef double result_type;
 concatenateClass(const Derived1 &A1, const Derived2 &A2) : Arg1(A1), Arg2(A2) {
    size1 = nimble_size_impl<Derived1>::getSize(Arg1);
    size2 = nimble_size_impl<Derived2>::getSize(Arg2);
    totalLength = size1 + size2;
  };
  
  result_type operator()(Index i) const //Eigen::DenseIndex //Assume Index1 type and Index2 type will always be the same, or cast-able.
  {
    //    std::cout<<"IN 1\n";
    if(i < size1)
      return nimble_eigen_coeff_impl< nimble_eigen_traits<Derived1>::nimbleUseLinearAccess, result_type, Derived1, Index >::getCoeff(Arg1, i); //generalization of Arg1(i) or Arg1.coeff(i) 
    else
      return nimble_eigen_coeff_impl< nimble_eigen_traits<Derived2>::nimbleUseLinearAccess, result_type, Derived2, Index >::getCoeff(Arg2, i - size1); //Arg2(i - size1);
  }

  result_type operator()(Index i, Index j) const // I don't think this should normally be called, but if it does, act like a vector
  {
    //std::cout<<"IN 2\n";
    return operator()(i);
  }
};

template<typename Index, typename Derived1, typename Derived2, typename Derived3, typename result_type = double>
class concatenate3Class {
 public:
  const Derived1 &Arg1;
  const Derived2 &Arg2;
  const Derived3 &Arg3;
  int size1, size2, size3, size12, totalLength;
 //  typedef double result_type;
 concatenate3Class(const Derived1 &A1, const Derived2 &A2, const Derived3 &A3) : Arg1(A1), Arg2(A2), Arg3(A3) {
    size1 = nimble_size_impl<Derived1>::getSize(Arg1);
    size2 = nimble_size_impl<Derived2>::getSize(Arg2);
    size12 = size1 + size2;
    size3 = nimble_size_impl<Derived3>::getSize(Arg3);
    totalLength = size12 + size3;
  };
  
  result_type operator()(Index i) const //Eigen::DenseIndex //Assume Index1 type and Index2 type will always be the same, or cast-able.
  {
    //std::cout<<"IN 1\n";
    if(i < size1)
      return nimble_eigen_coeff_impl< nimble_eigen_traits<Derived1>::nimbleUseLinearAccess, result_type, Derived1, Index >::getCoeff(Arg1, i); //generalization of Arg1(i) or Arg1.coeff(i) 
    else
      if(i < size12)
	return nimble_eigen_coeff_impl< nimble_eigen_traits<Derived2>::nimbleUseLinearAccess, result_type, Derived2, Index >::getCoeff(Arg2, i - size1); //Arg2(i - size1);
      else
	return nimble_eigen_coeff_impl< nimble_eigen_traits<Derived3>::nimbleUseLinearAccess, result_type, Derived3, Index >::getCoeff(Arg3, i - size12);
  }

  result_type operator()(Index i, Index j) const // I don't think this should normally be called, but if it does, act like a vector
  {
    //std::cout<<"IN 2\n";
    return operator()(i);
  }
};

template<typename Index, typename Derived1, typename Derived2, typename Derived3, typename Derived4, typename result_type = double>
class concatenate4Class {
 public:
  const Derived1 &Arg1;
  const Derived2 &Arg2;
  const Derived3 &Arg3;
  const Derived4 &Arg4;
  int size1, size2, size3, size4, size12, size123, totalLength;
 //  typedef double result_type;
 concatenate4Class(const Derived1 &A1, const Derived2 &A2, const Derived3 &A3, const Derived4 &A4) : Arg1(A1), Arg2(A2), Arg3(A3), Arg4(A4) {
    size1 = nimble_size_impl<Derived1>::getSize(Arg1);
    size2 = nimble_size_impl<Derived2>::getSize(Arg2);
    size12 = size1 + size2;
    size3 = nimble_size_impl<Derived3>::getSize(Arg3);
    size123 = size12 + size3;
    size4 = nimble_size_impl<Derived4>::getSize(Arg4);
    totalLength = size123 + size4;
  };
  
  result_type operator()(Index i) const //Eigen::DenseIndex //Assume Index1 type and Index2 type will always be the same, or cast-able.
  {
    //std::cout<<"IN 1\n";
    if(i < size1)
      return nimble_eigen_coeff_impl< nimble_eigen_traits<Derived1>::nimbleUseLinearAccess, result_type, Derived1, Index >::getCoeff(Arg1, i); //generalization of Arg1(i) or Arg1.coeff(i) 
    else
      if(i < size12)
	return nimble_eigen_coeff_impl< nimble_eigen_traits<Derived2>::nimbleUseLinearAccess, result_type, Derived2, Index >::getCoeff(Arg2, i - size1); //Arg2(i - size1);
      else
	if(i < size123)
	  return nimble_eigen_coeff_impl< nimble_eigen_traits<Derived3>::nimbleUseLinearAccess, result_type, Derived3, Index >::getCoeff(Arg3, i - size12);
	else
	  return nimble_eigen_coeff_impl< nimble_eigen_traits<Derived4>::nimbleUseLinearAccess, result_type, Derived4, Index >::getCoeff(Arg4, i - size123);
  }

  result_type operator()(Index i, Index j) const // I don't think this should normally be called, but if it does, act like a vector
  {
    //std::cout<<"IN 2\n";
    return operator()(i);
  }
};


namespace Eigen{
  namespace internal{
    template<typename Index, typename Derived1, typename Derived2, typename result_type>
    // 3 options for linear_access:
      //1. turn it off: struct functor_has_linear_access<gConcatenateClass<Derived1, Derived2> > { enum { ret = 0 }; }; 
      //2. turn it on: (and let it be resolved by nimble_eigen_coeff_impl>:
      struct functor_has_linear_access<concatenateClass<Index, Derived1, Derived2, result_type> > { enum { ret = 1}; }; 
      //3. Set it once according to arguments: (problem here is that if it is off for this expression, that may make expressions using this one forbidden from using linear access, which is a bit harsh: struct functor_has_linear_access<gConcatenateClass<Derived1, Derived2> > { enum { ret = traits<Derived1>::Flags & traits<Derived2>::Flags & LinearAccessBit }; };
    template<typename Index, typename Derived1, typename Derived2, typename result_type>
      struct functor_traits<concatenateClass<Index, Derived1, Derived2, result_type> >
      {
	enum
	{
	  Cost = 10, // there are templated costs available to pick up for this
	  PacketAccess = false, // happy to keep this false for now
	  IsRepeatable = true // default was false. 
	};
      };

    template<typename Index, typename Derived1, typename result_type>
      struct functor_has_linear_access<concatenate1Class<Index, Derived1, result_type> > { enum { ret = 1}; }; 
    template<typename Index, typename Derived1, typename result_type>
      struct functor_traits<concatenate1Class<Index, Derived1, result_type> > { enum { Cost = 10, PacketAccess = false, IsRepeatable = true }; };

    template<typename Index, typename Derived1, typename Derived2, typename Derived3, typename result_type>
      struct functor_has_linear_access<concatenate3Class<Index, Derived1, Derived2, Derived3, result_type> > { enum { ret = 1}; }; 
    template<typename Index, typename Derived1, typename Derived2, typename Derived3, typename result_type>
      struct functor_traits<concatenate3Class<Index, Derived1, Derived2, Derived3, result_type> > { enum { Cost = 10, PacketAccess = false, IsRepeatable = true }; };

    template<typename Index, typename Derived1, typename Derived2, typename Derived3, typename Derived4, typename result_type>
      struct functor_has_linear_access<concatenate4Class<Index, Derived1, Derived2, Derived3, Derived4, result_type> > { enum { ret = 1}; }; 
    template<typename Index, typename Derived1, typename Derived2, typename Derived3, typename Derived4, typename result_type>
      struct functor_traits<concatenate4Class<Index, Derived1, Derived2, Derived3, Derived4, result_type> > { enum { Cost = 10, PacketAccess = false, IsRepeatable = true }; };
  }
}

template<typename returnDerived, typename return_scalar = double>
struct concatenate_impl {
  typedef Eigen::Index Index;
  template<typename Derived1, typename Derived2>
    static CwiseNullaryOp<concatenateClass<Index, Derived1, Derived2, return_scalar>, returnDerived > concatenate(const Derived1 &A1, const Derived2 &A2) {
    concatenateClass<Index, Derived1, Derived2, return_scalar> c(A1, A2);
    return(CwiseNullaryOp<concatenateClass<Index, Derived1, Derived2, return_scalar>, returnDerived >(c.totalLength, 1, c));
  }
  template<typename Derived1>
    static CwiseNullaryOp<concatenate1Class<Index, Derived1, return_scalar>, returnDerived > concatenate(const Derived1 &A1) {
    concatenate1Class<Index, Derived1, return_scalar> c(A1);
    return(CwiseNullaryOp<concatenate1Class<Index, Derived1, return_scalar>, returnDerived >(c.totalLength, 1, c));
  }
  template<typename Derived1, typename Derived2, typename Derived3>
    static CwiseNullaryOp<concatenate3Class<Index, Derived1, Derived2, Derived3, return_scalar>, returnDerived > concatenate(const Derived1 &A1, const Derived2 &A2, const Derived3 &A3) {
    concatenate3Class<Index, Derived1, Derived2, Derived3, return_scalar> c(A1, A2, A3);
    return(CwiseNullaryOp<concatenate3Class<Index, Derived1, Derived2, Derived3, return_scalar>, returnDerived >(c.totalLength, 1, c));
  }
  template<typename Derived1, typename Derived2, typename Derived3, typename Derived4>
    static CwiseNullaryOp<concatenate4Class<Index, Derived1, Derived2, Derived3, Derived4, return_scalar>, returnDerived > concatenate(const Derived1 &A1, const Derived2 &A2, const Derived3 &A3, const Derived4 &A4) {
    concatenate4Class<Index, Derived1, Derived2, Derived3, Derived4, return_scalar> c(A1, A2, A3, A4);
    return(CwiseNullaryOp<concatenate4Class<Index, Derived1, Derived2, Derived3, Derived4, return_scalar>, returnDerived >(c.totalLength, 1, c));
  }
};

#define nimCd concatenate_impl<MatrixXd>::concatenate
#define nimCi concatenate_impl<MatrixXi>::concatenate
#define nimCb concatenate_impl<MatrixXb>::concatenate

// rep, rep(x, times, each)
template<typename Index1, typename Derived1> // times and each should be scalar but if non-scalar first value is used
class repClass {
public:
  const Derived1 &Arg1;
  int times, each;
  typedef enum{useEach, useTimes, useBoth} eachTimesCaseType;
  eachTimesCaseType eachTimesCase;
  int sizeArg1, outputLength;
  typedef double result_type;
 repClass(const Derived1 &A1, int timesIn, int eachIn) :
  Arg1(A1),
    times(timesIn),
    each(eachIn) {
    sizeArg1 = nimble_size_impl<Derived1>::getSize(Arg1);
    if(each > 1)
      if(times > 1)
	eachTimesCase = useBoth;
      else
	eachTimesCase = useEach;
    else
      eachTimesCase = useTimes;
    // may need to handle times = 0 or each = 0
    outputLength = sizeArg1 * timesIn * eachIn;
  };
 repClass(const Derived1 &A1, int timesIn, int eachIn, int outputLengthIn) : //ignore timesIn in this case
  Arg1(A1),
    each(eachIn),
    outputLength(outputLengthIn)
    {      
      sizeArg1 = nimble_size_impl<Derived1>::getSize(Arg1);
      times = ceil( outputLengthIn / static_cast<double>(sizeArg1 * each) ); // only used for purposes of setting eachTimesCase next
      if(each > 1)
	if(times > 1)
	  eachTimesCase = useBoth;
	else
	  eachTimesCase = useEach;
      else
	eachTimesCase = useTimes;
      // may need to handle times = 0 or each = 0
    };

  result_type operator()(Index1 i) const //Eigen::DenseIndex
  {
    Index1 iUse;
    //    std::cout<<"IN 1 rep\n";
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
    return nimble_eigen_coeff_impl< nimble_eigen_traits<Derived1>::nimbleUseLinearAccess, result_type, Derived1, Index1 >::getCoeff(Arg1, iUse);
  }

  result_type operator()(Index1 i, Index1 j) const
  {
    //std::cout<<"IN 2 rep\n";
    return operator()(i);
  }
};

namespace Eigen{
  namespace internal{
    template<typename Index, typename Derived1>
      struct functor_has_linear_access<repClass<Index, Derived1> > { enum { ret = 1}; }; 
    template<typename Index, typename Derived1>
      struct functor_traits<repClass<Index, Derived1> > { enum { Cost = 10, PacketAccess = false, IsRepeatable = true}; };
  }
}

template<typename returnDerived>
struct rep_impl {
  typedef Eigen::Index IndexReturn;

  template<typename Derived1, typename DerivedEach> // timesIn can always be int here because if it is non-int this code isn't used
  static CwiseNullaryOp<repClass<IndexReturn, Derived1>, returnDerived > rep(const Derived1 &A1, int timesIn, const DerivedEach &each) {
    repClass<IndexReturn, Derived1> repObj(A1,
					   timesIn,
					   nimble_eigen_coeff_impl< nimble_eigen_traits<DerivedEach>::nimbleUseLinearAccess, int, DerivedEach, unsigned int >::getCoeff(each, 0));
    return(CwiseNullaryOp<repClass<IndexReturn, Derived1>, returnDerived >(repObj.outputLength, 1, repObj));
  }
  
  template<typename Derived1, typename DerivedLengthOut, typename DerivedEach>
    static CwiseNullaryOp<repClass<IndexReturn, Derived1>, returnDerived > rep(const Derived1 &A1, int timesIn, const DerivedLengthOut &length_out, const DerivedEach &eachIn) {
    repClass<IndexReturn, Derived1> repObj(A1, timesIn,
					   nimble_eigen_coeff_impl< nimble_eigen_traits<DerivedEach>::nimbleUseLinearAccess, int, DerivedEach, unsigned int >::getCoeff(eachIn, 0),
					   nimble_eigen_coeff_impl< nimble_eigen_traits<DerivedLengthOut>::nimbleUseLinearAccess, int, DerivedLengthOut, unsigned int >::getCoeff(length_out, 0));
    return(CwiseNullaryOp<repClass<IndexReturn, Derived1>, returnDerived >(repObj.outputLength, 1, repObj));
  }
};

#define nimRepd rep_impl<MatrixXd>::rep
#define nimRepi rep_impl<MatrixXi>::rep
#define nimRepb rep_impl<MatrixXb>::rep

// sequences, seq(from, to, by, length.out)  not implementing along.with for now
// not implementing this for `:` because we already have special treatment `:` in for-loop context, so want to tread carefully there

typedef enum{useBy, useLength, useByAndLength} byOrLength;

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
      //    printf("Add some checking to seqClass constructor and deal with inconsistent scalar types\n");
      length_out = calcSeqLength(fromIn, toIn, byIn);
      //      length_out = 1 + static_cast<int>(floor(static_cast<double>(toIn) - static_cast<double>(from)) / static_cast<double>(byIn));
    };
  
  typedef Eigen::Index Index;
  result_type operator()(Index i) const //Eigen::DenseIndex
  {
    //std::cout<<"IN 1 seq\n";
    return( from + static_cast<int>(i) * by );
  }
  result_type operator()(Index i, Index j) const
  {
    //std::cout<<"IN 2 seq\n";
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
      //      printf("Add some checking to seqClass constructor and deal with inconsistent scalar types\n");
      if(length_out == 1) {
	//std::cout<<"setting by to 0\n";
	by = 0;
      } else {
	//std::cout<<static_cast<double>(toIn)<<" "<<static_cast<double>(from)<<" "<<static_cast<double>(length_out)<<"\n";
	by = (static_cast<double>(toIn) - static_cast<double>(from)) / (static_cast<double>(length_out) - 1.);
      }
    };

  
  typedef Eigen::Index Index;
  result_type operator()(Index i) const //Eigen::DenseIndex
  {
    //std::cout<<"IN 1 seq\n";
    return( from + static_cast<int>(i) * by );
  }
  result_type operator()(Index i, Index j) const
  {
    //std::cout<<"IN 2 seq\n";
    if(j != 0) printf("Problem calling seq in C++ with two indices\n");
    return( from + static_cast<int>(i) * by);
  }
};


template<typename DerivedOut, typename scalarFrom, typename scalarBy>
  class seqClass<DerivedOut, scalarFrom, int, scalarBy, useByAndLength> {
public:
  scalarFrom from;
  scalarBy by;
  unsigned int length_out;
  
  typedef double result_type; // need to pull this from DerivedOut
 seqClass(scalarFrom fromIn, scalarBy byIn, unsigned int length_outIn) :
  from(fromIn),
    by(byIn),
    length_out(length_outIn)
    {};

  
  typedef Eigen::Index Index;
  result_type operator()(Index i) const //Eigen::DenseIndex
  {
    return( from + static_cast<int>(i) * by );
  }
  result_type operator()(Index i, Index j) const
  {
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

  template<typename DerivedOut, typename scalarFrom, typename scalarTo, typename scalarBy>
    struct functor_has_linear_access<seqClass<DerivedOut, scalarFrom, scalarTo, scalarBy, useByAndLength> > { enum { ret = 1}; }; 
  template<typename DerivedOut, typename scalarFrom, typename scalarTo, typename scalarBy>
    struct functor_traits<seqClass<DerivedOut, scalarFrom, scalarTo, scalarBy, useByAndLength> >
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
  template<typename scalarFrom, typename scalarTo, typename scalarBy>
    static CwiseNullaryOp<seqClass<DerivedOut, scalarFrom,  scalarTo, scalarBy, useByAndLength>, DerivedOut > seqByLen(scalarFrom from, scalarTo to, scalarBy by, unsigned int len) {
    seqClass<DerivedOut, scalarFrom, scalarTo, scalarBy, useByAndLength> seqObj(from, by, len);
    return(CwiseNullaryOp<seqClass<DerivedOut, scalarFrom, scalarTo, scalarBy, useByAndLength> , DerivedOut >(seqObj.length_out, 1, seqObj));
  }

};

#define nimSeqByD seq_impl<MatrixXd>::seqBy
#define nimSeqLenD seq_impl<MatrixXd>::seqLen
#define nimSeqByLenD seq_impl<MatrixXd>::seqByLen

#define nimSeqByI seq_impl<MatrixXi>::seqBy
#define nimSeqLenI seq_impl<MatrixXi>::seqLen
#define nimSeqByLenI seq_impl<MatrixXi>::seqByLen

// coeffSetter, for cases like X[indices] <- Y
template<typename IndexType, typename DerivedTarget, typename DerivedIndex1, typename DerivedIndex2>
class coeffSetterClass {
public:
  DerivedTarget &target;
  const DerivedIndex1 &I1;
  const DerivedIndex2 &I2;
  int dim1, dim2, totSize;
  //  coeffSetterClass(DerivedTarget &targetIn, const DerivedIndex1 &I1in, const DerivedIndex2 &I2in) :
 coeffSetterClass(DerivedTarget &targetIn, const DerivedIndex1 &I1in, const DerivedIndex2 &I2in) :
    target(targetIn),
    I1(I1in),
    I2(I2in) {
    dim1 = nimble_size_impl<DerivedIndex1>::getSize(I1);
    dim2 = nimble_size_impl<DerivedIndex2>::getSize(I2);
    totSize = dim1 * dim2;
  }
  int size() const {return(totSize);}
  int rows() const {return(dim1);}
  int cols() const {return(dim2);}

  template<typename fromType>
  void operator=(const fromType &from) {
    //printf("In operator=\n");
    if(nimble_size_impl<fromType>::getSize(from) < totSize) {
      printf("PROBLEM\n");
      return;
    }
    for(int i = 0 ; i  < totSize; i++) {
      coeffRef(i) = nimble_eigen_coeff_impl< bool(nimble_eigen_traits<fromType>::nimbleUseLinearAccess),
      typename nimble_eigen_traits<fromType>::Scalar, fromType, IndexType >::getCoeff(from, i);
      //      from(i);
    }
  }
  template<typename fromType>
  void fill(const fromType &from) {
    Scalar val = nimble_eigen_coeff_impl< bool(nimble_eigen_traits<fromType>::nimbleUseLinearAccess), Scalar, fromType, IndexType >::getCoeff(from, 0);
    for(int i = 0 ; i  < totSize; i++) {
      coeffRef(i) = val;
    }
  }
  // this will only work for Eigen types 
  typedef typename Eigen::internal::traits<DerivedTarget>::Scalar Scalar;
  typedef typename nimble_eigen_traits<DerivedIndex1>::Scalar I1Scalar;
  typedef typename nimble_eigen_traits<DerivedIndex2>::Scalar I2Scalar;
  
  Scalar &coeffRef(IndexType i) const {
    std::div_t divRes = div(static_cast<int>(i), dim1);
    return target.coeffRef(nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedIndex1>::nimbleUseLinearAccess), I1Scalar, DerivedIndex1, IndexType >::getCoeff(I1, divRes.rem)-1,
			   nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedIndex2>::nimbleUseLinearAccess), I2Scalar, DerivedIndex2, IndexType >::getCoeff(I2, divRes.quot)-1);

    // use % to get the i-th total element
  }
  Scalar &coeffRef(IndexType i, IndexType j) const {
    return target.coeffRef( nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedIndex1>::nimbleUseLinearAccess), I1Scalar, DerivedIndex1, IndexType >::getCoeff(I1, i)-1,
			    nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedIndex2>::nimbleUseLinearAccess), I2Scalar, DerivedIndex2, IndexType >::getCoeff(I2, j)-1);
  }
  // accesses coeffRef of objects using provided indices.
};


template <typename DerivedTarget, typename DerivedIndex1, typename DerivedIndex2>
  coeffSetterClass<Eigen::Index, DerivedTarget, DerivedIndex1, DerivedIndex2> coeffSetter(DerivedTarget &targetIn, const DerivedIndex1 &I1in, const DerivedIndex2 &I2in ) {
  return coeffSetterClass<Eigen::Index, DerivedTarget, DerivedIndex1, DerivedIndex2>(targetIn, I1in, I2in);
}
// want coeffSetter(A, indices1, indices2) to work.  We pull out scalar type of A inside coeffSetterClass
  
// nonseqIndexed

template<typename IndexObj, typename DerivedObj, typename DerivedI1, typename DerivedI2, typename result_type = double>
class nonseqIndexedClass {
 public:
  const DerivedObj &obj;
  const DerivedI1 &index1;
  const DerivedI2 &index2;
  int dim1, dim2;
 //  typedef double result_type;
 nonseqIndexedClass(const DerivedObj &s, const DerivedI1 &i1, const DerivedI2 &i2) :
  obj(s),
    index1(i1),
    index2(i2) {
      dim1 = nimble_size_impl<DerivedI1>::getSize(i1);
      dim2 = nimble_size_impl<DerivedI2>::getSize(i2);
    }

  result_type operator()(IndexObj i) const //Eigen::DenseIndex
  {
    //std::cout<<"IN 1\n";
    std::div_t divRes = div(static_cast<int>(i), dim1);
    return obj.coeff(nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedI1>::nimbleUseLinearAccess), int, DerivedI1, IndexObj >::getCoeff(index1, divRes.rem) - 1,
		     nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedI2>::nimbleUseLinearAccess), int, DerivedI2, IndexObj >::getCoeff(index2, divRes.quot) - 1); // This type of the index argument is confusing.  What is being passed is a type from std::div_t, which ought to be castable to any Eigen Index type I hope.
    //index1(divRes.rem)-1, index2(floor(divRes.quot))-1);
  }
  result_type operator()(IndexObj i, IndexObj j) const
  {
    //std::cout<<"IN 2\n";
    return obj.coeff(nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedI1>::nimbleUseLinearAccess), int, DerivedI1, IndexObj >::getCoeff(index1, i) - 1,
		     nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedI2>::nimbleUseLinearAccess), int, DerivedI2, IndexObj >::getCoeff(index2, j) - 1);

    //return obj.coeff(index1(i)-1,
    //		     index2(j)-1);
  }
};

namespace Eigen{
  namespace internal{
    template<typename IndexObj, typename DerivedObj, typename Derived1, typename Derived2, typename result_type>
      struct functor_has_linear_access<nonseqIndexedClass<IndexObj, DerivedObj, Derived1, Derived2, result_type> > { enum { ret = 1}; }; 
    template<typename IndexObj, typename DerivedObj, typename Derived1, typename Derived2, typename result_type>
      struct functor_traits<nonseqIndexedClass<IndexObj, DerivedObj, Derived1, Derived2, result_type> > { enum {Cost = 10, PacketAccess = false, IsRepeatable = true }; }; 
  }
}

template<typename returnDerived, typename result_type = double>
struct nonseqIndexed_impl {
  typedef Eigen::Index IndexReturn;
  template<typename DerivedObj, typename DerivedI1, typename DerivedI2>
    static CwiseNullaryOp<nonseqIndexedClass<IndexReturn, DerivedObj, DerivedI1, DerivedI2, result_type >, returnDerived > nonseqIndexed(const DerivedObj &s, const DerivedI1 &i1, const DerivedI2 &i2) {
    nonseqIndexedClass<IndexReturn, DerivedObj, DerivedI1, DerivedI2, result_type > nonseqIndexedObj(s, i1, i2);
    return(CwiseNullaryOp<nonseqIndexedClass<IndexReturn, DerivedObj, DerivedI1, DerivedI2, result_type >, returnDerived >(nonseqIndexedObj.dim1, nonseqIndexedObj.dim2, nonseqIndexedObj));
  }
};

#define nimNonseqIndexedd nonseqIndexed_impl<MatrixXd>::nonseqIndexed
#define nimNonseqIndexedi nonseqIndexed_impl<MatrixXi>::nonseqIndexed
#define nimNonseqIndexedb nonseqIndexed_impl<MatrixXb>::nonseqIndexed

// get first element or length.  used for lengths of return values of recycling rule r functions needed for sizeExprs

template<typename NimArrOfSomeKind>
int rFunLength(const NimArrOfSomeKind &Arg) { //
  if(Arg.size() == 1) return Arg[0];
  return Arg.size();
}

// vectorization of any scalar function: R's so-called "Recycling Rule"
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
    std::div_t divRes = div(static_cast<int>(i % size), static_cast<int>(Arg.rows()));
    return Arg.coeff(divRes.rem, divRes.quot);
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
// double RRtest_add(double a1, double a2) {std::cout<<a1<<" "<<a2<<" "<<a1+a2<<"\n"; return a1 + a2;}

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
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2)); \
  }\
  RETURNSCALARTYPE operator()(Index i, Index j) const {\
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2)); \
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
  typedef Eigen::Index IndexReturn;\
  template<typename Derived1, typename Derived2>\
  static CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2>, DerivedReturn >\
  FUNNAME ## _RecyclingRule(const Derived1 &A1, const Derived2 &A2) {\
    FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2> obj(A1, A2);\
    return(CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2>, DerivedReturn >(obj.outputSize, 1, obj));\
  }\
};


// This is for "r" functions where the first argument, or the length of the first argument if not 1, gives the result size
#define MAKE_RECYCLING_RULE_CLASS_r3(FUNNAME, RETURNSCALARTYPE)		\
  template<typename Index, typename DerivedN, typename DerivedA1, typename DerivedA2, typename DerivedA3> \
    class FUNNAME ## RecyclingRuleClass {				\
public:									\
    const DerivedN &ArgN;						\
    const DerivedA1 &Arg1;						\
    const DerivedA2 &Arg2;						\
    const DerivedA2 &Arg3;						\
    std::vector<RETURNSCALARTYPE> values;				\
    unsigned int outputSize;						\
    FUNNAME ## RecyclingRuleClass(const DerivedN &AN, const DerivedA1 &A1, const DerivedA2 &A2, const DerivedA3 &A3 ) : \
    ArgN(AN), Arg1(A1), Arg2(A2), Arg3(A3)				\
    {									\
      int sizeN = nimble_size_impl<DerivedN>::getSize(ArgN);		\
      if(sizeN > 1) outputSize = sizeN;					\
      else outputSize = nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedN>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedN, Index >::getCoeff(ArgN, 0); \
      int size1 = nimble_size_impl<DerivedA1>::getSize(Arg1);		\
      int size2 = nimble_size_impl<DerivedA2>::getSize(Arg2);		\
      int size3 = nimble_size_impl<DerivedA3>::getSize(Arg3);		\
      values.reserve(outputSize);					\
      for(int i = 0; i < outputSize; i++) {				\
	values.push_back(FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
				 nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2), \
				 nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA3>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA3, Index >::getCoeff(Arg3, i, size3))); \
      }									\
    }									\
  RETURNSCALARTYPE operator()(Index i) const { \
    return values[i];			       \
  }\
  RETURNSCALARTYPE operator()(Index i, Index j) const {\
    return values[i];				       \
  }\
}; \
\
 namespace Eigen{\
  namespace internal{\
    template<typename Index, typename DerivedN, typename Derived1, typename Derived2, typename Derived3> \
      struct functor_has_linear_access<FUNNAME ## RecyclingRuleClass<Index, DerivedN, Derived1, Derived2, Derived3> > { enum { ret = 1}; }; \
    template<typename Index, typename DerivedN, typename Derived1, typename Derived2, typename Derived3> \
      struct functor_traits<FUNNAME ## RecyclingRuleClass<Index, DerivedN, Derived1, Derived2, Derived3> > { enum {Cost = 10, PacketAccess = false, IsRepeatable = true}; }; \
  }\
}\
\
template<typename DerivedReturn>\
struct FUNNAME ## _RR_impl {\
  typedef Eigen::Index IndexReturn;\
  template<typename DerivedN, typename Derived1, typename Derived2, typename Derived3>	\
    static CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, DerivedN, Derived1, Derived2, Derived3>, DerivedReturn > \
    FUNNAME ## _RecyclingRule(const DerivedN &AN, const Derived1 &A1, const Derived2 &A2, const Derived3 &A3) { \
    FUNNAME ## RecyclingRuleClass<IndexReturn, DerivedN, Derived1, Derived2, Derived3> obj(AN, A1, A2, A3); \
    return(CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, DerivedN, Derived1, Derived2, Derived3>, DerivedReturn >(obj.outputSize, 1, obj)); \
  }\
};

// r2
#define MAKE_RECYCLING_RULE_CLASS_r2(FUNNAME, RETURNSCALARTYPE)		\
  template<typename Index, typename DerivedN, typename DerivedA1, typename DerivedA2> \
    class FUNNAME ## RecyclingRuleClass {				\
public:									\
    const DerivedN &ArgN;						\
    const DerivedA1 &Arg1;						\
    const DerivedA2 &Arg2;						\
    std::vector<RETURNSCALARTYPE> values;				\
    unsigned int outputSize;						\
    FUNNAME ## RecyclingRuleClass(const DerivedN &AN, const DerivedA1 &A1, const DerivedA2 &A2 ) : \
    ArgN(AN), Arg1(A1), Arg2(A2)					\
    {									\
      int sizeN = nimble_size_impl<DerivedN>::getSize(ArgN);		\
      if(sizeN > 1) outputSize = sizeN;					\
      else outputSize = nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedN>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedN, Index >::getCoeff(ArgN, 0); \
      int size1 = nimble_size_impl<DerivedA1>::getSize(Arg1);		\
      int size2 = nimble_size_impl<DerivedA2>::getSize(Arg2);		\
      values.reserve(outputSize);					\
      for(int i = 0; i < outputSize; i++) {				\
	values.push_back(FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
				 nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2))); \
      }									\
    }									\
  RETURNSCALARTYPE operator()(Index i) const { \
    return values[i];			       \
  }\
  RETURNSCALARTYPE operator()(Index i, Index j) const {\
    return values[i];				       \
  }\
}; \
\
 namespace Eigen{\
  namespace internal{\
    template<typename Index, typename DerivedN, typename Derived1, typename Derived2> \
      struct functor_has_linear_access<FUNNAME ## RecyclingRuleClass<Index, DerivedN, Derived1, Derived2> > { enum { ret = 1}; }; \
    template<typename Index, typename DerivedN, typename Derived1, typename Derived2> \
      struct functor_traits<FUNNAME ## RecyclingRuleClass<Index, DerivedN, Derived1, Derived2> > { enum {Cost = 10, PacketAccess = false, IsRepeatable = true}; }; \
  }\
}\
\
template<typename DerivedReturn>\
struct FUNNAME ## _RR_impl {\
  typedef Eigen::Index IndexReturn;\
  template<typename DerivedN, typename Derived1, typename Derived2>	\
    static CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, DerivedN, Derived1, Derived2>, DerivedReturn > \
    FUNNAME ## _RecyclingRule(const DerivedN &AN, const Derived1 &A1, const Derived2 &A2) { \
    FUNNAME ## RecyclingRuleClass<IndexReturn, DerivedN, Derived1, Derived2> obj(AN, A1, A2); \
    return(CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, DerivedN, Derived1, Derived2>, DerivedReturn >(obj.outputSize, 1, obj)); \
  }\
};


// r1
#define MAKE_RECYCLING_RULE_CLASS_r1(FUNNAME, RETURNSCALARTYPE)		\
  template<typename Index, typename DerivedN, typename DerivedA1> \
    class FUNNAME ## RecyclingRuleClass {				\
public:									\
    const DerivedN &ArgN;						\
    const DerivedA1 &Arg1;						\
    std::vector<RETURNSCALARTYPE> values;				\
    unsigned int outputSize;						\
    FUNNAME ## RecyclingRuleClass(const DerivedN &AN, const DerivedA1 &A1 ) : \
    ArgN(AN), Arg1(A1)					\
    {									\
      int sizeN = nimble_size_impl<DerivedN>::getSize(ArgN);		\
      if(sizeN > 1) outputSize = sizeN;					\
      else outputSize = nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedN>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedN, Index >::getCoeff(ArgN, 0); \
      int size1 = nimble_size_impl<DerivedA1>::getSize(Arg1);		\
      values.reserve(outputSize);					\
      for(int i = 0; i < outputSize; i++) {				\
	values.push_back(FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1))); \
      }									\
    }									\
  RETURNSCALARTYPE operator()(Index i) const { \
    return values[i];			       \
  }\
  RETURNSCALARTYPE operator()(Index i, Index j) const {\
    return values[i];				       \
  }\
}; \
\
 namespace Eigen{\
  namespace internal{\
    template<typename Index, typename DerivedN, typename Derived1> \
      struct functor_has_linear_access<FUNNAME ## RecyclingRuleClass<Index, DerivedN, Derived1> > { enum { ret = 1}; }; \
    template<typename Index, typename DerivedN, typename Derived1> \
      struct functor_traits<FUNNAME ## RecyclingRuleClass<Index, DerivedN, Derived1> > { enum {Cost = 10, PacketAccess = false, IsRepeatable = true}; }; \
  }\
}\
\
template<typename DerivedReturn>\
struct FUNNAME ## _RR_impl {\
  typedef Eigen::Index IndexReturn;\
  template<typename DerivedN, typename Derived1>	\
    static CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, DerivedN, Derived1>, DerivedReturn > \
    FUNNAME ## _RecyclingRule(const DerivedN &AN, const Derived1 &A1) { \
    FUNNAME ## RecyclingRuleClass<IndexReturn, DerivedN, Derived1> obj(AN, A1); \
    return(CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, DerivedN, Derived1>, DerivedReturn >(obj.outputSize, 1, obj)); \
  }\
};


// This case doesn't even really recycle, but for consistency:
#define MAKE_RECYCLING_RULE_CLASS1_1scalar(FUNNAME, RETURNSCALARTYPE) \
template<typename Index, typename DerivedA1, typename DerivedA2> \
class FUNNAME ## RecyclingRuleClass { \
public: \
  const DerivedA1 &Arg1;\
  const DerivedA2 &Arg2;\
  unsigned int size1, outputSize;\
  FUNNAME ## RecyclingRuleClass(const DerivedA1 &A1, const DerivedA2 &A2 ) :\
  Arg1(A1), Arg2(A2)\
{\
  outputSize = size1 = nimble_size_impl<DerivedA1>::getSize(Arg1); \
  } \
  RETURNSCALARTYPE operator()(Index i) const { \
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedA2>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, 0)); \
  }\
  RETURNSCALARTYPE operator()(Index i, Index j) const {\
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedA2>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, 0)); \
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
  typedef Eigen::Index IndexReturn;\
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
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA3>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA3, Index >::getCoeff(Arg3, i, size3)); \
  }\
  RETURNSCALARTYPE operator()(Index i, Index j) const {\
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA3>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA3, Index >::getCoeff(Arg3, i, size3)); \
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
  typedef Eigen::Index IndexReturn;\
  template<typename Derived1, typename Derived2, typename Derived3>			\
    static CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3>, DerivedReturn > \
    FUNNAME ## _RecyclingRule(const Derived1 &A1, const Derived2 &A2, const Derived3 &A3) { \
    FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3> obj(A1, A2, A3); \
    return(CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3>, DerivedReturn >(obj.outputSize, 1, obj)); \
  }\
};

#define MAKE_RECYCLING_RULE_CLASS2_1scalar(FUNNAME, RETURNSCALARTYPE) \
  template<typename Index, typename DerivedA1, typename DerivedA2, typename DerivedA3>	\
class FUNNAME ## RecyclingRuleClass { \
public: \
  const DerivedA1 &Arg1;\
  const DerivedA2 &Arg2;\
  const DerivedA3 &Arg3;\
  unsigned int size1, size2, outputSize;				\
  FUNNAME ## RecyclingRuleClass(const DerivedA1 &A1, const DerivedA2 &A2, const DerivedA3 &A3 ) : \
  Arg1(A1), Arg2(A2), Arg3(A3) \
{\
  outputSize = size1 = nimble_size_impl<DerivedA1>::getSize(Arg1); \
  size2 = nimble_size_impl<DerivedA2>::getSize(Arg2); \
  if(size2 > outputSize) outputSize = size2;	      \
}						      \
  RETURNSCALARTYPE operator()(Index i) const { \
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2), \
		   nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedA3>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA3, Index >::getCoeff(Arg3, 0)); \
  }\
  RETURNSCALARTYPE operator()(Index i, Index j) const {\
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2), \
		   nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedA3>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA3, Index >::getCoeff(Arg3, 0)); \
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
  typedef Eigen::Index IndexReturn;\
  template<typename Derived1, typename Derived2, typename Derived3>			\
    static CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3>, DerivedReturn > \
    FUNNAME ## _RecyclingRule(const Derived1 &A1, const Derived2 &A2, const Derived3 &A3) { \
    FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3> obj(A1, A2, A3); \
    return(CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3>, DerivedReturn >(obj.outputSize, 1, obj)); \
  }\
};

#define MAKE_RECYCLING_RULE_CLASS2_2scalar(FUNNAME, RETURNSCALARTYPE) \
  template<typename Index, typename DerivedA1, typename DerivedA2, typename DerivedA3, typename DerivedA4>	\
class FUNNAME ## RecyclingRuleClass { \
public: \
  const DerivedA1 &Arg1;\
  const DerivedA2 &Arg2;\
  const DerivedA3 &Arg3;\
  const DerivedA4 &Arg4;						\
  unsigned int size1, size2, outputSize;				\
  FUNNAME ## RecyclingRuleClass(const DerivedA1 &A1, const DerivedA2 &A2, const DerivedA3 &A3, const DerivedA4 &A4 ) : \
  Arg1(A1), Arg2(A2), Arg3(A3), Arg4(A4)				\
{\
  outputSize = size1 = nimble_size_impl<DerivedA1>::getSize(Arg1); \
  size2 = nimble_size_impl<DerivedA2>::getSize(Arg2); \
  if(size2 > outputSize) outputSize = size2;	      \
}						      \
  RETURNSCALARTYPE operator()(Index i) const { \
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2), \
		   nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedA3>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA3, Index >::getCoeff(Arg3, 0), \
		   nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedA4>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA4, Index >::getCoeff(Arg4, 0)); \
  }\
  RETURNSCALARTYPE operator()(Index i, Index j) const {\
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2), \
		   nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedA3>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA3, Index >::getCoeff(Arg3, 0), \
    		   nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedA4>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA4, Index >::getCoeff(Arg4, 0)); \
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
  typedef Eigen::Index IndexReturn;\
  template<typename Derived1, typename Derived2, typename Derived3, typename Derived4>	\
    static CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3, Derived4>, DerivedReturn > \
    FUNNAME ## _RecyclingRule(const Derived1 &A1, const Derived2 &A2, const Derived3 &A3, const Derived4 &A4) { \
    FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3, Derived4> obj(A1, A2, A3, A4); \
    return(CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3, Derived4>, DerivedReturn >(obj.outputSize, 1, obj)); \
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
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA3>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA3, Index >::getCoeff(Arg3, i, size3), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA4>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA4, Index >::getCoeff(Arg4, i, size4)); \
  }\
  RETURNSCALARTYPE operator()(Index i, Index j) const {\
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA3>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA3, Index >::getCoeff(Arg3, i, size3), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA4>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA4, Index >::getCoeff(Arg4, i, size4)); \
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
  typedef Eigen::Index IndexReturn;\
  template<typename Derived1, typename Derived2, typename Derived3, typename Derived4>	\
    static CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3, Derived4>, DerivedReturn > \
    FUNNAME ## _RecyclingRule(const Derived1 &A1, const Derived2 &A2, const Derived3 &A3, const Derived4 &A4) { \
    FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3, Derived4> obj(A1, A2, A3, A4); \
    return(CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3, Derived4>, DerivedReturn >(obj.outputSize, 1, obj)); \
  }\
};


// 3 recycling rule arguments and 1 argument whose first value only will be used
#define MAKE_RECYCLING_RULE_CLASS3_1scalar(FUNNAME, RETURNSCALARTYPE) \
  template<typename Index, typename DerivedA1, typename DerivedA2, typename DerivedA3, typename DerivedA4>	\
class FUNNAME ## RecyclingRuleClass { \
public: \
  const DerivedA1 &Arg1;\
  const DerivedA2 &Arg2;\
  const DerivedA3 &Arg3;\
  const DerivedA4 &Arg4;\
  unsigned int size1, size2, size3, outputSize;			\
  FUNNAME ## RecyclingRuleClass(const DerivedA1 &A1, const DerivedA2 &A2, const DerivedA3 &A3, const DerivedA4 &A4 ) : \
  Arg1(A1), Arg2(A2), Arg3(A3), Arg4(A4) \
{\
  outputSize = size1 = nimble_size_impl<DerivedA1>::getSize(Arg1); \
  size2 = nimble_size_impl<DerivedA2>::getSize(Arg2); \
  size3 = nimble_size_impl<DerivedA3>::getSize(Arg3); \
  if(size2 > outputSize) outputSize = size2; \
  if(size3 > outputSize) outputSize = size3; \
  } \
  RETURNSCALARTYPE operator()(Index i) const { \
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA3>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA3, Index >::getCoeff(Arg3, i, size3), \
		   nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedA4>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA4, Index >::getCoeff(Arg4, 0)); \
  }\
  RETURNSCALARTYPE operator()(Index i, Index j) const {\
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA3>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA3, Index >::getCoeff(Arg3, i, size3), \
		   nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedA4>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA4, Index >::getCoeff(Arg4, 0)); \
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
  typedef Eigen::Index IndexReturn;\
  template<typename Derived1, typename Derived2, typename Derived3, typename Derived4>	\
    static CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3, Derived4>, DerivedReturn > \
    FUNNAME ## _RecyclingRule(const Derived1 &A1, const Derived2 &A2, const Derived3 &A3, const Derived4 &A4) { \
    FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3, Derived4> obj(A1, A2, A3, A4); \
    return(CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3, Derived4>, DerivedReturn >(obj.outputSize, 1, obj)); \
  }\
};


// 3 recycling rule arguments and 2 argument whose first value only will be used
#define MAKE_RECYCLING_RULE_CLASS3_2scalar(FUNNAME, RETURNSCALARTYPE) \
  template<typename Index, typename DerivedA1, typename DerivedA2, typename DerivedA3, typename DerivedA4, typename DerivedA5> \
class FUNNAME ## RecyclingRuleClass { \
public: \
  const DerivedA1 &Arg1;\
  const DerivedA2 &Arg2;\
  const DerivedA3 &Arg3;\
  const DerivedA4 &Arg4;\
  const DerivedA5 &Arg5;\
  unsigned int size1, size2, size3, outputSize;			\
  FUNNAME ## RecyclingRuleClass(const DerivedA1 &A1, const DerivedA2 &A2, const DerivedA3 &A3, const DerivedA4 &A4, const DerivedA5 &A5 ) : \
  Arg1(A1), Arg2(A2), Arg3(A3), Arg4(A4), Arg5(A5)			\
{\
  outputSize = size1 = nimble_size_impl<DerivedA1>::getSize(Arg1); \
  size2 = nimble_size_impl<DerivedA2>::getSize(Arg2); \
  size3 = nimble_size_impl<DerivedA3>::getSize(Arg3); \
  if(size2 > outputSize) outputSize = size2; \
  if(size3 > outputSize) outputSize = size3; \
  } \
  RETURNSCALARTYPE operator()(Index i) const { \
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA3>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA3, Index >::getCoeff(Arg3, i, size3), \
		   nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedA4>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA4, Index >::getCoeff(Arg4, 0), \
    		   nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedA5>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA5, Index >::getCoeff(Arg5, 0)); \
  }\
  RETURNSCALARTYPE operator()(Index i, Index j) const {\
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA3>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA3, Index >::getCoeff(Arg3, i, size3), \
		   nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedA4>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA4, Index >::getCoeff(Arg4, 0), \
    		   nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedA5>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA5, Index >::getCoeff(Arg5, 0)); \
  }\
}; \
\
 namespace Eigen{\
  namespace internal{\
    template<typename Index, typename Derived1, typename Derived2, typename Derived3, typename Derived4, typename Derived5> \
      struct functor_has_linear_access<FUNNAME ## RecyclingRuleClass<Index, Derived1, Derived2, Derived3, Derived4, Derived5> > { enum { ret = 1}; }; \
    template<typename Index, typename Derived1, typename Derived2, typename Derived3, typename Derived4, typename Derived5> \
      struct functor_traits<FUNNAME ## RecyclingRuleClass<Index, Derived1, Derived2, Derived3, Derived4, Derived5> > { enum {Cost = 10, PacketAccess = false, IsRepeatable = true}; }; \
  }\
}\
\
template<typename DerivedReturn>\
struct FUNNAME ## _RR_impl {\
  typedef Eigen::Index IndexReturn;\
  template<typename Derived1, typename Derived2, typename Derived3, typename Derived4, typename Derived5>	\
    static CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3, Derived4, Derived5>, DerivedReturn > \
    FUNNAME ## _RecyclingRule(const Derived1 &A1, const Derived2 &A2, const Derived3 &A3, const Derived4 &A4, const Derived5 &A5) { \
    FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3, Derived4, Derived5> obj(A1, A2, A3, A4, A5); \
    return(CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3, Derived4, Derived5>, DerivedReturn >(obj.outputSize, 1, obj)); \
  }\
};

// 4 recycling rule arguments and 1 argument whose first value only will be used
#define MAKE_RECYCLING_RULE_CLASS4_1scalar(FUNNAME, RETURNSCALARTYPE) \
  template<typename Index, typename DerivedA1, typename DerivedA2, typename DerivedA3, typename DerivedA4, typename DerivedA5> \
class FUNNAME ## RecyclingRuleClass { \
public: \
  const DerivedA1 &Arg1;\
  const DerivedA2 &Arg2;\
  const DerivedA3 &Arg3;\
  const DerivedA4 &Arg4;\
  const DerivedA5 &Arg5;\
  unsigned int size1, size2, size3, size4, outputSize;	\
  FUNNAME ## RecyclingRuleClass(const DerivedA1 &A1, const DerivedA2 &A2, const DerivedA3 &A3, const DerivedA4 &A4, const DerivedA5 &A5 ) : \
  Arg1(A1), Arg2(A2), Arg3(A3), Arg4(A4), Arg5(A5)	\
{\
  outputSize = size1 = nimble_size_impl<DerivedA1>::getSize(Arg1); \
  size2 = nimble_size_impl<DerivedA2>::getSize(Arg2); \
  size3 = nimble_size_impl<DerivedA3>::getSize(Arg3); \
  size4 = nimble_size_impl<DerivedA4>::getSize(Arg4); \
  if(size2 > outputSize) outputSize = size2;   \
  if(size3 > outputSize) outputSize = size3; \
  if(size4 > outputSize) outputSize = size4; \
  }					       \
  RETURNSCALARTYPE operator()(Index i) const { \
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA3>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA3, Index >::getCoeff(Arg3, i, size3), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA4>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA4, Index >::getCoeff(Arg4, i, size4), \
		   nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedA5>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA5, Index >::getCoeff(Arg5, 0)); \
  }\
  RETURNSCALARTYPE operator()(Index i, Index j) const {\
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA3>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA3, Index >::getCoeff(Arg3, i, size3), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA4>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA4, Index >::getCoeff(Arg4, i, size4), \
		   nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedA5>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA5, Index >::getCoeff(Arg5, 0)); \
  }\
}; \
\
 namespace Eigen{\
  namespace internal{\
    template<typename Index, typename Derived1, typename Derived2, typename Derived3, typename Derived4, typename Derived5> \
      struct functor_has_linear_access<FUNNAME ## RecyclingRuleClass<Index, Derived1, Derived2, Derived3, Derived4, Derived5> > { enum { ret = 1}; }; \
    template<typename Index, typename Derived1, typename Derived2, typename Derived3, typename Derived4, typename Derived5> \
      struct functor_traits<FUNNAME ## RecyclingRuleClass<Index, Derived1, Derived2, Derived3, Derived4, Derived5> > { enum {Cost = 10, PacketAccess = false, IsRepeatable = true}; }; \
  }\
}\
\
template<typename DerivedReturn>\
struct FUNNAME ## _RR_impl {\
  typedef Eigen::Index IndexReturn;\
  template<typename Derived1, typename Derived2, typename Derived3, typename Derived4, typename Derived5>	\
    static CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3, Derived4, Derived5>, DerivedReturn > \
    FUNNAME ## _RecyclingRule(const Derived1 &A1, const Derived2 &A2, const Derived3 &A3, const Derived4 &A4, const Derived5 &A5) { \
    FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3, Derived4, Derived5> obj(A1, A2, A3, A4, A5); \
    return(CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3, Derived4, Derived5>, DerivedReturn >(obj.outputSize, 1, obj)); \
  }\
};


// 4 recycling rule arguments and 2 argument whose first value only will be used
#define MAKE_RECYCLING_RULE_CLASS4_2scalar(FUNNAME, RETURNSCALARTYPE) \
  template<typename Index, typename DerivedA1, typename DerivedA2, typename DerivedA3, typename DerivedA4, typename DerivedA5, typename DerivedA6> \
class FUNNAME ## RecyclingRuleClass { \
public: \
  const DerivedA1 &Arg1;\
  const DerivedA2 &Arg2;\
  const DerivedA3 &Arg3;\
  const DerivedA4 &Arg4;\
  const DerivedA5 &Arg5;\
  const DerivedA6 &Arg6;\
  unsigned int size1, size2, size3, size4, outputSize;			\
  FUNNAME ## RecyclingRuleClass(const DerivedA1 &A1, const DerivedA2 &A2, const DerivedA3 &A3, const DerivedA4 &A4, const DerivedA5 &A5, const DerivedA6 &A6 ) : \
  Arg1(A1), Arg2(A2), Arg3(A3), Arg4(A4), Arg5(A5), Arg6(A6)		\
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
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA3>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA3, Index >::getCoeff(Arg3, i, size3), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA4>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA4, Index >::getCoeff(Arg4, i, size4), \
    		   nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedA5>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA5, Index >::getCoeff(Arg5, 0), \
                   nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedA6>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA6, Index >::getCoeff(Arg6, 0)); \
  }\
  RETURNSCALARTYPE operator()(Index i, Index j) const {\
    return FUNNAME(nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA1>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA1, Index >::getCoeff(Arg1, i, size1), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA2>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA2, Index >::getCoeff(Arg2, i, size2), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA3>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA3, Index >::getCoeff(Arg3, i, size3), \
		   nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedA4>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA4, Index >::getCoeff(Arg4, i, size4), \
    		   nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedA5>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA5, Index >::getCoeff(Arg5, 0), \
		   nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedA6>::nimbleUseLinearAccess), RETURNSCALARTYPE, DerivedA6, Index >::getCoeff(Arg6, 0)); \
  }\
}; \
\
 namespace Eigen{\
  namespace internal{\
    template<typename Index, typename Derived1, typename Derived2, typename Derived3, typename Derived4, typename Derived5, typename Derived6> \
      struct functor_has_linear_access<FUNNAME ## RecyclingRuleClass<Index, Derived1, Derived2, Derived3, Derived4, Derived5, Derived6> > { enum { ret = 1}; }; \
    template<typename Index, typename Derived1, typename Derived2, typename Derived3, typename Derived4, typename Derived5, typename Derived6> \
      struct functor_traits<FUNNAME ## RecyclingRuleClass<Index, Derived1, Derived2, Derived3, Derived4, Derived5, Derived6> > { enum {Cost = 10, PacketAccess = false, IsRepeatable = true}; }; \
  }\
}\
\
template<typename DerivedReturn>\
struct FUNNAME ## _RR_impl {\
  typedef Eigen::Index IndexReturn;\
  template<typename Derived1, typename Derived2, typename Derived3, typename Derived4, typename Derived5, typename Derived6>	\
    static CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3, Derived4, Derived5, Derived6>, DerivedReturn > \
    FUNNAME ## _RecyclingRule(const Derived1 &A1, const Derived2 &A2, const Derived3 &A3, const Derived4 &A4, const Derived5 &A5, const Derived6 &A6) { \
    FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3, Derived4, Derived5, Derived6> obj(A1, A2, A3, A4, A5, A6); \
    return(CwiseNullaryOp<FUNNAME ## RecyclingRuleClass<IndexReturn, Derived1, Derived2, Derived3, Derived4, Derived5, Derived6>, DerivedReturn >(obj.outputSize, 1, obj)); \
  }\
};


MAKE_RECYCLING_RULE_CLASS2(RRtest_add, double) // only return type is needed here, and correct number of arguments as last digit of macro name

// number of arguments is always 2 more than number of distribution parameters (1 for x and 1 for log)
MAKE_RECYCLING_RULE_CLASS3_1scalar(dbinom, double)
MAKE_RECYCLING_RULE_CLASS2_1scalar(dexp_nimble, double)
MAKE_RECYCLING_RULE_CLASS2_1scalar(dexp, double)
MAKE_RECYCLING_RULE_CLASS3_1scalar(dnbinom, double)
MAKE_RECYCLING_RULE_CLASS2_1scalar(dpois, double)
MAKE_RECYCLING_RULE_CLASS2_1scalar(dchisq, double)
MAKE_RECYCLING_RULE_CLASS3_1scalar(dbeta, double)
MAKE_RECYCLING_RULE_CLASS3_1scalar(dnorm, double)
MAKE_RECYCLING_RULE_CLASS3_1scalar(dgamma, double)
MAKE_RECYCLING_RULE_CLASS3_1scalar(dinvgamma, double)
MAKE_RECYCLING_RULE_CLASS3_1scalar(ddexp, double)
MAKE_RECYCLING_RULE_CLASS3_1scalar(dlnorm, double)
MAKE_RECYCLING_RULE_CLASS3_1scalar(dlogis, double)
MAKE_RECYCLING_RULE_CLASS3_1scalar(dunif, double)
MAKE_RECYCLING_RULE_CLASS3_1scalar(dweibull, double)
MAKE_RECYCLING_RULE_CLASS4_1scalar(dt_nonstandard, double)
MAKE_RECYCLING_RULE_CLASS2_1scalar(dt, double)

MAKE_RECYCLING_RULE_CLASS3_2scalar(pbinom, double)
MAKE_RECYCLING_RULE_CLASS3_2scalar(pnbinom, double)
MAKE_RECYCLING_RULE_CLASS2_2scalar(ppois, double)
MAKE_RECYCLING_RULE_CLASS2_2scalar(pexp_nimble, double)
MAKE_RECYCLING_RULE_CLASS2_2scalar(pexp, double)
MAKE_RECYCLING_RULE_CLASS2_2scalar(pchisq, double)
MAKE_RECYCLING_RULE_CLASS3_2scalar(pbeta, double)
MAKE_RECYCLING_RULE_CLASS3_2scalar(pnorm, double)
MAKE_RECYCLING_RULE_CLASS3_2scalar(pgamma, double)
MAKE_RECYCLING_RULE_CLASS3_2scalar(pinvgamma, double)
MAKE_RECYCLING_RULE_CLASS3_2scalar(pdexp, double)
MAKE_RECYCLING_RULE_CLASS3_2scalar(plnorm, double)
MAKE_RECYCLING_RULE_CLASS3_2scalar(plogis, double)
MAKE_RECYCLING_RULE_CLASS3_2scalar(punif, double)
MAKE_RECYCLING_RULE_CLASS3_2scalar(pweibull, double)
MAKE_RECYCLING_RULE_CLASS4_2scalar(pt_nonstandard, double)
MAKE_RECYCLING_RULE_CLASS2_2scalar(pt, double)

MAKE_RECYCLING_RULE_CLASS3_2scalar(qbinom, double)
MAKE_RECYCLING_RULE_CLASS2_2scalar(qexp_nimble, double)
MAKE_RECYCLING_RULE_CLASS2_2scalar(qexp, double)
MAKE_RECYCLING_RULE_CLASS3_2scalar(qnbinom, double)
MAKE_RECYCLING_RULE_CLASS2_2scalar(qpois, double)
MAKE_RECYCLING_RULE_CLASS2_2scalar(qchisq, double)
MAKE_RECYCLING_RULE_CLASS3_2scalar(qbeta, double)
MAKE_RECYCLING_RULE_CLASS3_2scalar(qnorm, double)
MAKE_RECYCLING_RULE_CLASS3_2scalar(qgamma, double)
MAKE_RECYCLING_RULE_CLASS3_2scalar(qinvgamma, double)
MAKE_RECYCLING_RULE_CLASS3_2scalar(qdexp, double)
MAKE_RECYCLING_RULE_CLASS3_2scalar(qlnorm, double)
MAKE_RECYCLING_RULE_CLASS3_2scalar(qlogis, double)
MAKE_RECYCLING_RULE_CLASS3_2scalar(qunif, double)
MAKE_RECYCLING_RULE_CLASS3_2scalar(qweibull, double)
MAKE_RECYCLING_RULE_CLASS4_2scalar(qt_nonstandard, double)
MAKE_RECYCLING_RULE_CLASS2_2scalar(qt, double)

MAKE_RECYCLING_RULE_CLASS_r2(rbinom, double)
MAKE_RECYCLING_RULE_CLASS_r1(rexp_nimble, double)
MAKE_RECYCLING_RULE_CLASS_r1(rexp, double)
MAKE_RECYCLING_RULE_CLASS_r2(rnbinom, double)
MAKE_RECYCLING_RULE_CLASS_r1(rpois, double)
MAKE_RECYCLING_RULE_CLASS_r1(rchisq, double)
MAKE_RECYCLING_RULE_CLASS_r2(rbeta, double)
MAKE_RECYCLING_RULE_CLASS_r2(rnorm, double)
MAKE_RECYCLING_RULE_CLASS_r2(rgamma, double)
MAKE_RECYCLING_RULE_CLASS_r2(rinvgamma, double)
MAKE_RECYCLING_RULE_CLASS_r2(rdexp, double)
MAKE_RECYCLING_RULE_CLASS_r2(rlnorm, double)
MAKE_RECYCLING_RULE_CLASS_r2(rlogis, double)
MAKE_RECYCLING_RULE_CLASS_r2(runif, double)
MAKE_RECYCLING_RULE_CLASS_r2(rweibull, double)
MAKE_RECYCLING_RULE_CLASS_r3(rt_nonstandard, double)
MAKE_RECYCLING_RULE_CLASS_r1(rt, double)

MAKE_RECYCLING_RULE_CLASS2_1scalar(bessel_k, double)
MAKE_RECYCLING_RULE_CLASS1_1scalar(pow_int, double)

// matrix, array, as.numeric, as.matrix, as.array

// need the additional parts below to make newMatrixClass work
// for array, will need to make something in nimbleEigenNimArr and always lift to a full copy operation
// can write as.numeric here
// as.matrix may be basically the same as matrix

template<typename Index, typename DerivedInput, typename result_type = double>
  class newMatrixClass {
 public:
  const DerivedInput &input;
  int dim1, dim2, totalLength, inputLength, inputRows;
  bool init; // would be a bit silly to call with init = FALSE, but it is allowed to simplify code generation
  bool recycle;
  //  typedef double result_type;
 newMatrixClass(const DerivedInput &inputIn, bool initIn, bool recycleIn, int rowsIn, int colsIn) :
  input(inputIn),
    init(initIn),
    recycle(recycleIn) {
      inputLength = nimble_size_impl<DerivedInput>::getSize(input);
      inputRows = nimble_size_impl<DerivedInput>::getRows(input);
      bool rowsProvided = rowsIn > 0;
      bool colsProvided = colsIn > 0;
      if(!rowsProvided) {
	if(!colsProvided) {
	  dim1 = inputLength;
	  dim2 = 1;
	} else {
	  dim2 = colsIn;
	  dim1 = floor((double(inputLength)-1) / double(colsIn)) + 1;
	}
      } else {
	if(!colsProvided) {
	  dim1 = rowsIn;
	  dim2 = floor((double(inputLength)-1) / double(rowsIn)) + 1;
	} else {
	  dim1 = rowsIn;
	  dim2 = colsIn;
	}
      }
      totalLength = dim1 * dim2;
    }
  result_type operator()(Index i) const 
  {
    if(init) {
      if(recycle) {
	return nimble_eigen_coeff_mod_impl< bool(nimble_eigen_traits<DerivedInput>::nimbleUseLinearAccess), result_type, DerivedInput, Index >::getCoeff(input, i, inputLength);
      } else {
	if(static_cast<int>(i) < inputLength) {
	  return nimble_eigen_coeff_impl< bool(nimble_eigen_traits<DerivedInput>::nimbleUseLinearAccess), result_type, DerivedInput, Index >::getCoeff(input, i);
	} else {
	  return 0;
	}
      }
    } else {
      return 0;
    }
  }

  result_type operator()(Index i, Index j) const // I don't think this should normally be called, but if it does, act like a vector
  {
    if(init) 
      return operator()(i + j*dim1);
    return 0;
  }
};

namespace Eigen{
  namespace internal{
    template<typename IndexObj, typename DerivedObj, typename result_type>
      struct functor_has_linear_access<newMatrixClass<IndexObj, DerivedObj, result_type> > { enum { ret = 1}; }; 
    template<typename IndexObj, typename DerivedObj, typename result_type>
      struct functor_traits<newMatrixClass<IndexObj, DerivedObj, result_type> > { enum {Cost = 10, PacketAccess = false, IsRepeatable = true }; }; 
  }
}

template<typename returnDerived, typename result_type = double>
struct newMatrix_impl {
  typedef Eigen::Index IndexReturn;
  template<typename DerivedObj>
  static CwiseNullaryOp<newMatrixClass<IndexReturn, DerivedObj, result_type >, returnDerived > newMatrix(const DerivedObj &s, bool initIn, bool recycle, int nRowIn, int nColIn) {
    newMatrixClass<IndexReturn, DerivedObj, result_type > obj(s, initIn, recycle, nRowIn, nColIn);
    return(CwiseNullaryOp<newMatrixClass<IndexReturn, DerivedObj, result_type >, returnDerived >(obj.dim1, obj.dim2, obj));
  }
};

#define nimNewMatrixD newMatrix_impl<MatrixXd>::newMatrix
#define nimNewMatrixI newMatrix_impl<MatrixXi>::newMatrix
#define nimNewMatrixB newMatrix_impl<MatrixXb>::newMatrix

#endif
