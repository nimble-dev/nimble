// This file is here in includes to make it distinct.
// It is not compiled into libnimble.a but rather copied
// to the on-the-fly working directly (normally tempdir())
// and compiled once-per-session there if AD is used.
#include <math.h>
#include <iostream>
#include <nimble/nimbleCppADbaseClass.h>
#include <nimble/NimArr.h>
#include <nimble/nimDerivs_dists.h>

#define USE_CPPAD_OPTIMIZE_FOR_MODEL_TAPES // comment this out to turn off atomics for nimDerivs(model$calculate(...),...)

void setOrdersFound(const NimArr<1, double> &derivOrders,
		    bool *ordersFound,
		    int &maxOrder) {
  /* ordersFound must have length (at least) 3.*/
  int orderSize = derivOrders.size();
  double const* array_derivOrders = derivOrders.getConstPtr();
  maxOrder = 0;
  maxOrder = *std::max_element(array_derivOrders, array_derivOrders + orderSize);
  std::fill(ordersFound, ordersFound + 3, false);
  for (int i = 0; i < orderSize; i++) {
    if ((array_derivOrders[i] > 2) | (array_derivOrders[i] < 0)) {
      printf("Error: Derivative orders must be between 0 and 2.\n");
    }
    ordersFound[static_cast<int>(array_derivOrders[i])] = true;
  }
}

bool check_inf_nan_gdi(double v) {
  if(((v == -std::numeric_limits<double>::infinity()) ||
      (v == std::numeric_limits<double>::infinity())) ||
     (std::isnan(v))) {
    return true;
  }
  return false;
}

bool check_inf_nan_gdi(CppAD::AD<double> v) {
  return false;
}

template<typename BASE>
inline double gd_getValue(const BASE& x) {
  return x;
}

template<>
inline double gd_getValue(const CppAD::AD<double>& x) {
  return CppAD::Value(x);
}

template<typename BASE, class TAPETYPE, class ADCLASS>
void getDerivs_internal(vector<BASE> &independentVars,
  TAPETYPE *ADtape,
  const NimArr<1, double> &derivOrders,
  const NimArr<1, double> &wrtVector,
  const NimArr<1, double> &outInds,
  const NimArr<1, BASE> &inDir,
  const NimArr<1, BASE> &outDir,
  nimSmartPtr<ADCLASS> &ansList) {
    // std::cout<<"entering getDerivs_internal"<<std::endl;
    // std::cout<<"independentVars: ";
    // for(size_t ijk = 0; ijk < independentVars.size(); ++ijk)
    //   std::cout<<independentVars[ijk]<<" ";
    // std::cout<<std::endl;
    // std::cout<<"derivOrders: ";
    //  for(size_t ijk = 0; ijk < derivOrders.size(); ++ijk)
    //    std::cout<<derivOrders[ijk]<<" ";
    //  std::cout<<std::endl;
    //  std::cout<<"wrtVector: ";
    //  for(size_t ijk = 0; ijk < wrtVector.size(); ++ijk)
    //    std::cout<<wrtVector[ijk]<<" ";
    //  std::cout<<std::endl;

    #ifdef _TIME_AD_GENERAL
    derivs_getDerivs_timer_start();
    derivs_tick_id();
    derivs_show_id();
    #endif
    std::size_t n = independentVars.size();  // dim of independent vars
    std::size_t m = ADtape->Range();

    int maxOrder;
    bool ordersFound[3];
    setOrdersFound(derivOrders, ordersFound, maxOrder);

    // std::cout<<"orders: "<<ordersFound[0]<<" "<<ordersFound[1]<<" "<<ordersFound[2]<<std::endl;
    // std::cout<<"maxOrder = "<<maxOrder<<std::endl;

    // std::cout<<"inDir info: size = "<<inDir.size()<<std::endl;
    // if(inDir.size() > 0) {
    //   std::cout<<"inDir[0] = "<<gd_getValue(inDir[0])<<std::endl;
    // } else {
    //   std::cout<<"inDir is empty."<<std::endl;
    // }

    bool outInds_provided = outInds.size() > 0 && outInds[0] != -1;
    bool inDir_provided = inDir.size() > 0 && !ISNA(gd_getValue(inDir[0]));
    bool outDir_provided = outDir.size() > 0 && !ISNA(gd_getValue(outDir[0]));
    bool wrt_provided = wrtVector.size() > 0 && wrtVector[0] != -1;
    //std::cout<<"outInds_provided: "<<outInds_provided<<" inDir_provided: "<<inDir_provided<<
    //  " outDir_provided: "<<outDir_provided<<" wrt_provided: "<<wrt_provided<<std::endl;
    /*
      Behavior of wrt, inDir, outInds and outDir (including some curious combinations):
      If inDir is provided, it is used for Forward order 1, if that is used.
        This results in the returned "x" dimension being 1.
      If outDir is provided, it is used for Reverse order 1 or 2, whichever is used (higher), if either.
           (If maxOrder==1, it may be done by Forward or Reverse, whereas maxOrder==2 always invokes Reverse.)
        This results in the returned "y" dimension being 1.
      If wrt is provided (conceptually, this could be called "inputInds") *without* inDir,
        it is used to filter result "x" dimensions for any order. This means
        for Forward mode only basis vectors in wrt directions will be used, and in
        Reverse mode, results will be filtered to only wrt directions.
        The returned "x" dimensions will be length(wrt).
      If wrt is provided *with* inDir,
        it is used to filter "x" dimension results from Reverse mode, whereas
        Forward mode will go in the inDir direction for "x". This will result in
        Reverse 1 returning "x" dimension 1 and Reverse 2 returning "x" dimensions (1, length(wrt)).
      If inDir *is* provided *without* wrt,
        if order includes 1, the Jacobian "x" dimension will be 1.
        If order includes 2, the Hessian "x" dimensions will be (1, n).
      If neither inDir nor wrt are provided,
        wrt defaults to 1:n.
      If outInds is provided *without* outDir,
        it is used to filter result "y" dimensions for any order (including 0). This means
        for Reverse mode only basis vectors in outInds directions will be used, and in
        Forward mode, results will be filtered to only outInds directions.
        The returned "y" dimensions will be length(outInds).
      If outInds is provided *with* outDir,
        it is used to filter "y" dimension results from Forward mode, whereas
        Reverse mode will go in the outDir direction for "y". This will result in
        Forward 1 returning "y" dimension length(outInds) and Reverse (1 or 2) returning "y" dimension 1.
      If outDir is provided *without* outInds,
        if maxOrder==1, the Jacobian "y" dimension will be 1.
        if maxOrder==2, the Jacobian "y" dimension (if requested) will be m.
        If maxOrder==2, the Hessian "y" dimension will be 1.
      N.B.: When maxOrder == 1, either Forward 1 or Reverse 1 (possibly with subgraph) method
        may be used for the Jacobian.
        If inDir is provided, it will force use of Forward 1 in the inDir direction.
        Otherwise, if outDir is provided, it will force use of Reverse 1 in the outDir direction.
        Otherwise, an automated decision will be made based on n and m.
    */

    if((inDir_provided || outDir_provided) && maxOrder == 0) {
      std::cout<<"Warning: inDir and/or outDir are provided, but only order 0 was requested, so inDir and outDir will be ignored."<<std::endl;
    }
    if(inDir_provided && outDir_provided && maxOrder == 1) {
      std::cout<<"Warning: Both inDir and outDir are provided, but order 2 was not requested, so outDir will be ignored."<<std::endl;
      outDir_provided = false;
    }
    if((wrt_provided && inDir_provided) && maxOrder != 2) {
      std::cout<<"Warning: Both wrt and inDir are provided, which is only meaningful for order 2 (not requested), so wrt will be ignored."<<std::endl;
      wrt_provided = false;
    }
    size_t num_inDirs;
    if(inDir_provided) {
      if( (inDir.size() % n) != 0)  { //inDir.size() != n) {
        std::cout << "length of inDir must be a multiple of length of independentVars." << std::endl;
      }
      num_inDirs = inDir.size() / n;
    }
    if(outDir_provided){
      if(outDir.size() != m) {
        std::cout << "outDir must have the same length as output." << std::endl;
      }
    }

    std::size_t length_wrt = wrtVector.size();            // dim of wrt vars
    if(length_wrt == 2){
      if(wrtVector[1] == -1){ // 2nd element -1 is a filler to ensure there is a vector out of compilation
        // The other input vectors were added at a later stage and don't need this kluge.
        length_wrt = 1;
      }
    }
    std::size_t length_outInds = outInds.size(); // dim of outInds

    // defaults are for all of wrt, outInds, inDir and outDir provided EMPTY
    size_t res_dimx_o1(n), res_dimx1_o2(n), res_dimx2_o2(n);
    size_t res_dimy_o0(m), res_dimy_o1(m), res_dimy_o2(m);
    // N.B. The "All" flags may not be checked when using inDir or outDir, but they are set for logical coherence anyway.
    bool wrtAllx_o1(true), wrtAllx1_o2(true), wrtAllx2_o2(true);
    bool outAlly_o0(true), outAlly_o1(true), outAlly_o2(true);
    if(!wrt_provided) {
      if(inDir_provided) {
        wrtAllx_o1 = false;  // not checked anyway
        wrtAllx1_o2 = false; // not checked anyway
        res_dimx_o1 = num_inDirs; //1
        res_dimx1_o2 = num_inDirs; //1;
      }
    } else { // wrt_provided
      wrtAllx_o1 = false;
      wrtAllx1_o2 = false;
      wrtAllx2_o2 = false;
      res_dimx2_o2 = length_wrt;
      if(inDir_provided) { // inDir_provided
        res_dimx_o1 = num_inDirs; //1;
        res_dimx1_o2 = num_inDirs; //1;
      } else {
        res_dimx_o1 = length_wrt;
        res_dimx1_o2 = length_wrt;
      }
    }
    if(outInds_provided) {
      outAlly_o0 = false;
      outAlly_o1 = false;
      outAlly_o2 = false;
      res_dimy_o0 = length_outInds;
      res_dimy_o1 = length_outInds;
      res_dimy_o2 = length_outInds;
    }
    if(outDir_provided) {
      if(maxOrder == 2) {
        outAlly_o2 = false;
        res_dimy_o2 = 1;
      }
      if(maxOrder == 1) {
        outAlly_o1 = false;
        res_dimy_o1 = 1;
      }
    }

    vector<BASE> value_ans;
    #ifdef _TIME_AD_GENERAL
    derivs_run_tape_timer_start();
    #endif

    int strategy_zero = 1; // 1 means do it, 0 means skip.
    int strategy_one = 0; // 0=skip. 1=forward regular. 2=reverse subgraph_rev_jac. 3=reverse normal
    int strategy_two= 0; // 0=skip. 1=reverse regular.
    if(maxOrder > 0) {
      if (ordersFound[1]) {
        // At least Jacobian is needed.
        if(!ordersFound[2]) {
          // Hessian is not needed.
          if((m <= n || outDir_provided) && !inDir_provided) { // possibly we should compare to wrt_n, but we have use cases where comparing to n works better
            // It's hard to really know, because reverse offers the subgraph_rev_jac option.
            // Use reverse mode b/c we have fewer (or equal) number of outputs as inputs, or because outDir is provided.
            if(m > 1 && !outDir_provided) {
              // If there are multiple outputs, use subgraph_rev_jac
              strategy_one = 2; // use subgraph_rev_jac if jacobian is the final order
              if(!ordersFound[0]) {
                // If value is not needed, we can skip Forward 0 because subgraph_rev_jac does it itself.
                strategy_zero = 0;
              }
            } else {
              // m == 1, or outDir_provided, so we can use regular reverse mode.
              strategy_one = 3;// simple reverse
            }
          } else {
            // Use forward mode b/c we have more outputs than inputs, or inDir was provided
            strategy_one = 1; // use forward
          }
        } else {
          // Use Forward mode because we need to continue to Hessian.
          // In this case, strategy_one is not actually checked below, but we set it here for consistency.
          strategy_one = 1;
        }
      } // done with ordersFound[1]==true
      if(ordersFound[2]) // Hessian is needed.
      strategy_two = 1; // not really used below, but set up for consistency.
    }
    // std::cout<<"n: "<<n<<" m: "<<m<<" length_wrt: "<<length_wrt<<std::endl;
    // std::cout<<"strategy flags: "<<strategy_zero<<" "<<strategy_one<<" "<<strategy_two<<std::endl;

    // to-do: get Forward(0) out of subgraph_rev_jac if both are needed/used.
    // For now we have to wastefully run Forward(0) even if it will be done again in subgraph_rev_jac.
    if(strategy_zero==1) {
      value_ans = ADtape->Forward(0, independentVars);
      //  std::cout<<"value_ans.size() = "<<value_ans.size()<<std::endl;
      #ifdef _TIME_AD_GENERAL
      derivs_run_tape_timer_stop();
      #endif
      if (ordersFound[0]) {
        ansList->value.setSize(res_dimy_o0, false, false);
        if(outAlly_o0) {
          std::copy(value_ans.begin(), value_ans.end(), ansList->value.getPtr());
        } else {
          BASE *LHS = ansList->value.getPtr();
          for(size_t iii=0;iii<res_dimy_o0;iii++,LHS++) {
            *LHS = value_ans[outInds[iii]-1];
          }
        }
      }
    } else {
      // Initialize value_ans with size m and all zeros when strategy_zero == 0
      value_ans.resize(m);
      std::fill(value_ans.begin(), value_ans.end(), 0.0);
    }
    if(maxOrder > 0){
      // std::cout<<"entering maxOrder> 0\n";
      vector<bool> infIndicators(m, false); // default values will be false
      for(size_t inf_ind = 0; inf_ind < m; inf_ind++){
        if(check_inf_nan_gdi(value_ans[inf_ind])) {
          infIndicators[inf_ind] = true;
        }
      }
      //std::cout<<"about to set jacobian and hessian sizes\n";
      if (ordersFound[1]) {
        ansList->jacobian.setSize(res_dimy_o1, res_dimx_o1, false, true); // true for fill_zero
      }
      //std::cout<<"done setting jacobian size\n";
      if (ordersFound[2]) {
        ansList->hessian.setSize(res_dimx1_o2, res_dimx2_o2, res_dimy_o2, false, true);
      }
      //std::cout<<"done setting hessian size\n";
      vector<BASE> cppad_derivOut;
      std::vector<BASE> w(m, 0);

      // begin replacement
      if (maxOrder == 1) {
        if(strategy_one == 2) { // subgraph_rev_jac
          // Reverse mode
          //   std::cout<<"subgraph_rev_jac version of jacobian\n";
          //   std::cout<<"select_domain: ";
          //   for (size_t i = 0; i < select_domain.size(); i++) {
          //     std::cout<<select_domain[i]<<" ";
          //   }
          //   std::cout<<std::endl;
          //   std::cout<<"select_range: ";
          //   for (size_t i = 0; i < select_range.size(); i++) {
          //     std::cout<<select_range[i]<<" ";
          //   }
          //   std::cout<<std::endl;
          // Use subgraph_rev_jac system
          std::vector<bool> select_domain(n, false);
          for (size_t i_wrt = 0; i_wrt < res_dimx_o1; i_wrt++) {
            int x_ind = wrtAllx_o1 ? i_wrt : wrtVector[i_wrt] - 1;
            select_domain[x_ind] = true;
          }
          std::vector<bool> select_range(m, false);
          for (size_t i_out = 0; i_out < res_dimy_o1; i_out++) {
            int y_ind = outAlly_o1 ? i_out : outInds[i_out] - 1;
            select_range[y_ind] = true;
          }

          std::vector<int> col_2_wrtVecm1;
          if(!wrtAllx_o1) {
            col_2_wrtVecm1.resize(n, -1);
            for(size_t i_wrt = 0; i_wrt < length_wrt; i_wrt++) {
              col_2_wrtVecm1[wrtVector[i_wrt] - 1] = i_wrt;
            }
          }
          std::vector<int> row_2_outIndm1;
          if(!outAlly_o1) {
            row_2_outIndm1.resize(m, -1);
            for(size_t i_out = 0; i_out < length_outInds; i_out++) {
              row_2_outIndm1[outInds[i_out] - 1] = i_out;
            }
          }
          // std::cout<<"done setting select_domain and select_range\n";
          CppAD::sparse_rcv<std::vector<size_t> , std::vector<BASE> > matrix_out;
          ADtape->subgraph_jac_rev(select_domain, select_range,
                                  independentVars, matrix_out);
            // std::cout<<"done with subgraph_jac_rev\n";
          size_t nnz = matrix_out.nnz();
          std::vector<size_t> row_major = matrix_out.row_major();
          size_t this_row, this_col, this_ind;
          BASE *LHS = ansList->jacobian.getPtr();
            // size_t max_row(0), max_col(0);
            // bool give_output = false;
          for(size_t k = 0; k < nnz; k++)
          {
            this_ind = row_major[k];
            this_row = matrix_out.row()[this_ind];
            this_col = matrix_out.col()[this_ind];
            int x1_ind = wrtAllx_o1 ? this_col : col_2_wrtVecm1[this_col];
            if(x1_ind < 0) {
              std::cout<<"Error: x1_ind is negative. This should not happen."<<std::endl;
              continue; // skip this entry
            }
            int y_ind = outAlly_o1 ? this_row : row_2_outIndm1[this_row];
            if(y_ind < 0) {
              std::cout<<"Error: y_ind is negative. This should not happen."<<std::endl;
              continue; // skip this entry
            }
            LHS[x1_ind*m + y_ind] = matrix_out.val()[this_ind];
              // if(this_row > max_row) {
              //   max_row = this_row;
              //   give_output = true;
              // }
              // if(this_col > max_col) {
              //   max_col = this_col;
              //   give_output = true;
              // }
              // if(give_output) {
              //   give_output = false;
              // //  std::cout<<"this_ind: "<<this_ind<<" this_row: "<<this_row<<" this_col: "<<this_col<<" dx1_ind: "<<dx1_ind<<std::endl;
              // }
              // std::cout<<"this_ind: "<<this_ind<<" this_row: "<<this_row<<
              // " this_col: "<<this_col<<" dx1_ind: "<<dx1_ind<<" value: "<<
              // matrix_out.val()[this_ind]<<std::endl;
          } // end k
            //std::cout<<"done filling jacobian\n";
        } // end strategy_one == 2
          else if(strategy_one == 3)
        { // regular reverse mode, used for m == 1 or outDir_provided
          // std::cout<<"regular reverse version of jacobian\n";
          for(size_t i_out = 0; i_out < res_dimy_o1; i_out++) {
            int y_ind(0);
            if(!outDir_provided) {
              y_ind = outAlly_o1 ? i_out : outInds[i_out] - 1;
              w[y_ind] = 1;
              if(!infIndicators[i_out]){
                  #ifdef _TIME_AD_GENERAL
                  derivs_run_tape_timer_start();
                  #endif
                cppad_derivOut = ADtape->Reverse(1, w);
                  #ifdef _TIME_AD_GENERAL
                  derivs_run_tape_timer_stop();
                  #endif
              } else {
                cppad_derivOut.resize(maxOrder*n, 0.0);
              }
            } else // end use_out
            {
              y_ind = 0; // outDir is used, so we only need one output.
              for(size_t iy = 0; iy < m; iy++) w[iy] = outDir[iy];
              cppad_derivOut = ADtape->Reverse(1, w);
            } // end use outDir.
            if (ordersFound[1]) { // will always be true if maxOrder == 1
              BASE *LHS = ansList->jacobian.getPtr() + i_out;
              if(!infIndicators[y_ind]){
                for (size_t i_wrt = 0; i_wrt < res_dimx_o1; ++i_wrt, LHS += res_dimy_o1) {
                  int x_ind = wrtAllx_o1 ? i_wrt : wrtVector[i_wrt] - 1;
                  *LHS = cppad_derivOut[y_ind + x_ind * res_dimy_o1];
                }
              } else {
                for (size_t i_wrt = 0; i_wrt < res_dimx_o1; i_wrt++) {
                  *LHS = CppAD::numeric_limits<BASE>::quiet_NaN();
                  LHS += outAlly_o1;
                }
              }
            }
            if(!outDir_provided) w[y_ind] = 0;
          } // end out_ind loop
          // end strategy_one == 3
        } else if(strategy_one == 1)
        { // m > n
            // Forward mode
            //std::cout<<"new version of jacobian\n";
          std::vector<BASE> x1(n, 0);
          for (size_t i_wrt = 0; i_wrt < res_dimx_o1; i_wrt++) {
            int x1_ind(0);
            if(!inDir_provided) {
              x1_ind = wrtAllx_o1 ? i_wrt : wrtVector[i_wrt] - 1;
              x1[x1_ind] = 1;
                #ifdef _TIME_AD_GENERAL
                derivs_run_tape_timer_start();
                #endif
              cppad_derivOut = ADtape->Forward(1, x1);
                #ifdef _TIME_AD_GENERAL
                derivs_run_tape_timer_stop();
                #endif
            } else { // end use_wrt
              x1_ind = 0; // inDir is used, so we only need one input.
              for(size_t ix = 0; ix < n; ix++) x1[ix] = inDir[n*i_wrt + ix];
              cppad_derivOut = ADtape->Forward(1, x1);
            } // end use inDir.
            if (ordersFound[1]) { // will always be true if maxOrder == 1
              BASE *LHS = ansList->jacobian.getPtr() + res_dimy_o1*i_wrt;

              for(size_t i_out = 0; i_out < res_dimy_o1; ++i_out, ++LHS) {
                int y_ind = outAlly_o1 ? i_out : outInds[i_out] - 1;
                *LHS = infIndicators[y_ind] ?
                  CppAD::numeric_limits<BASE>::quiet_NaN() :
                  cppad_derivOut[y_ind];
              }
            }
            if(!inDir_provided) x1[x1_ind] = 0;
          }// end vec_ind loop
           // end strategy_one == 1
        } else {
          std::cout<<"Error: in getDerivs: strategy_one is not set correctly for maxOrder == 1."<<std::endl;
        }
        // end maxOrder == 1
      } else {
       // std::cout<<"Entering maxOrder > 1\n";
        // maxOrder > 1: outer loop over vec_ind, inner loop over dy_ind
        // strategy_one and strategy_two are actually ignored here for now, b/c we wouldn't be here if not needed.
        std::vector<BASE> cppad_derivOut_F1;
        std::vector<BASE> x1(n, 0);
        for (size_t i1_wrt = 0; i1_wrt < res_dimx1_o2; i1_wrt++) {
          int x1_ind(0);
          if(!inDir_provided) {
            x1_ind = wrtAllx1_o2 ? i1_wrt : wrtVector[i1_wrt] - 1;
            x1[x1_ind] = 1;
             #ifdef _TIME_AD_GENERAL
              derivs_run_tape_timer_start();
              #endif
            cppad_derivOut_F1 = ADtape->Forward(1, x1);
              #ifdef _TIME_AD_GENERAL
              derivs_run_tape_timer_stop();
              #endif
          } else { // end use_wrt
            x1_ind = 0; // inDir is used, so we only need one input.
            for(size_t ix = 0; ix < n; ix++) x1[ix] = inDir[n*i1_wrt + ix];
            cppad_derivOut_F1 = ADtape->Forward(1, x1);
          } // end use inDir.
          // std::cout<<"contents of cppad_derivOut_F1: ";
          // for(size_t ijk = 0; ijk < cppad_derivOut_F1.size(); ++ijk)
          //   std::cout<<cppad_derivOut_F1[ijk]<<" ";
          // std::cout<<std::endl;

          for (size_t i_out = 0; i_out < res_dimy_o2; i_out++) {
            int y_ind(0);
            if(!outDir_provided) {
              y_ind = outAlly_o2 ? i_out : outInds[i_out] - 1;
              w[y_ind] = 1;
              if(!infIndicators[y_ind]){
                  #ifdef _TIME_AD_GENERAL
                  derivs_run_tape_timer_start();
                  #endif
                cppad_derivOut = ADtape->Reverse(2, w);
                  #ifdef _TIME_AD_GENERAL
                  derivs_run_tape_timer_stop();
                  #endif
              } else {
                cppad_derivOut.resize(maxOrder * n, 0.0);
              }
            } else { // outDir_provided is true
              y_ind = 0; // outDir is used, so we only need one output.
              for(size_t iy = 0; iy < m; iy++) w[iy] = outDir[iy];
                #ifdef _TIME_AD_GENERAL
                derivs_run_tape_timer_start();
                #endif
              cppad_derivOut = ADtape->Reverse(2, w);
                #ifdef _TIME_AD_GENERAL
                derivs_run_tape_timer_stop();
                #endif
            } // end outDir_provided
              // record Hessian outputs
            if(!infIndicators[y_ind]){
              for (size_t i2_wrt = 0; i2_wrt < res_dimx2_o2; i2_wrt++) {
                int x2_ind(0);
                x2_ind = wrtAllx2_o2 ? i2_wrt : wrtVector[i2_wrt] - 1;

                ansList->hessian[res_dimx1_o2 * res_dimx2_o2 * i_out + res_dimx1_o2 * i2_wrt + i1_wrt] =
                cppad_derivOut[x2_ind * 2 + 1];
              }
            } else {
              for (size_t i2_wrt = 0; i2_wrt < res_dimx2_o2; i2_wrt++) {
                ansList->hessian[res_dimx1_o2 * res_dimx2_o2 * i_out + res_dimx1_o2 * i2_wrt + i1_wrt] =
                CppAD::numeric_limits<BASE>::quiet_NaN();
              }
            }
              // end use_out2
            if(!outDir_provided) w[y_ind] = 0; // outDir is used, so we only need one output.
          } // end out_ind loop

          if (ordersFound[1]) {
            BASE *LHS = ansList->jacobian.getPtr() + res_dimy_o1 * i1_wrt;

            for(size_t i_out = 0; i_out < res_dimy_o1; ++i_out, ++LHS) {
              int y_ind = outAlly_o1 ? i_out : outInds[i_out] - 1;
              *LHS = infIndicators[y_ind] ?
               CppAD::numeric_limits<BASE>::quiet_NaN() :
               cppad_derivOut_F1[y_ind];
            }
          }
          if(!inDir_provided) x1[x1_ind] = 0;
        }
      }
    } // end else

      // end replacement
      #ifdef _TIME_AD_GENERAL
      derivs_getDerivs_timer_stop();
      #endif
  }

void nimbleFunctionCppADbase::getDerivs_meta(nimbleCppADinfoClass &ADinfo,
					     const NimArr<1, double> &derivOrders,
					     const NimArr<1, double> &wrtVector,
                const NimArr<1, double> &outInds,
                const NimArr<1, CppAD::AD<double> > &inDir,
                const NimArr<1, CppAD::AD<double> > &outDir,
					     const nimbleCppADrecordingInfoClass &nimRecInfo,
					     nimSmartPtr<NIMBLE_ADCLASS_META> &ansList) {
  //  std::cout<<"Entering getDerivs_meta"<<std::endl;
  //  std::cout<<"ADinfo is at :"<< &ADinfo <<"\n";

  //  if(!nimRecInfo.recording_cp()) return;

  bool orderIncludesZero(false);
  for(int i = 0; i < derivOrders.size(); ++i) {
    orderIncludesZero |= (derivOrders[i] == 0);
  }
  // std::cout << "orderIncludesZero = " << orderIncludesZero << std::endl;
  bool oldUpdateModel = ADinfo.updateModel();
  ADinfo.updateModel() = orderIncludesZero;

  // We need to use the tricks to have CppAD statics in (potentially) two compilation units
  // that were not linked together (one for model, one for the current algorithm)
  // be matching.  This makes it so the CppAD tape handle and atomic information are shared.
  // This is necessary even in this double taping step because of our atomic classes,
  // which might reside in the other compilation unit.  During double taping, an atomic
  // such as lgamma puts itself or other statics onto the new tape, and returns
  // CppAD::AD variables created in the other compilation unit.
  set_CppAD_tape_info_for_model my_tape_info_RAII_; // must be at function scope, not declared inside if(){}

  if(ADinfo.nodeFunPtrSet()) {
    //    std::cout<<"tape_id and handle:"<<  nimRecInfo.tape_id_cp() <<" "<< nimRecInfo.tape_handle_cp() <<"\n";
    //   std::cout<<"atomic info:"<<nimRecInfo.atomic_vec_ptr_cp()<<"\n";
    my_tape_info_RAII_.set_from_nodeFunPtr(ADinfo.nodeFunPtr(),
					   nimRecInfo.tape_id_cp(), //CppAD::AD<double>::get_tape_id_nimble(),
					   nimRecInfo.tape_handle_cp());//CppAD::AD<double>::get_tape_handle_nimble());
    set_CppAD_atomic_info_for_model(ADinfo.nodeFunPtr(),
				    nimRecInfo.atomic_vec_ptr_cp());
    //    std::cout<<"done setting nodeFunPtr\n";
  }

  CppAD::ADFun< CppAD::AD<double>, double > innerTape;
  innerTape = ADinfo.ADtape()->base2ad();
  innerTape.new_dynamic(ADinfo.dynamicVars_meta);

  //  std::cout<<" after making inner tape\n";
  //  std::cout<<"tape_id and handle:"<< CppAD::AD<double>::get_tape_id_nimble()<<" "<< CppAD::AD<double>::get_tape_handle_nimble()<<"\n";
  //  std::cout<<"atomic info:"<<CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage()<<"\n";
  getDerivs_internal< CppAD::AD<double>,
		      CppAD::ADFun< CppAD::AD<double>, double >,
		      NIMBLE_ADCLASS_META>(ADinfo.independentVars_meta,
					   &innerTape,
					   derivOrders,
					   wrtVector,
              outInds,
					   inDir,
					   outDir,
					   ansList);
  ADinfo.updateModel() = oldUpdateModel;
  //  std::cout<<"Exiting getDerivs_meta"<<std::endl;
}

void nimbleFunctionCppADbase::getDerivs(nimbleCppADinfoClass &ADinfo,
                                        const NimArr<1, double> &derivOrders,
                                        const NimArr<1, double> &wrtVector,
                                        const NimArr<1, double> &outInds,
                                        const NimArr<1, double> &inDir,
                                        const NimArr<1, double> &outDir,
                                        nimSmartPtr<NIMBLE_ADCLASS> &ansList) {
  // std::cout<<"Entering getDerivs"<<std::endl;
  bool orderIncludesZero(false);
  for(int i = 0; i < derivOrders.size(); ++ i) {
    orderIncludesZero |= (derivOrders[i] == 0);
  }
  //  std::cout << "orderIncludesZero = " << orderIncludesZero << std::endl;
  bool oldUpdateModel = ADinfo.updateModel();
  ADinfo.updateModel() = orderIncludesZero;
  getDerivs_internal<double,
                     CppAD::ADFun<double>,
                     NIMBLE_ADCLASS>(ADinfo.independentVars,
                                     ADinfo.ADtape(),
                                     derivOrders,
                                     wrtVector,
                                     outInds,
                                     inDir,
                                     outDir,
                                     ansList);
  ADinfo.updateModel() = oldUpdateModel;
  //  std::cout<<"Exiting getDerivs"<<std::endl;
}


CppAD::ADFun<double>* calculate_recordTape(NodeVectorClassNew_derivs &NV,
                                           bool includeExtraOutputs,
                                           nimbleCppADinfoClass &ADinfo) {
  vector< CppAD::AD<double> > dependentVars(1);
  NimArr<1, double> NimArrValues;
  NimArr<1, CppAD::AD<double> > NimArrValues_AD;

  // 1. Copy all constantNodes values from model -> model_AD
  int length_constant = NV.model_constant_accessor.getTotalLength();
  if(length_constant > 0) {
    NimArr<1, double> NimArrValues;
    NimArr<1, CppAD::AD<double> > NimArrValues_AD;
    NimArrValues.setSize(length_constant);
    NimArrValues_AD.setSize(length_constant);
    getValues(NimArrValues, NV.model_constant_accessor);
    std::copy( NimArrValues.getPtr(),
               NimArrValues.getPtr() + length_constant,
               NimArrValues_AD.getPtr());
    setValues_AD_AD(NimArrValues_AD, NV.model_AD_constant_accessor);
  }

  // 2. Copy all wrtNodes values from model -> model_AD, AND
  // 3. Copy all wrtNodes values from model -> independentVars, AND
  // 4. [Deleted]
  int length_wrt = NV.model_wrt_accessor.getTotalLength();
  int length_independent = length_wrt;
  // std::cout<<"recording with length "<<length_independent<<std::endl;
  // std::cout<<" length_wrt = "<<length_wrt<<std::endl;
  vector< CppAD::AD<double> > independentVars(length_independent);
  if(length_wrt > 0) {
    NimArrValues.setSize(length_wrt);
    getValues(NimArrValues, NV.model_wrt_accessor);
    // 2
    NimArrValues_AD.setSize(length_wrt);
    std::copy( NimArrValues.getPtr(),
               NimArrValues.getPtr() + length_wrt,
               NimArrValues_AD.getPtr());
    setValues_AD_AD(NimArrValues_AD, NV.model_AD_wrt_accessor);
    // 3
    std::copy(  NimArrValues.getPtr(),
                NimArrValues.getPtr() + length_wrt,
                independentVars.begin() );
  }
  // 5a. Copy all extraInputNodes values from model -> model_AD (ditto, may be redundant)
  int length_extraInput = NV.model_extraInput_accessor.getTotalLength();
  if(length_extraInput > 0) {
    NimArrValues.setSize(length_extraInput);
    NimArrValues_AD.setSize(length_extraInput);
    getValues(NimArrValues, NV.model_extraInput_accessor);
    std::copy( NimArrValues.getPtr(),
               NimArrValues.getPtr() + length_extraInput,
               NimArrValues_AD.getPtr());
    setValues_AD_AD(NimArrValues_AD, NV.model_AD_extraInput_accessor);
  }
  // 5b. Copy all extraInputNodes into dynamicVars
  // std::cout<<"Don't forget to set the CppAD statics as needed"<<std::endl;
  vector< CppAD::AD<double> > dynamicVars;
  dynamicVars.resize(length_extraInput);
  if(length_extraInput > 0) {
    std::copy( NimArrValues_AD.getPtr(),
               NimArrValues_AD.getPtr() + length_extraInput,
               dynamicVars.begin() );
  }

  // 6. Start taping
  size_t abort_op_index = 0;    // per CppAD examples, these match CppAD default values
  bool   record_compare = true; // but must be provided explicitly to get to the dynamic parameter (4th) argument
  // std::cout<<"recording with "<<dynamicVars.size()<<std::endl;
  // std::cout<<"Before independent: tape handle address = "<< CppAD::AD<double>::get_handle_address_nimble() <<std::endl;
  CppAD::Independent(independentVars, abort_op_index, record_compare, dynamicVars);
  ADinfo.set_internal_tape(CppAD::AD<double>::get_tape_handle_nimble());
  //  std::cout<<"After independent: tape handle address = "<< CppAD::AD<double>::get_handle_address_nimble() <<std::endl;
  {
    set_CppAD_tape_info_for_model my_tape_info_RAII_(NV,
                                                     CppAD::AD<double>::get_tape_id_nimble(),
                                                     CppAD::AD<double>::get_tape_handle_nimble());
    set_CppAD_atomic_info_for_model(NV, CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage());

    // 7. [deleted]
    // 8. [deleted]

    // 9. Copy all wrtNodes AD objects from independentVars -> model_AD.
    if(length_extraInput > 0) {
      NimArrValues_AD.setSize(length_extraInput);
      std::copy(dynamicVars.begin(),
                dynamicVars.begin() + length_extraInput,
                NimArrValues_AD.getPtr());
      setValues_AD_AD(NimArrValues_AD, NV.model_AD_extraInput_accessor);
    }
    if(length_wrt > 0) {
      NimArrValues_AD.setSize(length_wrt);
      std::copy(independentVars.begin(),
                independentVars.begin() + length_wrt,
                NimArrValues_AD.getPtr());
      setValues_AD_AD(NimArrValues_AD, NV.model_AD_wrt_accessor);
    }
    // 10. call calculate.  This also sets up the extraOutput step
    nimbleCppADrecordingInfoClass recordingInfo(true, &ADinfo);
    CppAD::AD<double> logProb = calculate_ADproxyModel(NV,
                                                       includeExtraOutputs, // if true, model will be updated from tape.
                                                       recordingInfo);
    dependentVars[0] = logProb;
    // 13. Finish taping, AND
    // 14. Call tape->optimize()
    // make it locally to get the right globals during recording and playback
    // DO NOT USE THE CONSTRUCTOR VERSION BECAUSE IT ALWAYS DOES .Forward(0)
    // INSTEAD MAKE THE BLANK OBJECT AND USE .Dependent(...)
    // TRY USING CppAD's vector type
  } // These {} ensure that the destructor for the my_tape_info_RAII_ is called before Dependent, which is necessary in some cases (depending on libnimble.a vs libnimble.so)
  CppAD::ADFun<double>* ansTape = new CppAD::ADFun<double>;
  ADinfo.sum_dummyOutputs_to_dependentVars(dependentVars);
  ansTape->Dependent(independentVars, dependentVars);
  //  std::cout<<"about to call optimize"<<std::endl;
  //  std::cout<<"tape handle address = "<< CppAD::AD<double>::get_handle_address_nimble() <<std::endl;
#ifdef USE_CPPAD_OPTIMIZE_FOR_MODEL_TAPES
  ansTape->optimize(); //("no_compare_op") makes almost no difference;
#endif
  // std::cout<<"done with optimize"<<std::endl;
  return ansTape;
}

void nimbleFunctionCppADbase::getDerivs_calculate_internal(nimbleCppADinfoClass &ADinfo,
                                                           // CppAD::ADFun<double>* &tapePtr,
                                                           NodeVectorClassNew_derivs &nodes,
                                                           const NimArr<1, double> &derivOrders,
                                                           const NimArr<1, double> &wrtVector,
                                                           bool do_update,
                                                           bool reset,
                                                           nimSmartPtr<NIMBLE_ADCLASS> ansList) {
  // If use_meta_tape is true, then double-recording will be used.
  // This means that a first tape will be recorded, then a second tape will be recorded of obtaining 1st order derivatives from the first tape.
  // When derivatives are requested, 0th order will be obtained from the regular model (not AD tape),
  // 1st order will be 0th order of the second tape, and 2nd order will be 1st order of the second tape.
  // If use_meta_tape is false, then the (single, first) tape will be used "directly" for 0th, 1st or 2nd order.
  using std::cout;
  using std::endl;
  bool use_meta_tape = true;
 // cout<<"in getDerivs_calculate_internal "<<ADinfo.ADtape_empty()<<" "<<reset<<endl;
  // Record tape(s) if this is the first time or if reset is true.
  NimArr<1, double > outInds, inDir, outDir; // being length 0 will make them not used.
  NimArr<1, CppAD::AD<double> > meta_inDir, meta_outDir;
  if(ADinfo.ADtape_empty() || reset) {
    // Delete previous tape if it exists.
    if(!ADinfo.ADtape_empty())
      ADinfo.ADtape_reset();
    if(!use_meta_tape) {
      ADinfo.ADtape() = calculate_recordTape(nodes, true, ADinfo); // sets internal tape for atomic tracking
    } else {
      CppAD::ADFun< double > *firstTape;
      firstTape = calculate_recordTape(nodes, false, ADinfo); // sets internal tape for atomic tracking
      CppAD::ADFun< CppAD::AD<double>, double > innerTape;
      // Make original tape use CppAD::AD<double> instead of double
      set_CppAD_atomic_info_for_model(nodes, CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage());
      innerTape = firstTape->base2ad();
      int length_wrt = nodes.model_wrt_accessor.getTotalLength();
      int length_independent = length_wrt;
      int length_extraInput = nodes.model_extraInput_accessor.getTotalLength();
      vector< CppAD::AD<double> > dependentVars(length_wrt); // This will be the jacobian from the first tape, i.e. value of the second tape
      vector< CppAD::AD<double> > independentVars(length_independent);
      vector< CppAD::AD<double> > dynamicVars(length_extraInput);
      size_t abort_op_index = 0;
      bool   record_compare = true;
      NimArr<1, double > NimArrVars;
      NimArrVars.setSize(length_wrt);
      // Initialize values of independentVars, before recording.
      getValues(NimArrVars, nodes.model_wrt_accessor);
      for(int iii = 0; iii < length_wrt; ++iii)
        independentVars[iii] = NimArrVars[iii];
      // std::copy(NimArrVars.getPtr(),
      // NimArrVars.getPtr() + length_wrt,
      // independentVars.begin());
      //
      // Initialize values of dynamicVars, before recording.
      if(length_extraInput > 0) {
        NimArr<1, double> NimArr_dynamicVars;
        NimArr_dynamicVars.setSize(length_extraInput);
        getValues(NimArr_dynamicVars, nodes.model_extraInput_accessor);
        for(int iii = 0; iii < length_extraInput; ++iii)
          dynamicVars[iii] = NimArr_dynamicVars[iii];
        // std::copy( NimArr_dynamicVars.getPtr(),
        // NimArr_dynamicVars.getPtr() + length_extraInput,
        // dynamicVars.begin() );
      }
      nimSmartPtr<NIMBLE_ADCLASS_META> ansList_meta = new NIMBLE_ADCLASS_META;
      // start recording new (second) tape
      CppAD::Independent(independentVars, abort_op_index, record_compare, dynamicVars);
      ADinfo.set_internal_tape(CppAD::AD<double>::get_tape_handle_nimble());
      // Trick CppAD statics to work across nimble compilation units
      {
        set_CppAD_tape_info_for_model my_tape_info_RAII_(nodes,
                                                         CppAD::AD<double>::get_tape_id_nimble(),
                                                         CppAD::AD<double>::get_tape_handle_nimble());
        set_CppAD_atomic_info_for_model(nodes, CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage());
        // Set up inputs to first tape (recorded in second tape)
        ADinfo.independentVars_meta.resize(length_wrt);
        ADinfo.dynamicVars_meta.resize(length_extraInput);
        if(length_extraInput > 0) {
          for(int iii = 0; iii < length_extraInput; ++iii)
            ADinfo.dynamicVars_meta[iii] = dynamicVars[iii]; // std::copy does not seem to work for CppAD recording
          // std::copy(dynamicVars.begin(),
          //  dynamicVars.end(),
          //  ADinfo.dynamicVars_meta.begin());
        }
        for(int iii = 0; iii < length_wrt; ++iii)
          ADinfo.independentVars_meta[iii] = independentVars[iii];
        // std::copy(independentVars.begin(),
        // independentVars.end(),
        // ADinfo.independentVars_meta.begin());
      NimArr<1, double> derivOrders_meta;
      derivOrders_meta.setSize(1);
      derivOrders_meta[0] = 1;
      // std::cout<<ADinfo.dynamicVars_meta.size()<<std::endl;
      innerTape.new_dynamic(ADinfo.dynamicVars_meta);
      getDerivs_internal< CppAD::AD<double>,
                          CppAD::ADFun< CppAD::AD<double>, double >,
                          NIMBLE_ADCLASS_META>(ADinfo.independentVars_meta,
                                               &innerTape,
                                               derivOrders_meta,
                                               wrtVector,
                                                outInds,
                                                meta_inDir,
                                                meta_outDir,
                                               ansList_meta);
      for(int iii = 0; iii < length_wrt; ++iii)
        dependentVars[iii] = ansList_meta->jacobian[iii];
      ADinfo.ADtape() = new CppAD::ADFun<double>;
      } // These {} ensure the RAII object's destructor is called before Dependent, which is important on OS's (linux) with libnimble.so instead of libnimble.a
      ADinfo.ADtape()->Dependent(dependentVars);
#ifdef USE_CPPAD_OPTIMIZE_FOR_MODEL_TAPES
      ADinfo.ADtape()->optimize();
#endif
      delete firstTape;
    }
  }

  // Recording, if needed, is done.
  // From here on is use of the tape(s).  This may be used much more often than recording section above.

  //std::cout<<"getDerivs_calculate_internal A"<<std::endl;
  // Copy values from the model into the independentVars
  int length_wrt = nodes.model_wrt_accessor.getTotalLength();
  int length_independent = length_wrt;
  ADinfo.independentVars.resize(length_independent);

  NimArr<1, double > NimArrVars;
  NimArrVars.setSize(length_wrt);
  getValues(NimArrVars, nodes.model_wrt_accessor);

  std::copy(NimArrVars.getPtr(),
            NimArrVars.getPtr() + length_wrt,
            ADinfo.independentVars.begin());
  /* set dynamic */
  // Copy extraInput (CppAD "dynamic") values from the model into the dynamicVars
  // *and* set them in the tape.
  size_t length_extraNodes_accessor = nodes.model_extraInput_accessor.getTotalLength();
  if(length_extraNodes_accessor > 0) {
    NimArr<1, double> NimArr_dynamicVars;
    NimArr_dynamicVars.setSize(length_extraNodes_accessor);
    getValues(NimArr_dynamicVars, nodes.model_extraInput_accessor);
    std::vector<double> dynamicVars(length_extraNodes_accessor);
    std::copy( NimArr_dynamicVars.getPtr(),
               NimArr_dynamicVars.getPtr() + length_extraNodes_accessor,
               dynamicVars.begin() );
    //std::cout<<"Setting new_dynamic to"<<std::endl;
    //for(int ijk = 0; ijk < length_extraNodes_accessor; ijk++)
    //  std::cout<<dynamicVars[ijk]<<" ";
    //std::cout<<std::endl;
    ADinfo.ADtape()->new_dynamic(dynamicVars);
  }
//  std::cout<<"getDerivs_calculate_internal B"<<std::endl;

  if(use_meta_tape) {
    // manage orders and use regular calculate for value
    // and tape for jacobian or hessian
  //std::cout<<"getDerivs_calculate_internal C"<<std::endl;
    int maxOrder;
    bool ordersFound[3];
    setOrdersFound(derivOrders, ordersFound, maxOrder);
 // std::cout<<"getDerivs_calculate_internal C2"<<std::endl;
    if(ordersFound[0]) {
      ansList->value.setSize(1, false, false);
      ansList->value[0] = calculate(nodes);
    }
  //    std::cout<<"getDerivs_calculate_internal C3"<<std::endl;

    NimArr<1, double> derivOrders_nested;
    int higherOrders = 0;
    if(ordersFound[1]) ++higherOrders;
    if(ordersFound[2]) ++higherOrders;
    if(higherOrders) {
      derivOrders_nested.setSize(higherOrders, false, false);
      higherOrders = 0;
      if(ordersFound[1]) derivOrders_nested[higherOrders++] = 0; // If Jacobian was requested, get value of meta tape
      if(ordersFound[2]) derivOrders_nested[higherOrders] = 1; // If Hessian was requested, get Jacobian of meta tape
      nimSmartPtr<NIMBLE_ADCLASS> ansList_nested = new NIMBLE_ADCLASS;
      //std::cout<<"about to call getDerivs_internal"<<std::endl;
//  std::cout<<"getDerivs_calculate_internal D"<<std::endl;
      getDerivs_internal<double,
                         CppAD::ADFun<double>,
                         NIMBLE_ADCLASS>(ADinfo.independentVars,
                                         ADinfo.ADtape(),
                                         derivOrders_nested,
                                         wrtVector, // NOTE: This will not behave fully correctly in non-default use without further thought.
                                          outInds,
                                          inDir,
                                          outDir,
                                         ansList_nested);
//  std::cout<<"getDerivs_calculate_internal E"<<std::endl;
      if(ordersFound[1]) {
        ansList->jacobian.setSize(1, length_wrt, false, false);
        for(int ii = 0; ii < length_wrt; ++ii) //We could rewrite this with better pointer magic
          ansList->jacobian[ ii ] = ansList_nested->value[ ii ];
      }
      if(ordersFound[2]) {
        ansList->hessian.setSize(length_wrt, length_wrt, 1, false, false);
        for(int ii = 0; ii < length_wrt; ++ii)
          for(int jj = 0; jj < length_wrt; ++jj)
            ansList->hessian[jj + ii*length_wrt ] = ansList_nested->jacobian[jj + ii*length_wrt]; //orientation shouldn't matter due to symmetry
      }
    }
  } else {
    /* run tape */
    //std::cout<<"running tape"<<std::endl;
    getDerivs_internal<double,
                       CppAD::ADFun<double>,
                       NIMBLE_ADCLASS>(ADinfo.independentVars,
                                       ADinfo.ADtape(),
                                       derivOrders,
                                       wrtVector,
                                        outInds,
                                       inDir,
                                       outDir,
                                       ansList);
  }
}
