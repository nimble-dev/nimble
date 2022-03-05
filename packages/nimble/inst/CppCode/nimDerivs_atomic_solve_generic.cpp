
/*
  This code is included in other files with several macros defined differently for each time it is included.
  This is not an aesthetically appealing implementation.  
  We tried a templated struct with static methods to dispatch different function calls in different cases.
However in the particular cases of Eigen code involved, this did not work well.  It gave crashes, unless 
the result was materialized locally, which defeats Eigen's efficiency.  Eigen documentation suggests such 
problems can occur if internal temporaries are created and then lost when the op is assigned and held on to 
as an unevaluated op.
 */

ATOMIC_SOLVE_CLASS::ATOMIC_SOLVE_CLASS(const std::string& name) :
  atomic_solve_base_class(), CppAD::atomic_three<double>(name)
{};
  
bool ATOMIC_SOLVE_CLASS::forward(
				     const CppAD::vector<double>&               parameter_x  ,
				     const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
				     size_t                              need_y       ,
				     size_t                              order_low    ,
				     size_t                              order_up     ,
				     const CppAD::vector<double>&               taylor_x     ,
				     CppAD::vector<double>&                     taylor_y     ) {
  //forward mode
  int nrow = order_up + 1;
#ifdef VERBOSE_ATOMIC_FS
  printf("In forwardsolve forward\n");
  std::cout<<"order_low = "<<order_low<<" order_up = "<<order_up<<" nrow = "<<nrow<<std::endl;
#endif
  int n = taylor_x.size()/nrow;
  int m = taylor_y.size() / nrow;
  int n1sq, n1, n2;
  if((!Aconstant()) && (!Bconstant())) n1sq = n-m;
  if(Aconstant()) n1sq = get_X_stored().size();
  if(Bconstant()) n1sq = n;
  if(Aconstant() && Bconstant())
    std::cout<<"atomic_forwardsolve is being used with both A and B constant.  This should not happen."<<std::endl;
  n1 = sqrt( static_cast<double>(n1sq) );
  n2 = m/n1;

  // Note on cache vs. constants:
  // The cache is used to save a value of taylor_y.
  // This is done because CppAD can optionally destroy it in its tape optimizations,
  // but then if later orders expect to use taylor_y, it is not there.  I have
  // talked about this with Brad Bell.  It is what it is for now.
  //
  // Constants are used for values that won't change, which are stored in X_stored (X as in taylor_x, the input).
  
  const double *Xptr;
  int row_mult;
  
  if(Aconstant()) 
    {Xptr = get_X_stored_ptr(); row_mult = 1;} else
    {Xptr = &taylor_x[0];       row_mult = nrow;}
  new (&Amap) EigenConstMap(Xptr, n1, n1, EigStrDyn(row_mult*n1, row_mult));
  //  EigenConstMap Amap(&taylor_x[0], n1, n1, EigStrDyn(nrow*n1, nrow) );
  
  if(order_low <= 0 & order_up >= 0) { // value
    //printf("In forward 0\n");
    // We could compile different cases depending on need for strides or not.
    //    EigenMap Ymap(&taylor_y[0], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    new (&Ymap) EigenMap(&taylor_y[0], n1, n2, EigStrDyn(nrow*n1, nrow ) );

    if(Bconstant())
      {Xptr = get_X_stored_ptr(); row_mult = 1; } else 
      { row_mult = nrow;
	if(Aconstant())
	  {Xptr = &taylor_x[0]; } else
	  {Xptr = &taylor_x[0 + n1sq*nrow];}
      }
    new (&Bmap) EigenConstMap(Xptr, n1, n2, EigStrDyn(row_mult*n1, row_mult ) );
    //    EigenConstMap Bmap(&taylor_x[0 + n1sq*nrow], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    // There is some indication in nimble EigenTypeDefs.h that we need to force a copy of B
    //    Ymap = Amap.template triangularView<Eigen::Lower>().solve(Bmap);
    Ymap = ND_SOLVE(Amap, Bmap).eval();
    double_cache.set_cache( 0, 0, order_up, taylor_x, taylor_y );
  }
  if(order_low <= 1 & order_up >= 1) {
    //printf("In forward >1\n");
    //      solve(A, dB - dA * Y)
    double_cache.check_and_set_cache(this,
				     parameter_x,
				     type_x,
				     0,
				     order_up,
				     taylor_x,
				     taylor_y.size());
    int cache_nrow = double_cache.nrow();
    // EigenMap Ymap(double_cache.taylor_y_ptr(), n1, n2, EigStrDyn(cache_nrow*n1, cache_nrow ) );
    new (&Ymap) EigenMap(double_cache.taylor_y_ptr(), n1, n2, EigStrDyn(cache_nrow*n1, cache_nrow ) );

    if(!Aconstant()) {
      new (&dA_map) EigenConstMap(&taylor_x[1], n1, n1, EigStrDyn(nrow*n1, nrow));
      // EigenConstMap dA_map(&taylor_x[1], n1, n1, EigStrDyn(nrow*n1, nrow) );
    }
    if(!Bconstant()) {
      if(Aconstant())
	{Xptr = &taylor_x[1];} else
	{Xptr = &taylor_x[1 + n1sq*nrow];}
      new (&dB_map) EigenConstMap(Xptr, n1, n2, EigStrDyn(nrow*n1, nrow));
      // EigenConstMap dB_map(&taylor_x[1 + n1sq*nrow], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    }

    new (&dY_map) EigenMap(&taylor_y[1], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    //EigenMap dY_map(&taylor_y[1], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    // if((!Aconstant()) && (!Bconstant())) {
    //   dY_map = ND_SOLVE(Amap, dB_map - ND_MAIN_TRI(dA_map) * Ymap).eval();
    //   //      dY_map = Amap.template triangularView<Eigen::Lower>().solve(dB_map - dA_map.template triangularView<Eigen::Lower>() * Ymap).eval();// This .eval() is necessary and I don't understand why.  Normally that would be for aliasing, but there should be no over-lapping points here.  There are inter-woven maps, and that's weird but should work.
    // } else if(Aconstant()) {
    //   dY_map = ND_SOLVE(Amap, dB_map).eval();
    //   //      dY_map = Amap.template triangularView<Eigen::Lower>().solve(dB_map).eval();
    // } else if(Bconstant()) {
    //   // std::cout<<"using Bconstant case"<<std::endl;
    //   dY_map = -ND_SOLVE(Amap, ND_MAIN_TRI(dA_map) * Ymap).eval();
    //   //      dY_map = -Amap.template triangularView<Eigen::Lower>().solve(dA_map.template triangularView<Eigen::Lower>() * Ymap).eval();
    // } else {
    //   // This case should never happen, and should have been warned above.
    // }

    //
    if(Avariable() && Bvariable()) {
      dY_map = ND_SOLVE(Amap, dB_map - ND_MAIN_TRI(dA_map) * Ymap).eval();
      //      dY_map = Amap.template triangularView<Eigen::Lower>().solve(dB_map - dA_map.template triangularView<Eigen::Lower>() * Ymap).eval();// This .eval() is necessary and I don't understand why.  Normally that would be for aliasing, but there should be no over-lapping points here.  There are inter-woven maps, and that's weird but should work.
    } else if(!Avariable()) {
      dY_map = ND_SOLVE(Amap, dB_map).eval();
      //      dY_map = Amap.template triangularView<Eigen::Lower>().solve(dB_map).eval();
    } else if(!Bvariable()) {
      dY_map = -ND_SOLVE(Amap, ND_MAIN_TRI(dA_map) * Ymap).eval();
      //      dY_map = -Amap.template triangularView<Eigen::Lower>().solve(dA_map.template triangularView<Eigen::Lower>() * Ymap).eval();
    } else {
      // Technically we could have A and B both dynamic (not variable, not constant).
      // I'm not sure order one will ever be called in this case.
    }
    //

    
    double_cache.set_cache( 1, 1, order_up, taylor_x, taylor_y );
  }
  return true;
}

bool ATOMIC_SOLVE_CLASS::forward(
				     const CppAD::vector<CppAD::AD<double> >&               parameter_x  ,
				     const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
				     size_t                              need_y       ,
				     size_t                              order_low    ,
				     size_t                              order_up     ,
				     const CppAD::vector<CppAD::AD<double> >&               taylor_x     ,
				     CppAD::vector<CppAD::AD<double> >&                     taylor_y     ) {
  // meta-forward mode
  int nrow = order_up + 1;
#ifdef VERBOSE_ATOMIC_FS
  printf("In forwardsolve meta-forward\n");
  std::cout<<"order_low = "<<order_low<<" order_up = "<<order_up<<" nrow = "<<nrow<<std::endl;
#endif
  int n = taylor_x.size()/nrow;
  int m = taylor_y.size() / nrow;
  int n1sq, n1, n2;
  if((!Aconstant()) && (!Bconstant())) n1sq = n-m;
  if(Aconstant()) n1sq = get_X_stored().size();
  if(Bconstant()) n1sq = n;
  if(Aconstant() && Bconstant())
    std::cout<<"atomic_forwardsolve is being used with both A and B constant.  This should not happen."<<std::endl;
  n1 = sqrt( static_cast<double>(n1sq) );
  n2 = m/n1;

  const CppAD::AD<double> *Xptr;
  int row_mult;
  
  if(Aconstant()) 
    {fill_X_AD_stored(); Xptr = get_X_AD_stored_ptr(); row_mult = 1;} else
    {                    Xptr = &taylor_x[0];       row_mult = nrow;}
  new (&mAmap) metaEigenConstMap(Xptr, n1, n1, EigStrDyn(row_mult*n1, row_mult));
  // metaEigenConstMap mAmap(&taylor_x[0], n1, n1, EigStrDyn(nrow*n1, nrow) );
  if(order_low <= 0 & order_up >= 0) { // value
      // We could compile different cases depending on need for strides or not.
    new (&mYmap) metaEigenMap(&taylor_y[0], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    //    metaEigenMap mYmap(&taylor_y[0], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    
    if(Bconstant())
      {fill_X_AD_stored(); Xptr = get_X_AD_stored_ptr(); row_mult = 1; } else 
      { row_mult = nrow;
	if(Aconstant())
	  {Xptr = &taylor_x[0]; } else
	  {Xptr = &taylor_x[0 + n1sq*nrow];}
      }
    new (&mBmap) metaEigenConstMap(Xptr, n1, n2, EigStrDyn(row_mult*n1, row_mult ) );
    // metaEigenConstMap mBmap(&taylor_x[0 + n1sq*nrow], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    // There is some indication in nimble EigenTypeDefs.h that we need to force a copy of B
    mYmap = ND_META_SOLVE(mAmap, mBmap);
    // Ymap = Amap.template triangularView<Eigen::Upper>().solve(Bmap);
    
    CppADdouble_cache.set_cache( 0, 0, order_up, taylor_x, taylor_y );
  }
  if(order_low <= 1 & order_up >= 1) {
    //printf("In forward >1\n");
    //      solve(A, dB - dA * Y)
    CppADdouble_cache.check_and_set_cache(this,
					  parameter_x,
					  type_x,
					  0,
					  order_up,
					  taylor_x,
					  taylor_y.size());
    int cache_nrow = CppADdouble_cache.nrow();
    new (&mYmap) metaEigenMap(CppADdouble_cache.taylor_y_ptr(), n1, n2, EigStrDyn(cache_nrow*n1, cache_nrow ) );
    // metaEigenMap Ymap(CppADdouble_cache.taylor_y_ptr(), n1, n2, EigStrDyn(cache_nrow*n1, cache_nrow ) );
    if(!Aconstant()) {
      new (&mdA_map) metaEigenConstMap(&taylor_x[1], n1, n1, EigStrDyn(nrow*n1, nrow));
      // metaEigenConstMap dA_map(&taylor_x[1], n1, n1, EigStrDyn(nrow*n1, nrow) );
    }
    
    if(!Bconstant()) {
      if(Aconstant())
	{Xptr = &taylor_x[1];} else
	{Xptr = &taylor_x[1 + n1sq*nrow];}
      new (&mdB_map) metaEigenConstMap(Xptr, n1, n2, EigStrDyn(nrow*n1, nrow));
      // metaEigenConstMap dB_map(&taylor_x[1 + n1sq*nrow], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    }

    new (&mdY_map) metaEigenMap(&taylor_y[1], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    //  metaEigenMap dY_map(&taylor_y[1], n1, n2, EigStrDyn(nrow*n1, nrow ) );

    if(Avariable() && Bvariable()) {
      mdY_map = ND_META_SOLVE(mAmap, mdB_map - nimDerivs_matmult( ND_MAIN_TRI(mdA_map), mYmap) );
      //      mdY_map = ND_META_SOLVE(mAmap, mdB_map - nimDerivs_matmult( ND_MAIN_TRI(mdA_map), mYmap));
      //mdY_map = nimDerivs_EIGEN_FS(mAmap, mdB_map - nimDerivs_matmult(mdA_map.template triangularView<Eigen::Lower>(), mYmap));

    //    dY_map = Amap.template triangularView<Eigen::Upper>().solve(dB_map - dA_map * Ymap).eval();// This .eval() is necessary and I don't understand why.  Normally that would be for aliasing, but there should be no over-lapping points here.  There are inter-woven maps, and that's weird but should work.
    } else if(!Avariable()) {
      mdY_map = ND_META_SOLVE(mAmap, mdB_map);
      // mdY_map = nimDerivs_EIGEN_FS(mAmap, mdB_map);
    } else if(!Bvariable()) {
      mdY_map = ND_META_SOLVE(mAmap, -nimDerivs_matmult(ND_MAIN_TRI(mdA_map), mYmap));
      // mdY_map = nimDerivs_EIGEN_FS(mAmap, -nimDerivs_matmult(mdA_map.template triangularView<Eigen::Lower>(), mYmap));
    } else {
      // Technically we could have A and B both dynamic (not variable, not constant).
      // I'm not sure order one will ever be called in this case.
    }
    CppADdouble_cache.set_cache( 1, 1, order_up, taylor_x, taylor_y );
  }
  return true;
}

bool ATOMIC_SOLVE_CLASS::reverse(
				     const CppAD::vector<double>&               parameter_x ,
				     const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				     size_t                              order_up    ,
				     const CppAD::vector<double>&               taylor_x    ,
				     const CppAD::vector<double>&               taylor_y    ,
				     CppAD::vector<double>&                     partial_x   ,
				     const CppAD::vector<double>&               partial_y   )
{
  //reverse mode
  int nrow = order_up + 1;
#ifdef VERBOSE_ATOMIC_FS
  printf("In forwardsolve reverse\n");
  std::cout<<" order_up = "<<order_up<<" nrow = "<<nrow<<std::endl;
#endif
  int n = taylor_x.size()/nrow;
  int m = taylor_y.size() / nrow;
  int n1sq, n1, n2;
  if((!Aconstant()) && (!Bconstant())) n1sq = n-m;
  if(Aconstant()) n1sq = get_X_stored().size();
  if(Bconstant()) n1sq = n;
  if(Aconstant() && Bconstant())
    std::cout<<"atomic_forwardsolve is being used with both A and B constant.  This should not happen."<<std::endl;
  n1 = sqrt( static_cast<double>(n1sq) );
  n2 = m/n1;

  double_cache.check_and_set_cache(this,
				   parameter_x,
				   type_x,
				   order_up >= 1 ? 1 : 0, // only use cached values up to order 1
				   order_up,
				   taylor_x,
				   taylor_y.size());
  int cache_nrow = double_cache.nrow();
  new (&YmapC) EigenConstMap(double_cache.taylor_y_ptr(), n1, n2, EigStrDyn(cache_nrow*n1, cache_nrow ) );
  //EigenConstMap YmapC(double_cache.taylor_y_ptr(), n1, n2, EigStrDyn(cache_nrow*n1, cache_nrow ) );
  //  EigenConstMap Ymap(&taylor_y[0], n1, n2, EigStrDyn(nrow*n1, nrow ) );

  const double *Xptr;
  double *pXptr;
  int row_mult;
  
  if(Aconstant()) 
    {Xptr = get_X_stored_ptr(); row_mult = 1;} else
    {Xptr = &taylor_x[0];       row_mult = nrow;
      new (&Aadjoint_map) EigenMap(&partial_x[0], n1, n1, EigStrDyn(nrow*n1, nrow));
      // EigenMap Aadjoint_map(&partial_x[0], n1, n1, EigStrDyn(nrow*n1, nrow) );
    }
  new (&Amap) EigenConstMap(Xptr, n1, n1, EigStrDyn(row_mult*n1, row_mult));
  //  EigenConstMap Amap(&taylor_x[0], n1, n1, EigStrDyn(nrow*n1, nrow) );
  
  if(Bconstant())
    { pXptr = get_Badjoint_memory(n1*n2); row_mult = 1;} else 
    { row_mult = nrow;
      if(Aconstant())
	{pXptr = &partial_x[0];            } else
	{pXptr = &partial_x[0 + n1sq*nrow];}
    }
  new (&Badjoint_map) EigenMap(pXptr, n1, n2, EigStrDyn(row_mult*n1, row_mult ) );
  // Note Badjoint_map points to internal class memory if B is constant because
  // Badjoint calculations are used in Aadjoint.  It is simpler to code them in this way below.
  // I don't think there would be much benefit to a separate Eigen expression in the constant-B
  // case because matrix multiplications are always materialized inside Eigen anyway, IIUC.
  //
  //EigenMap Badjoint_map(&partial_x[0 + n1sq*nrow], n1, n2, EigStrDyn(nrow*n1, nrow ) );

  // Bmap is not needed.
  
  if(order_up >= 0) {
    /* EigenTemplateTypes<double>::typeEigenConstMapStrd Bmap(&taylor_x[0 + n1sq*nrow], */
    /* 							 n1, n2, EigStrDyn(nrow*n1, nrow ) ); */
    new (&Yadjoint_map) EigenConstMap(&partial_y[0], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    //EigenConstMap Yadjoint_map(&partial_y[0], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    Badjoint_map = ND_TRANS_SOLVE(Amap, Yadjoint_map).eval();
    //    Badjoint_map = Amap.transpose().template triangularView<Eigen::Upper>().solve(Yadjoint_map).eval();
    if(order_up == 0) { // otherwise this gets included below
      if(Avariable())
	Aadjoint_map = ND_MAIN_TRI(-Badjoint_map * YmapC.transpose());
    }
  }
  if(order_up >= 1) {
    new (&Ydot_adjoint_map) EigenConstMap(&partial_y[1], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    // EigenConstMap Ydot_adjoint_map(&partial_y[1], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    new (&Ydot_map) EigenConstMap(double_cache.taylor_y_ptr() + 1, n1, n2, EigStrDyn(cache_nrow*n1, cache_nrow ) );
    // EigenConstMap Ydot_map(double_cache.taylor_y_ptr() + 1, n1, n2, EigStrDyn(cache_nrow*n1, cache_nrow ) );
    // EigenConstMap Ydot_map(&taylor_y[1], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    
    // Again we calculate a B quantity, Bdot_adjoint,
    // even if B is constant, because it goes into Adot_adjoint later.
    if(Bconstant()) 
      {pXptr = get_Bdot_adjoint_memory(n1*n2); row_mult = 1;} else
      {row_mult = nrow;
	if(Aconstant())
	  {pXptr = &partial_x[1];    } else
	  {pXptr = &partial_x[1 + n1sq*nrow];}
      }
    new (&Bdot_adjoint_map) EigenMap(pXptr, n1, n2, EigStrDyn(row_mult*n1, row_mult ) );
    Bdot_adjoint_map = ND_TRANS_SOLVE(Amap, Ydot_adjoint_map).eval();
    // Bdot_adjoint_map = Amap.transpose().template triangularView<Eigen::Upper>().solve(Ydot_adjoint_map).eval(); // This eval is necessary.  I am not clear when evals are necessary after solves or not.
      // EigenMap Bdot_adjoint_map(&partial_x[1 + n1sq*nrow], n1, n2, EigStrDyn(nrow*n1, nrow ) );
      /* Bdot_adjoint = A^-T Ydot_adjoint */

    if(Avariable()) {
      // Comes after Bdot_adjoint_map because that is used in Adjoint_map
      new (&Adot_map) EigenConstMap(&taylor_x[1], n1, n1, EigStrDyn(nrow*n1, nrow) );
      // EigenConstMap Adot_map(&taylor_x[1], n1, n1, EigStrDyn(nrow*n1, nrow) );
      new (&Adot_adjoint_map) EigenMap(&partial_x[1], n1, n1, EigStrDyn(nrow*n1, nrow) );
      // EigenMap Adot_adjoint_map(&partial_x[1], n1, n1, EigStrDyn(nrow*n1, nrow) );
 
      /* Adot_adjoint = -Bdot_adjoint Y^T */
      Adot_adjoint_map = ND_MAIN_TRI(-Bdot_adjoint_map * YmapC.transpose());
      // Adot_adjoint_map = (-Bdot_adjoint_map * YmapC.transpose()).template triangularView<Eigen::Lower>();

      Eigen::MatrixXd Adot_stuff = ND_TRANS_TRI(Adot_map); // doesn't work inside solve
      // Eigen::MatrixXd Adot_stuff = Adot_map.transpose().template triangularView<Eigen::Upper>(); // doesn't work inside solve
      /* Badjoint = A^-T * Yadjoint - (A^-1 Adot A^-1)^T  Ydot_adjoint */
      // Update Badjoint even if B is constant, pointing to internal memory, for purposes of Aadjoint
      // If A is constant, Adot is all 0, so this term is not needed.
      Badjoint_map -= ND_SOLVE(Amap,
			       ND_TRANS_SOLVE(Amap, Adot_stuff).eval().transpose()).eval().transpose() * Ydot_adjoint_map;
      // Badjoint_map -= Amap.template triangularView<Eigen::Lower>().solve( Amap.transpose().template triangularView<Eigen::Upper>().solve(Adot_stuff).eval().transpose() ).eval().transpose() * Ydot_adjoint_map;

      // The Badjoint and Bdot_adjoint in these terms is algebraic simplification.
      // The adjoints are non-zero if A is not constant.
      Aadjoint_map = ND_MAIN_TRI(-Badjoint_map * YmapC.transpose());
      // Aadjoint_map = (-Badjoint_map * YmapC.transpose()).template triangularView<Eigen::Lower>();
      Aadjoint_map -= ND_MAIN_TRI(Bdot_adjoint_map * Ydot_map.transpose());
      // Aadjoint_map -= (Bdot_adjoint_map * Ydot_map.transpose()).template triangularView<Eigen::Lower>(); 
	/* Aadjoint = -Badjoint * Y^T  - (A^-T Ydot_adjoint Ydot^T)= -Badjoint * Y^T  - ( Bdot_adjoint Ydot^T) , including both Badjoint terms*/
    }
    /* Note: A^-1 Z A^-1 = solve(A, solve(A^T, Z^T)^T) */
  }

  if(order_up >= 2) {
    printf("Unsupported reverse order requested\n");
    return false;
  }
  return true;
}

bool ATOMIC_SOLVE_CLASS::reverse(
				     const CppAD::vector<CppAD::AD<double> >&               parameter_x ,
				     const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				     size_t                              order_up    ,
				     const CppAD::vector<CppAD::AD<double> >&               taylor_x    ,
				     const CppAD::vector<CppAD::AD<double> >&               taylor_y    ,
				     CppAD::vector<CppAD::AD<double> >&                     partial_x   ,
				     const CppAD::vector<CppAD::AD<double> >&               partial_y   )
{
  //meta-reverse mode
  int nrow = order_up + 1;
#ifdef VERBOSE_ATOMIC_FS
  printf("In forwardsolve meta-reverse\n");
  std::cout<<" order_up = "<<order_up<<" nrow = "<<nrow<<std::endl;
#endif
  int n = taylor_x.size()/nrow;
  int m = taylor_y.size() / nrow;
  int n1sq, n1, n2;
  if((!Aconstant()) && (!Bconstant())) n1sq = n-m;
  if(Aconstant()) n1sq = get_X_stored().size();
  if(Bconstant()) n1sq = n;
  if(Aconstant() && Bconstant())
    std::cout<<"atomic_forwardsolve is being used with both A and B constant.  This should not happen."<<std::endl;
  n1 = sqrt( static_cast<double>(n1sq) );
  n2 = m/n1;

  CppAD::AD<double> *pXptr;
  const CppAD::AD<double> *Xptr;
  int row_mult;
  
  CppADdouble_cache.check_and_set_cache(this,
					parameter_x,
					type_x,
					order_up >= 1 ? 1 : 0, // only use cached values up to order 1
					order_up,
					taylor_x,
					taylor_y.size());
  int cache_nrow = CppADdouble_cache.nrow();
  new (&mYmapC) metaEigenConstMap(CppADdouble_cache.taylor_y_ptr(), n1, n2, EigStrDyn(cache_nrow*n1, cache_nrow ) );
  // metaEigenConstMap Ymap(CppADdouble_cache.taylor_y_ptr(), n1, n2, EigStrDyn(cache_nrow*n1, cache_nrow ) );
  // metaEigenConstMap Ymap(&taylor_y[0], n1, n2, EigStrDyn(nrow*n1, nrow ) );

  if(Aconstant()) 
    {fill_X_AD_stored(); Xptr = get_X_AD_stored_ptr(); row_mult = 1;} else
    {                    Xptr = &taylor_x[0];          row_mult = nrow;
      new (&mAadjoint_map) metaEigenMap(&partial_x[0], n1, n1, EigStrDyn(nrow*n1, nrow));
      // metaEigenMap Aadjoint_map(&partial_x[0], n1, n1, EigStrDyn(nrow*n1, nrow) );
    }
  new (&mAmap) metaEigenConstMap(Xptr, n1, n1, EigStrDyn(row_mult*n1, row_mult));
  // metaEigenConstMap Amap(&taylor_x[0], n1, n1, EigStrDyn(nrow*n1, nrow) );
  
  if(Bconstant())
    { pXptr = get_Badjoint_AD_memory(n1*n2); row_mult = 1; } else 
    { row_mult = nrow;
      if(Aconstant())
	{pXptr = &partial_x[0];            } else
	{pXptr = &partial_x[0 + n1sq*nrow];}
    }
  new (&mBadjoint_map) metaEigenMap(pXptr, n1, n2, EigStrDyn(row_mult*n1, row_mult ) );
  // metaEigenMap Badjoint_map(&partial_x[0 + n1sq*nrow], n1, n2, EigStrDyn(nrow*n1, nrow ) );
  if(order_up >= 0) {
    /* EigenTemplateTypes<double>::typeEigenConstMapStrd Bmap(&taylor_x[0 + n1sq*nrow], */
    /* 							 n1, n2, EigStrDyn(nrow*n1, nrow ) ); */
    new (&mYadjoint_map) metaEigenConstMap(&partial_y[0], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    // metaEigenConstMap Yadjoint_map(&partial_y[0], n1, n2, EigStrDyn(nrow*n1, nrow ) );

    mBadjoint_map = ND_META_TRANS_SOLVE( ND_TRANS_TRI(mAmap), mYadjoint_map);
    //mBadjoint_map = nimDerivs_EIGEN_BS(mAmap.transpose().template triangularView<Eigen::Upper>(), mYadjoint_map);
    // Badjoint_map = Amap.transpose().template triangularView<Eigen::Lower>().solve(Yadjoint_map);
    if(order_up == 0) {
      if(Avariable())
	 mAadjoint_map = ND_MAIN_TRI(nimDerivs_matmult(-mBadjoint_map, mYmapC.transpose()));
       // mAadjoint_map = nimDerivs_matmult(-mBadjoint_map, mYmapC.transpose()).template triangularView<Eigen::Lower>();//-Badjoint_map * Ymap.transpose(); // otherwise this gets included below
    }
  }
  if(order_up >= 1) {
    new (&mYdot_adjoint_map) metaEigenConstMap(&partial_y[1], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    // metaEigenConstMap Ydot_adjoint_map(&partial_y[1], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    new (&mYdot_map) metaEigenConstMap(CppADdouble_cache.taylor_y_ptr() + 1, n1, n2, EigStrDyn(cache_nrow*n1, cache_nrow ) );
   // metaEigenConstMap Ydot_map(CppADdouble_cache.taylor_y_ptr() + 1, n1, n2, EigStrDyn(cache_nrow*n1, cache_nrow ) );
    // metaEigenConstMap Ydot_map(&taylor_y[1], n1, n2, EigStrDyn(nrow*n1, nrow ) );


    // Again we calculate a B quantity, Bdot_adjoint,
    // even if B is constant, because it goes into Adot_adjoint later.
    if(Bconstant()) 
      {pXptr = get_Bdot_adjoint_AD_memory(n1*n2); row_mult = 1;} else 
      {row_mult = nrow;
	if(Aconstant())
	  {pXptr = &partial_x[1];    } else
	  {pXptr = &partial_x[1 + n1sq*nrow];}
      }
    new (&mBdot_adjoint_map) metaEigenMap(pXptr, n1, n2, EigStrDyn(row_mult*n1, row_mult ) );
    mBdot_adjoint_map = ND_META_TRANS_SOLVE(ND_TRANS_TRI(mAmap), mYdot_adjoint_map);
    // mBdot_adjoint_map = nimDerivs_EIGEN_BS(mAmap.transpose().template triangularView<Eigen::Upper>(), mYdot_adjoint_map);
    
    if(Avariable()) {
      new (&mAdot_map) metaEigenConstMap(&taylor_x[1], n1, n1, EigStrDyn(nrow*n1, nrow) );
      // metaEigenConstMap Adot_map(&taylor_x[1], n1, n1, EigStrDyn(nrow*n1, nrow) );
      new (&mAdot_adjoint_map) metaEigenMap(&partial_x[1], n1, n1, EigStrDyn(nrow*n1, nrow) );
      // metaEigenMap Adot_adjoint_map(&partial_x[1], n1, n1, EigStrDyn(nrow*n1, nrow) );
      /* Adot_adjoint = -Bdot_adjoint Y^T */
      mAdot_adjoint_map = ND_MAIN_TRI(nimDerivs_matmult(-mBdot_adjoint_map, mYmapC.transpose()));
      // mAdot_adjoint_map = nimDerivs_matmult(-mBdot_adjoint_map, mYmapC.transpose()).template triangularView<Eigen::Lower>();

    /* Note: A^-1 Z A^-1 = solve(A, solve(A^T, Z^T)^T) */
      
    /* Badjoint = A^-T * Yadjoint - (A^-1 Adot A^-1)^T  Ydot_adjoint */      

      // In meta case, this is materialized on pass-by-copy to nimDerivs_EIGEN_BS[FS] anyway,
      // so we don't need this "_stuff" temporary variable to keep Eigen happy.
      //      MatrixXd_CppAD mAdot_stuff = mAdot_map.transpose().template triangularView<Eigen::Upper>(); // doesn't work inside solve
      
      // Update Badjoint even if B is constant, pointing to internal memory, for purposes of Aadjoint
      //      mBadjoint_map -= nimDerivs_matmult( ND_META_TRANS_SOLVE( mAmap,
      //						      ND_META_TRANS_SOLVE(mAmap.transpose(), mAdot_stuff).transpose() ).transpose() ,
      //				  mYdot_adjoint_map);

      mBadjoint_map -= nimDerivs_matmult(
					 ND_META_SOLVE(mAmap,
						       ND_META_TRANS_SOLVE(ND_TRANS_TRI(mAmap),
									   ND_TRANS_TRI(mAdot_map)).transpose()).transpose() ,
					 mYdot_adjoint_map);

      // mBadjoint_map -= nimDerivs_matmult(
      // 					 nimDerivs_EIGEN_FS(mAmap,
      // 							    nimDerivs_EIGEN_BS(mAmap.transpose().template triangularView<Eigen::Upper>(),
      // 									       mAdot_map.transpose().template triangularView<Eigen::Upper>()).transpose()).transpose() ,
      // 					 mYdot_adjoint_map);

      /* Aadjoint = -Badjoint * Y^T  - (A^-T Ydot_adjoint Ydot^T)= -Badjoint * Y^T  - ( Bdot_adjoint Ydot^T) , including both Badjoint terms*/
      mAadjoint_map = ND_MAIN_TRI(nimDerivs_matmult(-mBadjoint_map, mYmapC.transpose()));
      //      mAadjoint_map = nimDerivs_matmult(-mBadjoint_map, mYmapC.transpose()).template triangularView<Eigen::Lower>();
      mAadjoint_map -= ND_MAIN_TRI(nimDerivs_matmult(mBdot_adjoint_map, mYdot_map.transpose()));
      //      mAadjoint_map -= nimDerivs_matmult(mBdot_adjoint_map, mYdot_map.transpose()).template triangularView<Eigen::Lower>();
      //    Aadjoint_map = -Badjoint_map * Ymap.transpose() - Bdot_adjoint_map * Ydot_map.transpose();
      
      //    Adot_adjoint_map = -Bdot_adjoint_map * Ymap.transpose();
    }
  }
  clear_Badjoint_AD_memory();
  clear_Bdot_adjoint_AD_memory();
  if(order_up >= 2) {
    printf("Unsupported reverse order requested\n");
      return false;
  }
  return true;
}


