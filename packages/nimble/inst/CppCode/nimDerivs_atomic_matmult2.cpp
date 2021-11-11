#include <nimble/nimDerivs_atomic_matmult.h>
#include <array>
/* This atomic does not need caching because adjoints in reverse mode do not use Y or Ydot. */

// Use of the "new in place", labeled here "NEW_DYNAMIC"  was a performance experiment.
// The idea was to have member data for every Eigen map to be used and then initialize them
// with new in place.  This would save instantiating new maps every time one is needed.
// The performance impact was negligible, so we generally operate with USE_NEW_DYNAMIC undefined (false, off, unused).
// But actually, the new in place works well to manage that some inputs might be constants.

// #define VERBOSE_ATOMIC_MATMULT
#define VERBOSE_ATOMIC_MATMULT_REGIONS

// #define VERBOSE_TRIANGULAR_CASES


atomic_matmult_class::atomic_matmult_class(const std::string& name) :
  CppAD::atomic_three<double>(name),
  n1_(0),
  X1cat_(unknown), X2cat_(unknown),
  x1_is_constant_(false), x2_is_constant_(false),
  x1_is_variable_(true), x2_is_variable_(true),
  X1mapC(0, 0, 0, EigStrDyn(0, 0)),
  X2mapC(0, 0, 0, EigStrDyn(0, 0)),
  dX1mapC(0, 0, 0, EigStrDyn(0, 0)), dX2mapC(0, 0, 0, EigStrDyn(0, 0)),
  Yadjoint_mapC(0, 0, 0, EigStrDyn(0, 0)), X1dot_mapC(0, 0, 0, EigStrDyn(0, 0)),
  X2dot_mapC(0, 0, 0, EigStrDyn(0, 0)), Ydot_adjoint_mapC(0, 0, 0, EigStrDyn(0, 0)),
  Ymap(0, 0, 0, EigStrDyn(0, 0)), dY_map(0, 0, 0, EigStrDyn(0, 0)),
  X1adjoint_map(0, 0, 0, EigStrDyn(0, 0)), X2adjoint_map(0, 0, 0, EigStrDyn(0, 0)),
  X1dot_adjoint_map(0, 0, 0, EigStrDyn(0, 0)),  X2dot_adjoint_map(0, 0, 0, EigStrDyn(0, 0)),
  mX1mapC(0, 0, 0, EigStrDyn(0, 0)), mX2mapC(0, 0, 0, EigStrDyn(0, 0)),
  mdX1mapC(0, 0, 0, EigStrDyn(0, 0)), mdX2mapC(0, 0, 0, EigStrDyn(0, 0)),
  mYadjoint_mapC(0, 0, 0, EigStrDyn(0, 0)), mX1dot_mapC(0, 0, 0, EigStrDyn(0, 0)),
  mX2dot_mapC(0, 0, 0, EigStrDyn(0, 0)), mYdot_adjoint_mapC(0, 0, 0, EigStrDyn(0, 0)),
  mYmap(0, 0, 0, EigStrDyn(0, 0)), mdY_map(0, 0, 0, EigStrDyn(0, 0)),
  mX1adjoint_map(0, 0, 0, EigStrDyn(0, 0)), mX2adjoint_map(0, 0, 0, EigStrDyn(0, 0)),
  mX1dot_adjoint_map(0, 0, 0, EigStrDyn(0, 0)),  mX2dot_adjoint_map(0, 0, 0, EigStrDyn(0, 0))  
{};

bool atomic_matmult_class::for_type(
				    const CppAD::vector<double>&               parameter_x ,
				    const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				    CppAD::vector<CppAD::ad_type_enum>&        type_y      )
{
#ifdef VERBOSE_ATOMIC_MATMULT
  printf("In matmult for_type\n");
#endif
  int n = static_cast<double>(type_x.size());
  int m = type_y.size();
  int n1 = get_n1(); //parameter_x[0];
  //std::cout<<"n1 = "<<n1<<std::endl;
  //std::cout<<"constant "<<CppAD::constant_enum<<std::endl;
  //std::cout<<"dynamic "<<CppAD::dynamic_enum<<std::endl;
  //std::cout<<"variable "<<CppAD::variable_enum<<std::endl;
  //    if(type_x[0] == CppAD::constant_enum) std::cout<<"constant"<<std::endl;
  //    if(type_x[0] == CppAD::dynamic_enum) std::cout<<"dynamic"<<std::endl;
  //    if(type_x[0] == CppAD::variable_enum) std::cout<<"variable"<<std::endl;
  int n3 = m/n1;
  int n2;
  if(X1constant()) n2 = n / n3;
  else if(X2constant()) n2 = n / n1; 
  else n2 = n/(n1 + n3);

  // HOT-WIRED FOR TESTING
  // for(size_t i = 0; i < type_y.size(); ++i)
  //   type_y[i] = CppAD::variable_enum;
  // return true;

#ifdef VERBOSE_ATOMIC_MATMULT
  std::cout << "n1 = "<< n1 <<" n2 = "<< n2 <<" n3 = "<< n3 <<" "<<X1constant()<<" "<<X2constant()<<std::endl;
#endif

  CppAD::vector<CppAD::ad_type_enum> x1RowTypes(n1);
  CppAD::vector<CppAD::ad_type_enum> x2ColTypes(n3);

  CppAD::ad_type_enum this_row_type;
  CppAD::ad_type_enum item_type;

  // Find the maximal type (constant, dynamic, variable) for each row of x1
  for(size_t i = 0; i < n1; ++i) {
    this_row_type = CppAD::constant_enum;
    if(!X1constant()) {
      for(size_t j = 0; j < n2; ++j) {
	item_type = type_x[i + j*n1];
	//  std::cout<<"x1("<<i<<","<< j<<") type = "<<item_type<<"\n";
	if(item_type == CppAD::variable_enum) {
	  this_row_type = CppAD::variable_enum;
	  break;
	} else {
	  if(item_type == CppAD::dynamic_enum)
	    this_row_type = CppAD::dynamic_enum;
	}
      }
    }
    x1RowTypes[i] = this_row_type;
  }
  int x2offset = X1constant() ? 0 : n1*n2;
  CppAD::ad_type_enum this_col_type;
  // Find the maximal type for each column
  for(size_t j = 0; j < n3; ++j) {
    this_col_type = CppAD::constant_enum;
    if(!X2constant()) {
      for(size_t i = 0; i < n2; ++i) {
	item_type = type_x[x2offset + i + j*n2];
	//std::cout<<"x2("<<i<<","<< j<<") type = "<<item_type<<"\n";
	if(item_type == CppAD::variable_enum) {
	  this_col_type = CppAD::variable_enum;
	  break;
	} else {
	  if(item_type == CppAD::dynamic_enum)
	    this_col_type = CppAD::dynamic_enum;
	}
      }
    }
    x2ColTypes[j] = this_col_type;
  }

  for(size_t i = 0; i < n1; ++i) {
    for(size_t j = 0; j < n3; ++j) {
      item_type = CppAD::constant_enum;
      if(x1RowTypes[i] == CppAD::variable_enum || x2ColTypes[j] == CppAD::variable_enum) {
	item_type = CppAD::variable_enum;
      } else {
	if(x1RowTypes[i] == CppAD::dynamic_enum || x2ColTypes[j] == CppAD::dynamic_enum) {
	  item_type = CppAD::dynamic_enum;
	}
      }
      type_y[i + j*n1] = item_type;
      //std::cout<<"y("<<i<<","<<j<<") type = "<<item_type<<"\n";
    }
  }
  //    size_t n = type_y.size();
  //    for(size_t i = 0; i < n; ++i) type_y[i] = CppAD::variable_enum;
  return true;
}

bool atomic_matmult_class::rev_depend(
				      const CppAD::vector<double>&          parameter_x ,
				      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				      CppAD::vector<bool>&                depend_x    ,
				      const CppAD::vector<bool>&          depend_y
				      ) {
#ifdef VERBOSE_ATOMIC_MATMULT
  printf("In matmult reverse_depend\n");
#endif
  // Here we put true in an element of depend_x if there are true depend_y elements that are functions of the depend_x element
  int n = static_cast<double>(type_x.size());
  int m = depend_y.size();
  int n1 = get_n1();// parameter_x[0];
  // std::cout<<"n1 = "<<n1<<std::endl;
  //    if(type_x[0] == CppAD::constant_enum) std::cout<<"constant"<<std::endl;
  //    if(type_x[0] == CppAD::dynamic_enum) std::cout<<"dynamic"<<std::endl;
  //    if(type_x[0] == CppAD::variable_enum) std::cout<<"variable"<<std::endl;
  int n3 = m/n1;
  int n2;
  if(X1constant()) n2 = n / n3;
  else if(X2constant()) n2 = n / n1; 
  else n2 = n/(n1 + n3);
  

  // HOT-WIRED FOR TESTING
  // for(size_t i = 0; i < depend_x.size(); ++i)
  //   depend_x[i] = true;
  // return true;
  
  CppAD::vector<bool> depend_x1Row(n1);
  CppAD::vector<bool> depend_x2Col(n3);
  for(size_t i = 0; i < n1; ++i) depend_x1Row[i] = false;
  for(size_t j = 0; j < n3; ++j) depend_x2Col[j] = false;
  bool this_depend;
  // If y[i, j] has depend == true, then
  // the row x1[i,] and the column x2[,j] must return depend = true.
  for(size_t i = 0; i < n1; ++i) {
    for(size_t j = 0; j < n3; ++j) {
      this_depend = depend_y[i + j*n1];
      depend_x1Row[i] |= this_depend;
      depend_x2Col[j] |= this_depend;
    }
  }
  //
  if(!X1constant()) {
    for(size_t i = 0; i < n1; ++i) {
      this_depend = depend_x1Row[i];
      for(size_t j = 0; j < n2; ++j) {
	depend_x[i + j*n1] = this_depend;
      }
    }
  }
  if(!X2constant()) {
    int x2offset = X1constant() ? 0 : n1*n2;
    for(size_t j = 0; j < n3; ++j) {
      this_depend = depend_x2Col[j];
      for(size_t i = 0; i < n2; ++i) {
	depend_x[x2offset + i + j*n2] = this_depend;
      }
    }
  }
  return true;
}

bool atomic_matmult_class::forward(
				   const CppAD::vector<double>&               parameter_x  ,
				   const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
				   size_t                              need_y       ,
				   size_t                              order_low    ,
				   size_t                              order_up     ,
				   const CppAD::vector<double>&               taylor_x     ,
				   CppAD::vector<double>&                     taylor_y     ) {
  //forward mode
  int nrow = order_up + 1;

  printf("In matmult forward\n");
  std::cout<<"need_y = "<<need_y<<std::endl;
  for( int i = 0; i < type_x.size(); ++i) std::cout<<type_x[i]<<"\t";
  std::cout<<std::endl;
  for( int i = 0; i < parameter_x.size(); ++i) std::cout<<parameter_x[i]<<"\t";
  std::cout<<std::endl;
  for( int i = 0; i < taylor_x.size(); ++i) std::cout<<parameter_x[i]<<"\t";
  std::cout<<std::endl;
  
#ifdef VERBOSE_ATOMIC_MATMULT
  printf("In matmult forward\n");
  std::cout<<"order_low = "<<order_low<<" order_up = "<<order_up<<" nrow = "<<nrow<<std::endl;
#endif
  int n = static_cast<double>(taylor_x.size()/nrow); // size of input
  int m = taylor_y.size() / nrow;                    // size of output
  int n1 = get_n1();                                 // rows of X1
  int n3 = m/n1;                                     // cols of X2
  int n2;
  if(X1constant()) n2 = n / n3;
  else if(X2constant()) n2 = n / n1; 
  else n2 = n/(n1 + n3);

  typedef EigenTemplateTypes<double>::typeEigenConstMapStrd EMapC;
  typedef EigenTemplateTypes<double>::typeEigenMapStrd      EMap;

  const double *Xptr;
  int row_mult;
  
  if(X1constant())
    {Xptr = get_X_stored_ptr(); row_mult = 1;} else
    {Xptr = &taylor_x[0];       row_mult = nrow;}
  new (&X1mapC) EMapC(Xptr, n1, n2, EigStrDyn(row_mult*n1, row_mult) );

#ifdef VERBOSE_ATOMIC_MATMULT
  //  cout<<"X1\n"<<X1mapC<<endl;
#endif
  if(X2constant()) 
    {Xptr = get_X_stored_ptr(); row_mult = 1; } else 
    { row_mult = nrow;
      if(X1constant())
	{Xptr = &taylor_x[0]; } else
	{Xptr = &taylor_x[0 + (n1*n2)*nrow]; }
    }
  new (&X2mapC) EMapC(Xptr, n2, n3, EigStrDyn(row_mult*n2, row_mult ) );
#ifdef VERBOSE_ATOMIC_MATMULT
  //  cout<<"X2\n"<<X2mapC<<endl;
#endif
  
  if(order_low <= 0 & order_up >= 0) { // value
    // We could compile different cases depending on need for strides or not.
    new (&Ymap) EMap(&taylor_y[0], n1, n3, EigStrDyn(nrow*n1, nrow ) );
    // Calculate value for any cases of dynamic or parameter CppAD types
    matmult_internal_respecting_upper_lower(X1mapC, X2mapC, Ymap, X1cat(), X2cat());
#ifdef VERBOSE_ATOMIC_MATMULT
    //    cout<<"Y\n"<<Ymap<<endl;
#endif
  }
  if(order_low <= 1 & order_up >= 1) { // forward 1
    // printf("In forward >1\n");
    new (&dY_map) EMap(&taylor_y[1], n1, n3, EigStrDyn(nrow*n1, nrow ) );
    if(!X1constant()) {
      new (&dX1mapC) EMapC(&taylor_x[1], n1, n2, EigStrDyn(nrow*n1, nrow) );
      if(X1variable())
	matmult_internal_respecting_upper_lower(dX1mapC, X2mapC, dY_map, X1cat(), X2cat());
      //dY_map = dX1mapC * X2mapC;
    }
    if(X2variable()) {
      if(!X1constant()) {
	new (&dX2mapC) EMapC(&taylor_x[1 + (n1*n2)*nrow], n2, n3, EigStrDyn(nrow*n2, nrow ) );
      } else {
	new (&dX2mapC) EMapC(&taylor_x[1 ], n2, n3, EigStrDyn(nrow*n2, nrow ) );
      }
      if(X1variable()) {
	matmult_internal_respecting_upper_lower_add(X1mapC, dX2mapC, dY_map, X1cat(), X2cat());
	// dY_map += X1mapC * dX2mapC;
      } else {
	matmult_internal_respecting_upper_lower(X1mapC, dX2mapC, dY_map, X1cat(), X2cat());
	// dY_map = X1mapC * dX2mapC;
      }
    }
    if(!(X1variable() || X2variable())) {
      dY_map.fill(0);
    }
#ifdef VERBOSE_ATOMIC_MATMULT
    //    cout<<"dY_map\n"<<dY_map<<endl;
#endif
    // dY_map = dX1mapC * X2mapC + X1mapC * dX2mapC;
  }
#ifdef VERBOSE_ATOMIC_MATMULT
  printf("Leaving matmult forward\n");
#endif
  return true;
}

bool atomic_matmult_class::forward(
				   const CppAD::vector<CppAD::AD<double> >&               parameter_x  ,
				   const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
				   size_t                              need_y       ,
				   size_t                              order_low    ,
				   size_t                              order_up     ,
				   const CppAD::vector<CppAD::AD<double> >&               taylor_x     ,
				   CppAD::vector<CppAD::AD<double> >&                     taylor_y     ) {
  //forward mode
  int nrow = order_up + 1;
#ifdef VERBOSE_ATOMIC_MATMULT
  printf("In matmult meta-forward\n");
  std::cout<<"order_low = "<<order_low<<" order_up = "<<order_up<<" nrow = "<<nrow<<std::endl;
#endif
  int n = static_cast<double>(taylor_x.size()/nrow);
  int m = taylor_y.size() / nrow;
  int n1 = get_n1(); //CppAD::Value(taylor_x[0]);
  int n3 = m/n1;
  int n2;
  if(X1constant()) n2 = n / n3;
  else if(X2constant()) n2 = n / n1; 
  else n2 = n/(n1 + n3);

#ifdef VERBOSE_ATOMIC_MATMULT
  std::cout << "n1 = "<< n1 <<" n2 = "<< n2 <<" n3 = "<< n3 <<" "<<X1constant()<<" "<<X2constant()<<std::endl;
#endif

  typedef EigenTemplateTypes<CppAD::AD<double> >::typeEigenConstMapStrd EMapC;
  typedef EigenTemplateTypes<CppAD::AD<double> >::typeEigenMapStrd      EMap;
  
  const CppAD::AD<double> *Xptr;
  int row_mult;

  if(X1constant())
    {fill_X_AD_stored(); Xptr = get_X_AD_stored_ptr(); row_mult = 1;} else
    {                    Xptr = &taylor_x[0];       row_mult = nrow;}
  new (&mX1mapC) EMapC(Xptr, n1, n2, EigStrDyn(row_mult*n1, row_mult) );
  
  if(X2constant()) 
    {fill_X_AD_stored(); Xptr = get_X_AD_stored_ptr(); row_mult = 1; } else 
    { row_mult = nrow;
      if(X1constant())
	{Xptr = &taylor_x[0]; } else
	{Xptr = &taylor_x[0 + (n1*n2)*nrow]; }
    }
  new (&mX2mapC) EMapC(Xptr, n2, n3, EigStrDyn(row_mult*n2, row_mult ) );  
  if(order_low <= 0 & order_up >= 0) { // value
    new (&mYmap) EMap(&taylor_y[0], n1, n3, EigStrDyn(nrow*n1, nrow ) );
    // This and other terms do not use X1cat() and X2cat() (triangular cases),
    // so the assumption is the zeros really are there are will be re-assessed as zeros
    // when the values pass through atomic_matmult again.
#ifdef VERBOSE_ATOMIC_MATMULT
    std::cout<<"recording mX1mapC %*% mX2mapC"<<std::endl;
#endif
    atomic_matmult(mX1mapC, mX2mapC, mTerm1); //Ymap = X1map * X2map;
#ifdef VERBOSE_ATOMIC_MATMULT
    std::cout<<"done recording mX1mapC %*% mX2mapC"<<std::endl;
#endif
    mYmap = mTerm1;
  }

  if(order_low <= 1 & order_up >= 1) { // forward 1
    //    printf("In forward >1\n");
    new (&mdY_map) EMap(&taylor_y[1], n1, n3, EigStrDyn(nrow*n1, nrow ) );
    if(!X1constant()) {
      new (&mdX1mapC) EMapC(&taylor_x[1], n1, n2, EigStrDyn(nrow*n1, nrow) );
#ifdef VERBOSE_ATOMIC_MATMULT
    std::cout<<"recording mdX1mapC %*% mX2mapC"<<std::endl;
#endif
    if(X1variable())
      atomic_matmult(mdX1mapC, mX2mapC, mTerm1);
#ifdef VERBOSE_ATOMIC_MATMULT
    std::cout<<"done recording mdX1mapC %*% mX2mapC"<<std::endl;
#endif
      mdY_map = mTerm1;
    }
    if(X2variable()) {
#ifdef VERBOSE_ATOMIC_MATMULT
    std::cout<<"recording mX1mapC %*% mdX2mapC"<<std::endl;
#endif
      if(!X1constant()) {
	new (&mdX2mapC) EMapC(&taylor_x[1 + (n1*n2)*nrow], n2, n3, EigStrDyn(nrow*n2, nrow ) );
      } else {
	new (&mdX2mapC) EMapC(&taylor_x[1 ], n2, n3, EigStrDyn(nrow*n2, nrow ) );
      }
      atomic_matmult(mX1mapC, mdX2mapC, mTerm2);
      if(X1variable()) {
	mdY_map += mTerm2;
	// dY_map += X1mapC * dX2mapC;
      } else {
	mdY_map = mTerm2;
	// dY_map = X1mapC * dX2mapC;
      }
      if(!(X1variable() || X2variable())) {
	dY_map.fill(0);
      }
#ifdef VERBOSE_ATOMIC_MATMULT
      std::cout<<"done recording mX1mapC %*% mdX2mapC"<<std::endl;
#endif
    }
    // dY_map = dX1mapC * X2mapC + X1mapC * dX2mapC;
  }
#ifdef VERBOSE_ATOMIC_MATMULT
  printf("Leaving meta-forward\n");
#endif
  return true;
}

bool atomic_matmult_class::reverse(
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
#ifdef VERBOSE_ATOMIC_MATMULT
  printf("In matmult reverse\n");
  std::cout<<"nrow = "<<nrow<<" order_up = "<<order_up<<std::endl;
#endif
  int n = static_cast<double>(taylor_x.size()/nrow);
  int m = taylor_y.size() / nrow;
  int n1 = get_n1();
  int n3 = m/n1;
  int n2;
  if(X1constant()) n2 = n / n3;
  else if(X2constant()) n2 = n / n1; 
  else n2 = n/(n1 + n3);

#ifdef VERBOSE_ATOMIC_MATMULT
  std::cout << "n1 = "<< n1 <<" n2 = "<< n2 <<" n3 = "<< n3 <<" "<<X1constant()<<" "<<X2constant()<<std::endl;
#endif

  typedef EigenTemplateTypes<double>::typeEigenConstMapStrd EMapC;
  typedef EigenTemplateTypes<double>::typeEigenMapStrd      EMap;

  const double *Xptr;
  double *pXptr;
  int row_mult;
  
  if(X1constant())
    {Xptr = get_X_stored_ptr(); row_mult = 1;} else
    {Xptr = &taylor_x[0];       row_mult = nrow;
      if(X1variable())
	new (&X1adjoint_map) EMap(&partial_x[0], n1, n2, EigStrDyn(nrow*n1, nrow) );}
  new (&X1mapC) EMapC(Xptr, n1, n2, EigStrDyn(row_mult*n1, row_mult) );
  
  if(X2constant()) 
    {Xptr = get_X_stored_ptr(); row_mult = 1; } else 
    { row_mult = nrow;
      if(X1constant())
	{Xptr = &taylor_x[0];                pXptr = &partial_x[0];} else
	{Xptr = &taylor_x[0 + (n1*n2)*nrow]; pXptr = &partial_x[0 + (n1*n2)*nrow];}
      if(X2variable())
	new (&X2adjoint_map) EMap(pXptr, n2, n3, EigStrDyn(nrow*n2, nrow ) );
    }
  new (&X2mapC) EMapC(Xptr, n2, n3, EigStrDyn(row_mult*n2, row_mult ) );
    
  new (&Yadjoint_mapC) EMapC(&partial_y[0], n1, n3, EigStrDyn(nrow*n1, nrow ) );

#ifdef VERBOSE_ATOMIC_MATMULT
  //  cout<<"X1\n"<<X1mapC<<endl;
  // cout<<"X2\n"<<X2mapC<<endl;
  // cout<<"Yadjoint\n"<<Yadjoint_mapC<<endl;
#endif
  
  if(order_up >= 0) {
    // reverse 1
    //    partial_x[0] = 0;
    if(X1variable()) {
      matmult_internal_respecting_upper_lower(Yadjoint_mapC, X2mapC.transpose(), X1adjoint_map,
					      unknown, transpose_X2cat());
      //    X1adjoint_map = Yadjoint_mapC * X2mapC.transpose();
#ifdef VERBOSE_ATOMIC_MATMULT
      //      std::cout<<"X1adjoint\n"<<X1adjoint_map<<std::endl;
#endif
    }
    if(X2variable()) {
      matmult_internal_respecting_upper_lower( X1mapC.transpose(), Yadjoint_mapC, X2adjoint_map,
					       transpose_X1cat(), unknown );
      // X2adjoint_map = X1mapC.transpose() * Yadjoint_mapC;
#ifdef VERBOSE_ATOMIC_MATMULT
      //      std::cout<<"X2adjoint\n"<<X2adjoint_map<<std::endl;
#endif
    }
  }
  if(order_up >= 1) {
    // reverse 2
    new (&Ydot_adjoint_mapC) EMapC(&partial_y[1], n1, n3, EigStrDyn(nrow*n1, nrow ) );
#ifdef VERBOSE_ATOMIC_MATMULT
    //    cout<<"Ydot_adjoint\n"<<Ydot_adjoint_mapC<<endl;
#endif
    
    if(X1variable()) {
      new (&X1dot_mapC) EMapC(&taylor_x[1], n1, n2, EigStrDyn(nrow*n1, nrow) );
      new (&X1dot_adjoint_map) EMap(&partial_x[1], n1, n2, EigStrDyn(nrow*n1, nrow) );
#ifdef VERBOSE_ATOMIC_MATMULT
      //      cout<<"X1dot\n"<<X1dot_mapC<<endl;
#endif
    }
    if(X2variable()) {
      if(X1variable()) {
	Xptr = &taylor_x[1 + (n1*n2)*nrow]; pXptr = &partial_x[1 + (n1*n2)*nrow];
      } else {
	Xptr = &taylor_x[1]; pXptr = &partial_x[1];
      }
      new (&X2dot_mapC) EMapC(Xptr, n2, n3, EigStrDyn(nrow*n2, nrow ) );
      new (&X2dot_adjoint_map) EMap(pXptr, n2, n3, EigStrDyn(nrow*n2, nrow ));
#ifdef VERBOSE_ATOMIC_MATMULT
      //      cout<<"X2dot\n"<<X2dot_mapC<<endl;
#endif
    }

    if((X1variable()) && (X2variable())) {
      matmult_internal_respecting_upper_lower_add(Ydot_adjoint_mapC, X2dot_mapC.transpose(), X1adjoint_map,
						  unknown, transpose_X2cat());
      //   X1adjoint_map += Ydot_adjoint_mapC * X2dot_mapC.transpose();
#ifdef VERBOSE_ATOMIC_MATMULT
      //      std::cout<<"X1adjoint\n"<<X1adjoint_map<<std::endl;
#endif
      
      matmult_internal_respecting_upper_lower_add(X1dot_mapC.transpose(), Ydot_adjoint_mapC, X2adjoint_map,
						   transpose_X1cat(), unknown);
#ifdef VERBOSE_ATOMIC_MATMULT
      //      std::cout<<"X2adjoint\n"<<X2adjoint_map<<std::endl;
#endif
      //   X2adjoint_map += X1dot_mapC.transpose() * Ydot_adjoint_mapC;
    }
    if(X1variable()) {
      matmult_internal_respecting_upper_lower(Ydot_adjoint_mapC, X2mapC.transpose(), X1dot_adjoint_map,
					      unknown, transpose_X2cat());
      //      X1dot_adjoint_map = Ydot_adjoint_mapC * X2mapC.transpose();
#ifdef VERBOSE_ATOMIC_MATMULT
      //      cout<<"X1dot_adjoint\n"<<X1dot_adjoint_map<<endl;
#endif
    }
    if(X2variable()) {
      matmult_internal_respecting_upper_lower(X1mapC.transpose(), Ydot_adjoint_mapC, X2dot_adjoint_map,
						  transpose_X1cat(), unknown);
#ifdef VERBOSE_ATOMIC_MATMULT
      //      cout<<"X2dot_adjoint\n"<<X2dot_adjoint_map<<endl;
#endif
      //      X2dot_adjoint_map = X1mapC.transpose() * Ydot_adjoint_mapC;
    }
  }
#ifdef VERBOSE_ATOMIC_MATMULT
  std::cout<<"Leaving matmult reverse"<<std::endl;
#endif

  if(order_up >= 2) {
    printf("Unsupported reverse order requested\n");
    return false;
  }
  return true;
}

bool atomic_matmult_class::reverse(
				   const CppAD::vector<CppAD::AD<double> >&               parameter_x ,
				   const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				   size_t                              order_up    ,
				   const CppAD::vector<CppAD::AD<double> >&               taylor_x    ,
				   const CppAD::vector<CppAD::AD<double> >&               taylor_y    ,
				   CppAD::vector<CppAD::AD<double> >&                     partial_x   ,
				   const CppAD::vector<CppAD::AD<double> >&               partial_y   )
{
  //reverse mode
  int nrow = order_up + 1;
#ifdef VERBOSE_ATOMIC_MATMULT
  printf("In matmult meta-reverse\n");
  std::cout<<"nrow = "<<nrow<<" order_up = "<<order_up<<std::endl;
#endif
  int n = static_cast<double>(taylor_x.size()/nrow);
  int m = taylor_y.size() / nrow;
  int n1 = get_n1(); //CppAD::Value(taylor_x[0]);
  int n3 = m/n1;
  int n2;
  if(X1constant()) n2 = n / n3;
  else if(X2constant()) n2 = n / n1; 
  else n2 = n/(n1 + n3);

#ifdef VERBOSE_ATOMIC_MATMULT
  std::cout << "n1 = "<< n1 <<" n2 = "<< n2 <<" n3 = "<< n3 <<" "<<X1constant()<<" "<<X2constant()<<std::endl;
#endif

  typedef EigenTemplateTypes<CppAD::AD<double> >::typeEigenConstMapStrd EMapC;
  typedef EigenTemplateTypes<CppAD::AD<double> >::typeEigenMapStrd      EMap;

  const CppAD::AD<double> *Xptr;
  CppAD::AD<double> *pXptr;
  int row_mult;

  
  // new (&mX1mapC) EigenTemplateTypes<CppAD::AD<double> >::typeEigenConstMapStrd(&taylor_x[0],
  //									       n1, n2, EigStrDyn(nrow*n1, nrow) );

  if(X1constant())
    {fill_X_AD_stored(); Xptr = get_X_AD_stored_ptr(); row_mult = 1;} else
    {                    Xptr = &taylor_x[0];          row_mult = nrow;
      if(X1variable())
	new (&mX1adjoint_map) EMap(&partial_x[0], n1, n2, EigStrDyn(nrow*n1, nrow) );}
  new (&mX1mapC) EMapC(Xptr, n1, n2, EigStrDyn(row_mult*n1, row_mult) );
  
  if(X2constant()) 
    {fill_X_AD_stored(); Xptr = get_X_AD_stored_ptr(); row_mult = 1; } else 
    { row_mult = nrow;
      if(X1constant())
	{Xptr = &taylor_x[0];                pXptr = &partial_x[0];} else
	{Xptr = &taylor_x[0 + (n1*n2)*nrow]; pXptr = &partial_x[0 + (n1*n2)*nrow];}
      if(X2variable())
	new (&mX2adjoint_map) EMap(pXptr, n2, n3, EigStrDyn(nrow*n2, nrow ) );
    }
  new (&mX2mapC) EMapC(Xptr, n2, n3, EigStrDyn(row_mult*n2, row_mult ) );
  new (&mYadjoint_mapC) EMapC(&partial_y[0], n1, n3, EigStrDyn(nrow*n1, nrow ) );

  if(order_up >= 0) {
    // reverse 1
    //    partial_x[0] = 0;
    if(X1variable()) {
      atomic_matmult(mYadjoint_mapC, mX2mapC.transpose(), mTerm1);
      mX1adjoint_map = mTerm1;
      //    X1adjoint_map = Yadjoint_mapC * X2mapC.transpose();
    }
    if(X2variable()) {
      atomic_matmult(mX1mapC.transpose(), mYadjoint_mapC,  mTerm2);
      mX2adjoint_map = mTerm2;
      // X2adjoint_map = X1mapC.transpose() * Yadjoint_mapC;
    }
  }

  if(order_up >= 1) {
    // reverse 2
    new (&mYdot_adjoint_mapC) EMapC(&partial_y[1], n1, n3, EigStrDyn(nrow*n1, nrow ) );
    
    if(X1variable()) {
      new (&mX1dot_mapC) EMapC(&taylor_x[1], n1, n2, EigStrDyn(nrow*n1, nrow) );
      new (&mX1dot_adjoint_map) EMap(&partial_x[1], n1, n2, EigStrDyn(nrow*n1, nrow) );
    }
    if(X2variable()) {
      if(X1variable()) {
	Xptr = &taylor_x[1 + (n1*n2)*nrow]; pXptr = &partial_x[1 + (n1*n2)*nrow];
      } else {
	Xptr = &taylor_x[1]; pXptr = &partial_x[1];
      }
      new (&mX2dot_mapC) EMapC(Xptr, n2, n3, EigStrDyn(nrow*n2, nrow ) );
      new (&mX2dot_adjoint_map) EMap(pXptr, n2, n3, EigStrDyn(nrow*n2, nrow ));
    }

    if((X1variable()) && (X2variable())) {
      atomic_matmult(mYdot_adjoint_mapC, mX2dot_mapC.transpose(), mTerm1);
      atomic_matmult(mX1dot_mapC.transpose() , mYdot_adjoint_mapC, mTerm2);
      mX1adjoint_map += mTerm1;
      mX2adjoint_map += mTerm2;
      //   X1adjoint_map += Ydot_adjoint_mapC * X2dot_mapC.transpose();
      //   X2adjoint_map += X1dot_mapC.transpose() * Ydot_adjoint_mapC;
    }
    if(X1variable()) {
      atomic_matmult(mYdot_adjoint_mapC, mX2mapC.transpose(), mTerm1);
      mX1dot_adjoint_map = mTerm1;
      //      X1dot_adjoint_map = Ydot_adjoint_mapC * X2mapC.transpose();
    }
    if(X2variable()) {
      atomic_matmult(mX1mapC.transpose(), mYdot_adjoint_mapC,  mTerm2);
      mX2dot_adjoint_map = mTerm2;
    //      X2dot_adjoint_map = X1mapC.transpose() * Ydot_adjoint_mapC;
    }
  }
#ifdef VERBOSE_ATOMIC_MATMULT
  std::cout<<"Leaving matmult meta-reverse"<<std::endl;
#endif
  if(order_up >= 2) {
    printf("Unsupported reverse order requested\n");
    return false;
  }
  return true;
}

void atomic_matmult_class::set_X_stored(const MatrixXd_CppAD &X) {
  int n1 = X.rows();
  int n2 = X.cols();
  X_stored.resize(n1 * n2);
  mat2vec_v(X, X_stored, 0);
}

// matmult atomic objects created by this factory are not
// cleaned up until the library is unloaded.
class atomic_matmult_factory {
public:
  std::vector<atomic_matmult_class *> created_objects;
  atomic_matmult_factory() {};
  ~atomic_matmult_factory() {
    for(size_t i = 0; i < created_objects.size(); ++i)
      delete created_objects[i];
  }
  atomic_matmult_class *create() {
    atomic_matmult_class *res = new atomic_matmult_class("atomic_matmult");
    created_objects.push_back(res);
    return res;
  }
};

atomic_matmult_factory global_matmult_factory;

matrix_category decide_matrix_category(const MatrixXd_CppAD &x) {
  int nRow = x.rows();
  int nCol = x.cols();
  if(nRow != nCol) return non_square;

  bool is_upper_diag(true);
  for(int i = 0; i < nRow; ++i) {
    for(int j = 0; j < i; ++j) {
      if(! (is_upper_diag &= CppAD::IdenticalZero(x(i, j) ) ) ) break;
    }
    if(!is_upper_diag) break;
  }
  if(is_upper_diag) return upper_diagonal;
  
  bool is_lower_diag(true);
  for(int i = 0; i < nRow; ++i) {
    for(int j = i+1; j < nCol; ++j) {
      if(! (is_lower_diag &= CppAD::IdenticalZero(x(i, j) ) ) ) break;
    }
    if(!is_lower_diag) break;
  }
  if(is_lower_diag) return lower_diagonal;

  return square_full;
}


void atomic_matmult_internal(const MatrixXd_CppAD &x1,
			     const MatrixXd_CppAD &x2,
			     MatrixXd_CppAD &y,
			     bool x1_is_constant, bool x2_is_constant) {
  int n1 = x1.rows(); // may not be general to all Eigen types
  int n2 = x1.cols();
  int n3 = x2.cols();
  
  if(n2 != x2.rows())
    cout<<"incommensurate matrices in atomic_matmult"<<endl;
  
  y.resize(n1, n3);

  if((n1 == 0) || (n2 == 0) || (n3 == 0))
    return;
  
#ifdef VERBOSE_ATOMIC_MATMULT
  cout<<"setting up an atomic call of dims ("<<n1<<", "<<n2<<", "<<x1_is_constant<<") %*% ("<<n2<<", "<<n3<<", "<<x2_is_constant<<") = ("<<n1<<", "<<n3<<")"<<endl;
#endif
  
  atomic_matmult_class * atomic_matmult = global_matmult_factory.create();
  atomic_matmult->X1cat() = decide_matrix_category(x1);
  atomic_matmult->X2cat() = decide_matrix_category(x2);
  atomic_matmult->X1constant() = x1_is_constant;
  atomic_matmult->X2constant() = x2_is_constant;

  auto dyn_cond = [](const CppAD::AD<double> &x)->bool {return CppAD::Dynamic(x);};
  bool x1_is_dynamic(false);
  int a, b, c, d; // dummies
  if(!x1_is_constant)
    x1_is_dynamic = delineate_condition_region(dyn_cont, x1, a, b, c, d);
  bool x2_is_dynamic(false);
  if(!x2_is_constant)
    x2_is_dynamic = delineate_condition_region(dyn_cont, x2, a, b, c, d);

  atomic_matmult->X1variable() = !(x1_is_constant || x1_is_dynamic);
  atomic_matmult->X2variable() = !(x2_is_constant || x2_is_dynamic);
  
  int xVecSize = 0;
  if(!x1_is_constant) xVecSize += n1*n2;
  if(!x2_is_constant) xVecSize += n2*n3;
  
  std::vector<CppAD::AD<double> > xVec(xVecSize);
  atomic_matmult->set_n1(n1);
  if(x1_is_constant) {
    atomic_matmult->set_X_stored(x1);
  } else {
    mat2vec(x1, xVec, 0);
  }
  if(x2_is_constant) {
    atomic_matmult->set_X_stored(x2);
  } else {
    if(x1_is_constant) {
      mat2vec(x2, xVec, 0);
    } else {
      mat2vec(x2, xVec, n1*n2);
    }
  }
  std::vector<CppAD::AD<double> > yVec(n1*n3);
  (*atomic_matmult)(xVec, yVec);
  vec2mat(yVec, y);
}

void atomic_matmult_diagnostic_message(int x1i, int x1j, int x1rows, int x1cols,
				       int x2i, int x2j, int x2rows, int x2cols,
				       int yi, int yj, int yrows, int ycols) {
  std::cout<<"y["<<yi<<":"<<yi + yrows-1<<", "<<yj<<":"<<yj+ycols-1<<"] += ";
  std::cout<<"x1["<<x1i<<":"<<x1i + x1rows-1<<", "<<x1j<<":"<<x1j+x1cols-1<<"] %*%";
  std::cout<<"x2["<<x2i<<":"<<x2i + x2rows-1<<", "<<x2j<<":"<<x2j+x2cols-1<<"]"<<std::endl;
}

void atomic_matmult_internal_const_var(const MatrixXd_CppAD &x1,
				       const std::array<int, 4> &x1i,
				       const MatrixXd_CppAD &x2,
				       const std::array<int, 4> &x2i,
				       MatrixXd_CppAD &y) {
  // The [x1|x2|y]i variables are in the order of block entries.
  // i.e. row start, col start, row extent, col extent

  // NZ = non-zero
  // Get x1startNZrow, x1endNZrow, x1startNZcol, x1endNZcol
  // If there are no all-zero rows or columns, these should equal x1i[0], x1i[1], x1i[0]+x1i[2], x1i[1]+x1i[3]
  // Then set x1NZrowExtent, x1NZcolExtent
  int x1startNZrow(x1i[0]);
  int x1endNZrow(x1i[0] + x1i[2]);
  int x1startNZcol(x1i[1]);
  int x1endNZcol(x1i[1] + x1i[3]);
  // std::cout<<"entering const_var delineate condition_region"<<std::endl;
  auto cond = [](const CppAD::AD<double> &x)->bool {return CppAD::IdenticalZero(x);};
  bool all_zero = delineate_condition_region(cond,
					     x1,
					     x1startNZrow, x1endNZrow, x1startNZcol, x1endNZcol,
					     false);
  // std::cout<<x1startNZrow<<" "<< x1endNZrow<<" "<< x1startNZcol<<" "<< x1endNZcol<<std::endl;
  if(all_zero) return;
  int x1NZrowExtent = x1endNZrow - x1startNZrow;
  int x1NZcolExtent = x1endNZcol - x1startNZcol;
  MatrixXd_CppAD y_const_var1;
  //  y_const_var1.resize(x1rowStart, x2colEnd - x2colStart);
  y_const_var1.resize(x1NZrowExtent, x2i[3]);
  // atomic_matmult_internal( x1.block(0, x2rowStart, x1rowStart, x2rowEnd - x2rowStart),
  // 			   x2.block(x2rowStart, x2colStart, x2rowEnd - x2rowStart, x2colEnd - x2colStart),
  // 			   y_const_var1, 1, 0);
  atomic_matmult_internal( x1.block(x1startNZrow, x1startNZcol, x1NZrowExtent, x1NZcolExtent),
			   x2.block(x2i[0] + (x1startNZcol - x1i[1]), x2i[1], x1NZcolExtent, x2i[3]),
			   y_const_var1, 1, 0);
  // y.block(0, x2colStart, x1rowStart, x2colEnd - x2colStart) += y_const_var1;
  y.block(x1startNZrow, x2i[1], x1NZrowExtent, x2i[3]) += y_const_var1;
}

void atomic_matmult_internal_var_const(const MatrixXd_CppAD &x1,
				       const std::array<int, 4> &x1i,
				       const MatrixXd_CppAD &x2,
				       const std::array<int, 4> &x2i,
				       MatrixXd_CppAD &y) {
  // The [x1|x2|y]i variables are in the order of block entries.
  // i.e. row start, col start, row extent, col extent

  // NZ = non-zero
  int x2startNZrow(x2i[0]);
  int x2endNZrow(x2i[0] + x2i[2]);
  int x2startNZcol(x2i[1]);
  int x2endNZcol(x2i[1] + x2i[3]);
  //  std::cout<<"entering var_const delineate condition_region"<<std::endl;
  auto cond = [](const CppAD::AD<double> &x)->bool {return CppAD::IdenticalZero(x);};
  bool all_zero = delineate_condition_region(cond,
					     x2,
					     x2startNZrow, x2endNZrow, x2startNZcol, x2endNZcol,
					     false);
  // std::cout<<x2startNZrow<<" "<< x2endNZrow<<" "<< x2startNZcol<<" "<< x2endNZcol<<std::endl;
  if(all_zero) return;
  int x2NZrowExtent = x2endNZrow - x2startNZrow;
  int x2NZcolExtent = x2endNZcol - x2startNZcol;
  MatrixXd_CppAD y_var_const1;
  y_var_const1.resize(x1i[2], x2NZcolExtent);
  atomic_matmult_internal( x1.block(x1i[0], x1i[1] + (x2startNZrow - x2i[0]), x1i[2], x2NZrowExtent),
			   x2.block(x2startNZrow, x2startNZcol, x2NZrowExtent, x2NZcolExtent),
			   y_var_const1, 0, 1);
  y.block(x1i[0], x2startNZcol, x1i[2], x2NZcolExtent) += y_var_const1;
}

void atomic_matmult(const MatrixXd_CppAD &x1,
		    const MatrixXd_CppAD &x2,
		    MatrixXd_CppAD &y,
		    bool debug) {
  int n1 = x1.rows(); // may not be general to all Eigen types
  int n2 = x1.cols();
  int n3 = x2.cols();

  if(n2 != x2.rows())
    cout<<"incommensurate matrices in atomic_matmult"<<endl;

  y.resize(n1, n3);

  if(debug) {
    cout<<"n1 = "<<n1<<" n2 = "<<n2<<" n3 = "<<n3<<endl;
    y.fill(0);
    //    return;
  }
  
  if( ( n1 == 1 ) && ( n2 == 1 ) && ( n3 == 1 ) ) {
    y(0,0) = x1(0,0) * x2(0,0);
    return;
  }
    
  int x1rowStart, x1rowEnd, x1colStart, x1colEnd;
  int x2rowStart, x2rowEnd, x2colStart, x2colEnd;
  bool x1_is_constant, x2_is_constant;

  // CppAD has three roles: constant, dynamic parameter, variable
  // Names of checking functions are confusing
  // CppAD::Constant(x) checks if x is constant (defined as not being part of the tape)
  // CppAD::Dynamic(x) if x is dynamic parameter or variable (This is the confusing one.  It is effectively !Constant(x), I believe.
  // CppAD::Parameter(x) checks if x is dynamic parameter
  // CppAD::Variable(x) checks if x is a variable
  //
  // Make this const_or_param
  auto const_cond = [](const CppAD::AD<double> &x)->bool {return CppAD::Constant(x);};
  // Make these x1_has_no_variables
  x1_is_constant = delineate_condition_region(const_cond, x1, x1rowStart, x1rowEnd, x1colStart, x1colEnd);
  x2_is_constant = delineate_condition_region(const_cond, x2, x2rowStart, x2rowEnd, x2colStart, x2colEnd);

  using std::cout;
#ifdef VERBOSE_ATOMIC_MATMULT_REGIONS
  cout<<"x1 region: "<<x1rowStart<<":"<<x1rowEnd<<", "<<x1colStart<<":"<<x1colEnd<<endl;
  cout<<"x2 region: "<<x2rowStart<<":"<<x2rowEnd<<", "<<x2colStart<<":"<<x2colEnd<<endl;
#endif
  if(debug) {
    cout<<"n1 = "<<n1<<" n2 = "<<n2<<" n3 = "<<n3<<endl;
    cout<<"x1\n";
    for(int iii = 0; iii < n1; ++iii) {
      for(int jjj = 0; jjj < n2; ++jjj)
	cout<<x1(iii,jjj)<<"\t";
      cout<<endl;
    }
      
    cout<<"x2\n";
    for(int iii = 0; iii < n2; ++iii) {
      for(int jjj = 0; jjj < n3; ++jjj)
	cout<<x2(iii,jjj)<<"\t";
      cout<<endl;
    }
    cout<<"x1 region: "<<x1rowStart<<":"<<x1rowEnd<<", "<<x1colStart<<":"<<x1colEnd<<endl;
    cout<<"x2 region: "<<x2rowStart<<":"<<x2rowEnd<<", "<<x2colStart<<":"<<x2colEnd<<endl;
    y.fill(0);
    //    return;
  }
  // We will go through the matrix multiplication inefficiently and do
  // every constant operation.  Then we will record multiple atomics for various blocks that have variables

  // y(i, k) = sum_j x1(i, j) x2(j, k)
  int B1B2s1End = x1colStart < x2rowStart ? x1colStart : x2rowStart; // first start. s1 is 0:first start
  int B1B2s3Sta = x1colEnd   > x2rowEnd   ? x1colEnd   : x2rowEnd;   // last end.    s3 is last end:n3
  int B1B2s2Sta = x1colEnd   < x2rowEnd   ? x1colEnd   : x2rowEnd;   // first end.   s2 is first end:last start. This is the region between non-overlapping variable regions and might often be empty
  int B1B2s2End = x1colStart > x2rowStart ? x1colStart : x2rowStart; // last start
  //                      A1*A2     , A1*B2 seg1, A1*B2 seg2, A1*C2     , B1*A2 seg1, B1*A2 seg2, B1*B2 seg1, B1*B2 seg2, B1*B2 seg3, B1*C2 seg1, B1*C2 seg2, C1*A2     , C1*B2 seg1, C1*B2 seg2, C1*C2
  const std::vector<int> iStart={0         , 0         , 0         , 0         , x1rowStart, x1rowStart, x1rowStart, x1rowStart, x1rowStart, x1rowStart, x1rowStart, x1rowEnd  , x1rowEnd  , x1rowEnd  , x1rowEnd };
  const std::vector<int> iEnd  ={x1rowStart, x1rowStart, x1rowStart, x1rowStart, x1rowEnd  , x1rowEnd  , x1rowEnd  , x1rowEnd  , x1rowEnd  , x1rowEnd  , x1rowEnd  , n1        , n1        , n1        , n1       };
  const std::vector<int> kStart={0         , x2colStart, x2colStart, x2colEnd  , 0         , 0         , x2colStart, x2colStart, x2colStart, x2colEnd  , x2colEnd  , 0         , x2colStart, x2colStart, x2colEnd };
  const std::vector<int> kEnd  ={x2colStart, x2colEnd  , x2colEnd  , n3        , x2colStart, x2colStart, x2colEnd  , x2colEnd  , x2colEnd  , n3        , n3        , x2colStart, x2colEnd  , x2colEnd  , n3       };
  const std::vector<int> jStart={0         , 0         , x2rowEnd  , 0         , 0         , x1colEnd  , 0         , B1B2s2Sta , B1B2s3Sta , 0         , x1colEnd  , 0         , 0         , x2rowEnd  , 0        };
  const std::vector<int> jEnd  ={n2        , x2rowStart, n2        , n2        , x1colStart, n2        , B1B2s1End , B1B2s2End , n2        , x1colStart, n2        , n2        , x2rowStart, n2        , n2       };

  y.fill(0);
  
  for(int seg = 0; seg < iStart.size(); ++seg) {
    for(int i = iStart[seg]; i < iEnd[seg]; ++i) {
      for(int k = kStart[seg]; k < kEnd[seg]; ++k) {
	CppAD::AD<double> oneVal = 0;
	for(int j = jStart[seg]; j < jEnd[seg]; ++j) {
	  if(!((CppAD::Constant(x1(i, j))) && (CppAD::Constant(x2(j, k))))) std::cout<<"There is some problem with a matrix multiplication for derivs."<<std::endl;
	  oneVal += x1(i, j) * x2(j, k);
	}
	y(i, k) += oneVal;
      }
    }
  }
  //
  if(x1rowStart > 0) {       // there are some constant first rows in x1, i.e. region A1.
    if(x2colStart < x2colEnd) { // there are some variable cols in x2, i.e. region B2
      // A1 * B2
      // populate up to 3 contributions to y[0:(x1rowStart-1), x1colStart:(x1colEnd-1)]
      // get constant * variable contribution to y[0:(x1rowStart-1), x1colStart:(x1colEnd-1)]
#ifdef VERBOSE_ATOMIC_MATMULT
      cout<<"A1 * B2"<<endl;
      atomic_matmult_diagnostic_message(0, x2rowStart, x1rowStart, x2rowEnd - x2rowStart,
					x2rowStart, x2colStart, x2rowEnd - x2rowStart, x2colEnd - x2colStart,
					0, x2colStart, x1rowStart, x2colEnd - x2colStart);
#endif
      atomic_matmult_internal_const_var(x1, std::array<int, 4>({0, x2rowStart, x1rowStart, x2rowEnd - x2rowStart}) ,
      					x2, std::array<int, 4>({x2rowStart, x2colStart, x2rowEnd - x2rowStart, x2colEnd - x2colStart}),
      					y);

      // MatrixXd_CppAD y_const_var1;
      // y_const_var1.resize(x1rowStart, x2colEnd - x2colStart);
      // atomic_matmult_internal( x1.block(0, x2rowStart, x1rowStart, x2rowEnd - x2rowStart),
      // 			       x2.block(x2rowStart, x2colStart, x2rowEnd - x2rowStart, x2colEnd - x2colStart),
      // 			       y_const_var1, 1, 0);//1, 0);
      // y.block(0, x2colStart, x1rowStart, x2colEnd - x2colStart) += y_const_var1;
    }
  } // End of A1 components

  if(x1rowStart < x1rowEnd) { // there are some variable rows in x1, i.e region B1
    if(x2colStart > 0) {      // there are some constant first cols in x2, i.e. region A2
      // B1 * A2
#ifdef VERBOSE_ATOMIC_MATMULT
      cout<<"B1 * A2"<<endl;
      // 3 cases here
      // set y[x1rowStart:(x1rowEnd-1), 0:(x2colStart-1)]
      atomic_matmult_diagnostic_message(x1rowStart, x1colStart, x1rowEnd-x1rowStart, x1colEnd - x1colStart,
					x1colStart, 0, x1colEnd - x1colStart, x2colStart,
					x1rowStart, 0, x1rowEnd - x1rowStart, x2colStart);
#endif

      atomic_matmult_internal_var_const(x1, std::array<int, 4>({x1rowStart, x1colStart, x1rowEnd-x1rowStart, x1colEnd - x1colStart}),
					x2, std::array<int, 4>({x1colStart, 0, x1colEnd - x1colStart, x2colStart}),
					y);

      // MatrixXd_CppAD y_var_const1;
      // y_var_const1.resize( x1rowEnd - x1rowStart, x2colStart );
      // atomic_matmult_internal( x1.block(x1rowStart, x1colStart, x1rowEnd-x1rowStart, x1colEnd - x1colStart),
      // 			       x2.block(x1colStart, 0, x1colEnd - x1colStart, x2colStart),
      // 			       y_var_const1, 0, 1);
      // y.block(x1rowStart, 0, x1rowEnd - x1rowStart, x2colStart ) += y_var_const1;
    }
    if(x2colStart < x2colEnd) { // B1 * B2, the most complication region to handle
#ifdef VERBOSE_ATOMIC_MATMULT
      cout<<"B1 * B2"<<endl;
#endif
      // There are up to 5 cases here!
      // constant * constant; constant * var OR var * constant; var * var OR constant * constant; var*constant OR constant*var; constant*constant;
      // The constant * constant pieces would have been dealt with above.
      // Start with the var * var region if there is is one, because it is the one region of interest if all elements of x1 and x2 are variables; fill to 0 otherwise
      if(((x1colStart < x2rowEnd) && (x1colEnd > x2rowStart)) || ((x2rowStart < x1colEnd) && (x2rowEnd > x1colStart))) {
#ifdef VERBOSE_ATOMIC_MATMULT
	cout<<"\t segment 1"<<endl;
#endif
	
	int largerStart = x1colStart > x2rowStart ? x1colStart : x2rowStart;
	int smallerEnd = x1colEnd < x2rowEnd ? x1colEnd : x2rowEnd;
	if(smallerEnd <= largerStart) cout<<"Some logic appears to be wrong in the var*var part of nimDerivs matrix mult."<<endl;
#ifdef VERBOSE_ATOMIC_MATMULT
	atomic_matmult_diagnostic_message(x1rowStart, largerStart, x1rowEnd-x1rowStart, smallerEnd - largerStart,
					  largerStart, x2colStart, smallerEnd - largerStart,  x2colEnd - x2colStart,
					  x1rowStart, x2colStart, x1rowEnd - x1rowStart, x2colEnd - x2colStart );
#endif
	MatrixXd_CppAD y_var_var;
	y_var_var.resize(x1rowEnd - x1rowStart, x2colEnd - x2colStart);
	atomic_matmult_internal( x1.block(x1rowStart, largerStart, x1rowEnd-x1rowStart, smallerEnd - largerStart),
				 x2.block(largerStart, x2colStart, smallerEnd - largerStart,  x2colEnd - x2colStart),
				 y_var_var, 0, 0);
	y.block(x1rowStart, x2colStart, x1rowEnd - x1rowStart, x2colEnd - x2colStart) += y_var_var;
      }
      if(x1colStart < x2rowStart) { // there is a var * constant region
#ifdef VERBOSE_ATOMIC_MATMULT
	cout<<"\t segment 2"<<endl;
#endif
	int start = x1colStart;
	int end = x2rowStart < x1colEnd ? x2rowStart : x1colEnd;
#ifdef VERBOSE_ATOMIC_MATMULT
	atomic_matmult_diagnostic_message(x1rowStart, start, x1rowEnd-x1rowStart, end-start,
					  start, x2colStart, end-start,  x2colEnd - x2colStart,
					  x1rowStart, x2colStart, x1rowEnd - x1rowStart, x2colEnd - x2colStart );
#endif
	atomic_matmult_internal_var_const(x1, std::array<int, 4>({x1rowStart , start , (x1rowEnd - x1rowStart) , end-start}),
					  x2, std::array<int, 4>({start , x2colStart , end-start , (x2colEnd - x2colStart)}),
					  y);

	// MatrixXd_CppAD y_var_const1;
	// y_var_const1.resize(  x1rowEnd - x1rowStart, x2colEnd - x2colStart );
	// atomic_matmult_internal( x1.block(x1rowStart , start , (x1rowEnd - x1rowStart) , end-start ),
	// 			 x2.block(start , x2colStart , end-start , (x2colEnd - x2colStart)),
	// 			y_var_const1, 0, 1);
	// y.block(x1rowStart, x2colStart, x1rowEnd - x1rowStart, x2colEnd - x2colStart) 
	// 		       += y_var_const1;
      }
      if( x2rowStart < x1colStart ) { // there is a constant * var region. mutually exclusive with previous case
#ifdef VERBOSE_ATOMIC_MATMULT
	cout<<"\t segment 3"<<endl;
#endif
	int start = x2rowStart; 
	int end = x1colStart < x2rowEnd ? x1colStart : x2rowEnd;
#ifdef VERBOSE_ATOMIC_MATMULT
	atomic_matmult_diagnostic_message(x1rowStart, start, x1rowEnd-x1rowStart, end-start,
					  start, x2colStart, end-start,  x2colEnd - x2colStart,
					  x1rowStart, x2colStart, x1rowEnd - x1rowStart, x2colEnd - x2colStart );
#endif

	atomic_matmult_internal_const_var(x1, std::array<int, 4>({x1rowStart , start , (x1rowEnd - x1rowStart) , end-start}) ,
					  x2, std::array<int, 4>({start , x2colStart , end-start , (x2colEnd - x2colStart)}),
					  y);
	// MatrixXd_CppAD y_const_var1;
	// y_const_var1.resize(  x1rowEnd - x1rowStart, x2colEnd - x2colStart );
	// atomic_matmult_internal( x1.block(x1rowStart , start , (x1rowEnd - x1rowStart) , end-start ),
	// 			 x2.block(start , x2colStart , end-start , (x2colEnd - x2colStart)),
	// 			 y_const_var1, 1, 0);
	// y.block(x1rowStart, x2colStart, x1rowEnd - x1rowStart, x2colEnd - x2colStart) 
	// 		       += y_const_var1;
      }
      if( x1colEnd > x2rowEnd ) { // there is a var * constant region
#ifdef VERBOSE_ATOMIC_MATMULT
	cout<<"\t segment 4"<<endl;
#endif
	
	int end = x1colEnd;
	int start = x1colStart > x2rowEnd ? x1colStart : x2rowEnd;
#ifdef VERBOSE_ATOMIC_MATMULT
	atomic_matmult_diagnostic_message(x1rowStart, start, x1rowEnd-x1rowStart, end-start,
					  start, x2colStart, end-start,  x2colEnd - x2colStart,
					  x1rowStart, x2colStart, x1rowEnd - x1rowStart, x2colEnd - x2colStart );
#endif

      atomic_matmult_internal_var_const(x1, std::array<int, 4>({x1rowStart , start , (x1rowEnd - x1rowStart) , end-start}),
					x2, std::array<int, 4>({start , x2colStart , end-start , (x2colEnd - x2colStart)}),
					y);

	// MatrixXd_CppAD y_var_const1;
	// y_var_const1.resize(  x1rowEnd - x1rowStart, x2colEnd - x2colStart );
	// atomic_matmult_internal( x1.block(x1rowStart , start , (x1rowEnd - x1rowStart) , end-start ),
	// 			 x2.block(start , x2colStart , end-start , (x2colEnd - x2colStart)),
	// 			 y_var_const1, 0, 1);
	// y.block(x1rowStart, x2colStart, x1rowEnd - x1rowStart, x2colEnd - x2colStart) 
	// 		       += y_var_const1;
      }
      if( x2rowEnd > x1colEnd ) { // there is a constant * var region. mutually exclusive with previous case
#ifdef VERBOSE_ATOMIC_MATMULT
	cout<<"\t segment 5"<<endl;
#endif
	int end = x2rowEnd;
	int start = x2rowStart > x1colEnd ? x2rowStart : x1colEnd;
#ifdef VERBOSE_ATOMIC_MATMULT
	atomic_matmult_diagnostic_message(x1rowStart, start, x1rowEnd-x1rowStart, end-start,
					  start, x2colStart, end-start,  x2colEnd - x2colStart,
					  x1rowStart, x2colStart, x1rowEnd - x1rowStart, x2colEnd - x2colStart );
#endif
	atomic_matmult_internal_const_var(x1, std::array<int, 4>({x1rowStart , start , (x1rowEnd - x1rowStart) , end-start}) ,
					  x2, std::array<int, 4>({start , x2colStart , end-start , (x2colEnd - x2colStart)}),
					  y);
	// MatrixXd_CppAD y_const_var1;
	// y_const_var1.resize(  x1rowEnd - x1rowStart, x2colEnd - x2colStart );

	// atomic_matmult_internal( x1.block(x1rowStart , start , (x1rowEnd - x1rowStart) , end-start ),
	// 			 x2.block(start , x2colStart , end-start , (x2colEnd - x2colStart)),
	// 			 y_const_var1, 1, 0);
	// y.block(x1rowStart, x2colStart, x1rowEnd - x1rowStart, x2colEnd - x2colStart) 
	// 		       += y_const_var1;
      }
      // B1 * B2 components are done.
      // 
      // B1 * C2, in this section because C2 only exists if B2 exists
      if(x2colEnd < n3) { // there are some constant final cols in x2.
	// B1 * C2, which is only relevant if there is a B2 region.  Otherwise A2 covers all of x2
	// 3 cases here
#ifdef VERBOSE_ATOMIC_MATMULT
	cout<<"B1 * C2"<<endl;	
	// set y[x1rowStart:(x1rowEnd-1), x2colEnd:(n3-1)]
	atomic_matmult_diagnostic_message(x1rowStart, x1colStart, x1rowEnd-x1rowStart, x1colEnd - x1colStart,
					  x1colStart, x2colEnd, x1colEnd - x1colStart, n3 - x2colEnd,
					  x1rowStart, x2colEnd, x1rowEnd - x1rowStart, n3 - x2colEnd );
#endif

      atomic_matmult_internal_var_const(x1, std::array<int, 4>({x1rowStart, x1colStart, x1rowEnd-x1rowStart, x1colEnd - x1colStart}),
					x2, std::array<int, 4>({x1colStart, x2colEnd, x1colEnd - x1colStart, n3 - x2colEnd}),
					y);
       
	// MatrixXd_CppAD y_var_const1;
	// y_var_const1.resize( x1rowEnd - x1rowStart, n3 - x2colEnd );
	// atomic_matmult_internal( x1.block(x1rowStart, x1colStart, x1rowEnd-x1rowStart, x1colEnd - x1colStart),
	// 			 x2.block(x1colStart, x2colEnd, x1colEnd - x1colStart, n3 - x2colEnd),
	// 			 y_var_const1, 0, 1);
	// y.block(x1rowStart, x2colEnd, x1rowEnd - x1rowStart, n3 - x2colEnd ) += y_var_const1;
      }
      //
      //
      // C1 * B2, in this section because C1 exists only if B1 exists
      if(x1rowEnd < n1 ) {
#ifdef VERBOSE_ATOMIC_MATMULT
	cout<<"C1 * B2"<<endl;	
	// populate up to 3 contributions to y[x1rowEnd:(n1-1), x1colStart:(x1colEnd-1)]
	// get constant * variable contribution to y[x1rowEnd:(n1-1), x1colStart:(x1colEnd-1)]
	atomic_matmult_diagnostic_message(x1rowEnd, x2rowStart, n1-x1rowEnd, x2rowEnd - x2rowStart,
					  x2rowStart, x2colStart, x2rowEnd - x2rowStart, x2colEnd - x2colStart,
					  x1rowEnd, x2colStart, n1-x1rowEnd, x2colEnd - x2colStart );
#endif
	atomic_matmult_internal_const_var(x1, std::array<int, 4>({x1rowEnd, x2rowStart, n1-x1rowEnd, x2rowEnd - x2rowStart}) ,
					  x2, std::array<int, 4>({x2rowStart, x2colStart, x2rowEnd - x2rowStart, x2colEnd - x2colStart}),
					  y);
	// MatrixXd_CppAD y_const_var1;
	// y_const_var1.resize(n1-x1rowEnd, x2colEnd - x2colStart);
	// atomic_matmult_internal( x1.block(x1rowEnd, x2rowStart, n1-x1rowEnd, x2rowEnd - x2rowStart),
	// 			 x2.block(x2rowStart, x2colStart, x2rowEnd - x2rowStart, x2colEnd - x2colStart),
	// 			 y_const_var1, 1, 0);
	// y.block(x1rowEnd, x2colStart, n1-x1rowEnd, x2colEnd - x2colStart) += y_const_var1;
      }
    } // end of B2 within B1
  } // end of B1
}    

MatrixXd_CppAD nimDerivs_matmult(const MatrixXd_CppAD &x1,
				 const MatrixXd_CppAD &x2,
				 bool debug) {
  MatrixXd_CppAD ans;
  
#ifdef VERBOSE_ATOMIC_MATMULT
  cout<<"Entering nimDerivs_matmult"<<endl;
#endif
  atomic_matmult(x1, x2, ans, debug);
#ifdef VERBOSE_ATOMIC_MATMULT
  cout<<"Leaving nimDerivs_matmult"<<endl;
#endif
  return ans;
}
