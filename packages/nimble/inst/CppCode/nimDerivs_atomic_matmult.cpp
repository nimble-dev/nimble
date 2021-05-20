#include <nimble/nimDerivs_atomic_matmult.h>

/* This atomic does not need caching because adjoints in reverse mode do not use Y or Ydot. */

// Use of the "new in place", labeled here "NEW_DYNAMIC"  was a performance experiment.
// The idea was to have member data for every Eigen map to be used and then initialize them
// with new in place.  This would save instantiating new maps every time one is needed.
// The performance impact was negligible, so we generally operate with USE_NEW_DYNAMIC undefined (false, off, unused).
// But actually, the new in place works well to manage that some inputs might be constants.

#ifdef USE_NEW_DYNAMIC
atomic_matmult_class::atomic_matmult_class(const std::string& name) :
  CppAD::atomic_three<double>(name),
  n1_(0),
  X1cat_(unknown), X2cat_(unknown),
  x1_is_constant_(false), x2_is_constant_(false),
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
  // printf("In matmult for_type\n");
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
  int n2 = n/(n1 + n3);


  // HOT-WIRED FOR TESTING
  for(size_t i = 0; i < type_y.size(); ++i)
    type_y[i] = CppAD::variable_enum;
  return true;

  // std::cout << "n1 = "<< n1 <<" n2 = "<< n2 <<" n3 = "<< n3 <<std::endl;
  CppAD::vector<CppAD::ad_type_enum> x1RowTypes(n1);
  CppAD::vector<CppAD::ad_type_enum> x2ColTypes(n3);

  CppAD::ad_type_enum this_row_type;
  CppAD::ad_type_enum item_type;
  for(size_t i = 0; i < n1; ++i) {
    this_row_type = CppAD::constant_enum;
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
    x1RowTypes[i] = this_row_type;
  }
  int n12 = n1*n2;
  CppAD::ad_type_enum this_col_type;
  for(size_t j = 0; j < n3; ++j) {
    this_col_type = CppAD::constant_enum;
    for(size_t i = 0; i < n2; ++i) {
      item_type = type_x[n12 + i + j*n2];
      //std::cout<<"x2("<<i<<","<< j<<") type = "<<item_type<<"\n";
      if(item_type == CppAD::variable_enum) {
	this_col_type = CppAD::variable_enum;
	break;
      } else {
	if(item_type == CppAD::dynamic_enum)
	  this_col_type = CppAD::dynamic_enum;
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
  // printf("In matmult reverse_depend\n");
  // Here we put true in an element of depend_x if there are true depend_y elements that are functions of the depend_x element
  int n = static_cast<double>(type_x.size());
  int m = depend_y.size();
  int n1 = get_n1();// parameter_x[0];
  // std::cout<<"n1 = "<<n1<<std::endl;
  //    if(type_x[0] == CppAD::constant_enum) std::cout<<"constant"<<std::endl;
  //    if(type_x[0] == CppAD::dynamic_enum) std::cout<<"dynamic"<<std::endl;
  //    if(type_x[0] == CppAD::variable_enum) std::cout<<"variable"<<std::endl;
  int n3 = m/n1;
  int n2 = n/(n1 + n3);
  int n12 = n1*n2;

  // HOT-WIRED FOR TESTING
  for(size_t i = 0; i < depend_x.size(); ++i)
    depend_x[i] = true;
  return true;
  
  CppAD::vector<bool> depend_x1Row(n1);
  CppAD::vector<bool> depend_x2Col(n3);
  for(size_t i = 0; i < n1; ++i) depend_x1Row[i] = false;
  for(size_t j = 0; j < n3; ++j) depend_x2Col[j] = false;
  bool this_depend;
  for(size_t i = 0; i < n1; ++i) {
    for(size_t j = 0; j < n3; ++j) {
      this_depend = depend_y[i + j*n1];
      depend_x1Row[i] |= this_depend;
      depend_x2Col[j] |= this_depend;
    }
  }
  for(size_t i = 0; i < n1; ++i) {
    this_depend = depend_x1Row[i];
    for(size_t j = 0; j < n2; ++j) {
      depend_x[i + j*n1] = this_depend;
    }
  }
  for(size_t j = 0; j < n3; ++j) {
    this_depend = depend_x2Col[j];
    for(size_t i = 0; i < n2; ++i) {
      depend_x[n12 + i + j*n2] = this_depend;
    }
  }
  return true;
}

template<class X1c, class X2c, class Yc >
void matmult_internal_respecting_upper_lower(const X1c &X1, const X2c &X2, Yc &Y,
					     const matrix_category &X1cat,
					     const matrix_category &X2cat) {
  if(X1cat == upper_diagonal) {
    if(X2cat == upper_diagonal) {
      MatrixXd temp = X2.template triangularView<Eigen::Upper>(); // This MatrixXd needs to be extracted from the template types in general
      Y = X1.template triangularView<Eigen::Upper>() * temp;
    } else if(X2cat == lower_diagonal) {
      MatrixXd temp = X2.template triangularView<Eigen::Lower>();
      Y = X1.template triangularView<Eigen::Upper>() * temp;
    } else {
      Y = X1.template triangularView<Eigen::Upper>() * X2;
    }      
  } else if(X1cat == lower_diagonal) {
    if(X2cat == upper_diagonal) {
      MatrixXd temp = X2.template triangularView<Eigen::Upper>();
      Y = X1.template triangularView<Eigen::Lower>() * temp;
    } else if(X2cat == lower_diagonal) {
      MatrixXd temp = X2.template triangularView<Eigen::Lower>();
      Y = X1.template triangularView<Eigen::Lower>() * temp;
    } else {
      Y = X1.template triangularView<Eigen::Lower>() * X2;
    }
  } else {
    if(X2cat == upper_diagonal) {
      Y = X1
	* X2.template triangularView<Eigen::Upper>();
    } else if(X1cat == lower_diagonal) {
      Y = X1
	* X2.template triangularView<Eigen::Lower>();
    } else {
      cout<<"general case"<<endl;
      Y = X1 * X2; // general case is here
    }
  }
}

template<class X1c, class X2c, class Yc >
void matmult_internal_respecting_upper_lower_add(const X1c &X1, const X2c &X2, Yc &Y,
						 const matrix_category &X1cat,
						 const matrix_category &X2cat) {
  if(X1cat == upper_diagonal) {
    if(X2cat == upper_diagonal) {
      MatrixXd temp = X2.template triangularView<Eigen::Upper>();
      Y += X1.template triangularView<Eigen::Upper>() * temp;
    } else if(X2cat == lower_diagonal) {
      MatrixXd temp = X2.template triangularView<Eigen::Lower>();
      Y += X1.template triangularView<Eigen::Upper>() * temp;
    } else {
      Y += X1.template triangularView<Eigen::Upper>() * X2;
    }      
  } else if(X1cat == lower_diagonal) {
    if(X2cat == upper_diagonal) {
      MatrixXd temp = X2.template triangularView<Eigen::Upper>();
      Y += X1.template triangularView<Eigen::Lower>() * temp;
    } else if(X2cat == lower_diagonal) {
      MatrixXd temp = X2.template triangularView<Eigen::Lower>();
      Y += X1.template triangularView<Eigen::Lower>() * temp;
    } else {
      Y += X1.template triangularView<Eigen::Lower>() * X2;
    }
  } else {
    if(X2cat == upper_diagonal) {
      Y += X1
	* X2.template triangularView<Eigen::Upper>();
    } else if(X1cat == lower_diagonal) {
      Y += X1
	* X2.template triangularView<Eigen::Lower>();
    } else {
      Y += X1 * X2; // general case is here
    }
  }
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
  printf("In matmult forward\n");
  int nrow = order_up + 1;
  // std::cout<<"nrow = "<<nrow<<std::endl;
  int n = static_cast<double>(taylor_x.size()/nrow); // size of input
  int m = taylor_y.size() / nrow;                    // size of output
  int n1 = get_n1();                                 // rows of X1
  int n3 = m/n1;                                     // cols of X2
  int n2;
  if(X1constant()) n2 = n / n3;
  else if(X2constant()) n2 = n / n1; 
  else n2 = n/(n1 + n3);

  std::cout << "n1 = "<< n1 <<" n2 = "<< n2 <<" n3 = "<< n3 <<std::endl;
  if(X1constant()) {
    new (&X1mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(get_X_stored_ptr(),
								    n1, n2, EigStrDyn(nrow*n1, nrow) );
  } else {
    new (&X1mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(&taylor_x[0],
								    n1, n2, EigStrDyn(nrow*n1, nrow) );
  }
  cout<<X1mapC<<endl;
  //   EigenTemplateTypes<double>::typeEigenConstMapStrd X1map(&taylor_x[0 + 1*nrow],
  //							   n1, n2, EigStrDyn(nrow*n1, nrow) );
  if(X2constant()) {
    new (&X2mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(get_X_stored_ptr(),
								    n2, n3, EigStrDyn(nrow*n2, nrow ) );
  } else {
    if(X1constant())
      new (&X2mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(&taylor_x[0],
								      n2, n3, EigStrDyn(nrow*n2, nrow ) );
    else
      new (&X2mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(&taylor_x[0 + (n1*n2)*nrow],
								      n2, n3, EigStrDyn(nrow*n2, nrow ) );
  }
  //  cout"X2\n"<<X2mapC<<endl;

  if(order_low <= 0 & order_up >= 0) { // value
    // We could compile different cases depending on need for strides or not.
    // We could use const maps.
    
    new (&Ymap) EigenTemplateTypes<double>::typeEigenMapStrd(&taylor_y[0],
							     n1, n3, EigStrDyn(nrow*n1, nrow ) );
    matmult_internal_respecting_upper_lower(X1mapC, X2mapC, Ymap, X1cat(), X2cat());
    //    Ymap = X1mapC * X2mapC;
    // for(int i = 0; i < X1mapC.rows(); ++i) {
    //   for(int j = 0; j < X1mapC.cols(); ++j) {
    // 	std::cout<<X1mapC(i,j)<<" ";
    //   }
    //   std::cout<<std::endl;
    // }      
  }
  if(order_low <= 1 & order_up >= 1) { // forward 1
    // printf("In forward >1\n");
    new (&dY_map) EigenTemplateTypes<double>::typeEigenMapStrd(&taylor_y[1],
							       n1, n3, EigStrDyn(nrow*n1, nrow ) );
    if(!X1constant()) {
      new (&dX1mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(&taylor_x[1],
								       n1, n2, EigStrDyn(nrow*n1, nrow) );
      
      matmult_internal_respecting_upper_lower(dX1mapC, X2mapC, dY_map, X1cat(), X2cat());
      //dY_map = dX1mapC * X2mapC;
    }
    if(!X2constant()) {
      if(!X1constant()) {
	new (&dX2mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(&taylor_x[1 + (n1*n2)*nrow],
									 n2, n3, EigStrDyn(nrow*n2, nrow ) );
	matmult_internal_respecting_upper_lower_add(X1mapC, dX2mapC, dY_map, X1cat(), X2cat());
	// dY_map += X1mapC * dX2mapC;
      } else {
	new (&dX2mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(&taylor_x[1 ],
									 n2, n3, EigStrDyn(nrow*n2, nrow ) );
	matmult_internal_respecting_upper_lower(X1mapC, dX2mapC, dY_map, X1cat(), X2cat());
	// dY_map = X1mapC * dX2mapC;
      }
    }
    // X1constant() and X2constant() should not both be true, or this atomic wouldn't exist.
    // dY_map = dX1mapC * X2mapC + X1mapC * dX2mapC;
  }
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
  //  printf("In matmult meta-forward\n");
  int nrow = order_up + 1;
  // std::cout<<"nrow = "<<nrow<<std::endl;
  int n = static_cast<double>(taylor_x.size()/nrow);
  int m = taylor_y.size() / nrow;
  int n1 = get_n1(); //CppAD::Value(taylor_x[0]);
  int n3 = m/n1;
  int n2 = n/(n1 + n3);
  new (&mX1mapC) EigenTemplateTypes<CppAD::AD<double> >::typeEigenConstMapStrd(&taylor_x[0],
									       n1, n2, EigStrDyn(nrow*n1, nrow) );
  new (&mX2mapC) EigenTemplateTypes<CppAD::AD<double> >::typeEigenConstMapStrd(&taylor_x[0 + (n1*n2)*nrow],
									       n2, n3, EigStrDyn(nrow*n2, nrow ) );
  if(order_low <= 0 & order_up >= 0) { // value
    new (&mYmap) EigenTemplateTypes<CppAD::AD<double> >::typeEigenMapStrd(&taylor_y[0],
									  n1, n3, EigStrDyn(nrow*n1, nrow ) );
    atomic_matmult(mX1mapC, mX2mapC, mTerm1); //Ymap = X1map * X2map;
    mYmap = mTerm1;
  }
  if(order_low <= 1 & order_up >= 1) { // forward 1
    // printf("In forward >1\n");
    new (&mdX1mapC) EigenTemplateTypes<CppAD::AD<double>>::typeEigenConstMapStrd(&taylor_x[1],
										 n1, n2, EigStrDyn(nrow*n1, nrow) );
    new (&mdX2mapC) EigenTemplateTypes<CppAD::AD<double>>::typeEigenConstMapStrd(&taylor_x[1 + (n1*n2)*nrow],
										 n2, n3, EigStrDyn(nrow*n2, nrow ) );
    new (&mdY_map) EigenTemplateTypes<CppAD::AD<double>>::typeEigenMapStrd(&taylor_y[1],
									   n1, n3, EigStrDyn(nrow*n1, nrow ) );

    atomic_matmult(mdX1mapC, mX2mapC, mTerm1);
    atomic_matmult(mX1mapC, mdX2mapC, mTerm2);
    mdY_map = mTerm1 + mTerm2;
  }
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
  // printf("In matmult reverse\n");
  int nrow = order_up + 1;
  // std::cout<<"nrow = "<<nrow<<std::endl;
  int n = static_cast<double>(taylor_x.size()/nrow);
  int m = taylor_y.size() / nrow;
  int n1 = get_n1();
  int n3 = m/n1;
  int n2;
  if(X1constant()) n2 = n / n3;
  else if(X2constant()) n2 = n / n1; 
  else n2 = n/(n1 + n3);
  
  if(X1constant()) {
    new (&X1mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(get_X_stored_ptr(),
								    n1, n2, EigStrDyn(nrow*n1, nrow) );
  } else {
    new (&X1mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(&taylor_x[0],
								    n1, n2, EigStrDyn(nrow*n1, nrow) );
    new (&X1adjoint_map) EigenTemplateTypes<double>::typeEigenMapStrd(&partial_x[0],
								      n1, n2, EigStrDyn(nrow*n1, nrow) );
  }
  if(X2constant()) {
    new (&X2mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(get_X_stored_ptr(),
								    n2, n3, EigStrDyn(nrow*n2, nrow ) );
  } else { // X2 is not constant
    if(X1constant()) {
      new (&X2mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(&taylor_x[0],
								      n2, n3, EigStrDyn(nrow*n2, nrow ) );
      new (&X2adjoint_map) EigenTemplateTypes<double>::typeEigenMapStrd(&partial_x[0],
									n2, n3, EigStrDyn(nrow*n2, nrow ) );
      
    } else { // neither is constant
      new (&X2mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(&taylor_x[0 + (n1*n2)*nrow],
								      n2, n3, EigStrDyn(nrow*n2, nrow ) );

      new (&X2adjoint_map) EigenTemplateTypes<double>::typeEigenMapStrd(&partial_x[0 + (n1*n2)*nrow],
									n2, n3, EigStrDyn(nrow*n2, nrow ) );
    }
  }
  new (&Yadjoint_mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(&partial_y[0],
									 n1, n3, EigStrDyn(nrow*n1, nrow ) );
  if(order_up >= 0) {
    // reverse 1
    //    partial_x[0] = 0;
    if(!X1constant()) {
      matmult_internal_respecting_upper_lower(Yadjoint_mapC, X2mapC.transpose(), X1adjoint_map,
					      unknown, transpose_X2cat());
      //    X1adjoint_map = Yadjoint_mapC * X2mapC.transpose();
    }
    if(!X2constant()) {
      matmult_internal_respecting_upper_lower( X1mapC.transpose(), Yadjoint_mapC, X2adjoint_map,
					       transpose_X1cat(), unknown );
      // X2adjoint_map = X1mapC.transpose() * Yadjoint_mapC;
    }
  }
  if(order_up >= 1) {
    // reverse 2
    new (&Ydot_adjoint_mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(&partial_y[1],
									       n1, n3, EigStrDyn(nrow*n1, nrow ) );

    if(!X1constant()) {
      new (&X1dot_mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(&taylor_x[1],
									  n1, n2, EigStrDyn(nrow*n1, nrow) );
      new (&X1dot_adjoint_map) EigenTemplateTypes<double>::typeEigenMapStrd(&partial_x[1],
									    n1, n2, EigStrDyn(nrow*n1, nrow) );
    }
    if(!X2constant()) {
      if(!X1constant()) {
	new (&X2dot_mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(&taylor_x[1 + (n1*n2)*nrow],
									    n2, n3, EigStrDyn(nrow*n2, nrow ) );
	new (&X2dot_adjoint_map) EigenTemplateTypes<double>::typeEigenMapStrd(&partial_x[1 + (n1*n2)*nrow],
									      n2, n3, EigStrDyn(nrow*n2, nrow ));
      } else {
	new (&X2dot_mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(&taylor_x[1],
									    n2, n3, EigStrDyn(nrow*n2, nrow ) );
	new (&X2dot_adjoint_map) EigenTemplateTypes<double>::typeEigenMapStrd(&partial_x[1],
									      n2, n3, EigStrDyn(nrow*n2, nrow ));
      }
    }

    if((!X1constant()) && (!X2constant())) {
      matmult_internal_respecting_upper_lower_add(Ydot_adjoint_mapC, X2dot_mapC.transpose(), X1adjoint_map,
						  unknown, transpose_X2cat());
      //   X1adjoint_map += Ydot_adjoint_mapC * X2dot_mapC.transpose();
      matmult_internal_respecting_upper_lower_add(X1dot_mapC.transpose(), Ydot_adjoint_mapC, X2adjoint_map,
						   transpose_X1cat(), unknown);
      //   X2adjoint_map += X1dot_mapC.transpose() * Ydot_adjoint_mapC;
    }
    if(!X1constant()) {
      matmult_internal_respecting_upper_lower(Ydot_adjoint_mapC, X2mapC.transpose(), X1dot_adjoint_map,
					      unknown, transpose_X2cat());
      //      X1dot_adjoint_map = Ydot_adjoint_mapC * X2mapC.transpose();
    }
    if(!X2constant()) {
      matmult_internal_respecting_upper_lower_add(X1mapC.transpose(), Ydot_adjoint_mapC, X2dot_adjoint_map,
						  transpose_X1cat(), unknown);
      //      X2dot_adjoint_map = X1mapC.transpose() * Ydot_adjoint_mapC;
    }
  }
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
  // printf("In matmult reverse\n");
  int nrow = order_up + 1;
  // std::cout<<"nrow = "<<nrow<<std::endl;
  int n = static_cast<double>(taylor_x.size()/nrow);
  int m = taylor_y.size() / nrow;
  int n1 = get_n1(); //CppAD::Value(taylor_x[0]);
  int n3 = m/n1;
  int n2 = n/(n1 + n3);
  new (&mX1mapC) EigenTemplateTypes<CppAD::AD<double> >::typeEigenConstMapStrd(&taylor_x[0],
									       n1, n2, EigStrDyn(nrow*n1, nrow) );
  new (&mX2mapC) EigenTemplateTypes<CppAD::AD<double> >::typeEigenConstMapStrd(&taylor_x[0 + (n1*n2)*nrow],
									       n2, n3, EigStrDyn(nrow*n2, nrow ) );
  new (&mYadjoint_mapC) EigenTemplateTypes<CppAD::AD<double> >::typeEigenConstMapStrd(&partial_y[0],
										      n1, n3, EigStrDyn(nrow*n1, nrow ) );
  new (&mX1adjoint_map) EigenTemplateTypes<CppAD::AD<double> >::typeEigenMapStrd(&partial_x[0],
										 n1, n2, EigStrDyn(nrow*n1, nrow) );
  new (&mX2adjoint_map) EigenTemplateTypes<CppAD::AD<double> >::typeEigenMapStrd(&partial_x[0 + (n1*n2)*nrow],
										 n2, n3, EigStrDyn(nrow*n2, nrow ) );
  if(order_up >= 0) {
    // reverse 1
    partial_x[0] = 0;
    atomic_matmult(mYadjoint_mapC, mX2mapC.transpose(), mTerm1);
    mX1adjoint_map = mTerm1;
    atomic_matmult(mX1mapC.transpose(), mYadjoint_mapC,  mTerm2);
    mX2adjoint_map = mTerm2;
  }
  if(order_up >= 1) {
    // reverse 2
    new (&mX1dot_mapC) EigenTemplateTypes<CppAD::AD<double> >::typeEigenConstMapStrd(&taylor_x[1],
										     n1, n2, EigStrDyn(nrow*n1, nrow) );
    new (&mX2dot_mapC) EigenTemplateTypes<CppAD::AD<double> >::typeEigenConstMapStrd(&taylor_x[1 + (n1*n2)*nrow],
										     n2, n3, EigStrDyn(nrow*n2, nrow ) );
    new (&mYdot_adjoint_mapC) EigenTemplateTypes<CppAD::AD<double> >::typeEigenConstMapStrd(&partial_y[1],
											    n1, n3, EigStrDyn(nrow*n1, nrow ) );
    new (&mX1dot_adjoint_map) EigenTemplateTypes<CppAD::AD<double> >::typeEigenMapStrd(&partial_x[1],
										       n1, n2, EigStrDyn(nrow*n1, nrow) );
    new (&mX2dot_adjoint_map) EigenTemplateTypes<CppAD::AD<double> >::typeEigenMapStrd(&partial_x[1 + (n1*n2)*nrow],
										       n2, n3, EigStrDyn(nrow*n2, nrow ) );

    atomic_matmult(mYdot_adjoint_mapC, mX2dot_mapC.transpose(), mTerm1);
    atomic_matmult(mX1dot_mapC.transpose() , mYdot_adjoint_mapC, mTerm2);
    mX1adjoint_map += mTerm1;
    mX2adjoint_map += mTerm2;

    atomic_matmult(mYdot_adjoint_mapC, mX2mapC.transpose(), mTerm1);
    mX1dot_adjoint_map = mTerm1;
    atomic_matmult(mX1mapC.transpose(), mYdot_adjoint_mapC,  mTerm2);
    mX2dot_adjoint_map = mTerm2;
  }
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

#else // USE_
// Use local Eigen maps
// This is now deprecated because of use of the matrix_category system above

atomic_matmult_class::atomic_matmult_class(const std::string& name) :
  CppAD::atomic_three<double>(name),
  n1_(0)
{};
bool atomic_matmult_class::for_type(
				    const CppAD::vector<double>&               parameter_x ,
				    const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				    CppAD::vector<CppAD::ad_type_enum>&        type_y      )
{
  // printf("In matmult for_type\n");
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
  int n2 = n/(n1 + n3);
  // std::cout << "n1 = "<< n1 <<" n2 = "<< n2 <<" n3 = "<< n3 <<std::endl;
  CppAD::vector<CppAD::ad_type_enum> x1RowTypes(n1);
  CppAD::vector<CppAD::ad_type_enum> x2ColTypes(n3);

  CppAD::ad_type_enum this_row_type;
  CppAD::ad_type_enum item_type;
  for(size_t i = 0; i < n1; ++i) {
    this_row_type = CppAD::constant_enum;
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
    x1RowTypes[i] = this_row_type;
  }
  int n12 = n1*n2;
  CppAD::ad_type_enum this_col_type;
  for(size_t j = 0; j < n3; ++j) {
    this_col_type = CppAD::constant_enum;
    for(size_t i = 0; i < n2; ++i) {
      item_type = type_x[n12 + i + j*n2];
      //std::cout<<"x2("<<i<<","<< j<<") type = "<<item_type<<"\n";
      if(item_type == CppAD::variable_enum) {
	this_col_type = CppAD::variable_enum;
	break;
      } else {
	if(item_type == CppAD::dynamic_enum)
	  this_col_type = CppAD::dynamic_enum;
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
  // printf("In matmult reverse_depend\n");
  // Here we put true in an element of depend_x if there are true depend_y elements that are functions of the depend_x element
  int n = static_cast<double>(type_x.size());
  int m = depend_y.size();
  int n1 = get_n1();// parameter_x[0];
  // std::cout<<"n1 = "<<n1<<std::endl;
  //    if(type_x[0] == CppAD::constant_enum) std::cout<<"constant"<<std::endl;
  //    if(type_x[0] == CppAD::dynamic_enum) std::cout<<"dynamic"<<std::endl;
  //    if(type_x[0] == CppAD::variable_enum) std::cout<<"variable"<<std::endl;
  int n3 = m/n1;
  int n2 = n/(n1 + n3);

  CppAD::vector<bool> depend_x1Row(n1);
  CppAD::vector<bool> depend_x2Col(n3);
  for(size_t i = 0; i < n1; ++i) depend_x1Row[i] = false;
  for(size_t j = 0; j < n3; ++j) depend_x2Col[j] = false;
  bool this_depend;
  for(size_t i = 0; i < n1; ++i) {
    for(size_t j = 0; j < n3; ++j) {
      this_depend = depend_y[i + j*n1];
      depend_x1Row[i] |= this_depend;
      depend_x2Col[j] |= this_depend;
    }
  }
  for(size_t i = 0; i < n1; ++i) {
    this_depend = depend_x1Row[i];
    for(size_t j = 0; j < n2; ++j) {
      depend_x[i + j*n1] = this_depend;
    }
  }
  int n12 = n1*n2;
  for(size_t j = 0; j < n3; ++j) {
    this_depend = depend_x2Col[j];
    for(size_t i = 0; i < n2; ++i) {
      depend_x[n12 + i + j*n2] = this_depend;
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
  //  printf("In matmult forward\n");
  int nrow = order_up + 1;
  // std::cout<<"nrow = "<<nrow<<std::endl;
  int n = static_cast<double>(taylor_x.size()/nrow);
  int m = taylor_y.size() / nrow;
  int n1 = get_n1(); //taylor_x[0];
  int n3 = m/n1;
  int n2 = n/(n1 + n3);

  // std::cout << "n1 = "<< n1 <<" n2 = "<< n2 <<" n3 = "<< n3 <<std::endl;
   
  EigenTemplateTypes<double>::typeEigenConstMapStrd X1mapC(&taylor_x[0],
							   n1, n2, EigStrDyn(nrow*n1, nrow) );
  EigenTemplateTypes<double>::typeEigenConstMapStrd X2mapC(&taylor_x[0 + (n1*n2)*nrow],
							   n2, n3, EigStrDyn(nrow*n2, nrow ) );
  if(order_low <= 0 & order_up >= 0) { // value
    // We could compile different cases depending on need for strides or not.
    // We could use const maps.
    EigenTemplateTypes<double>::typeEigenMapStrd Ymap(&taylor_y[0],
						      n1, n3, EigStrDyn(nrow*n1, nrow ) );
    Ymap = X1mapC * X2mapC;
    // for(int i = 0; i < X1mapC.rows(); ++i) {
    //   for(int j = 0; j < X1mapC.cols(); ++j) {
    // 	std::cout<<X1mapC(i,j)<<" ";
    //   }
    //   std::cout<<std::endl;
    // }      
  }
  if(order_low <= 1 & order_up >= 1) { // forward 1
    // printf("In forward >1\n");
    EigenTemplateTypes<double>::typeEigenConstMapStrd dX1mapC(&taylor_x[1],
							      n1, n2, EigStrDyn(nrow*n1, nrow) );
    EigenTemplateTypes<double>::typeEigenConstMapStrd dX2mapC(&taylor_x[1 + (n1*n2)*nrow],
							      n2, n3, EigStrDyn(nrow*n2, nrow ) );
    EigenTemplateTypes<double>::typeEigenMapStrd dY_map(&taylor_y[1],
							n1, n3, EigStrDyn(nrow*n1, nrow ) );
    dY_map = dX1mapC * X2mapC + X1mapC * dX2mapC;
  }
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
  //  printf("In matmult meta-forward\n");
  int nrow = order_up + 1;
  // std::cout<<"nrow = "<<nrow<<std::endl;
  int n = static_cast<double>(taylor_x.size()/nrow);
  int m = taylor_y.size() / nrow;
  int n1 = get_n1(); //CppAD::Value(taylor_x[0]);
  int n3 = m/n1;
  int n2 = n/(n1 + n3);
  EigenTemplateTypes<CppAD::AD<double>>::typeMatrixXd mTerm1, mTerm2;
  EigenTemplateTypes<CppAD::AD<double> >::typeEigenConstMapStrd mX1mapC(&taylor_x[0],
									n1, n2, EigStrDyn(nrow*n1, nrow) );
  EigenTemplateTypes<CppAD::AD<double> >::typeEigenConstMapStrd mX2mapC(&taylor_x[0 + (n1*n2)*nrow],
									n2, n3, EigStrDyn(nrow*n2, nrow ) );
  if(order_low <= 0 & order_up >= 0) { // value
    EigenTemplateTypes<CppAD::AD<double> >::typeEigenMapStrd mYmap(&taylor_y[0],
								   n1, n3, EigStrDyn(nrow*n1, nrow ) );
    atomic_matmult(mX1mapC, mX2mapC, mTerm1); //Ymap = X1map * X2map;
    mYmap = mTerm1;
  }
  if(order_low <= 1 & order_up >= 1) { // forward 1
    // printf("In forward >1\n");
    EigenTemplateTypes<CppAD::AD<double>>::typeEigenConstMapStrd mdX1mapC(&taylor_x[1],
									  n1, n2, EigStrDyn(nrow*n1, nrow) );
    EigenTemplateTypes<CppAD::AD<double>>::typeEigenConstMapStrd mdX2mapC(&taylor_x[1 + (n1*n2)*nrow],
									  n2, n3, EigStrDyn(nrow*n2, nrow ) );
    EigenTemplateTypes<CppAD::AD<double>>::typeEigenMapStrd mdY_map(&taylor_y[1],
								    n1, n3, EigStrDyn(nrow*n1, nrow ) );

    atomic_matmult(mdX1mapC, mX2mapC, mTerm1);
    atomic_matmult(mX1mapC, mdX2mapC, mTerm2);
    mdY_map = mTerm1 + mTerm2;
  }
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
  // printf("In matmult reverse\n");
  int nrow = order_up + 1;
  // std::cout<<"nrow = "<<nrow<<std::endl;
  int n = static_cast<double>(taylor_x.size()/nrow);
  int m = taylor_y.size() / nrow;
  int n1 = get_n1();//taylor_x[0];
  int n3 = m/n1;
  int n2 = n/(n1 + n3);
  EigenTemplateTypes<double>::typeEigenConstMapStrd X1mapC(&taylor_x[0],
							   n1, n2, EigStrDyn(nrow*n1, nrow) );
  EigenTemplateTypes<double>::typeEigenConstMapStrd X2mapC(&taylor_x[0 + (n1*n2)*nrow],
							   n2, n3, EigStrDyn(nrow*n2, nrow ) );
  EigenTemplateTypes<double>::typeEigenConstMapStrd Yadjoint_mapC(&partial_y[0],
								  n1, n3, EigStrDyn(nrow*n1, nrow ) );
  EigenTemplateTypes<double>::typeEigenMapStrd X1adjoint_map(&partial_x[0],
							     n1, n2, EigStrDyn(nrow*n1, nrow) );
  EigenTemplateTypes<double>::typeEigenMapStrd X2adjoint_map(&partial_x[0 + (n1*n2)*nrow],
							     n2, n3, EigStrDyn(nrow*n2, nrow ) );
  if(order_up >= 0) {
    // reverse 1
    partial_x[0] = 0;
    X1adjoint_map = Yadjoint_mapC *  X2mapC.transpose();
    X2adjoint_map = X1mapC.transpose() * Yadjoint_mapC;
  }
  if(order_up >= 1) {
    // reverse 2
    EigenTemplateTypes<double>::typeEigenConstMapStrd X1dot_mapC(&taylor_x[1],
								 n1, n2, EigStrDyn(nrow*n1, nrow) );
    EigenTemplateTypes<double>::typeEigenConstMapStrd X2dot_mapC(&taylor_x[1 + (n1*n2)*nrow],
								 n2, n3, EigStrDyn(nrow*n2, nrow ) );
    EigenTemplateTypes<double>::typeEigenConstMapStrd Ydot_adjoint_mapC(&partial_y[1],
									n1, n3, EigStrDyn(nrow*n1, nrow ) );
    EigenTemplateTypes<double>::typeEigenMapStrd X1dot_adjoint_map(&partial_x[1],
								   n1, n2, EigStrDyn(nrow*n1, nrow) );
    EigenTemplateTypes<double>::typeEigenMapStrd X2dot_adjoint_map(&partial_x[1 + (n1*n2)*nrow],
								   n2, n3, EigStrDyn(nrow*n2, nrow ) );
    X1adjoint_map += Ydot_adjoint_mapC * X2dot_mapC.transpose();
    X2adjoint_map += X1dot_mapC.transpose() * Ydot_adjoint_mapC;
    X1dot_adjoint_map = Ydot_adjoint_mapC * X2mapC.transpose();
    X2dot_adjoint_map = X1mapC.transpose() * Ydot_adjoint_mapC;
  }
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
  // printf("In matmult reverse\n");
  int nrow = order_up + 1;
  // std::cout<<"nrow = "<<nrow<<std::endl;
  int n = static_cast<double>(taylor_x.size()/nrow);
  int m = taylor_y.size() / nrow;
  int n1 = get_n1(); //CppAD::Value(taylor_x[0]);
  int n3 = m/n1;
  int n2 = n/(n1 + n3);
  EigenTemplateTypes<CppAD::AD<double>>::typeMatrixXd mTerm1, mTerm2;
  EigenTemplateTypes<CppAD::AD<double> >::typeEigenConstMapStrd mX1mapC(&taylor_x[0],
									n1, n2, EigStrDyn(nrow*n1, nrow) );
  EigenTemplateTypes<CppAD::AD<double> >::typeEigenConstMapStrd mX2mapC(&taylor_x[0 + (n1*n2)*nrow],
									n2, n3, EigStrDyn(nrow*n2, nrow ) );
  EigenTemplateTypes<CppAD::AD<double> >::typeEigenConstMapStrd mYadjoint_mapC(&partial_y[0],
									       n1, n3, EigStrDyn(nrow*n1, nrow ) );
  EigenTemplateTypes<CppAD::AD<double> >::typeEigenMapStrd mX1adjoint_map(&partial_x[0],
									  n1, n2, EigStrDyn(nrow*n1, nrow) );
  EigenTemplateTypes<CppAD::AD<double> >::typeEigenMapStrd mX2adjoint_map(&partial_x[0 + (n1*n2)*nrow],
									  n2, n3, EigStrDyn(nrow*n2, nrow ) );
  if(order_up >= 0) {
    // reverse 1
    partial_x[0] = 0;
    atomic_matmult(mYadjoint_mapC, mX2mapC.transpose(), mTerm1);
    mX1adjoint_map = mTerm1;
    atomic_matmult(mX1mapC.transpose(), mYadjoint_mapC,  mTerm2);
    mX2adjoint_map = mTerm2;
  }
  if(order_up >= 1) {
    // reverse 2
    EigenTemplateTypes<CppAD::AD<double> >::typeEigenConstMapStrd mX1dot_mapC(&taylor_x[1],
									      n1, n2, EigStrDyn(nrow*n1, nrow) );
    EigenTemplateTypes<CppAD::AD<double> >::typeEigenConstMapStrd mX2dot_mapC(&taylor_x[1 + (n1*n2)*nrow],
									      n2, n3, EigStrDyn(nrow*n2, nrow ) );
    EigenTemplateTypes<CppAD::AD<double> >::typeEigenConstMapStrd mYdot_adjoint_mapC(&partial_y[1],
										     n1, n3, EigStrDyn(nrow*n1, nrow ) );
    EigenTemplateTypes<CppAD::AD<double> >::typeEigenMapStrd mX1dot_adjoint_map(&partial_x[1],
										n1, n2, EigStrDyn(nrow*n1, nrow) );
    EigenTemplateTypes<CppAD::AD<double> >::typeEigenMapStrd mX2dot_adjoint_map(&partial_x[1 + (n1*n2)*nrow],
										n2, n3, EigStrDyn(nrow*n2, nrow ) );

    atomic_matmult(mYdot_adjoint_mapC, mX2dot_mapC.transpose(), mTerm1);
    atomic_matmult(mX1dot_mapC.transpose() , mYdot_adjoint_mapC, mTerm2);
    mX1adjoint_map += mTerm1;
    mX2adjoint_map += mTerm2;

    atomic_matmult(mYdot_adjoint_mapC, mX2mapC.transpose(), mTerm1);
    mX1dot_adjoint_map = mTerm1;
    atomic_matmult(mX1mapC.transpose(), mYdot_adjoint_mapC,  mTerm2);
    mX2dot_adjoint_map = mTerm2;
  }
  if(order_up >= 2) {
    printf("Unsupported reverse order requested\n");
    return false;
  }
  return true;
}

#endif

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
    for(int j = 0; j <= i; ++j) {
      if(! (is_upper_diag &= CppAD::IdenticalZero(x(i, j) ) ) ) break;
    }
    if(!is_upper_diag) break;
  }
  if(is_upper_diag) return upper_diagonal;
  
  bool is_lower_diag(true);
  for(int i = 0; i < nRow; ++i) {
    for(int j = i; j < nCol; ++j) {
      if(! (is_lower_diag &= CppAD::IdenticalZero(x(i, j) ) ) ) break;
    }
    if(!is_lower_diag) break;
  }
  if(is_lower_diag) return lower_diagonal;

  return square_full;
}

bool delineate_constant_region(const MatrixXd_CppAD &x,
			       int &rowStart, int &rowEnd,
			       int &colStart, int &colEnd) {
  // All "end" values are C-style, that is one past the last valid index.
  int nRow = x.rows();
  int nCol = x.cols();
  rowStart = nRow; // first row that is not all constant
  rowEnd = 0;      // one past last row that is not all constant
  colStart = nCol; // first col that is not all constant
  colEnd = 0;      // one past last col that is not all constant

  for(int i = 0; i < nRow; ++i) {
    bool this_row_constant(true);
    for(int j = 0; j < nCol; ++j) {
      if(!(this_row_constant &= CppAD::Constant(x(i, j )))) break;
    }
    if(!this_row_constant) {
      rowStart = i;
      break;
    }
  }

  for(int i = nRow-1; i >= 0; --i) {
    bool this_row_constant(true);
    for(int j = 0; j < nCol; ++j) {
      if(!(this_row_constant &= CppAD::Constant(x(i, j )))) break;
    }
    if(!this_row_constant) {
      rowEnd = i + 1;
      break;
    }
  }

  for(int j = 0; j < nCol; ++j) {
    bool this_col_constant(true);
    for(int i = 0; i < nRow; ++i) {
      if(!(this_col_constant &= CppAD::Constant(x(i, j )))) break;
    }
    if(!this_col_constant) {
      colStart = j;
      break;
    }
  }
  
  for(int j = nCol-1; j >= 0; --j) {
    bool this_col_constant(true);
    for(int i = 0; i < nRow; ++i) {
      if(!(this_col_constant &= CppAD::Constant(x(i, j )))) break;
    }
    if(!this_col_constant) {
      colEnd = j + 1;
      break;
    }
  }
  if((rowStart == nRow) && (rowEnd == 0) && (colStart == nCol) && (colEnd == 0))
    return true;
  
  return false; // not all constant
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
  cout<<"setting up an atomic call of dims ("<<n1<<", "<<n2<<", "<<x1_is_constant<<") %*% ("<<n2<<", "<<n3<<", "<<x2_is_constant<<") = ("<<n1<<", "<<n3<<")"<<endl;
  
  atomic_matmult_class * atomic_matmult = global_matmult_factory.create();
  atomic_matmult->X1cat() = decide_matrix_category(x1);
  atomic_matmult->X2cat() = decide_matrix_category(x2);
  atomic_matmult->X1constant() = x1_is_constant;
  atomic_matmult->X2constant() = x2_is_constant;
  
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

void atomic_matmult(const MatrixXd_CppAD &x1,
		    const MatrixXd_CppAD &x2,
		    MatrixXd_CppAD &y) {
  int n1 = x1.rows(); // may not be general to all Eigen types
  int n2 = x1.cols();
  int n3 = x2.cols();

  if(n2 != x2.rows())
    cout<<"incommensurate matrices in atomic_matmult"<<endl;

  y.resize(n1, n3);

  int x1rowStart, x1rowEnd, x1colStart, x1colEnd;
  int x2rowStart, x2rowEnd, x2colStart, x2colEnd;
  bool x1_is_constant, x2_is_constant;

  x1_is_constant = delineate_constant_region(x1, x1rowStart, x1rowEnd, x1colStart, x1colEnd);
  x2_is_constant = delineate_constant_region(x2, x2rowStart, x2rowEnd, x2colStart, x2colEnd);

  using std::cout;
  cout<<"x1 region: "<<x1rowStart<<":"<<x1rowEnd<<", "<<x1colStart<<":"<<x1colEnd<<endl;
  cout<<"x2 region: "<<x2rowStart<<":"<<x2rowEnd<<", "<<x2colStart<<":"<<x2colEnd<<endl;
  
  // We will go through the matrix multiplication inefficiently and do
  // every constant operation.  Then we will record multiple atomics for various blocks that have variables
  for(int i = 0; i < n1; ++i) {
    for(int k = 0; k < n3; ++k) {
      CppAD::AD<double> oneVal = 0;
      for(int j = 0; j < n2; ++j) {
	if((CppAD::Constant(x1(i, j))) && (CppAD::Constant(x2(j, k)))) {
	  oneVal += x1(i, j) * x2(j, k);
	}
      }
      y(i, k) = oneVal;
    }
  }
  //
  if(x1rowStart > 0) {       // there are some constant first rows in x1, i.e. region A1.
    if(x2colStart < x2colEnd) { // there are some variable cols in x2, i.e. region B2
      // A1 * B2
      // populate up to 3 contributions to y[0:(x1rowStart-1), x1colStart:(x1colEnd-1)]
      // get constant * variable contribution to y[0:(x1rowStart-1), x1colStart:(x1colEnd-1)]
      cout<<"A1 * B2"<<endl;
      MatrixXd_CppAD y_const_var1;
      y_const_var1.resize(x1rowStart, x2colEnd - x2colStart);
      atomic_matmult_internal( x1.block(0, x2rowStart, x1rowStart, x2rowEnd - x2rowStart),
			       x2.block(x2rowStart, x2colStart, x2rowEnd - x2rowStart, x2colEnd - x2colStart),
			       y_const_var1, 1, 0);//1, 0);
      y.block(0, x2colStart, x1rowStart, x2colEnd - x2colStart) += y_const_var1;
    }
  } // End of A1 components
  if(x1rowStart < x1rowEnd) { // there are some variable rows in x1, i.e region B1
    if(x2colStart > 0) {      // there are some constant first cols in x2, i.e. region A2
      // B1 * A2
      // 3 cases here
      // set y[x1rowStart:(x1rowEnd-1), 0:(x2colStart-1)]
      MatrixXd_CppAD y_var_const1;
      y_var_const1.resize( x1rowEnd - x1rowStart, x2colStart );
      atomic_matmult_internal( x1.block(x1rowStart, x1colStart, x1rowEnd-x1rowStart, x1colEnd - x1colStart),
			       x2.block(x1colStart, 0, x1colEnd - x1colStart, x2colStart),
			       y_var_const1, 0, 1);
      y.block(x1rowStart, 0, x1rowEnd - x1rowStart, x2colStart ) += y_var_const1;
    }
    if(x2colStart < x2colEnd) { // B1 * B2, the most complication region to handle
      // There are up to 5 cases here!
      // constant * constant; constant * var OR var * constant; var * var OR constant * constant; var*constant OR constant*var; constant*constant;
      // The constant * constant pieces would have been dealt with above.
      // Start with the var * var region if there is is one, because it is the one region of interest if all elements of x1 and x2 are variables; fill to 0 otherwise
      if(((x1colStart < x2rowEnd) && (x1colEnd > x2rowStart)) || ((x2rowStart < x1colEnd) && (x2rowEnd > x1colStart))) { 
	int largerStart = x1colStart > x2rowStart ? x1colStart : x2rowStart;
	int smallerEnd = x1colEnd < x2rowEnd ? x1colEnd : x2rowEnd;
	if(smallerEnd <= largerStart) cout<<"Some logic appears to be wrong in the var*var part of nimDerivs matrix mult."<<endl;
	MatrixXd_CppAD y_var_var;
	y_var_var.resize(x1rowEnd - x1rowStart, x2colEnd - x2colStart);
	atomic_matmult_internal( x1.block(x1rowStart, largerStart, x1rowEnd-x1rowStart, smallerEnd - largerStart),
				 x2.block(largerStart, x2colStart, smallerEnd - largerStart,  x2colEnd - x2colStart),
				 y_var_var, 0, 0);
	y.block(x1rowStart, x2colStart, x1rowEnd - x1rowStart, x2colEnd - x2colStart) += y_var_var;
      }
      if(x1colStart < x2rowStart) { // there is a var * constant region
	int start = x1colStart;
	int end = x2rowStart < x1colEnd ? x2rowStart : x1colEnd;
	MatrixXd_CppAD y_var_const1;
	atomic_matmult_internal( x1.block(x1rowStart , start , (x1rowEnd - x1rowStart) , end-start ),
				 x2.block(start , x2colStart , end-start , (x2colEnd - x2colStart)),
				y_var_const1, 0, 1);
	y.block(x1rowStart, x2colStart, x1rowEnd - x1rowStart, x2colEnd - x2colStart) 
			       += y_var_const1;
      }
      if( x2rowStart < x1colStart ) { // there is a constant * var region. mutually exclusive with previous case
	int start = x2rowStart; 
	int end = x1colStart < x2rowEnd ? x1colStart : x2rowEnd;
	MatrixXd_CppAD y_const_var1;
	atomic_matmult_internal( x1.block(x1rowStart , start , (x1rowEnd - x1rowStart) , end-start ),
				 x2.block(start , x2colStart , end-start , (x2colEnd - x2colStart)),
				 y_const_var1, 1, 0);
	y.block(x1rowStart, x2colStart, x1rowEnd - x1rowStart, x2colEnd - x2colStart) 
			       += y_const_var1;
      }
      if( x1colEnd > x2rowEnd ) { // there is a var * constant region
	int end = x1colEnd;
	int start = x1colStart > x2rowEnd ? x1colStart : x2rowEnd;
	MatrixXd_CppAD y_var_const1;
	atomic_matmult_internal( x1.block(x1rowStart , start , (x1rowEnd - x1rowStart) , end-start ),
				 x2.block(start , x2colStart , end-start , (x2colEnd - x2colStart)),
				 y_var_const1, 0, 1);
	y.block(x1rowStart, x2colStart, x1rowEnd - x1rowStart, x2colEnd - x2colStart) 
			       += y_var_const1;
      }
      if( x2rowEnd > x1colEnd ) { // there is a constant * var region. mutually exclusive with previous case
	int end = x2rowEnd;
	int start = x2rowStart > x1colEnd ? x2rowStart : x1colEnd;
	MatrixXd_CppAD y_const_var1;
	atomic_matmult_internal( x1.block(x1rowStart , start , (x1rowEnd - x1rowStart) , end-start ),
				 x2.block(start , x2colStart , end-start , (x2colEnd - x2colStart)),
				 y_const_var1, 1, 0);
	y.block(x1rowStart, x2colStart, x1rowEnd - x1rowStart, x2colEnd - x2colStart) 
			       += y_const_var1;
      }
      // B1 * B2 components are done.
      // 
      // B1 * C2, in this section because C2 only exists if B2 exists
      if(x2colEnd < n3) { // there are some constant final cols in x2.
	// B1 * C2, which is only relevant if there is a B2 region.  Otherwise A2 covers all of x2
	// 3 cases here
	// set y[x1rowStart:(x1rowEnd-1), x2colEnd:(n3-1)]
	MatrixXd_CppAD y_var_const1;
	y_var_const1.resize( x1rowEnd - x1rowStart, n3 - x2colEnd );
	atomic_matmult_internal( x1.block(x1rowStart, x1colStart, x1rowEnd-x1rowStart, x1colEnd - x1colStart),
				 x2.block(x1colStart, x2colEnd, x1colEnd - x1colStart, n3 - x2colEnd),
				 y_var_const1, 0, 1);
	y.block(x1rowStart, x2colEnd, x1rowEnd - x1rowStart, n3 - x2colEnd ) += y_var_const1;
      }
      //
      //
      // C1 * B2, in this section because C1 exists only if B1 exists
      if(x1rowEnd < n1 ) {
	// populate up to 3 contributions to y[x1rowEnd:(n1-1), x1colStart:(x1colEnd-1)]
	// get constant * variable contribution to y[x1rowEnd:(n1-1), x1colStart:(x1colEnd-1)]
	MatrixXd_CppAD y_const_var1;
	y_const_var1.resize(n1-x1rowEnd, x2colEnd - x2colStart);
	atomic_matmult_internal( x1.block(x1rowEnd, x2rowStart, n1-x1rowEnd, x2rowEnd - x2rowStart),
				 x2.block(x2rowStart, x2colStart, x2rowEnd - x2rowStart, x2colEnd - x2colStart),
				 y_const_var1, 1, 0);
	y.block(x1rowEnd, x2colStart, n1-x1rowEnd, x2colEnd - x2colStart) += y_const_var1;
      }
    } // end of B2 within B1
  } // end of B1
}    

// void atomic_matmult_old(const MatrixXd_CppAD &x1,
// 		    const MatrixXd_CppAD &x2,
// 		    MatrixXd_CppAD &y) {
//   // Cases we want to handle include:
//   // 1. x1 and/or x2 are upper or lower diagonal
//   // 2. x1 or x2 have structural (constant) zeros.
//   // 3. x1 and/or x2 are constant, not parameter or variable (CppAD terms)
//   // 4. x1 and/or x2 are parameter.  That means they can change but not have derivatives tracked.
//   // These can be handled with some mix of identification here and/or inside the atomic
//   //
//   // 1. Diagonal: Identify that here because then we can reduce the number of elements
//   //              being copied around on CppAD tapes.  Use matmult flags to handle via Eigen.
//   // 2. Structural zeros: Identify that here because then we can reduce number of elements.
//   //                       and construct the return object partly in and partly outside of the atomic.
//   // 3. Constant (non-zero): Identify that here and use a different version of the atomics that holds the constant values.
//   // 4. Parameter: Identify that here and set matmult flags.
//   atomic_matmult_class * atomic_matmult = global_matmult_factory.create();  
//   int n1 = x1.rows(); // may not be general to all Eigen types
//   int n2 = x1.cols();
//   int n3 = x2.cols();

//   if(n2 != x2.rows())
//     cout<<"incommensurate matrices in atomic_matmult"<<endl;

//   y.resize(n1, n3);

//   // Look for entirely constant rows and or columns at start and/or end of regions
//   int x1rowStart, x1rowEnd, x1colStart, x1colEnd;
//   int x2rowStart, x2rowEnd, x2colStart, x2colEnd;
//   bool x1_is_constant, x2_is_constant;
//   x1_is_constant = delineate_constant_region(x1, x1rowStart, x1rowEnd, x1colStart, x1colEnd);
//   x2_is_constant = delineate_constant_region(x2, x2rowStart, x2rowEnd, x2colStart, x2colEnd);


//   // This gets somewhat complicated.  The goal is to break the matrix multiplication into
//   // blocks that are either both constant, one constant, or neither constant.
//   // Blocks that are both constant can be calculated *once* here, without the calculations
//   // needing to be recorded in AD.
//   // Blocks that have one side constant need to be recorded, using our atomic
//   //   but with flags so that only relevant parts of forward and reverse steps are done.
//   // Blocks with neight side constant need the atomic with both components of forward and reverse steps.
//   //
//   // x1 is split by rows into       x1 = [ constant A1 ; variable B1 ; constant C1]
//   // x2 is split by columns into    x2 = [ constant A2 | variable B2 | constant C2]
//   // Within x1, we also know the columns for which the variable block in B1 starts and ends.
//   // Within x2, we also know the rows for which the variable block in B2 starts and ends.
//   // A1 * A2, A1 * C2, C1 * A2 and C1 * C2 are easy.
//   // For A1 * B2 and C1 * B2, we separate into constant * constant, constant * variable, and constant * constant terms.
//   // For B1 * A2 and B1 * C2, we separate into constant * constant, variable * constant, and constant * constant terms.
//   // For B1 * B2, we separate into potentially five terms.
//   //
//   // Any of the regions can be potentially empty.  In practice many often will be empty.
//   // 
//   // If x1 is all constant, it will be all A1
//   // If x2 is all constant, it will be all A2
//   //
//   // An important use case is double taping of second-order reverse mode (or multiple first-order reverse to get to second derivatives).
//   // In those cases, there can be matrix multiplications using inputs that may have a single element non-zero.  It
//   // is important for efficiency that we record only what is needed for that single element.  This whole scheme should support that.
  
//   // Now we go through regions
//   if(x1rowStart > 0) {       // there are some constant first rows in x1, i.e. region A1.
//     if(x2colStart > 0) {     // there are some constant first cols in x2, i.e. region A2
//       // A1 * A2
//       // there is a constant * constant region
//       // set y[0:(x1rowStart-1), 0:(x2colStart-1)] as constant
//       y.block(0, 0, x1rowStart, x2colStart )
// 	= x1.block(0, 0, x1rowStart, n2) * x2.block(0, 0, n2, x2colStart);
//     }
//     if(x2colStart < x2colEnd) { // there are some variable cols in x2, i.e. region B2
//       // A1 * B2
//       // populate up to 3 contributions to y[0:(x1rowStart-1), x1colStart:(x1colEnd-1)]
//       // get constant * variable contribution to y[0:(x1rowStart-1), x1colStart:(x1colEnd-1)]
//       MatrixXd_CppAD y_const_var1;
//       y_const_var1.resize(x1rowStart, x2colEnd - x2colStart);
//       atomic_matmult_internal( x1.block(0, x2rowStart, x1rowStart, x2rowEnd - x2rowStart),
// 			       x2.block(x2rowStart, x2colStart, x2rowEnd - x2rowStart, x2colEnd - x2colStart),
// 			       y_const_var1, 1, 0);
//       y.block(0, x2colStart, x1rowStart, x2colEnd - x2colStart) =  y_const_var1;
      
//       if(x2RowStart > 0) { // there are some constant rows in x2
// 	// there is a constant * constant contribution to y[0:(x1rowStart-1), x1colStart:(x1colEnd-1)] from before the constant * variable region
// 	 y.block(0, x2colStart, x1rowStart, x2colEnd - x2colStart)
// 	  += x1.block( 0, 0, x1rowStart, x2RowStart) * x2.block(0, x2colStart, x2RowStart, x2colEnd - x2colStart);
//       }
//       if(x2RowEnd < n2) {
// 	// there is a constant * constant contribution to y[0:(x1rowStart-1), x1colStart:(x1colEnd-1)] from after the constant * variable region
// 	y.block(0, x2colStart, x1rowStart, x2colEnd - x2colStart)
// 	  += x1.block( 0, x2rowEnd, x1rowStart, n2 - x2rowEnd  ) * x2.block(x2rowEnd, x2colStart, n2 - x2rowEnd , x2colEnd - x2colStart);
//       }
//       //
//       if(x2colEnd < n3) { // there are some constant final cols in x2, i.e. C2
// 	// This is only relevant if there was a variable region.  Otherwise it will have been inluded as A2.
// 	// Hence this is inside the  if(x2colStart < x2colEnd) {} block
// 	//
// 	// A1 * C2
// 	// set y[0:(x1rowStart-1), x2colEnd:(n3-1)] as constant
// 	y.block(0, x2colEnd, x1rowStart, n3 - x2colEnd )
// 	  = x1.block(0, 0, x1rowStart, n2) * x2.block(0, x2colEnd, n2, n3 - x2colEnd);
//       }
//     }
//   }
//   //
//   if(x1rowStart < x1rowEnd) { // there are some variable rows in x1, i.e region B1
//     if(x2colStart > 0) {      // there are some constant first cols in x2, i.e. region A2
//       // B1 * A2
//       // 3 cases here
//       // set y[x1rowStart:(x1rowEnd-1), 0:(x2colStart-1)]
//       MatrixXd_CppAD y_var_const1;
//       y_var_const1.resize( x1rowEnd - x1rowStart, x2colStart );
//       atomic_matmult_internal( x1.block(x1rowStart, x1colStart, x1rowEnd-x1rowStart, x1colEnd - x1colStart),
// 			       x2.block(x1rowStart, 0, x1colEnd - x1colStart, x2colStart),
// 			       y_var_const1, 0, 1);
//       y.block(x1rowStart, 0, x1rowEnd - x1rowStart, x2colStart ) = y_var_const1;
//       if(x1colStart > 0) {
// 	y.block(x1rowStart, 0, x1rowEnd - x1rowStart, x2colStart )
// 	  += x1.block(x1rowStart, 0 , x1rowEnd - x1rowStart, x1colStart) * x2.block(0, 0, x1colStart, x2colStart);
//       }
//       if(x1colEnd < n2) {
// 	y.block(x1rowStart, 0, x1rowEnd - x1rowStart, x2colStart )
// 	  += x1.block(x1rowStart, x1colEnd, x1rowEnd - x1rowStart, n2 - x1colEnd) * x2.block(x1colEnd , 0 ,  n2 - x1colEnd  , x2colStart );
//       }
//     }

//     if(x2colStart < x2colEnd) { // B1 * B2, the most complication region to handle
//       // There are up to 5 cases here!
//       // constant * constant; constant * var OR var * constant; var * var OR constant * constant; var*constant OR constant*var; constant*constant;
//       // All x1 rows go from x1rowStart for size (x1rowEnd - x1rowStart)
//       // All x2 cols go from x2colStart for size (x2colEnd - x2colStart)
//       //
//       // Start with the var * var region if there is is one, because it is the one region of interest if all elements of x1 and x2 are variables; fill to 0 otherwise
//       if(((x1colStart < x2rowEnd) && (x1colEnd > x2rowStart)) || ((x2rowStart < x1colEnd) && (x2rowEnd > x1colStart))) { 
// 	int largerStart = x1colStart > x2rowStart ? x1colStart : x2rowStart;
// 	int smallerEnd = x1colEnd < x2rowEnd ? x1colEnd : x2rowEnd;
// 	if(smallerEnd <= largerStart) cout<<"Some logic appears to be wrong in the var*var part of nimDerivs matrix mult."<<endl;
// 	MatrixXd_CppAD y_var_var;
// 	y_var_var.resize(x1rowEnd - x1rowStart, x2colEnd - x2colStart);
// 	atomic_matmult_internal( x1.block(x1rowStart, largerStart, x1rowEnd-x1rowStart, smallerEnd - largerStart);
// 				 x2.block(largerStart, x2colStart, smallerEnd - largerStart,  x2colEnd - x2colStart),
// 				 y_var_const1, 1, 1);
// 	y.block(x1rowStart, x2colStart, x1rowEnd - x1rowStart, x2colEnd - x2colStart) = y_var_var;
//       } else {
// 	y.block(x1rowStart, x2colStart, x1rowEnd - x1rowStart, x2colEnd - x2colStart).fill(0);
//       }
//       //
//       if(x1colStart > 0 && x2rowStart > 0) { // there is a constant * constant region to start
// 	int smallerStart = x1colStart < x2rowStart ? x1colStart : x2rowStart;
// 	y.block(x1rowStart, x2colStart, x1rowEnd - x1rowStart, x2colEnd - x2colStart)
// 					+= x1.block(x1rowStart, 0 , x1rowEnd - x1rowStart, smallerStart ) * x2.block(0, x2colStart, smallerStart, x2colEnd - x2colStart);
//       }
//       if(x1colStart < x2rowStart) { // there is a constant * var region
// 	int start = x1colStart;
// 	int end = x2rowStart < x1colEnd ? x2rowStart : x1colEnd;
// 	MatrixXd_CppAD y_const_var1;
// 	atomic_matmult_internal( x1.block(x1rowStart , start , (x1rowEnd - x1rowStart) , end-start ),
// 				 x2.block(start , x2colStart , end-start , (x2colEnd - x2colStart)),
// 				 y_const_var1, 1, 0);
// 	y.block(x1rowStart, x2colStart, x1rowEnd - x1rowStart, x2colEnd - x2colStart) 
// 			       += y_const_var1;
//       }
//       if( x2rowStart < x1colStart ) { // there is a var * constant region. mutually exclusive with previous case
// 	int start = x1colStart;
// 	int end = x2rowStart < x1colEnd ? x2rowStart : x1colEnd;

//       }
// 	// THERE COULD BE A CONSTANT*CONSTANT REGION BETWEEN TWO CONSTANT * VAR OR VAR * CONSTANT REGIONS
    
//       //
//       if(x2colEnd < n3) { // there are some constant final cols in x2.
// 	// B1 * C2, which is only relevant if there is a B2 region.  Otherwise A2 covers all of x2
// 	// 3 cases here
// 	// set y[x1rowStart:(x1rowEnd-1), x2colEnd:(n3-1)]
// 	MatrixXd_CppAD y_var_const1;
// 	y_var_const1.resize( x1rowEnd - x1rowStart, n3 - x2colEnd );
// 	atomic_matmult_internal( x1.block(x1rowStart, x1colStart, x1rowEnd-x1rowStart, x1colEnd - x1colStart),
// 				 x2.block(x1rowStart, x2colEnd, x1colEnd - x1colStart, n3 - x2colEnd),
// 				 y_var_const1, 0, 1);
// 	y.block(x1rowStart, x2colEnd, x1rowEnd - x1rowStart, n3 - x2colEnd ) = y_var_const1;
// 	if(x1colStart > 0) {
// 	  y.block(x1rowStart, x2colEnd, x1rowEnd - x1rowStart, n3 - x2colEnd )
// 	    += x1.block(x1rowStart, 0 , x1rowEnd - x1rowStart, x1colStart) * x2.block(0, x2colEnd, x1colStart,  n3 - x2colEnd );
// 	}
// 	if(x1colEnd < n2) {
// 	  y.block(x1rowStart, x2colEnd, x1rowEnd - x1rowStart, n3 - x2colEnd )
// 	    += x1.block(x1rowStart, x1colEnd, x1rowEnd - x1rowStart, n2 - x1colEnd) * x2.block(x1colEnd , x2colEnd ,  n2 - x1colEnd  , n3 - x2colEnd );
// 	}
//       }
//     }
//   }
//   //
//   if(x2colStart < x2colEnd && x1rowEnd < n1) { // there are some constant last rows in x1, i.e. region C1, only relevant if there is a B1 region.
//     if(x2colStart > 0) {     // there are some constant first cols in x2, i.e. region A2
//       // C1 * A2
//       // there is a constant * constant region
//       // set y[x1rowEnd:(n1-1), 0:(x2colStart-1)] as constant
//       y.block(x1rowEnd, n1, n1-x1rowEnd, x2colStart )
// 	= x1.block(x1rowEnd, 0, n1-x1rowEnd, n2) * x2.block(0, 0, n2, x2colStart);
//     }
//     if(x2colStart < x2colEnd) { // there are some non-constant cols in x2, i.e. region B2
//       // C1 * B2
//       // populate up to 3 contributions to y[x1rowEnd:(n1-1), x1colStart:(x1colEnd-1)]
//       // get constant * variable contribution to y[x1rowEnd:(n1-1), x1colStart:(x1colEnd-1)]
//       MatrixXd_CppAD y_const_var1;
//       y_const_var1.resize(n1-x1rowEnd, x2colEnd - x2colStart);
//       atomic_matmult_internal( x1.block(x1rowEnd, x2rowStart, n1-x1rowEnd, x2rowEnd - x2rowStart),
// 			       x2.block(x2rowStart, x2colStart, x2rowEnd - x2rowStart, x2colEnd - x2colStart),
// 			       y_const_var1, 1, 0);
//       y.block(x1rowEnd, x2colStart, n1-x1rowEnd, x2colEnd - x2colStart) =  y_const_var1;
      
//       if(x2RowStart > 0) { // there are some constant rows in x2
// 	// there is a constant * constant contribution to y[x1rowEnd:(n1-1), x1colStart:(x1colEnd-1)] from before the constant * variable region
// 	y.block(x1rowEnd, x2colStart, n1-x1rowEnd, x2colEnd - x2colStart) 
// 	  += x1.block(x1rowEnd, 0, n1-x1rowEnd, x2RowStart) * x2.block(0, x2colStart, x2RowStart, x2colEnd - x2colStart);
//       }
//       if(x2RowEnd < n2) {
// 	// there is a constant * constant contribution to y[x1rowEnd:(n1-1), x1colStart:(x1colEnd-1)] from after the constant * variable region
// 	y.block(x1rowEnd, x2colStart, n1-x1rowEnd, x2colEnd - x2colStart) 
// 	  += x1.block( x1rowEnd, x2rowEnd, n1-x1rowEnd, n2 - x2rowEnd  ) * x2.block(x2rowEnd, x2colStart, n2 - x2rowEnd , x2colEnd - x2colStart);
//       }
//       //
//       if(x2colEnd < n3) { // there are some constant final cols in x2, i.e. region C2
// 	// C1 * C2
// 	// set y[x1rowEnd:(n1-1), x2colEnd:(n3-1)] as constant
// 	y.block(x1rowEnd, x2colEnd, n1-x1rowEnd, n3 - x2colEnd )
// 	  = x1.block(x1rowEnd, 0, n1-x1rowEnd, n2) * x2.block(0, x2colEnd, n2, n3 - x2colEnd);
//       }
//     }
//   }
  
//   std::vector<CppAD::AD<double> > xVec(n1*n2 + n2*n3);
//   atomic_matmult->set_n1(n1);
//   mat2vec(x1, xVec, 0);
//   mat2vec(x2, xVec, n1*n2);
//   std::vector<CppAD::AD<double> > yVec(n1*n3);
//   (*atomic_matmult)(xVec, yVec);
  
//   vec2mat(yVec, y);
// }

MatrixXd_CppAD nimDerivs_matmult(const MatrixXd_CppAD &x1,
				 const MatrixXd_CppAD &x2) {
  MatrixXd_CppAD ans;
  cout<<"Entering nimDerivs_matmult"<<endl;
  atomic_matmult(x1, x2, ans);
  cout<<"Leaving nimDerivs_matmult"<<endl;
  return ans;
}
