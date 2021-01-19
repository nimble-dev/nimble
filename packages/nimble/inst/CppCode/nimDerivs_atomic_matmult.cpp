#include <nimble/nimDerivs_atomic_matmult.h>

#ifdef USE_NEW_DYNAMIC
atomic_matmult_class::atomic_matmult_class(const std::string& name) :
  CppAD::atomic_three<double>(name),
  n1_(0),
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
   
  new (&X1mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(&taylor_x[0],
								  n1, n2, EigStrDyn(nrow*n1, nrow) );
  //   EigenTemplateTypes<double>::typeEigenConstMapStrd X1map(&taylor_x[0 + 1*nrow],
  //							   n1, n2, EigStrDyn(nrow*n1, nrow) );
  new (&X2mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(&taylor_x[0 + (n1*n2)*nrow],
								  n2, n3, EigStrDyn(nrow*n2, nrow ) );
  if(order_low <= 0 & order_up >= 0) { // value
    // We could compile different cases depending on need for strides or not.
    // We could use const maps.
    new (&Ymap) EigenTemplateTypes<double>::typeEigenMapStrd(&taylor_y[0],
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
    new (&dX1mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(&taylor_x[1],
								     n1, n2, EigStrDyn(nrow*n1, nrow) );
    new (&dX2mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(&taylor_x[1 + (n1*n2)*nrow],
								     n2, n3, EigStrDyn(nrow*n2, nrow ) );
    new (&dY_map) EigenTemplateTypes<double>::typeEigenMapStrd(&taylor_y[1],
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
  int n1 = get_n1();//taylor_x[0];
  int n3 = m/n1;
  int n2 = n/(n1 + n3);
  new (&X1mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(&taylor_x[0],
								  n1, n2, EigStrDyn(nrow*n1, nrow) );
  new (&X2mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(&taylor_x[0 + (n1*n2)*nrow],
								  n2, n3, EigStrDyn(nrow*n2, nrow ) );
  new (&Yadjoint_mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(&partial_y[0],
									 n1, n3, EigStrDyn(nrow*n1, nrow ) );
  new (&X1adjoint_map) EigenTemplateTypes<double>::typeEigenMapStrd(&partial_x[0],
								    n1, n2, EigStrDyn(nrow*n1, nrow) );
  new (&X2adjoint_map) EigenTemplateTypes<double>::typeEigenMapStrd(&partial_x[0 + (n1*n2)*nrow],
								    n2, n3, EigStrDyn(nrow*n2, nrow ) );
  if(order_up >= 0) {
    // reverse 1
    partial_x[0] = 0;
    X1adjoint_map = Yadjoint_mapC *  X2mapC.transpose();
    X2adjoint_map = X1mapC.transpose() * Yadjoint_mapC;
  }
  if(order_up >= 1) {
    // reverse 2
    new (&X1dot_mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(&taylor_x[1],
									n1, n2, EigStrDyn(nrow*n1, nrow) );
    new (&X2dot_mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(&taylor_x[1 + (n1*n2)*nrow],
									n2, n3, EigStrDyn(nrow*n2, nrow ) );
    new (&Ydot_adjoint_mapC) EigenTemplateTypes<double>::typeEigenConstMapStrd(&partial_y[1],
									       n1, n3, EigStrDyn(nrow*n1, nrow ) );
    new (&X1dot_adjoint_map) EigenTemplateTypes<double>::typeEigenMapStrd(&partial_x[1],
									  n1, n2, EigStrDyn(nrow*n1, nrow) );
    new (&X2dot_adjoint_map) EigenTemplateTypes<double>::typeEigenMapStrd(&partial_x[1 + (n1*n2)*nrow],
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
#else
// Use local Eigen maps

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
  //   EigenTemplateTypes<double>::typeEigenConstMapStrd X1map(&taylor_x[0 + 1*nrow],
  //							   n1, n2, EigStrDyn(nrow*n1, nrow) );
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

void atomic_matmult(const MatrixXd_CppAD &x1,
		    const MatrixXd_CppAD &x2,
		    MatrixXd_CppAD &y) {
  atomic_matmult_class * atomic_matmult = global_matmult_factory.create();
  int n1 = x1.rows(); // may not be general to all Eigen types
  int n2 = x1.cols();
  int n3 = x2.cols();
  std::vector<CppAD::AD<double> > xVec(n1*n2 + n2*n3);
  atomic_matmult->set_n1(n1);
  mat2vec(x1, xVec, 0);
  mat2vec(x2, xVec, n1*n2);
  std::vector<CppAD::AD<double> > yVec(n1*n3);
  (*atomic_matmult)(xVec, yVec);
  y.resize(n1, n3);
  vec2mat(yVec, y);
}

MatrixXd_CppAD nimDerivs_matmult(const MatrixXd_CppAD &x1,
				 const MatrixXd_CppAD &x2) {
  MatrixXd_CppAD ans;
  atomic_matmult(x1, x2, ans);
  return ans;
}
