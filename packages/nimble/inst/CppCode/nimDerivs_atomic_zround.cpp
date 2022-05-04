/*********************************/
/*****     log_zround        ****/
/*********************************/
/*
 y = f(x) = round(x)

 We thought it should work to use CPPAD_DISCRETE_FUNCTION(double, nimround)
 where nimround is a function that returns round(x).

 However, a CPPAD_DISCRETE_FUNCTION fails through double-taping.
*/
#include <nimble/nimDerivs_atomic_zround.h>


atomic_zround_class::atomic_zround_class(const std::string& name) : 
  CppAD::atomic_three<double>(name)
{ }

bool atomic_zround_class::for_type(
				   const CppAD::vector<double>&               parameter_x ,
				   const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				   CppAD::vector<CppAD::ad_type_enum>&        type_y      )
{
  type_y[0] = type_x[0];
  return true;
}
// Not sure this is ever needed.
bool atomic_zround_class::for_type(
				   const CppAD::vector<CppAD::AD<double> >&               parameter_x ,
				   const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				   CppAD::vector<CppAD::ad_type_enum>&        type_y      )
{
  type_y[0] = type_x[0];
  return true;
}

bool atomic_zround_class::rev_depend(
				     const CppAD::vector<double>&          parameter_x ,
				     const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				     CppAD::vector<bool>&                depend_x    ,
				     const CppAD::vector<bool>&          depend_y
				     ) {
  depend_x[0] = depend_y[0];
  return true;
}

bool atomic_zround_class::forward(
     const CppAD::vector<double>&               parameter_x  ,
     const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
     size_t                              need_y       ,
     size_t                              order_low    ,
     size_t                              order_up     ,
     const CppAD::vector<double>&               taylor_x     ,
     CppAD::vector<double>&                     taylor_y     )
{
  if(order_low <= 0 & order_up >= 0) {
    taylor_y[0] = round(taylor_x[0]);    
  }
  if(order_low <= 1 & order_up >= 1) {
    taylor_y[1] = 0.;
  }
  return true;
}

bool atomic_zround_class::forward(
				   const CppAD::vector< CppAD::AD<double> >&               parameter_x  ,
				   const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
				   size_t                              need_y       ,
				   size_t                              order_low    ,
				   size_t                              order_up     ,
				   const CppAD::vector< CppAD::AD<double> >&               taylor_x     ,
				   CppAD::vector< CppAD::AD<double> >&                     taylor_y     ) {
  if(order_low <= 0 & order_up >= 0) {
    taylor_y[0] = nimDerivs_zround(taylor_x[0]);
  }
  if(order_low <= 1 & order_up >= 1) {
    taylor_y[1] = nimDerivs_zround(CppAD::AD<double>(0.));
  }
  return true;
}

bool atomic_zround_class::reverse(
     const CppAD::vector<double>&               parameter_x ,
     const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
     size_t                              order_up    ,
     const CppAD::vector<double>&               taylor_x    ,
     const CppAD::vector<double>&               taylor_y    ,
     CppAD::vector<double>&                     partial_x   ,
     const CppAD::vector<double>&               partial_y   )
{
  if(order_up >= 0) {
    partial_x[0] = 0.;
  }
  if(order_up >= 1) {
    partial_x[1] = 0.;
  }
  return true;
}

bool atomic_zround_class::reverse(
				   const CppAD::vector< CppAD::AD<double> >&               parameter_x ,
				   const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				   size_t                              order_up    ,
				   const CppAD::vector< CppAD::AD<double> >&               taylor_x    ,
				   const CppAD::vector< CppAD::AD<double> >&               taylor_y    ,
				   CppAD::vector< CppAD::AD<double> >&                     partial_x   ,
				   const CppAD::vector< CppAD::AD<double> >&               partial_y   ) {
  if(order_up >= 0) {
    partial_x[0] = nimDerivs_zround(CppAD::AD<double>(0.));
  }
  if(order_up >= 1) {
    partial_x[1] = nimDerivs_zround(CppAD::AD<double>(0.));
  }
  return true;
};

CppAD::AD<double> nimDerivs_zround(const CppAD::AD<double> x) {
  atomic_zround_class* atomic_zround;
  bool recording = CppAD::AD<double>::get_tape_handle_nimble() != nullptr;
  if(!recording) {
    atomic_zround = new atomic_zround_class("atomic_zround");
  } else {
    atomic_zround = track_atomic_zround(CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr(),
					  CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage() );
  }
  CppAD::vector< CppAD::AD<double> > in(1);
  in[0] = x;
  CppAD::vector< CppAD::AD<double> > out(1);
  (*atomic_zround)(in, out);
  if(!recording) {
    delete atomic_zround;
  }
  return out[0];
}

/**************************/

double nimDerivs_floor_class::operator()(double x) {return floor(x);}
double nimDerivs_floor_class::operator()(int x) {return floor(x);}
CppAD::AD<double> nimDerivs_floor_class::operator()(const CppAD::AD<double> &x) {return nimDerivs_floor(x);}

atomic_floor_class::atomic_floor_class(const std::string& name) : 
  atomic_discrete_class<nimDerivs_floor_class>(name)
{ }

CppAD::AD<double> nimDerivs_floor(const CppAD::AD<double> x) {
    atomic_floor_class* atomic_floor;
  bool recording = CppAD::AD<double>::get_tape_handle_nimble() != nullptr;
  if(!recording) {
    atomic_floor = new atomic_floor_class("atomic_floor");
  } else {
    atomic_floor = track_atomic_floor(CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr(),
					  CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage() );
  }
  CppAD::vector< CppAD::AD<double> > in(1);
  in[0] = x;
  CppAD::vector< CppAD::AD<double> > out(1);
  (*atomic_floor)(in, out);
  if(!recording) {
    delete atomic_floor;
  }
  return out[0];
}

/**************************/

double nimDerivs_ceil_class::operator()(double x) {return ceil(x);}
double nimDerivs_ceil_class::operator()(int x) {return ceil(x);}
CppAD::AD<double> nimDerivs_ceil_class::operator()(const CppAD::AD<double> &x) {return nimDerivs_ceil(x);}

atomic_ceil_class::atomic_ceil_class(const std::string& name) : 
  atomic_discrete_class<nimDerivs_ceil_class>(name)
{ }

CppAD::AD<double> nimDerivs_ceil(const CppAD::AD<double> x) {
    atomic_ceil_class* atomic_ceil;
  bool recording = CppAD::AD<double>::get_tape_handle_nimble() != nullptr;
  if(!recording) {
    atomic_ceil = new atomic_ceil_class("atomic_ceil");
  } else {
    atomic_ceil = track_atomic_ceil(CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr(),
					  CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage() );
  }
  CppAD::vector< CppAD::AD<double> > in(1);
  in[0] = x;
  CppAD::vector< CppAD::AD<double> > out(1);
  (*atomic_ceil)(in, out);
  if(!recording) {
    delete atomic_ceil;
  }
  return out[0];
}


/**************************/

double nimDerivs_ftrunc_class::operator()(double x) {return ftrunc(x);}
double nimDerivs_ftrunc_class::operator()(int x) {return ftrunc(x);}
CppAD::AD<double> nimDerivs_ftrunc_class::operator()(const CppAD::AD<double> &x) {return nimDerivs_ftrunc(x);}

atomic_ftrunc_class::atomic_ftrunc_class(const std::string& name) : 
  atomic_discrete_class<nimDerivs_ftrunc_class>(name)
{ }

CppAD::AD<double> nimDerivs_ftrunc(const CppAD::AD<double> x) {
    atomic_ftrunc_class* atomic_ftrunc;
  bool recording = CppAD::AD<double>::get_tape_handle_nimble() != nullptr;
  if(!recording) {
    atomic_ftrunc = new atomic_ftrunc_class("atomic_ftrunc");
  } else {
    atomic_ftrunc = track_atomic_ftrunc(CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr(),
					  CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage() );
  }
  CppAD::vector< CppAD::AD<double> > in(1);
  in[0] = x;
  CppAD::vector< CppAD::AD<double> > out(1);
  (*atomic_ftrunc)(in, out);
  if(!recording) {
    delete atomic_ftrunc;
  }
  return out[0];
}

/**************************/

double nimDerivs_nimRound_class::operator()(double x) {return nimRound(x);}
double nimDerivs_nimRound_class::operator()(int x) {return nimRound(x);}
CppAD::AD<double> nimDerivs_nimRound_class::operator()(const CppAD::AD<double> &x) {return nimDerivs_nimRound(x);}

atomic_nimRound_class::atomic_nimRound_class(const std::string& name) : 
  atomic_discrete_class<nimDerivs_nimRound_class>(name)
{ }

CppAD::AD<double> nimDerivs_nimRound(const CppAD::AD<double> x) {
  atomic_nimRound_class* atomic_nimRound;
  bool recording = CppAD::AD<double>::get_tape_handle_nimble() != nullptr;
  if(!recording) {
    atomic_nimRound = new atomic_nimRound_class("atomic_nimRound");
  } else {
    atomic_nimRound = track_atomic_nimRound(CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr(),
					  CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage() );
  }
  CppAD::vector< CppAD::AD<double> > in(1);
  in[0] = x;
  CppAD::vector< CppAD::AD<double> > out(1);
  (*atomic_nimRound)(in, out);
  if(!recording) {
    delete atomic_nimRound;
  }
  return out[0];
}

