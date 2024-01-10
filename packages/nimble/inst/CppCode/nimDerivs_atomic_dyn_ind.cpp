#include <nimble/nimDerivs_atomic_dyn_ind.h>

atomic_dyn_ind_get_class::atomic_dyn_ind_get_class(const std::string& name) :
  CppAD::atomic_three<double>(name) {};

bool atomic_dyn_ind_get_class::for_type(
                                        const CppAD::vector<double>&               parameter_x ,
                                        const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
                                        CppAD::vector<CppAD::ad_type_enum>&        type_y      )
{
  // for_type and rev_depend are calling once when recording.
  // not during each playback.
  // Therefore they can't depend on the value of index.
  CppAD::ad_type_enum max_type = type_x[0];
  size_t nx = type_x.size();
  if(nx > 1) {
    for(size_t i = 1; i < nx; ++i)
      max_type = std::max(max_type, type_x[i]);
  }
  type_y[0] = max_type;
  return true;
}


bool atomic_dyn_ind_get_class::rev_depend(
                                          const CppAD::vector<double>&          parameter_x ,
                                          const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
                                          CppAD::vector<bool>&                depend_x    ,
                                          const CppAD::vector<bool>&          depend_y
                                          ) {
  size_t n = type_x.size() - 1;
  for(size_t i = 0; i < n; ++i) depend_x[i] = depend_y[0];
  depend_x[n] = depend_y[0];
  return true;
}

bool atomic_dyn_ind_get_class::forward(
                                       const CppAD::vector<double>&               parameter_x  ,
                                       const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
                                       size_t                              need_y       ,
                                       size_t                              order_low    ,
                                       size_t                              order_up     ,
                                       const CppAD::vector<double>&               taylor_x     ,
                                       CppAD::vector<double>&                     taylor_y     ) {
  size_t nrow = order_up + 1;
  size_t ncol = taylor_x.size()/nrow;
  //size_t nx = ncol-1; // i.e. length of the actual "x"
  size_t i_index = ncol-1;
  //  std::cout<<"in get::forward "<<order_low<<" "<<order_up<<" "<<taylor_x[0 + i_index*nrow]<<" "<<parameter_x[0 + i_index*nrow]<<" ";
  int index = static_cast<int>(taylor_x[0 + i_index*nrow]);
  if(order_low == 0)
    taylor_y[0] = taylor_x[0 + index * nrow];
  if(order_low <= 1 & order_up >= 1) {
    taylor_y[1] = taylor_x[1 + index * nrow];
  }
  if(order_low <= 2 & order_up >= 2) {
    taylor_y[2] = taylor_x[2 + index * nrow];
  }
  //std::cout<<index<<" "<<taylor_y[0]<<std::endl;
  //                  for(size_t i = 0; i < taylor_x.size();++i) std::cout<<taylor_x[i]<<" ";
  //                                                                        std::cout<<std::endl;
  return true;
}

bool atomic_dyn_ind_get_class::forward(
                                       const CppAD::vector<CppAD::AD<double> >&               parameter_x  ,
                                       const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
                                       size_t                              need_y       ,
                                       size_t                              order_low    ,
                                       size_t                              order_up     ,
                                       const CppAD::vector<CppAD::AD<double> >&               taylor_x     ,
                                       CppAD::vector<CppAD::AD<double> >&                     taylor_y     ) {
  size_t nrow = order_up + 1;
  size_t ncol = taylor_x.size()/nrow;
  size_t nx = ncol-1; // i.e. length of the actual "x"
  size_t i_index = ncol-1;
  CppAD::AD<double> index = taylor_x[0 + i_index*nrow];
  if(order_low == 0) {
    taylor_y[0] = dyn_ind_get(taylor_x, index, 0, nrow, nx);
  }
  if(order_low <= 1 & order_up >= 1) {
    taylor_y[1] = dyn_ind_get(taylor_x, index, 1, nrow, nx);
  }
  if(order_low <= 1 & order_up >= 1) {
    taylor_y[2] = dyn_ind_get(taylor_x, index, 2, nrow, nx);
  }
  return true;
}

bool atomic_dyn_ind_get_class::reverse(
                                       const CppAD::vector<double>&               parameter_x ,
                                       const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
                                       size_t                              order_up    ,
                                       const CppAD::vector<double>&               taylor_x    ,
                                       const CppAD::vector<double>&               taylor_y    ,
                                       CppAD::vector<double>&                     partial_x   ,
                                       const CppAD::vector<double>&               partial_y   )
{
  size_t nrow = order_up + 1;
  size_t ncol = taylor_x.size()/nrow;
  size_t nx = ncol-1;
  size_t i_index = ncol-1;
  int index = static_cast<int>(taylor_x[0 + i_index*nrow]);
  // std::cout<<"in get reverse "<<order_up<<" "<<index<<std::endl;
  for(size_t order = 0; order <= order_up; ++order) {
    for(size_t i = 0; i < nx; ++i) partial_x[order + i * nrow] = 0;
    partial_x[order + index * nrow] = partial_y[order];
    partial_x[order + i_index *nrow] = 0;
  }
  return true;
}


bool atomic_dyn_ind_get_class::reverse(
                                       const CppAD::vector<CppAD::AD<double> >&               parameter_x ,
                                       const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
                                       size_t                              order_up    ,
                                       const CppAD::vector<CppAD::AD<double> >&               taylor_x    ,
                                       const CppAD::vector<CppAD::AD<double> >&               taylor_y    ,
                                       CppAD::vector<CppAD::AD<double> >&                     partial_x   ,
                                       const CppAD::vector<CppAD::AD<double> >&               partial_y   )
{
  size_t nrow = order_up + 1;
  size_t ncol = taylor_x.size()/nrow;
  size_t nx = ncol-1;
  size_t i_index = ncol-1;
  CppAD::vector<CppAD::AD<double> > new_x(nx);
  CppAD::AD<double> index = taylor_x[0 + i_index*nrow];

  for(size_t order = 0; order <= order_up; ++order) {
    for(size_t i = 0; i < nx; ++i) partial_x[order + i*nrow] = CppAD::AD<double>(0);
    dyn_ind_set(partial_x, index, partial_y[order], order, nrow, nx );
    partial_x[order + i_index * nrow] = CppAD::AD<double>(0);
  }
  return true;
}

atomic_dyn_ind_set_class::atomic_dyn_ind_set_class(const std::string& name) :
  CppAD::atomic_three<double>(name) {};


bool atomic_dyn_ind_set_class::for_type(
                      const CppAD::vector<double>&               parameter_x ,
                      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
                      CppAD::vector<CppAD::ad_type_enum>&        type_y      )
{
  CppAD::ad_type_enum max_type = type_x[0];
  size_t length_x = type_x.size();
  size_t nx = length_x-2;
  size_t ny = type_y.size();
  assert(nx == ny);
  if(length_x > 1) { // should always be true
    for(size_t i = 1; i < length_x; ++i)
      max_type = std::max(max_type, type_x[i]);
  }
//  max_type = std::max(max_type, type_x[length_x-1]);
  for(size_t i = 0; i < ny; ++i) type_y[i] = max_type;
  //std::cout<<"exiting dyn_ind_set for_type with max_type = "<<max_type<<std::endl;
  //           std::cout<<"reference types: "<<CppAD::constant_enum <<" "<<CppAD::dynamic_enum <<" "<<CppAD::variable_enum<<std::endl;
  return true;
}

bool atomic_dyn_ind_set_class::rev_depend(
                        const CppAD::vector<double>&          parameter_x ,
                        const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
                        CppAD::vector<bool>&                depend_x    ,
                        const CppAD::vector<bool>&          depend_y
                        ) {
  // To-do:
  // Here is a place-holder
  bool any_depend_y = depend_y[0];
  size_t ny = depend_y.size();
  if(ny > 1) {
    for(size_t i = 1; i < ny; ++i)
      any_depend_y = any_depend_y || depend_y[i];
  }
  size_t length_x = depend_x.size();
  size_t nx = length_x-2;
  assert(nx == ny);
  for(size_t i = 0; i < length_x; ++i) depend_x[i] = any_depend_y;
  //std::cout<<"exiting dyn_ind_set rev_depend with any_depend_y = "<<any_depend_y<<std::endl;
  return true;
}

bool atomic_dyn_ind_set_class::forward(
                     const CppAD::vector<double>&               parameter_x  ,
                     const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
                     size_t                              need_y       ,
                     size_t                              order_low    ,
                     size_t                              order_up     ,
                     const CppAD::vector<double>&               taylor_x     ,
                     CppAD::vector<double>&                     taylor_y     ) {
  size_t nrow = order_up + 1;
  size_t ny = taylor_y.size()/nrow;
  //size_t length_x = ny + 2;
  size_t i_index = ny;
  size_t i_x = ny+1;
  int index = static_cast<int>(taylor_x[0 + i_index*nrow]);
  //std::cout<<"dyn_ind_set forward order_low = "<<order_low <<" order_up = "<<order_up <<" index = "<<index <<std::endl;
  for(size_t order = order_low; order <= order_up; order++) {
    for(size_t i = 0; i < ny; ++i) taylor_y[order + i*nrow] = taylor_x[order + i*nrow];
    taylor_y[order + index*nrow] = taylor_x[order + i_x*nrow];
  }
  return true;
}

bool atomic_dyn_ind_set_class::forward(
                     const CppAD::vector<CppAD::AD<double> >&               parameter_x  ,
                     const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
                     size_t                              need_y       ,
                     size_t                              order_low    ,
                     size_t                              order_up     ,
                     const CppAD::vector<CppAD::AD<double> >&               taylor_x     ,
                     CppAD::vector<CppAD::AD<double> >&                     taylor_y     ) {
  size_t nrow = order_up + 1;
  size_t ny = taylor_y.size()/nrow;
  //size_t length_x = ny + 2;
  size_t i_index = ny;
  size_t i_x = ny+1;
  CppAD::AD<double> index = taylor_x[0 + i_index*nrow];
  for(size_t order = order_low; order <= order_up; order++) {
    for(size_t i = 0; i < ny; ++i) taylor_y[order + i * nrow] = taylor_x[order + i * nrow];
    dyn_ind_set(taylor_y, index, taylor_x[order + i_x * nrow], order, nrow);
  }
  return true;
}

bool atomic_dyn_ind_set_class::reverse(
                     const CppAD::vector<double>&               parameter_x ,
                     const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
                     size_t                              order_up    ,
                     const CppAD::vector<double>&               taylor_x    ,
                     const CppAD::vector<double>&               taylor_y    ,
                     CppAD::vector<double>&                     partial_x   ,
                     const CppAD::vector<double>&               partial_y   )
{
  size_t nrow = order_up + 1;
  size_t ny = taylor_y.size()/nrow;
  //size_t length_x = ny + 2;
  size_t i_index = ny;
  size_t i_x = ny+1;
  int index = static_cast<int>(taylor_x[0 + i_index*nrow]);
  //std::cout<<"dyn_ind_set reverse order_up = "<<order_up <<" index = "<<index <<std::endl;

  for(size_t order = 0; order <= order_up; ++order) {
    for(size_t i = 0; i < ny; ++i) partial_x[order + i*nrow] = partial_y[order + i*nrow];
    partial_x[order + i_index*nrow] = 0;
    partial_x[order + index*nrow] = 0;
    partial_x[order + i_x*nrow] = partial_y[order + index*nrow];
  }
  return true;
}

bool atomic_dyn_ind_set_class::reverse(
                     const CppAD::vector<CppAD::AD<double> >&               parameter_x ,
                     const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
                     size_t                              order_up    ,
                     const CppAD::vector<CppAD::AD<double> >&               taylor_x    ,
                     const CppAD::vector<CppAD::AD<double> >&               taylor_y    ,
                     CppAD::vector<CppAD::AD<double> >&                     partial_x   ,
                     const CppAD::vector<CppAD::AD<double> >&               partial_y   )
{
  size_t nrow = order_up + 1;
  size_t ny = taylor_y.size()/nrow;
  //size_t length_x = ny + 2;
  size_t i_index = ny;
  size_t i_x = ny+1;

  CppAD::AD<double> index = taylor_x[0 + i_index*nrow];

  for(size_t order = 0; order <= order_up; ++order) {
    for(size_t i = 0; i < ny; ++i) partial_x[order + i*nrow] = partial_y[order + i*nrow];
    partial_x[order + i_index*nrow] = CppAD::AD<double>(0);
    dyn_ind_set(partial_x, index, CppAD::AD<double>(0), order, nrow, ny);
    partial_x[order + i_x * nrow]   = dyn_ind_get(partial_y, index, order, nrow);
  }
  return true;
}

CppAD::AD<double> dyn_ind_get(const CppAD::vector<CppAD::AD<double> > &x,
                              const CppAD::AD<double> &index,
                              const size_t offset,
                              const size_t nrow,
                              const int nx_) {
  if(CppAD::Constant(index)) return x[offset + CppAD::Value(index)*nrow];
  atomic_dyn_ind_get_class *atomic_dyn_ind_get;
  //static atomic_dyn_ind_get_class dyn_get("dyn_get");
  CppAD::vector<CppAD::AD<double> > y(1);
  size_t nx;
  if(nx_ < 0) {
    nx = x.size() / nrow;
  } else {
    nx = nx_;
  }
  CppAD::vector<CppAD::AD<double> > packed_x(nx+1);
  for(size_t i = 0; i < nx; ++i) packed_x[i] = x[offset + i*nrow];
  packed_x[nx] = index;

  bool recording = CppAD::AD<double>::get_tape_handle_nimble() != nullptr;

//  std::cout<<"calling dyn_ind_get "<<recording<<" "<<CppAD::Value(index)<<std::endl;

  if(!recording) {
    atomic_dyn_ind_get = new atomic_dyn_ind_get_class("atomic_dyn_ind_get");
  } else {
    void *tape_mgr = CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr();
    atomic_dyn_ind_get = new_atomic_dyn_ind_get(tape_mgr, "atomic_dyn_ind_get");
  }

  (*atomic_dyn_ind_get)(packed_x, y);

  if(!recording) {
    delete atomic_dyn_ind_get;
  } else {
    track_nimble_atomic(atomic_dyn_ind_get,
                        CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr(),
                        CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage() );
  }


  return y[0];
}

// y[i] = x;
// done as y2[] = y[]; y2[i] = x; y[] = y2[];
void dyn_ind_set(CppAD::vector<CppAD::AD<double> > &y,
                 const CppAD::AD<double> &index,
                 const CppAD::AD<double> &x,
                 const size_t offset,
                 const size_t nrow,
                 const int ny_ ) {
  if(CppAD::Constant(index)) {
    y[offset + CppAD::Value(index)*nrow] = x;
    return;
  }
  atomic_dyn_ind_set_class *atomic_dyn_ind_set;
  //static atomic_dyn_ind_set_class dyn_set("dyn_set");
  size_t ny;
  if(ny_ < 0) {
    ny = y.size() / nrow;
  } else {
    ny = ny_;
  }
  CppAD::vector<CppAD::AD<double> > packed_x(ny + 2);
  for(size_t i = 0; i < ny; ++i) packed_x[i] = y[offset + i*nrow];
  packed_x[ny] = index;
  packed_x[ny+1] = x;

  bool recording = CppAD::AD<double>::get_tape_handle_nimble() != nullptr;

  if(!recording) {
    atomic_dyn_ind_set = new atomic_dyn_ind_set_class("atomic_dyn_ind_set");
  } else {
    void *tape_mgr = CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr();
    atomic_dyn_ind_set = new_atomic_dyn_ind_set(tape_mgr, "atomic_dyn_ind_set");
  }

  if(offset == 0 && nrow == 1 && ny_ < 0) {
    (*atomic_dyn_ind_set)(packed_x, y);
  } else {
    CppAD::vector<CppAD::AD<double> > temp_y(ny);
    (*atomic_dyn_ind_set)(packed_x, temp_y);
    for(size_t i = 0; i < ny; ++i) y[offset + i*nrow] = temp_y[i];
  }

  if(!recording) {
    delete atomic_dyn_ind_set;
  } else {
    track_nimble_atomic(atomic_dyn_ind_set,
                        CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr(),
                        CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage() );
  }

}
