#ifndef NIMDERIVS_ATOMIC_CACHE__
#define NIMDERIVS_ATOMIC_CACHE__

#include <nimble/nimDerivs_atomic_cache.h>

// #define DEBUG_ATOMIC_CACHE_

template<typename T>
atomic_cache_class<T>::atomic_cache_class() :
  taylor_x_cache_nrow(0),
  taylor_y_cache_nrow(0),
  current_cache_order(-1)
{}

template<typename T>
T* atomic_cache_class<T>::taylor_y_ptr() {
  return &(taylor_y_cache[0]);
}

template<typename T>
int atomic_cache_class<T>::nrow() const {
  return taylor_y_cache_nrow;
}

template<typename T>
void atomic_cache_class<T>::check_and_set_cache_size(size_t order_up_set,
						     size_t order_up,
						     size_t v_size,
						     CppAD::vector<T>& v_cache,
						     int &cache_nrow) {
  int nrow = order_up + 1; // pertains to v
  int nrow_set = order_up_set + 1;
  if(nrow_set > cache_nrow) {
  int n = static_cast<int>(static_cast<double>(v_size)/nrow);
#ifdef DEBUG_ATOMIC_CACHE_
    std::cout<<"resizing from "<<cache_nrow<<" to "<<v_size<<" for nrow = "<<nrow<<" n = "<<n<<std::endl;
#endif
    /* We need to preserve contents because lower orders might already be cached. */
    if(cache_nrow > 0) {
      CppAD::vector<T> temp(v_cache);
      v_cache.resize( n * nrow_set );
      for(size_t irow = 0; irow < cache_nrow; ++irow) { // irow is order
	for(size_t i = 0; i < n; ++i) {
	  v_cache[irow + i*nrow] = temp[irow + i*cache_nrow]; // restore contents, which had cache_nrow, to new cache, which has nrow
	}
      }
    } else {
      v_cache.resize( n * nrow_set );
    }
    cache_nrow = nrow_set;
  }
}

/*
A reason this is not quite as simple as it might seem is the following.

On a 0th order call to forward, taylor_x will have nrow = 1.
On a 1st order call to forward, taylor_x will have nrow = 2.
We expand the cache to maximum seen on previous calls.
Hence, the cache may be bigger than the input on a 0th-order call.
*/

template<typename T>
void atomic_cache_class<T>::show_cache_generic(const CppAD::vector< T >& v_cache,
					       int &cache_nrow) {
  std::cout<<"cache_nrow = "<<cache_nrow<<std::endl;
  int n = static_cast<int>(static_cast<double>(v_cache.size())/cache_nrow);
  for(size_t irow = 0; irow < cache_nrow; ++irow) {
    for(size_t i = 0; i < n; ++i) {
      std::cout<<v_cache[irow + i * cache_nrow]<<"\t";
    }
    std::cout<<std::endl;
  }
}

template<typename T>
void atomic_cache_class<T>::show_taylor_y() {
  std::cout<<"taylor_y_cache:"<<std::endl;
  show_cache_generic(taylor_y_cache, taylor_y_cache_nrow);
}

template<typename T>
void atomic_cache_class<T>::set_cache_generic(size_t order_low_set,
					      size_t order_up_set,
					      size_t order_up,
					      const CppAD::vector< T >& v,
					      CppAD::vector< T >& v_cache,
					      int &cache_nrow
					      ) {
  int nrow = order_up + 1; // pertains to v
  if(nrow <= 0) return;
  int nrow_set = order_up_set + 1;
  /* 
    Expanding cache is destructive, but we don't need the old values.
    Caches will be expanded but not contracted.  That is because we 
    assume that whatever order is needed will typically be used repeatedly. 
  */
  check_and_set_cache_size(order_up_set, order_up, v.size(), v_cache, cache_nrow);
  /*
    Copy up to order_up.
    It is guaranteed that cache_nrow >= nrow.
  */
  typename CppAD::vector< T >::const_iterator v_iter;
  typename CppAD::vector< T >::iterator v_cache_iter;
  int n = static_cast<int>(static_cast<double>(v.size())/nrow);
  
  for(size_t irow = order_low_set; irow < nrow_set; ++irow) { // irow is order
    // v[i + j*nrow] is the j-th element of the i-th order
#ifdef DEBUG_ATOMIC_CACHE_
    std::cout<<"setting cache for order = "<<irow<<": ";
#endif
    v_iter = v.begin() + irow;
    v_cache_iter = v_cache.begin() + irow;
    const typename CppAD::vector< T >::const_iterator v_end(v_iter + n*nrow);
    for( ; v_iter != v_end; v_iter += nrow, v_cache_iter += cache_nrow) {
      *v_cache_iter = *v_iter;
#ifdef DEBUG_ATOMIC_CACHE_
      std::cout<<*v_iter<<"\t";
#endif
    }
#ifdef DEBUG_ATOMIC_CACHE_
    std::cout<<std::endl;
#endif
  }
}

template<typename T>
void atomic_cache_class<T>::set_Xcache(size_t  order_low_set,
				       size_t  order_up_set,
				       size_t order_up,
				       const CppAD::vector< T >& taylor_x) {
#ifdef DEBUG_ATOMIC_CACHE_
  std::cout<<"setting Xcache"<<std::endl;
#endif
  set_cache_generic(order_low_set, order_up_set, order_up, taylor_x, taylor_x_cache, taylor_x_cache_nrow);
}

template<typename T>
bool atomic_cache_class<T>::check_Xcache(size_t order_up_check,
					 size_t order_up,
					 const CppAD::vector< T >& taylor_x) {
#ifdef DEBUG_ATOMIC_CACHE_
  std::cout<<"checking Xcache (outcome is false unless \"all true\" is shown)..."<<std::endl;
#endif
  if(order_up_check > current_cache_order) return(false);
  int nrow = order_up + 1; // for taylor_x
  typename CppAD::vector< T >::const_iterator v_iter;
  typename CppAD::vector< T >::iterator v_cache_iter;
  int n = static_cast<int>(static_cast<double>(taylor_x.size())/nrow);
  int nrow_check = order_up_check + 1;
  
  for(size_t irow = 0; irow < nrow_check; ++irow) { // irow is order
    // v[i + j*nrow] is the j-th element of the i-th order
    v_iter = taylor_x.begin() + irow;
    v_cache_iter = taylor_x_cache.begin() + irow;
    const typename CppAD::vector< T >::const_iterator v_end(v_iter + n*nrow);
    for( ; v_iter != v_end; v_iter += nrow, v_cache_iter += taylor_x_cache_nrow) {
      if(*v_cache_iter != *v_iter) return(false);
    }
#ifdef DEBUG_ATOMIC_CACHE_
    std::cout<<"true for order "<<irow<<"..."<<std::endl;
#endif
  }
#ifdef DEBUG_ATOMIC_CACHE_
  std::cout<<"all true"<<std::endl;
#endif    
  return true;
}

template<typename T>
void atomic_cache_class<T>::set_Ycache( size_t order_low_set,
					size_t order_up_set,
					size_t order_up,
					const CppAD::vector< T >& taylor_y) {
#ifdef DEBUG_ATOMIC_CACHE_
  std::cout<<"setting Ycache"<<std::endl;
#endif
  set_cache_generic(order_low_set, order_up_set, order_up, taylor_y, taylor_y_cache, taylor_y_cache_nrow);
}

template<typename T>
void atomic_cache_class<T>::set_cache(size_t                   order_low_set,
				      size_t                   order_up_set,
				      size_t                   order_up, 
				   const CppAD::vector< T >&   taylor_x,
				   CppAD::vector< T >&         taylor_y
				   ) {
  /*
    It should be guaranteed that taylor_x and taylor_y are well-defined and 
    up to date with each other (that is Y = f(X), Y' = g(X, X'), etc., up to order_up.
  */
  if(&taylor_y == &taylor_y_cache) return; // This means forward was called from check_and_set_cache.
  
#ifdef DEBUG_ATOMIC_CACHE_
  std::cout<<"setting X and Y caches with order_up = "<<order_up<<std::endl;
#endif  
  set_Xcache(order_low_set, order_up_set, order_up, taylor_x);
  set_Ycache(order_low_set, order_up_set, order_up, taylor_y);
  current_cache_order = order_up_set;
}
 
template<typename T>
template<typename S>
void atomic_cache_class<T>::check_and_set_cache(S *owner,
						const CppAD::vector< T >&        parameter_x,
						const CppAD::vector<CppAD::ad_type_enum>&  type_x,
						size_t                           order_up_check,
						size_t                           order_up,
						const CppAD::vector< T >&        taylor_x,
						size_t taylor_y_size) {
  // 1. Check if all taylor_x orders relevant to the sweep status match
  //    what is in the cache.
  if(order_up_check > order_up)
    std::cout<<"Something is wrong in check_and_set_caches"<<std::endl;
  
  bool cache_is_current = check_Xcache(order_up_check, order_up, taylor_x);
#ifdef DEBUG_ATOMIC_CACHE_
  std::cout<<"check and set caches with cache_is_current = "<<cache_is_current<<std::endl;
#endif
  
  // 2. Otherwise, calculate Ymap and Ydot_map into the cache by calling forward.
  if(!cache_is_current) {
    check_and_set_cache_size(order_up_check, order_up, taylor_y_size, taylor_y_cache, taylor_y_cache_nrow);
    int nrow = order_up + 1;
    int n = static_cast<int>(static_cast<double>(taylor_x.size())/nrow);
    int m = static_cast<int>(static_cast<double>(taylor_y_size)/nrow);
    int temp_nrow = order_up_check + 1;
    
    CppAD::vector< T > taylor_x_temp( n * temp_nrow ); // These can be necessary because order_up might not match order_up_check and hence migt not match sizes of taylor_x and taylor_y
    CppAD::vector< T > taylor_y_temp( m * temp_nrow );
    // We need to copy from taylor_x to taylor_x_temp, possibly not for all orders.
    // This is equivalent to caching into taylor_x_temp.
    set_cache_generic(0, order_up_check, order_up, taylor_x, taylor_x_temp, temp_nrow);

    owner->forward( parameter_x,
		    type_x,
		    CppAD::variable_enum,
		    0, order_up_check,
		    taylor_x_temp,
		    taylor_y_temp);
    //    set_Xcache(0, order_up, taylor_x);
    //    current_cache_order = order_up;
  }
}

// template class atomic_cache_class<double>;
// template class atomic_cache_class<CppAD::AD<double> >;

#endif
