#ifndef NIMDERIVS_ATOMIC_CACHE_H__
#define NIMDERIVS_ATOMIC_CACHE_H__

template<typename T>
class atomic_cache_class {
public:
  atomic_cache_class();
  T* taylor_y_ptr();
  int nrow() const;
  void show_taylor_y();
  void set_cache(size_t order_low_set,
		 size_t order_up_set,
		 size_t order_up,
		 const CppAD::vector< T >&   taylor_x,
		 CppAD::vector< T >&         taylor_y
		 );
  template<typename S>
    void check_and_set_cache(S *owner,
			   const CppAD::vector< T >&        parameter_x,
			   const CppAD::vector<CppAD::ad_type_enum>&  type_x,
			   size_t                              order_up_check,
			     size_t order_up,
			   const CppAD::vector< T >&      taylor_x,
			   size_t taylor_y_size);

private:
  CppAD::vector< T > taylor_x_cache;
  CppAD::vector< T > taylor_y_cache;

  int taylor_x_cache_nrow = -1;
  int taylor_y_cache_nrow = -1;
  int current_cache_order = -1;
  
  void show_cache_generic(const CppAD::vector< T >& v,
			  int &cache_nrow);
  void set_cache_generic(size_t order_low_set,
			 size_t order_up_set,
			 size_t order_up,
		 const CppAD::vector< T >& v,
		 CppAD::vector< T >& v_cache,
		 int &cache_nrow
		 );
  void check_and_set_cache_size(size_t order_up_set,
				size_t order_up,
				size_t v_size,
				CppAD::vector<T>& v_cache,
				int &cache_nrow);
  void set_Xcache(size_t order_low_set,
		  size_t  order_up_set,
		  size_t order_up,
		  const CppAD::vector< T >& taylor_x);
  bool check_Xcache(size_t  order_up_check,
		    size_t order_up,
		    const CppAD::vector< T >& taylor_x);
  void set_Ycache( size_t order_low_set,
		   size_t  order_up_set,
		   size_t order_up,
		   const CppAD::vector< T >& taylor_y);    
};

#endif
