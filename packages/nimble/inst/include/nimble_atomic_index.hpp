#define USING_CPPAD_IN_NIMBLE

/* The following class exists only for its static method. */
/* It is important that this not cross compilation units (e.g. via alternative design using inline) */
/* because the point is to be able to have one compilation unit set the vec_ptr of another. */
template <class Base>
class atomic_index_info_vec_manager_nimble {
public:
  static std::vector<atomic_index_info>* manage(int set = 0, // 0=return stored ptr, 1 = set stored ptr to provided value and return, 2 = set stored ptr to internval object
						std::vector<atomic_index_info>* input_vec_ptr = 0
						) {
    static std::vector<atomic_index_info> vec;
    static std::vector<atomic_index_info>* vec_ptr;
    static bool first=true; // initialization happens only once
    if(first) {             // This has the effect that on first call vec_ptr is set to &vec.
      vec_ptr = &vec;
      first = false;
    }
    //    std::cout<<" [&vec_ptr:"<< &vec_ptr<<"] ";
    //    std::cout<<"entering atomic info manage. vec = "<<&vec<<" set = "<<set<<" vec_ptr = "<<vec_ptr<<" input_vec_ptr = "<<input_vec_ptr<<std::endl;
    switch(set) {
    case 0 :
      return vec_ptr;
    case 1 :
      vec_ptr = input_vec_ptr;
      return vec_ptr;
    case 2:
      vec_ptr = &vec;
      return vec_ptr;
    default:
      std::cout<<"Any invalid value was provided to manage-atomic_index_info_vec"<<std::endl;
      return &vec;
    }
  }
};
