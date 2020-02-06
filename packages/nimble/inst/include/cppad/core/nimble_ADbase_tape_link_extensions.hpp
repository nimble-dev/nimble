  // Special method added for nimble.  See ad.hpp.
template<class Base>
void AD<Base>::set_tape_info_nimble(tape_id_t tape_id, local::ADTape<Base>* tape_handle_, bool recover) {
  // Modified from AD<Base>::tape_ptr
  size_t thread = size_t( tape_id % CPPAD_MAX_NUM_THREADS );
  CPPAD_ASSERT_KNOWN(
    thread == thread_alloc::thread_num(),
    "In set_tape_info_nimble: thread from tape_id does not match thread from thread_num."
    );
  static tape_id_t saved_tape_id;
  static local::ADTape<Base>* saved_tape_handle;
  local::ADTape<Base>** tape_handle_ptr = tape_handle(thread);
  if(recover) {
    *tape_id_ptr(thread) = saved_tape_id;
    *tape_handle(thread) = saved_tape_handle;
  } else {
    saved_tape_id = *tape_id_ptr(thread);
    saved_tape_handle = *tape_handle(thread);
    *tape_id_ptr(thread) = tape_id;
    *tape_handle(thread) = tape_handle_;
  }
}

template<class Base>
void* AD<Base>::get_handle_address_nimble() {
  void* ans = static_cast<void *>(tape_ptr());
  return ans;
}
