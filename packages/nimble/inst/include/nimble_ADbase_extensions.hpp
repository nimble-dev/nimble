// Added for use by nimble when different shared libraries are being managed from R.
// In such a case, a need arises to set fake tape_id and tape ptr values.
// When CppAD version is updated, this will need to be inserted into the new version.
// The definition for set_tape_info_nimble is in tape_link.hpp, where related CppAD methods appear.
static void set_tape_info_nimble(tape_id_t tape_id, local::ADTape<Base>* tape_handle_, bool recover = false);
static void* get_handle_address_nimble(); // used only for diagnostic purposes
static tape_id_t get_tape_id_nimble();
static local::ADTape<Base>* get_tape_handle_nimble();
