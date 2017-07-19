#include <base64/base64.h>
#include <nimble/tensorflow.h>

#define TF_CHECK_OK(status)                                           \
  {                                                                   \
    if (TF_GetCode(status) != TF_OK) {                                \
      NIMERROR("Tensorflow error at " __FILE__ ":%d\n  %s", __LINE__, \
               TF_Message(status));                                   \
    }                                                                 \
  }

// This is passed as the deallocator of memory that is not owned.
void fake_deallocator(void* data, size_t len, void* arg) {}  // Does nothing.

NimTf_Runner::NimTf_Runner(const std::string& graphDefBase64,
                           const std::string& configBase64,
                           const std::vector<std::string>& inputNames,
                           const std::vector<std::string>& outputNames)
    : status_(TF_NewStatus()),
      graph_(TF_NewGraph()),
      session_(NULL),
      inputs_(inputNames.size()),
      input_values_(inputNames.size()),
      outputs_(outputNames.size()),
      output_values_(outputNames.size()),
      input_pos_(0),
      output_pos_(outputNames.size()) {
  // Initialize graph.
  {
    std::string graphDef;
    Base64::Decode(graphDefBase64, &graphDef);
    TF_Buffer buffer;
    buffer.data = graphDef.data();
    buffer.length = graphDef.size();
    TF_ImportGraphDefOptions* options = TF_NewImportGraphDefOptions();
    TF_GraphImportGraphDef(graph_, &buffer, options, status_);
    TF_CHECK_OK(status_);
    TF_DeleteImportGraphDefOptions(options);
  }

  // Initialize session.
  {
    std::string config;
    Base64::Decode(configBase64, &config);
    TF_SessionOptions* options = TF_NewSessionOptions();
    TF_SetConfig(options, config.data(), config.size(), status_);
    TF_CHECK_OK(status_);
    session_ = TF_NewSession(graph_, options, status_);
    TF_CHECK_OK(status_);
    TF_DeleteSessionOptions(options);
  }

  // Set inputs
  for (int i = 0; i < inputNames.size(); ++i) {
    inputs_[i].oper = TF_GraphOperationByName(graph_, inputNames[i].c_str());
    inputs_[i].index = 0;
  }

  // Set outputs.
  for (int i = 0; i < outputNames.size(); ++i) {
    outputs_[i].oper = TF_GraphOperationByName(graph_, outputNames[i].c_str());
    outputs_[i].index = 0;
  }
}

NimTf_Runner::~NimTf_Runner() {
  NIM_ASSERT_EQ(input_pos_, 0);
  NIM_ASSERT_EQ(output_pos_, outputs_.size());

  // Clean up.
  TF_CloseSession(session_, status_);
  TF_CHECK_OK(status_);
  TF_DeleteSession(session_, status_);
  TF_CHECK_OK(status_);
  TF_DeleteGraph(graph_);
  TF_DeleteStatus(status_);
}

void NimTf_Runner::NimTf_setInput(double& scalar) {
  NIM_ASSERT_LT(input_pos_, inputs_.size());
  NIM_ASSERT_EQ(output_pos_, outputs_.size());
  input_values_[input_pos_] = TF_NewTensor(
      TF_DOUBLE, NULL, 0, &scalar, sizeof(double), fake_deallocator, NULL);
  input_pos_ += 1;
}

void NimTf_Runner::NimTf_setInput(NimArrBase<double>& nimArr) {
  NIM_ASSERT1(!nimArr.isMap(), "Cannot handle mapped array");
  NIM_ASSERT_LT(input_pos_, inputs_.size());
  NIM_ASSERT_EQ(output_pos_, outputs_.size());
  static std::vector<int64_t> dims;
  const int n = nimArr.numDims();
  dims.resize(n);
  size_t byteSize = sizeof(double);
  for (int i = 0; i < n; ++i) {
    const int dim = nimArr.dimSize(i);
    dims[n - i - 1] = dim;  // Note the transpose.
    byteSize *= dim;
  }
  input_values_[input_pos_] =
      TF_NewTensor(TF_DOUBLE, dims.data(), dims.size(), nimArr.getPtr(),
                   byteSize, fake_deallocator, NULL);
  input_pos_ += 1;
}

void NimTf_Runner::NimTf_run() {
  NIM_ASSERT_EQ(input_pos_, inputs_.size());
  NIM_ASSERT_EQ(output_pos_, outputs_.size());

  // Run the graph.
  {
    const TF_Buffer* run_options = NULL;
    const int ntargets = 0;
    const TF_Operation** targets = NULL;
    TF_Buffer* run_metadata = NULL;
    TF_SessionRun(session_, run_options,                                    //
                  inputs_.data(), input_values_.data(), inputs_.size(),     //
                  outputs_.data(), output_values_.data(), outputs_.size(),  //
                  targets, ntargets,                                        //
                  run_metadata, status_);
    TF_CHECK_OK(status_);
  }

  // Clean up input values.
  for (int i = 0; i < input_values_.size(); ++i) {
    TF_DeleteTensor(input_values_[i]);
  }
  input_pos_ = 0;
  output_pos_ = 0;
}

void NimTf_Runner::NimTf_getOutput(double& scalar) {
  NIM_ASSERT_EQ(input_pos_, 0);
  NIM_ASSERT_LT(output_pos_, outputs_.size());
  TF_Tensor* tensor = output_values_[output_pos_];
  scalar = *static_cast<double*>(TF_TensorData(tensor));

  // Clean up output values.
  TF_DeleteTensor(tensor);
  output_pos_ += 1;
}

void NimTf_Runner::NimTf_getOutput(NimArrBase<double>& nimArr) {
  NIM_ASSERT1(!nimArr.isMap(), "Cannot handle mapped array");
  NIM_ASSERT_EQ(input_pos_, 0);
  NIM_ASSERT_LT(output_pos_, outputs_.size());
  TF_Tensor* tensor = output_values_[output_pos_];
  const int n = TF_NumDims(tensor);
  NIM_ASSERT_EQ(nimArr.numDims(), TF_NumDims(tensor));
  for (int i = 0; i < n; ++i) {
    // Note the transpose:
    NIM_ASSERT_EQ(nimArr.dimSize(i), TF_Dim(tensor, n - i - 1));
  }
  memcpy(nimArr.getPtr(), TF_TensorData(tensor), TF_TensorByteSize(tensor));

  // Clean up output values.
  TF_DeleteTensor(tensor);
  output_pos_ += 1;
}
