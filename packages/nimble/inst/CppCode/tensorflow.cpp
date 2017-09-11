/*
 * NIMBLE: an R package for programming with BUGS models.
 * Copyright (C) 2014-2017 Perry de Valpine, Christopher Paciorek,
 * Daniel Turek, Clifford Anderson-Bergman, Nick Michaud, Fritz Obermeyer,
 * Duncan Temple Lang.
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, a copy is available at
 * https://www.R-project.org/Licenses/
 */

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
      gradients_(inputNames.size()),
      gradient_values_(inputNames.size()),
      input_pos_(0),
      output_pos_(outputNames.size()),
      gradient_pos_(inputNames.size()),
      has_gradients_(true) {
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

  // Set gradients.
  for (int i = 0; i < inputNames.size(); ++i) {
    const std::string name = inputNames[i] + "_gradient";
    gradients_[i].oper = TF_GraphOperationByName(graph_, name.c_str());
    gradients_[i].index = 0;
    if (gradients_[i].oper == NULL) {
      has_gradients_ = false;
      break;
    }
  }
}

NimTf_Runner::~NimTf_Runner() {
  NIM_ASSERT_EQ(input_pos_, 0);
  NIM_ASSERT_EQ(output_pos_, num_outputs());
  NIM_ASSERT_EQ(gradient_pos_, num_inputs());

  // Clean up.
  TF_CloseSession(session_, status_);
  TF_CHECK_OK(status_);
  TF_DeleteSession(session_, status_);
  TF_CHECK_OK(status_);
  TF_DeleteGraph(graph_);
  TF_DeleteStatus(status_);
}

void NimTf_Runner::NimTf_setInput(double& scalar) {
  NIM_ASSERT_LT(input_pos_, num_inputs());
  NIM_ASSERT_EQ(output_pos_, num_outputs());
  NIM_ASSERT_EQ(gradient_pos_, num_inputs());

  input_values_[input_pos_] = TF_NewTensor(
      TF_DOUBLE, NULL, 0, &scalar, sizeof(double), fake_deallocator, NULL);
  input_pos_ += 1;
}

void NimTf_Runner::NimTf_setInput(NimArrBase<double>& nimArr) {
  NIM_ASSERT1(!nimArr.isMap(), "Cannot handle mapped array");
  NIM_ASSERT_LT(input_pos_, num_inputs());
  NIM_ASSERT_EQ(output_pos_, num_outputs());
  NIM_ASSERT_EQ(gradient_pos_, num_inputs());

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

void NimTf_Runner::NimTf_setInput(const std::vector<int64_t>& dims,
                                  double* data) {
  NIM_ASSERT_LT(input_pos_, num_inputs());
  NIM_ASSERT_EQ(output_pos_, num_outputs());
  NIM_ASSERT_EQ(gradient_pos_, num_inputs());

  size_t byteSize = sizeof(double);
  for (int i = 0; i < dims.size(); ++i) {
    byteSize *= dims[i];
  }
  input_values_[input_pos_] =
      TF_NewTensor(TF_DOUBLE, dims.data(), dims.size(), data, byteSize,
                   fake_deallocator, NULL);
  input_pos_ += 1;
}

void NimTf_Runner::NimTf_run() {
  NIM_ASSERT_EQ(input_pos_, num_inputs());
  NIM_ASSERT_EQ(output_pos_, num_outputs());
  NIM_ASSERT_EQ(gradient_pos_, num_inputs());

  // Run the graph.
  {
    const TF_Buffer* run_options = NULL;
    const int ntargets = 0;
    const TF_Operation** targets = NULL;
    TF_Buffer* run_metadata = NULL;
    TF_SessionRun(session_, run_options,                                  //
                  inputs_.data(), input_values_.data(), num_inputs(),     //
                  outputs_.data(), output_values_.data(), num_outputs(),  //
                  targets, ntargets,                                      //
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

void NimTf_Runner::NimTf_runGradient() {
  NIM_ASSERT1(has_gradients_, "Tensorflow graph is missing gradients");
  NIM_ASSERT_EQ(input_pos_, num_inputs());
  NIM_ASSERT_EQ(output_pos_, num_outputs());
  NIM_ASSERT_EQ(gradient_pos_, num_inputs());

  // Run the graph.
  {
    const TF_Buffer* run_options = NULL;
    const int ntargets = 0;
    const TF_Operation** targets = NULL;
    TF_Buffer* run_metadata = NULL;
    TF_SessionRun(session_, run_options,                                     //
                  inputs_.data(), input_values_.data(), num_inputs(),        //
                  gradients_.data(), gradient_values_.data(), num_inputs(),  //
                  targets, ntargets,                                         //
                  run_metadata, status_);
    TF_CHECK_OK(status_);
  }

  // Clean up input values.
  for (int i = 0; i < input_values_.size(); ++i) {
    TF_DeleteTensor(input_values_[i]);
  }
  input_pos_ = 0;
  gradient_pos_ = 0;
}

void NimTf_Runner::NimTf_getOutput(double& scalar) {
  NIM_ASSERT_EQ(input_pos_, 0);
  NIM_ASSERT_LT(output_pos_, num_outputs());
  NIM_ASSERT_EQ(gradient_pos_, num_inputs());

  TF_Tensor* tensor = output_values_[output_pos_];
  scalar = *static_cast<double*>(TF_TensorData(tensor));

  // Clean up output values.
  TF_DeleteTensor(tensor);
  output_pos_ += 1;
}

void NimTf_Runner::NimTf_getOutput(NimArrBase<double>& nimArr) {
  NIM_ASSERT1(!nimArr.isMap(), "Cannot handle mapped array");
  NIM_ASSERT_EQ(input_pos_, 0);
  NIM_ASSERT_LT(output_pos_, num_outputs());
  NIM_ASSERT_EQ(gradient_pos_, num_inputs());

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

void NimTf_Runner::NimTf_getGradient(const std::vector<int64_t>& dims,
                                     double* data) {
  NIM_ASSERT_EQ(input_pos_, 0);
  NIM_ASSERT_EQ(output_pos_, num_outputs());
  NIM_ASSERT_LT(gradient_pos_, num_inputs());

  TF_Tensor* tensor = gradient_values_[gradient_pos_];
  const int n = TF_NumDims(tensor);
  NIM_ASSERT_EQ(dims.size(), TF_NumDims(tensor));
  for (int i = 0; i < n; ++i) {
    NIM_ASSERT_EQ(dims[i], TF_Dim(tensor, n - i - 1));  // Note the transpose.
  }
  memcpy(data, TF_TensorData(tensor), TF_TensorByteSize(tensor));

  // Clean up gradient values.
  TF_DeleteTensor(tensor);
  gradient_pos_ += 1;
}

#if NIMBLE_HAVE_CPPAD

NimTf_Op::NimTf_Op(NimTf_Runner& runner)
    : CppAD::atomic_base<double>("tensorflow"), runner_(runner) {
  NIM_ASSERT_EQ(runner.num_outputs(), 1);
}

void NimTf_Op::NimTf_setInput(ad_double& scalar) {
  packed_arg_.push_back(0);
  packed_arg_.push_back(scalar);
}

void NimTf_Op::NimTf_setInput(NimArrBase<ad_double>& nimArr) {
  const int n = nimArr.numDims();
  packed_arg_.push_back(n);
  int64_t size = 1;
  for (int i = 0; i < n; ++i) {
    int64_t dim = nimArr.dimSize(i);
    packed_arg_.push_back(n);
    size *= dim;
  }
  for (int64_t i = 0; i < size; ++i) {
    packed_arg_.push_back(i);
  }
}

void NimTf_Op::NimTf_run() {
  packed_result_.resize(1);
  (*this)(packed_arg_, packed_result_);
  packed_arg_.clear();
}

void NimTf_Op::NimTf_getOutput(ad_double& scalar) {
  scalar = packed_result_[0];
}

void NimTf_Op::NimTf_getOutput(NimArrBase<ad_double>& nimArr) {
  const int n = nimArr.numDims();
  for (int i = 0; i < n; ++i) {
    NIM_ASSERT_EQ(nimArr.dimSize(i), 1);
  }
  nimArr.getPtr()[0] = packed_result_[0];
}

bool NimTf_Op::forward(size_t p, size_t q, const CppAD::vector<bool>& vx,
                       CppAD::vector<bool>& vy, const CppAD::vector<double>& tx,
                       CppAD::vector<double>& ty) {
  // This implements only computation of a single scalar value (no derivatives).
  if (p != 0) return false;
  if (q != 0) return false;
  NIM_ASSERT_EQ(ty.size(), 1);

  // TODO Is it safe to allow tensorflow to mutate this input data?
  double* mutable_tx = const_cast<double*>(tx.data());

  // Set inputs.
  size_t pos = 0;
  for (int i = 0; i < runner_.num_inputs(); ++i) {
    // Set array numDims.
    const int numDims = static_cast<int>(tx[pos++]);

    // Set array dims.
    std::vector<int64_t> dims;
    int64_t size = 1;
    for (int j = 0; j < numDims; ++j) {
      const int dim = static_cast<int>(tx[pos++]);
      dims.push_back(dim);
      size *= dim;
    }
    std::reverse(dims.begin(), dims.end());  // Note the transpose.

    // Set array data.
    runner_.NimTf_setInput(dims, mutable_tx + pos);
    pos += size;
  }

  // Run the tensorflow graph.
  runner_.NimTf_run();

  // Get outputs.
  vy[0] = true;
  runner_.NimTf_getOutput(ty[0]);

  return true;
}

bool NimTf_Op::reverse(size_t q, const CppAD::vector<double>& tx,
                       const CppAD::vector<double>& ty,
                       CppAD::vector<double>& px,
                       const CppAD::vector<double>& py) {
  // This implements only computation of a single scalar gradient.
  if (q != 1) return false;
  NIM_ASSERT_EQ(ty.size(), 1);

  // TODO Is it safe to allow tensorflow to mutate this input data?
  double* mutable_tx = const_cast<double*>(tx.data());

  // Set inputs.
  size_t pos = 0;
  for (int i = 0; i < runner_.num_inputs(); ++i) {
    // Set array numDims.
    const int numDims = static_cast<int>(tx[pos++]);

    // Set array dims.
    std::vector<int64_t> dims;
    int size = 1;
    for (int j = 0; j < numDims; ++j) {
      const int dim = static_cast<int>(tx[pos++]);
      dims.push_back(dim);
      size *= dim;
    }
    std::reverse(dims.begin(), dims.end());  // Note the transpose.

    // Set array data.
    runner_.NimTf_setInput(dims, mutable_tx + pos);
    pos += size;
  }

  // Run the tensorflow graph.
  runner_.NimTf_runGradient();

  // Get outputs.
  px.resize(tx.size());
  pos = 0;
  for (int i = 0; i < runner_.num_inputs(); ++i) {
    // Get array numDims.
    px[pos] = 0;
    const int numDims = static_cast<int>(tx[pos]);
    pos += 1;

    // Get array dims.
    std::vector<int64_t> dims;
    int64_t size = 1;
    for (int j = 0; j < numDims; ++j) {
      px[pos] = 0;
      const int dim = static_cast<int>(tx[pos]);
      dims.push_back(dim);
      size *= dim;
      pos += 1;
    }
    std::reverse(dims.begin(), dims.end());  // Note the transpose.

    // Get array data.
    runner_.NimTf_getGradient(dims, px.data() + pos);
    for (int j = 0; j < size; ++j) {
      px[pos + j] *= py[0];
    }
    pos += size;
  }

  return true;
}

#endif  // NIMBLE_HAVE_CPPAD
