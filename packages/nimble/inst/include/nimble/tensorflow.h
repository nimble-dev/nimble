#ifndef __NIMBLE_TENSORFLOW_H
#define __NIMBLE_TENSORFLOW_H

#include <nimble/NimArr.h>
#include <nimble/Utils.h>
#include <tensorflow/c/c_api.h>
#include <string>
#include <vector>

// This class runs a single tensorflow graph. Instances should be created as
// static local variables of a nimbleFunction.
// Example usage:
//  NimArr<double, 1> axpy(NimArr<double, 1> ARG1_a_,
//                         NimArr<double, 1> ARG2_x_,
//                         NimArr<double, 1> ARG3_y_) {
//    static NimTfRunner* runner = 
//      NimTfBuilder("0123456789abcdef")  // The string is a graphDefBase64.
//        .withInput("ARG1_a_")
//        .withInput("ARG2_x_")
//        .withInput("ARG3_y_")
//        .withOutput("RESULT_")
//        .build();
//
//    NimArr<double, 1> result;
//    result->initialize(0.0, false, ARG1_x_.dimSize(0));
//    runner->setInput(ARG1_a_);
//    runner->setInput(ARG2_x_);
//    runner->setInput(ARG3_y_);
//    runner->run();
//    runner->getOutput(result);
//    return result;
//  }
class NimTfRunner {
  TF_Status* status_;
  TF_Graph* graph_;
  TF_Session* session_;

  std::vector<TF_Output> inputs_;
  std::vector<TF_Tensor*> input_values_;

  std::vector<TF_Output> outputs_;
  std::vector<TF_Tensor*> output_values_;

  int input_pos_;
  int output_pos_;

 public:
  NimTfRunner(const std::string& graphDefBase64,
              const std::vector<std::string>& inputNames,
              const std::vector<std::string>& outputName);
  ~NimTfRunner();

  // These must be called in strict order.
  void setInput(NimArrBase<double>& nimArr);
  void run();
  void getOutput(NimArrBase<double>& nimArr);

 private:
  // Safety measures:
  NimTfRunner(const NimTfRunner&) { NIMERROR("Illegal"); }
  NimTfRunner& operator=(const NimTfRunner&) { NIMERROR("Illegal"); }
};

// Helper to construct static instances of NimTfRunner.
// Example usage:
//  static NimTfRunner* nimTf =
//     NimTfBuilder("0123456789abcdef")  // The string is a graphDefBase64.
//      .withInput("ARG1_a_")
//      .withInput("ARG2_x_")
//      .withInput("ARG3_y_")
//      .withOutput("RESULT_")
//      .build();
class NimTfBuilder {
  std::string graphDefBase64_;
  std::vector<std::string> inputNames_;
  std::vector<std::string> outputNames_;

 public:
  NimTfBuilder(const std::string& graphDefBase64)
      : graphDefBase64_(graphDefBase64) {}

  NimTfBuilder& withInput(const std::string& name) {
    inputNames_.push_back(name);
    return *this;
  }
  NimTfBuilder& withOutput(const std::string& name) {
    outputNames_.push_back(name);
    return *this;
  }

  NimTfRunner* build() {
    return new NimTfRunner(graphDefBase64_, inputNames_, outputNames_);
  }
};

#endif  // __NIMBLE_TENSORFLOW_H
