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
//    static NimTf_Runner& runner =
//      NimTf_Builder("0123456789abcdef",  // This string is a graphDefBase64.
//                    "0123456789abcdef")  // This string is a configBase64.
//        .NimTf_withInput("ARG1_a_")
//        .NimTf_withInput("ARG2_x_")
//        .NimTf_withInput("ARG3_y_")
//        .NimTf_withOutput("RESULT_")
//        .NimTf_build();
//
//    NimArr<1, double> result;
//    result.initialize(0.0, false, ARG2_x_.dimSize(0));
//    runner->NimTf_setInput(ARG1_a_);
//    runner->NimTf_setInput(ARG2_x_);
//    runner->NimTf_setInput(ARG3_y_);
//    runner->NimTf_run();
//    runner->NimTf_getOutput(result);
//    return result;
//  }
class NimTf_Runner {
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
  NimTf_Runner(const std::string& graphDefBase64,
               const std::string& configBase64,
               const std::vector<std::string>& inputNames,
               const std::vector<std::string>& outputName);
  ~NimTf_Runner();

  // These must be called in strict order.
  void NimTf_setInput(double& scalar);
  void NimTf_setInput(NimArrBase<double>& nimArr);
  void NimTf_run();
  void NimTf_getOutput(double& scalar);
  void NimTf_getOutput(NimArrBase<double>& nimArr);

 private:
  // Safety measures:
  NimTf_Runner(const NimTf_Runner&) { NIMERROR("Illegal"); }
  NimTf_Runner& operator=(const NimTf_Runner&) { NIMERROR("Illegal"); }
};

// Helper to construct static instances of NimTf_Runner.
// Example usage:
//  static NimTf_Runner& nimTf =
//     *NimTf_Builder("0123456789abcdef",  // The string is a graphDefBase64.
//      .NimTF_withConfig("0123456789abcdef")  // The string is a configBase64.
//      .NimTf_withInput("ARG1_a_")
//      .NimTf_withInput("ARG2_x_")
//      .NimTf_withInput("ARG3_y_")
//      .NimTf_withOutput("RESULT_")
//      .NimTf_build();
// Note that this uses the static pointer trick to avoid errors due to
// out-of-order destructor calls at program shutdown; doing so leaks memory but
// is common practice.
class NimTf_Builder {
  std::string graphDefBase64_;
  std::string configBase64_;
  std::vector<std::string> inputNames_;
  std::vector<std::string> outputNames_;

 public:
  NimTf_Builder(const std::string& graphDefBase64)
      : graphDefBase64_(graphDefBase64) {}

  NimTf_Builder& NimTf_withConfig(const std::string& configBase64) {
    configBase64_ = configBase64;
    return *this;
  }
  NimTf_Builder& NimTf_withInput(const std::string& name) {
    inputNames_.push_back(name);
    return *this;
  }
  NimTf_Builder& NimTf_withOutput(const std::string& name) {
    outputNames_.push_back(name);
    return *this;
  }

  NimTf_Runner* NimTf_build() {
    return new NimTf_Runner(graphDefBase64_, configBase64_, inputNames_,
                            outputNames_);
  }
};

#endif  // __NIMBLE_TENSORFLOW_H
