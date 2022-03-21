/******************************************************************************
Copyright (c) 2022, Farbod Farshidian. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

 * Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************/

#pragma once

#include <ocs2_mpcnet/MpcnetDefinitionBase.h>

namespace ocs2 {
namespace legged_robot {

/**
 * MPC-Net definitions for legged robot.
 */
class LeggedRobotMpcnetDefinition : public ocs2::mpcnet::MpcnetDefinitionBase {
 public:
  /**
   * Constructor.
   * @param [in] defaultState : Default state.
   */
  LeggedRobotMpcnetDefinition(const vector_t& defaultState) : defaultState_(defaultState) {}

  /**
   * Default destructor.
   */
  ~LeggedRobotMpcnetDefinition() override = default;

  /**
   * @see MpcnetDefinitionBase::getGeneralizedTime
   */
  vector_t getGeneralizedTime(scalar_t t, const ModeSchedule& modeSchedule) override;

  /**
   * @see MpcnetDefinitionBase::getRelativeState
   */
  vector_t getRelativeState(scalar_t t, const vector_t& x, const TargetTrajectories& targetTrajectories) override;

  /**
   * @see MpcnetDefinitionBase::getInputTransformation
   */
  matrix_t getInputTransformation(scalar_t t, const vector_t& x) override;

  /**
   * @see MpcnetDefinitionBase::validState
   */
  bool validState(const vector_t& x) override;

 private:
  vector_t defaultState_;
};

}  // namespace legged_robot
}  // namespace ocs2
