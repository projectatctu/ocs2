/******************************************************************************
Copyright (c) 2017, Farbod Farshidian. All rights reserved.

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

#include <utility>
#include <vector>

#include <ocs2_core/Types.h>
#include <ocs2_core/control/ControllerBase.h>

#include "Rollout_Settings.h"

namespace ocs2 {

/**
 * This class is an interface class for forward rollout of the system dynamics.
 */
class RolloutBase {
 public:
  using time_interval_t = std::pair<scalar_t, scalar_t>;
  using time_interval_array_t = std::vector<time_interval_t>;

  /**
   * Default constructor.
   *
   * @param [in] rolloutSettings: The rollout settings.
   */
  RolloutBase(Rollout_Settings rolloutSettings = Rollout_Settings()) : rolloutSettings_(std::move(rolloutSettings)) {}

  /**
   * Default destructor.
   */
  virtual ~RolloutBase() = default;

  /**
   * Returns the rollout settings.
   *
   * @return The rollout settings.
   */
  Rollout_Settings& settings() { return rolloutSettings_; }

  /**
   * Returns the rollout settings.
   *
   * @return The rollout settings.
   */
  const Rollout_Settings& settings() const { return rolloutSettings_; }

  /**
   * The kills the integrator inside the rollout.
   */
  virtual void abortRollout() {}

  /**
   * The enables the integrator inside the rollout to start again.
   */
  virtual void reactivateRollout() {}

  /**
   * Returns pointer to the class.
   *
   * @return A raw pointer to the class.
   */
  virtual RolloutBase* clone() const = 0;

  /**
   * Forward integrate the system dynamics with given controller. It uses the given control policies and initial state,
   * to integrate the system dynamics in time period [initTime, finalTime].
   *
   * @param [in] initTime: The initial time.
   * @param [in] initState: The initial state.
   * @param [in] finalTime: The final time.
   * @param [in] controller: control policy.
   * @param [in] eventTimes: The sorted event times array which can cover time beyond initTime and finalTime {s_0 < ... < s_{n-1}}.
   * @param [out] timeTrajectory: The time trajectory stamp.
   * @param [out] postEventIndicesStock: Indices containing past-the-end index of events trigger.
   * @param [out] stateTrajectory: The state trajectory.
   * @param [out] inputTrajectory: The control input trajectory.
   *
   * @return The final state (state jump is considered if it took place)
   */
  vector_t run(scalar_t initTime, const vector_t& initState, scalar_t finalTime, ControllerBase* controller,
               const scalar_array_t& eventTimes, scalar_array_t& timeTrajectory, size_array_t& postEventIndicesStock,
               vector_array_t& stateTrajectory, vector_array_t& inputTrajectory);

  /**
   * Prints out the rollout.
   *
   * @param [in] timeTrajectory: The time trajectory stamp.
   * @param [in] postEventIndicesStock: An array of the post-event indices.
   * @param [in] stateTrajectory: The state trajectory.
   * @param [in] inputTrajectory: The control input trajectory.
   */
  static void display(const scalar_array_t& timeTrajectory, const size_array_t& postEventIndicesStock,
                      const vector_array_t& stateTrajectory, const vector_array_t* const inputTrajectory);

 protected:
  /**
   * Forward integrate the system dynamics with given controller. It uses the given control policies and initial state,
   * to integrate the system dynamics in time period [initTime, finalTime].
   *
   * @param [in] timeIntervalArray: An array of the rollout's start and final times.
   * @param [in] initState: The initial state.
   * @param [in] controller: control policy.
   * @param [out] timeTrajectory: The time trajectory stamp.
   * @param [out] postEventIndicesStock: Indices containing past-the-end index of events trigger.
   * @param [out] stateTrajectory: The state trajectory.
   * @param [out] inputTrajectory: The control input trajectory.
   *
   * @return The final state (state jump is considered if it took place)
   */
  virtual vector_t runImpl(time_interval_array_t timeIntervalArray, const vector_t& initState, ControllerBase* controller,
                           scalar_array_t& timeTrajectory, size_array_t& postEventIndicesStock, vector_array_t& stateTrajectory,
                           vector_array_t& inputTrajectory) = 0;

  /**
   * Checks for the numerical stability if Rollout_Settings::checkNumericalStability_ is true.
   *
   * @param [in] timeTrajectory: The time trajectory stamp.
   * @param [in] postEventIndicesStock: Indices containing past-the-end index of events trigger.
   * @param [in] stateTrajectory: The state trajectory.
   * @param [in] inputTrajectory: The control input trajectory.
   */
  void checkNumericalStability(ControllerBase* controller, const scalar_array_t& timeTrajectory, const size_array_t& postEventIndicesStock,
                               const vector_array_t& stateTrajectory, const vector_array_t& inputTrajectory) const;

 private:
  Rollout_Settings rolloutSettings_;
};

}  // namespace ocs2
