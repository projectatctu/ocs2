//
// Created by rgrandia on 17.02.20.
//

#include "ocs2_quadruped_interface/QuadrupedInterface.h"

#include <ocs2_core/misc/LoadData.h>
#include <ocs2_switched_model_interface/core/Rotations.h>
#include <ocs2_switched_model_interface/core/SwitchedModelStateEstimator.h>
#include <ocs2_switched_model_interface/foot_planner/FeetZDirectionPlanner.h>
#include <ocs2_switched_model_interface/foot_planner/cpg/SplineCPG.h>

namespace switched_model {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
QuadrupedInterface::QuadrupedInterface(const kinematic_model_t& kinematicModel, const ad_kinematic_model_t& adKinematicModel,
                                       const com_model_t& comModel, const ad_com_model_t& adComModel, const std::string& pathToConfigFolder)

    : kinematicModelPtr_(kinematicModel.clone()),
      adKinematicModelPtr_(adKinematicModel.clone()),
      comModelPtr_(comModel.clone()),
      adComModelPtr_(adComModel.clone()) {
  loadSettings(pathToConfigFolder + "/task.info");

  dynamicsPtr_.reset(new system_dynamics_t(*adKinematicModelPtr_, *adComModelPtr_, modelSettings_.recompileLibraries_));
  dynamicsDerivativesPtr_.reset(dynamicsPtr_->clone());
  constraintsPtr_.reset(new constraint_t(*adKinematicModelPtr_, *adComModelPtr_, logicRulesPtr_, modelSettings_));
  costFunctionPtr_.reset(new cost_function_t(*comModelPtr_, logicRulesPtr_, Q_, R_, QFinal_));
  operatingPointsPtr_.reset(new operating_point_t(*comModelPtr_, logicRulesPtr_));
  timeTriggeredRolloutPtr_.reset(new time_triggered_rollout_t(*dynamicsPtr_, rolloutSettings_));
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
void QuadrupedInterface::loadSettings(const std::string& pathToConfigFile) {
  slqSettings_.loadSettings(pathToConfigFile, "slq", true);
  rolloutSettings_.loadSettings(pathToConfigFile, "slq.rollout");
  mpcSettings_.loadSettings(pathToConfigFile, true);
  modelSettings_.loadSettings(pathToConfigFile, true);

  std::cerr << std::endl;

  // partitioning times
  scalar_t timeHorizon;
  size_t numPartitions;
  ocs2::loadData::loadPartitioningTimes(pathToConfigFile, timeHorizon, numPartitions, partitioningTimes_, true);

  // display
  std::cerr << "Time Partition: {";
  for (const auto& timePartition : partitioningTimes_) {
    std::cerr << timePartition << ", ";
  }
  std::cerr << "\b\b}" << std::endl;

  // initial state of the switched system
  Eigen::Matrix<scalar_t, RBD_STATE_DIM, 1> initRbdState;
  ocs2::loadData::loadEigenMatrix(pathToConfigFile, "initialRobotState", initRbdState);
  SwitchedModelStateEstimator switchedModelStateEstimator(*comModelPtr_);
  initialState_ = switchedModelStateEstimator.estimateComkinoModelState(initRbdState);

  // cost function components
  ocs2::loadData::loadEigenMatrix(pathToConfigFile, "Q", Q_);
  ocs2::loadData::loadEigenMatrix(pathToConfigFile, "R", R_);
  ocs2::loadData::loadEigenMatrix(pathToConfigFile, "Q_final", QFinal_);

  // costs over Cartesian velocities
  Eigen::Matrix<double, 12, 12> J_allFeet;
  for (int leg = 0; leg < 4; ++leg) {
    Eigen::Matrix<double, 6, 12> J_thisfoot = kinematicModelPtr_->baseToFootJacobianInBaseFrame(leg, getJointPositions(initRbdState));
    J_allFeet.block<3, 12>(3 * leg, 0) = J_thisfoot.bottomRows<3>();
  }
  R_.template block<12, 12>(12, 12) = (J_allFeet.transpose() * R_.template block<12, 12>(12, 12) * J_allFeet).eval();

  // load the mode sequence template
  std::cerr << std::endl;
  mode_sequence_template_t initialModeSequenceTemplate;
  loadModeSequenceTemplate(pathToConfigFile, "initialModeSequenceTemplate", initialModeSequenceTemplate, false);
  std::cerr << std::endl;
  if (initialModeSequenceTemplate.templateSubsystemsSequence_.empty()) {
    throw std::runtime_error("initialModeSequenceTemplate.templateSubsystemsSequence should have at least one entry.");
  }
  if (initialModeSequenceTemplate.templateSwitchingTimes_.size() != initialModeSequenceTemplate.templateSubsystemsSequence_.size() + 1) {
    throw std::runtime_error(
        "initialModeSequenceTemplate.templateSwitchingTimes size should be 1 + "
        "size_of(initialModeSequenceTemplate.templateSubsystemsSequence).");
  }

  size_array_t initSwitchingModes = initialModeSequenceTemplate.templateSubsystemsSequence_;
  initSwitchingModes.push_back(string2ModeNumber("STANCE"));

  auto& templateSwitchingTimes = initialModeSequenceTemplate.templateSwitchingTimes_;
  scalar_array_t initEventTimes = scalar_array_t(templateSwitchingTimes.begin() + 1, templateSwitchingTimes.end());

  // display
  std::cerr << "Initial Switching Modes: {";
  for (const auto& switchingMode : initSwitchingModes) {
    std::cerr << switchingMode << ", ";
  }
  std::cerr << "\b\b}" << std::endl;
  std::cerr << "Initial Event Times:     {";
  for (const auto& switchingtime : initEventTimes) {
    std::cerr << switchingtime << ", ";
  }
  if (!initEventTimes.empty()) {
    std::cerr << "\b\b}" << std::endl;
  } else {
    std::cerr << "}" << std::endl;
  }

  // load the mode sequence template
  std::cerr << std::endl;
  loadModeSequenceTemplate(pathToConfigFile, "defaultModeSequenceTemplate", defaultModeSequenceTemplate_, true);
  std::cerr << std::endl;

  // logic rule
  using cpg_t = SplineCPG<scalar_t>;
  using feet_z_planner_t = FeetZDirectionPlanner<scalar_t, cpg_t>;
  using feet_z_planner_ptr_t = std::shared_ptr<feet_z_planner_t>;
  feet_z_planner_ptr_t feetZPlannerPtr(new feet_z_planner_t(modelSettings_.swingLegLiftOff_, 1.0 /*swingTimeScale*/,
                                                            modelSettings_.liftOffVelocity_, modelSettings_.touchDownVelocity_));

  logicRulesPtr_ = std::shared_ptr<logic_rules_t>(new logic_rules_t(feetZPlannerPtr, modelSettings_.phaseTransitionStanceTime_));
  logicRulesPtr_->setModeSequence(initSwitchingModes, initEventTimes);
}

}  // namespace switched_model
