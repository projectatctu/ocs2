#include "ocs2_double_slit_example/DoubleSlitInterface.h"
#include "ocs2_double_slit_example/ros_comm/MRT_ROS_Double_Slit.h"
#include "ocs2_double_slit_example/ros_comm/MRT_ROS_Dummy_Double_Slit.h"
#include "ocs2_double_slit_example/definitions.h"

using namespace ocs2;
using namespace double_slit;

int main(int argc, char **argv)
{
	// task file
	if (argc <= 1) throw std::runtime_error("No task file specified. Aborting.");
	std::string taskFileFolderName = std::string(argv[1]);

	// double_slitInterface
	DoubleSlitInterface double_slitInterface(taskFileFolderName);
	double rollout_dt;
	loadScalar(double_slitInterface.taskFile_, "pathIntegral.rollout_settings.minTimeStep", rollout_dt);

	typedef MRT_ROS_Double_Slit mrt_t;
	typedef mrt_t::BASE::Ptr mrt_base_ptr_t;
	typedef mrt_t::scalar_t scalar_t;
	typedef mrt_t::system_observation_t system_observation_t;

	mrt_base_ptr_t mrtPtr(new mrt_t("double_slit"));

	// Dummy double_slit
	MRT_ROS_Dummy_Linear_System dummyDoubleSlit(
			mrtPtr,
			double_slitInterface.mpcSettings().mrtDesiredFrequency_,
			double_slitInterface.mpcSettings().mpcDesiredFrequency_,
			double_slitInterface.getDynamicsPtr().get(),
			Rollout_Settings(1e-9, 1e-6, 5000, rollout_dt, IntegratorType::EULER, false, true));

	dummyDoubleSlit.launchNodes(argc, argv);

	// Run dummy
	MRT_ROS_Dummy_Linear_System::system_observation_t initObservation;
	double_slitInterface.getInitialState(initObservation.state());

	MRT_ROS_Dummy_Linear_System::cost_desired_trajectories_t initCostDesiredTraj;
	initCostDesiredTraj.desiredTimeTrajectory().push_back(0.0);
	initCostDesiredTraj.desiredInputTrajectory().push_back(DoubleSlitInterface::input_vector_t::Zero());
	initCostDesiredTraj.desiredStateTrajectory().push_back(DoubleSlitInterface::state_vector_t::Zero());

	dummyDoubleSlit.run(initObservation, initCostDesiredTraj);

	// Successful exit
	return 0;
}


