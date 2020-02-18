//
// Created by rgrandia on 17.02.20.
//

#pragma once

#include <ros/node_handle.h>

#include "QuadrupedLoopshapingInterface.h"

namespace switched_model_loopshaping {

void quadrupedLoopshapingMpcNode(ros::NodeHandle& nodeHandle, const QuadrupedLoopshapingInterface& quadrupedInterface);

}  // namespace switched_model_loopshaping