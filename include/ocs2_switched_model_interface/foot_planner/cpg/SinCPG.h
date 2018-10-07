/*
 * SinCPG.h
 *
 *  Created on: Jul 5, 2016
 *      Author: farbod
 */

#ifndef SINCPG_H_
#define SINCPG_H_

#include <cmath>
#include <memory>

#include "ocs2_switched_model_interface/foot_planner/cpg/CPG_BASE.h"

namespace switched_model {

template <typename scalar_t = double>
class SinCPG : public CPG_BASE<scalar_t>
{
public:
	typedef std::shared_ptr<SinCPG<scalar_t>> Ptr;

	typedef CPG_BASE<scalar_t> Base;

	SinCPG(const scalar_t& swingLegLiftOff = 0.15, const scalar_t& swingTimeScale = 1.0)
	: Base(swingLegLiftOff, swingTimeScale)
	{}

	~SinCPG() {}

	void setConstant() final {
		Base::set(0.0, 1.0, 0.0);
		Base::maxHight_ = 0.0;
	}

	scalar_t calculatePosition(const scalar_t& time) const final {
		const scalar_t omega = M_PI/(Base::finalTime_-Base::startTime_);
		scalar_t position = Base::maxHight_ * pow(sin( omega*(time-Base::startTime_) ), 2);
		return position;
	}

	scalar_t calculateVelocity(const scalar_t& time) const final {
		const scalar_t omega = M_PI/(Base::finalTime_-Base::startTime_);
		scalar_t velocity = 2*omega * Base::maxHight_ * cos( omega*(time-Base::startTime_) ) *
				sin( omega*(time-Base::startTime_) );
		return velocity;
	}

	scalar_t calculateStartTimeDerivative(const scalar_t& time) const final {
		throw std::runtime_error("calculateStartTimeDerivative method is not implemented.");
	}

	scalar_t calculateFinalTimeDerivative(const scalar_t& time) const final {
		throw std::runtime_error("calculateFinalTimeDerivative method is not implemented.");
	}

private:

};


}  // end of switched_model namespace


#endif /* SINCPG_H_ */
