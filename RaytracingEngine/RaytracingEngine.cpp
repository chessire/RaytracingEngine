#include "RaytracingEngine.h"
#include "GlobalConst.h"

#include <Eigen/Dense>

using namespace Eigen;

RaytracingEngine* RaytracingEngine::Get()
{
	static RaytracingEngine deepRT;
	return &deepRT;
}

void RaytracingEngine::Initialize()
{
	Eigen::initParallel();
}

RaytracingEngine::RaytracingEngine()
{
}

RaytracingEngine::~RaytracingEngine()
{
}

void RaytracingEngine::operator=(const RaytracingEngine& raytracingEngine)
{
	UNUSED(raytracingEngine);
}
