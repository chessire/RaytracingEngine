#include "DeepRT.h"

#include <Eigen/Core>
#include "Rendering/Renderer.h"
#include "Rendering/BVH/BVH.h"
#include "Rendering/Scene/Scene.h"

using namespace Eigen;

DeepRT* DeepRT::Get()
{
	static DeepRT deepRT;
	return &deepRT;
}

void DeepRT::Initialize()
{
	Eigen::initParallel();

	Renderer::Get()->Intialize();
	BVH::Get()->Intialize();
	Scene::Get()->Intialize();
}

DeepRT::DeepRT()
{
}

DeepRT::~DeepRT()
{
}

void DeepRT::operator=(const DeepRT& deeprt)
{
	UNUSED(deeprt);
}
