#include "Scene.h"

Scene* Scene::Get()
{
	static Scene scene;
	return &scene;
}

void Scene::Intialize()
{
}

Scene::Scene()
{
}

Scene::~Scene()
{
}

void Scene::operator=(const Scene& scene)
{
	UNUSED(scene);
}
