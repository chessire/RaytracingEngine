#ifndef __SCENE__
#define __SCENE__

#include "GlobalConst.h"

class Scene
{
public:
	static Scene* Get();

	void Intialize();

private:
	Scene();
	~Scene();
	void operator=(const Scene& scene);
};

#endif // __SCENE__