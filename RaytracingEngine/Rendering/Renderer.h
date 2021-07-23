#ifndef __RENDERER__
#define __RENDERER__

#include "GlobalConst.h"

class Renderer
{
public:
	static Renderer* Get();

	void Intialize();

	void Render();

private:
	Renderer();
	~Renderer();
	void operator=(const Renderer& renderer);
};

#endif //__RENDERER__