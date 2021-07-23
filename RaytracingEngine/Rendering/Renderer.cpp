#include "Renderer.h"

Renderer* Renderer::Get()
{
	static Renderer renderer;

	return &renderer;
}

void Renderer::Intialize()
{
}

void Renderer::Render()
{
}

Renderer::Renderer()
{
}

Renderer::~Renderer()
{
}

void Renderer::operator=(const Renderer& renderer)
{
	UNUSED(renderer);
}
