#include "BVH.h"

BVH* BVH::Get()
{
	static BVH bvh;
	return &bvh;
}

void BVH::Intialize()
{
}

BVH::BVH()
{
}

BVH::~BVH()
{
}

void BVH::operator=(const BVH& bvh)
{
	UNUSED(bvh);
}
