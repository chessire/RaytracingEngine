#ifndef __BVH__
#define __BVH__

#include "GlobalConst.h"

class BVH
{
public:
	static BVH* Get();

	void Intialize();

private:
	BVH();
	~BVH();
	void operator=(const BVH& renderer);
};

#endif // __BVH__