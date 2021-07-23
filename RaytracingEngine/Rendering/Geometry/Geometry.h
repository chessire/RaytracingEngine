#ifndef __GEOMETRY__

#include "../Material/Material.h"

class Geometry
{
public:
	Geometry(const Material& m) : material(m) {}
	virtual ~Geometry();

	const Material& GetMaterial() const { return material; }

private:
	Material material;
};

#endif // __GEOMETRY__
