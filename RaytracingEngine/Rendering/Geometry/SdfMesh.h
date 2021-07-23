#ifndef __SDF_MESH__

#include "Geometry.h"

class SdfMesh : public Geometry
{
public:
	SdfModel(const Material& m) : Geometry(m) { }
	virtual ~SdfModel() {}

	virtual float Sdf(const Vector3f& point) const = 0;
	virtual bool TryGetNormal(const Vector3f& point, Vector3f& n) const = 0;
};

#endif // __SDF_MESH__
