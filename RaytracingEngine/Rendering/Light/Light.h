#ifndef __LIGHT__
#define __LIGHT__

class Light {
public:
	Light(const Vector3f &p, const float i) : position(p), intensity(i) {}
	Vector3f position;
	float intensity;
};

#endif // __LIGHT__