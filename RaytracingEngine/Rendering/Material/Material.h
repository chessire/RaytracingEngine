#ifndef __MATERIAL__

class Material
{
public:
	Material(const float r, const Vector4f &a, const Vector3f &color, const float spec) : refractive_index(r), albedo(a), diffuse_color(color), specular_exponent(spec) {}
	Material() : refractive_index(1), albedo(1, 0, 0, 0), diffuse_color(), specular_exponent() {}

	float refractive_index;
	Vector4f albedo;
	Vector3f diffuse_color;
	float specular_exponent;
};

#endif // __MATERIAL__
