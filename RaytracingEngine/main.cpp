#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <Eigen/Core>

#include "GlobalConst.h"
#include "DeepRT.h"

// test

using namespace Eigen;

struct Light {
    Light(const Vector3f &p, const float i) : position(p), intensity(i) {}
    Vector3f position;
    float intensity;
};

struct Material {
    Material(const float r, const Vector4f &a, const Vector3f &color, const float spec) : refractive_index(r), albedo(a), diffuse_color(color), specular_exponent(spec) {}
    Material() : refractive_index(1), albedo(1,0,0,0), diffuse_color(), specular_exponent() {}
    float refractive_index;
    Vector4f albedo;
    Vector3f diffuse_color;
    float specular_exponent;
};

class SdfModel
{
	Material material;
public:
	SdfModel(const Material& m) : material(m) { }
	virtual ~SdfModel() {}

	virtual float sdf(const Vector3f& point) const = 0;
	virtual bool try_get_normal(const Vector3f& point, Vector3f& n) const = 0;
	const Material& get_material() const { return material; }
};

struct Sphere : public SdfModel {
	Vector3f center;
	float radius;

	Sphere(const Vector3f &c, const float r, const Material &m) : SdfModel(m), center(c), radius(r) {}

	float sdf(const Vector3f& point) const final
	{
		return (point - center).norm() - radius;
	}

	bool try_get_normal(const Vector3f& point, Vector3f& n) const final
	{
		n = (point - center).normalized();

		return true;
	}
};

Vector4f sin_vec(Vector4f v)
{
	return Vector4f(sinf(v[0]), sinf(v[1]), sinf(v[2]), sinf(v[3]));
}

struct StanfordBunny : public SdfModel {
	Vector3f center;

	StanfordBunny(const Vector3f &c, const Material &m) : SdfModel(m), center(c) {}

	float sdf(const Vector3f& point) const final
	{
		Vector3f p = point - center;

		//sdf is undefined outside the unit sphere, uncomment to witness the abominations
		float p_norm = p.norm();
		if (p_norm > 1.) {
			return p_norm - 0.8f;
		}

		//neural networks can be really compact... when they want to be
		Vector4f f00 = sin_vec(p[1] * Vector4f(-3.02, 1.95, -3.42, -.60) + p[2] * Vector4f(3.08, .85, -2.25, -.24) - p[0] * Vector4f(-.29, 1.16, -3.74, 2.89) + Vector4f(-.71, 4.50, -3.24, -3.50));
		Vector4f f01 = sin_vec(p[1] * Vector4f(-.40, -3.61, 3.23, -.14) + p[2] * Vector4f(-.36, 3.64, -3.91, 2.66) - p[0] * Vector4f(2.90, -.54, -2.75, 2.71) + Vector4f(7.02, -5.41, -1.12, -7.41));
		Vector4f f02 = sin_vec(p[1] * Vector4f(-1.77, -1.28, -4.29, -3.20) + p[2] * Vector4f(-3.49, -2.81, -.64, 2.79) - p[0] * Vector4f(3.15, 2.14, -3.85, 1.83) + Vector4f(-2.07, 4.49, 5.33, -2.17));
		Vector4f f03 = sin_vec(p[1] * Vector4f(-.49, .68, 3.05, .42) + p[2] * Vector4f(-2.87, .78, 3.78, -3.41) - p[0] * Vector4f(-2.65, .33, .07, -.64) + Vector4f(-3.24, -5.90, 1.14, -4.71));
		Vector4f f10 = sin_vec((Matrix4f()<<-.34, .06, -.59, -.76, .10, -.19, -.12, .44, .64, -.02, -.26, .15, -.16, .21, .91, .15).finished().transpose()*f00 +
			(Matrix4f()<<.01, .54, -.77, .11, .06, -.14, .43, .51, -.18, .08, .39, .20, .33, -.49, -.10, .19).finished().transpose()*f01 +
			(Matrix4f()<<.27, .22, .43, .53, .18, -.17, .23, -.64, -.14, .02, -.10, .16, -.13, -.06, -.04, -.36).finished().transpose()*f02 +
			(Matrix4f()<<-.13, .29, -.29, .08, 1.13, .02, -.83, .32, -.32, .04, -.31, -.16, .14, -.03, -.20, .39).finished().transpose()*f03 +
			Vector4f(.73, -4.28, -1.56, -1.80)) / 1.0 + f00;
		Vector4f f11 = sin_vec((Matrix4f()<<-1.11, .55, -.12, -1.00, .16, .15, -.30, .31, -.01, .01, .31, -.42, -.29, .38, -.04, .71).finished().transpose()*f00 +
			(Matrix4f()<<.96, -.02, .86, .52, -.14, .60, .44, .43, .02, -.15, -.49, -.05, -.06, -.25, -.03, -.22).finished().transpose()*f01 +
			(Matrix4f()<<.52, .44, -.05, -.11, -.56, -.10, -.61, -.40, -.04, .55, .32, -.07, -.02, .28, .26, -.49).finished().transpose()*f02 +
			(Matrix4f()<<.02, -.32, .06, -.17, -.59, .00, -.24, .60, -.06, .13, -.21, -.27, -.12, -.14, .58, -.55).finished().transpose()*f03 +
			Vector4f(-2.24, -3.48, -.80, 1.41)) / 1.0 + f01;
		Vector4f f12 = sin_vec((Matrix4f()<<.44, -.06, -.79, -.46, .05, -.60, .30, .36, .35, .12, .02, .12, .40, -.26, .63, -.21).finished().transpose()*f00 +
			(Matrix4f()<<-.48, .43, -.73, -.40, .11, -.01, .71, .05, -.25, .25, -.28, -.20, .32, -.02, -.84, .16).finished().transpose()*f01 +
			(Matrix4f()<<.39, -.07, .90, .36, -.38, -.27, -1.86, -.39, .48, -.20, -.05, .10, -.00, -.21, .29, .63).finished().transpose()*f02 +
			(Matrix4f()<<.46, -.32, .06, .09, .72, -.47, .81, .78, .90, .02, -.21, .08, -.16, .22, .32, -.13).finished().transpose()*f03 +
			Vector4f(3.38, 1.20, .84, 1.41)) / 1.0 + f02;
		Vector4f f13 = sin_vec((Matrix4f()<<-.41, -.24, -.71, -.25, -.24, -.75, -.09, .02, -.27, -.42, .02, .03, -.01, .51, -.12, -1.24).finished().transpose()*f00 +
			(Matrix4f()<<.64, .31, -1.36, .61, -.34, .11, .14, .79, .22, -.16, -.29, -.70, .02, -.37, .49, .39).finished().transpose()*f01 +
			(Matrix4f()<<.79, .47, .54, -.47, -1.13, -.35, -1.03, -.22, -.67, -.26, .10, .21, -.07, -.73, -.11, .72).finished().transpose()*f02 +
			(Matrix4f()<<.43, -.23, .13, .09, 1.38, -.63, 1.57, -.20, .39, -.14, .42, .13, -.57, -.08, -.21, .21).finished().transpose()*f03 +
			Vector4f(-.34, -3.28, .43, -.52)) / 1.0 + f03;
		f00 = sin_vec((Matrix4f()<<-.72, .23, -.89, .52, .38, .19, -.16, -.88, .26, -.37, .09, .63, .29, -.72, .30, -.95).finished().transpose()*f10 +
			(Matrix4f()<<-.22, -.51, -.42, -.73, -.32, .00, -1.03, 1.17, -.20, -.03, -.13, -.16, -.41, .09, .36, -.84).finished().transpose()*f11 +
			(Matrix4f()<<-.21, .01, .33, .47, .05, .20, -.44, -1.04, .13, .12, -.13, .31, .01, -.34, .41, -.34).finished().transpose()*f12 +
			(Matrix4f()<<-.13, -.06, -.39, -.22, .48, .25, .24, -.97, -.34, .14, .42, -.00, -.44, .05, .09, -.95).finished().transpose()*f13 +
			Vector4f(.48, .87, -.87, -2.06)) / 1.4 + f10;
		f01 = sin_vec((Matrix4f()<<-.27, .29, -.21, .15, .34, -.23, .85, -.09, -1.15, -.24, -.05, -.25, -.12, -.73, -.17, -.37).finished().transpose()*f10 +
			(Matrix4f()<<-1.11, .35, -.93, -.06, -.79, -.03, -.46, -.37, .60, -.37, -.14, .45, -.03, -.21, .02, .59).finished().transpose()*f11 +
			(Matrix4f()<<-.92, -.17, -.58, -.18, .58, .60, .83, -1.04, -.80, -.16, .23, -.11, .08, .16, .76, .61).finished().transpose()*f12 +
			(Matrix4f()<<.29, .45, .30, .39, -.91, .66, -.35, -.35, .21, .16, -.54, -.63, 1.10, -.38, .20, .15).finished().transpose()*f13 +
			Vector4f(-1.72, -.14, 1.92, 2.08)) / 1.4 + f11;
		f02 = sin_vec((Matrix4f()<<1.00, .66, 1.30, -.51, .88, .25, -.67, .03, -.68, -.08, -.12, -.14, .46, 1.15, .38, -.10).finished().transpose()*f10 +
			(Matrix4f()<<.51, -.57, .41, -.09, .68, -.50, -.04, -1.01, .20, .44, -.60, .46, -.09, -.37, -1.30, .04).finished().transpose()*f11 +
			(Matrix4f()<<.14, .29, -.45, -.06, -.65, .33, -.37, -.95, .71, -.07, 1.00, -.60, -1.68, -.20, -.00, -.70).finished().transpose()*f12 +
			(Matrix4f()<<-.31, .69, .56, .13, .95, .36, .56, .59, -.63, .52, -.30, .17, 1.23, .72, .95, .75).finished().transpose()*f13 +
			Vector4f(-.90, -3.26, -.44, -3.11)) / 1.4 + f12;
		f03 = sin_vec((Matrix4f()<<.51, -.98, -.28, .16, -.22, -.17, -1.03, .22, .70, -.15, .12, .43, .78, .67, -.85, -.25).finished().transpose()*f10 +
			(Matrix4f()<<.81, .60, -.89, .61, -1.03, -.33, .60, -.11, -.06, .01, -.02, -.44, .73, .69, 1.02, .62).finished().transpose()*f11 +
			(Matrix4f()<<-.10, .52, .80, -.65, .40, -.75, .47, 1.56, .03, .05, .08, .31, -.03, .22, -1.63, .07).finished().transpose()*f12 +
			(Matrix4f()<<-.18, -.07, -1.22, .48, -.01, .56, .07, .15, .24, .25, -.09, -.54, .23, -.08, .20, .36).finished().transpose()*f13 +
			Vector4f(-1.11, -4.28, 1.02, -.23)) / 1.4 + f13;

		return f00.dot(Vector4f(.09, .12, -.07, -.03)) + f01.dot(Vector4f(-.04, .07, -.08, .05)) +
			f02.dot(Vector4f(-.01, .06, -.02, .07)) + f03.dot(Vector4f(-.05, .07, .03, .04));
	}

	bool try_get_normal(const Vector3f& point, Vector3f& n) const final
	{
		Vector3f other = point - Vector3f(0.001f, 0.f, 0.f);

		float sdf_point = sdf(point);
		float sdf_other = sdf(other);
		n = Vector3f(sdf_point - sdf_other, sdf_point - sdf_other, sdf_point - sdf_other).normalized();
		return true;
	}
};

void fresnel(const Vector3f &I, const Vector3f &N, const float &ior, float &kr)
{
	float cosi = std::clamp(I.dot(N), -1.f, 1.f);
	float etai = 1, etat = ior;
	if (cosi > 0) { std::swap(etai, etat); }
	// Compute sini using Snell's law
	float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
	// Total internal reflection
	if (sint >= 1) {
		kr = 1;
	}
	else {
		float cost = sqrtf(std::max(0.f, 1 - sint * sint));
		cosi = fabsf(cosi);
		float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
		float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
		kr = (Rs * Rs + Rp * Rp) / 2;
	}
	// As a consequence of the conservation of energy, transmittance is given by:
	// kt = 1 - kr;
}

Vector3f reflect(const Vector3f &I, const Vector3f &N) {
	return I - N * 2.f*(I.dot(N));
}

Vector3f refract(const Vector3f &I, const Vector3f &N, const float eta_t, const float eta_i=1.f) { // Snell's law
    float cosi = - std::max(-1.f, std::min(1.f, I.dot(N)));
    if (cosi<0) return refract(I, -N, eta_i, eta_t); // if the ray comes from the inside the object, swap the air and the media
    float eta = eta_i / eta_t;
    float k = 1 - eta*eta*(1 - cosi*cosi);
    return k<0 ? Vector3f(1,0,0) : I*eta + N*(eta*cosi - sqrtf(k)); // k<0 = total reflection, no ray to refract. I refract it anyways, this has no physical meaning
}

float scene_sdf(const Vector3f& point, const std::vector<const SdfModel*> &models, const SdfModel*& hit_model)
{
	float min_dist = MAX_DISTANCE;
	hit_model = nullptr;
	for (const auto& model : models)
	{
		float dist = model->sdf(point);
		if (dist < EPSILON)
			continue;

		if (min_dist > dist)
		{
			min_dist = dist;
			hit_model = model;
		}
	}
	return min_dist;
}

bool ray_marching(const Vector3f &orig, const Vector3f &dir, const std::vector<const SdfModel*> &models, Vector3f &hit, Vector3f &N, Material &material)
{
	float depth = EPSILON;
	for (int i = 0; i < MAX_MARCHING_STEPS; ++i)
	{
		const SdfModel* hit_model = nullptr;
		float dist = scene_sdf(orig + dir * depth, models, hit_model);

		if (hit_model == nullptr)
			break;

		depth += dist;
		if (dist < EPSILON)
		{
			hit = orig + dir * depth;
			if (hit_model->try_get_normal(hit, N) == false)
				std::cout << "normal bug!" << std::endl;
			material = hit_model->get_material();
			return true;
		}

		// for checkboard(temp comment)
		if (depth >= MAX_DISTANCE)
			break;
	}

	float checkerboard_dist = std::numeric_limits<float>::max();
	if (fabs(dir.y()) > EPSILON) {
		float d = -(orig.y() + 4) / dir.y(); // the checkerboard plane has equation y = -4
		Vector3f pt = orig + dir * d;
		if (d > 0 && fabs(pt.x()) < 10 && pt.z() < -10 && pt.z() > -30) {
			checkerboard_dist = d;
			hit = pt;
			N = Vector3f(0, 1, 0);
			material.diffuse_color = (int(.5*hit.x() + 1000) + int(.5*hit.z())) & 1 ? Vector3f(.3, .3, .3) : Vector3f(.3, .2, .1);
		}
	}

	return std::min(MAX_DISTANCE, checkerboard_dist)<1000;
}

Vector3f cast_ray(const Vector3f &orig, const Vector3f &dir, const std::vector<const SdfModel*> &models, const std::vector<Light> &lights, size_t depth=0) {
    Vector3f point, N;
    Material material;

    if (depth>4 || !ray_marching(orig, dir, models, point, N, material)) {
        return Vector3f(0.2, 0.7, 0.8); // background color
	}

	Vector3f refract_color(0.f, 0.f, 0.f);
	// compute fresnelt
	float kr;
	fresnel(dir, N, material.refractive_index, kr);
	// compute refraction if it is not a case of total internal reflection
	if (kr < 1) {
		Vector3f refract_dir = refract(dir, N, material.refractive_index).normalized();
		Vector3f refract_orig = refract_dir.dot(N) < 0.f ? (point - N * EPSILON).eval() : (point + N * EPSILON).eval();
		refract_color = cast_ray(refract_orig, refract_dir, models, lights, depth + 1);
	}

	Vector3f reflect_dir = reflect(dir, N).normalized();
	Vector3f reflect_orig = reflect_dir.dot(N) < 0.f ? (point - N * EPSILON).eval() : (point + N * EPSILON).eval(); // offset the original point to avoid occlusion by the object itself
	Vector3f reflect_color = cast_ray(reflect_orig, reflect_dir, models, lights, depth + 1);

    float diffuse_light_intensity = 0, specular_light_intensity = 0;
    for (size_t i=0; i<lights.size(); i++) {
        Vector3f light_dir      = (lights[i].position - point).normalized();
        float light_distance = (lights[i].position - point).norm();

		Vector3f shadow_orig = light_dir.dot(N) < 0.f ? (point - N * EPSILON).eval() : (point + N * EPSILON).eval(); // checking if the point lies in the shadow of the lights[i]
        Vector3f shadow_pt, shadow_N;
        Material tmpmaterial;
        if (ray_marching(shadow_orig, light_dir, models, shadow_pt, shadow_N, tmpmaterial) && (shadow_pt-shadow_orig).norm() < light_distance)
            continue;

        diffuse_light_intensity  += lights[i].intensity * std::max(0.f, light_dir.dot(N));
        specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir, N).dot(dir)), material.specular_exponent)*lights[i].intensity;
    }

    return material.diffuse_color * diffuse_light_intensity * material.albedo[0] + Vector3f(1., 1., 1.)*specular_light_intensity * material.albedo[1] + reflect_color * material.albedo[2] * kr + refract_color * material.albedo[3] * (1 - kr);
}

void render(const std::vector<const SdfModel*> &models, const std::vector<Light> &lights) {
    const int   width    = 1024;
    const int   height   = 768;
    const float fov      = M_PI/3.;
    std::vector<Vector3f> framebuffer(width*height);

//#pragma omp parallel for
	for (size_t j = 0; j < height; j++) { // actual rendering loop
        for (size_t i = 0; i<width; i++) {
            float dir_x =  (i + 0.5) -  width/2.;
            float dir_y = -(j + 0.5) + height/2.;    // this flips the image at the same time
            float dir_z = -height/(2.*tan(fov/2.));
            framebuffer[i+j*width] = cast_ray(Vector3f(0.f, 0.f, 0.f), Vector3f(dir_x, dir_y, dir_z).normalized(), models, lights);
        }
    }

    std::ofstream ofs; // save the framebuffer to file
    ofs.open("./out.ppm",std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height*width; ++i) {
        Vector3f &c = framebuffer[i];
		float max = std::max(c[0], std::max(c[1], c[2]));
		if (max > 1) c = c * (1. / max);
        for (size_t j = 0; j<3; j++) {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
}

int main() {
	DeepRT::Get()->Initialize();

    Material      ivory(0.0, Vector4f(0.6,  0.3, 0.1, 0.0), Vector3f(0.4, 0.4, 0.3),   50.);
    Material      glass(1.5, Vector4f(0.0,  0.5, 0.1, 0.8), Vector3f(0.6, 0.7, 0.8),  125.);
    Material red_rubber(0.0, Vector4f(0.9,  0.1, 0.0, 0.0), Vector3f(0.3, 0.1, 0.1),   10.);
    Material     mirror(0.0, Vector4f(0.0, 10.0, 0.8, 0.0), Vector3f(1.0, 1.0, 1.0), 1425.);

    std::vector<const SdfModel*> models;
    models.push_back(new Sphere(Vector3f(-3,    0,   -16), 2,      ivory));
	models.push_back(new Sphere(Vector3f(-1.0, -1.5, -12), 2,      glass));
	models.push_back(new Sphere(Vector3f( 1.5, -0.5, -18), 3, red_rubber));
	models.push_back(new Sphere(Vector3f(7, 5, -18), 4, mirror));
	models.push_back(new StanfordBunny(Vector3f(0, 3, -10), red_rubber));

    std::vector<Light>  lights;
    lights.push_back(Light(Vector3f(-20, 20,  20), 1.5));
    lights.push_back(Light(Vector3f( 30, 50, -25), 1.8));
	lights.push_back(Light(Vector3f(30, 20, 30), 1.7));

	render(models, lights);

	for (auto model : models)
	{
		delete model;
	}
	models.clear();

    return 0;
}

