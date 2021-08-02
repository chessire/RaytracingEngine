#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "Math/Math.h"
#include "RaytracingEngine.h"

struct Light {
    Light(const Vector3f& p, const float i) : position(p), intensity(i) {}
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

class Model
{
	Material material;
public:
	Model(const Material& m) : material(m) { }
	virtual ~Model() {}

	virtual bool RayIntersect(const Vector3f& orig, const Vector3f& dir, float& t0) const = 0;
	virtual Vector3f GetNormal(Vector3f point) const = 0;
	const Material& GetMaterial() const { return material; }
};

struct Sphere : public Model {
	Vector3f center;
	float radius;

	Sphere(const Vector3f &c, const float r, const Material &m) : Model(m), center(c), radius(r) {}

	bool RayIntersect(const Vector3f& orig, const Vector3f& dir, float& t0) const final
	{
		Vector3f L = center - orig;
		float tca = L.dot(dir);
		float d2 = L.dot(L) - tca * tca;
		if (d2 > radius * radius) return false;
		float thc = sqrtf(radius * radius - d2);
		t0 = tca - thc;
		float t1 = tca + thc;
		if (t0 < 0) t0 = t1;
		if (t0 < 0) return false;
		return true;
	}

	Vector3f GetNormal(Vector3f point) const final
	{
		return (point - center).normalized();
	}
};

void Fresnel(const Vector3f &I, const Vector3f &N, const float &ior, float &kr)
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

Vector3f Reflect(const Vector3f &I, const Vector3f &N) {
	return I - N * 2.f*(I.dot(N));
}

Vector3f Refract(const Vector3f &I, const Vector3f &N, const float eta_t, const float eta_i=1.f) { // Snell's law
    float cosi = - std::max(-1.f, std::min(1.f, I.dot(N)));
    if (cosi<0) return Refract(I, -N, eta_i, eta_t); // if the ray comes from the inside the object, swap the air and the media
    float eta = eta_i / eta_t;
    float k = 1 - eta*eta*(1 - cosi*cosi);
    return k<0 ? Vector3f(1,0,0) : I*eta + N*(eta*cosi - sqrtf(k)); // k<0 = total reflection, no ray to refract. I refract it anyways, this has no physical meaning
}

float SceneIntersect(const Vector3f& orig, const Vector3f& dir, const std::vector<const Model*>& models, Vector3f& hit, Vector3f& N, Material& material)
{
	float spheres_dist = std::numeric_limits<float>::max();
	for (size_t i = 0; i < models.size(); i++) {
		float dist_i;
		if (models[i]->RayIntersect(orig, dir, dist_i) && dist_i < spheres_dist) {
			spheres_dist = dist_i;
			hit = orig + dir * dist_i;
			N = models[i]->GetNormal(hit);
			material = models[i]->GetMaterial();
		}
	}

	float checkerboard_dist = std::numeric_limits<float>::max();
	if (fabs(dir[1]) > 1e-3) {
		float d = -(orig[1] + 4) / dir[1]; // the checkerboard plane has equation y = -4
		Vector3f pt = orig + dir * d;
		if (d > 0 && fabs(pt[0]) < 10 && pt[2]<-10 && pt[2]>-30 && d < spheres_dist) {
			checkerboard_dist = d;
			hit = pt;
			N = Vector3f(0, 1, 0);
			material.diffuse_color = (int(.5f * hit[1] + 1000) + int(.5f * hit[2])) & 1 ? Vector3f(.3f, .3f, .3f) : Vector3f(.3f, .2f, .1f);
		}
	}
	return std::min(spheres_dist, checkerboard_dist) < 1000;
}

Vector3f CastRay(const Vector3f &orig, const Vector3f &dir, const std::vector<const Model*> &models, const std::vector<Light> &lights, size_t depth = 0) {
    Vector3f point, N;
    Material material;

    if (depth>4 || !SceneIntersect(orig, dir, models, point, N, material)) {
        return Vector3f(0.2f, 0.7f, 0.8f); // background color
	}

	Vector3f refract_color(0.f, 0.f, 0.f);
	// compute fresnelt
	float kr;
	Fresnel(dir, N, material.refractive_index, kr);
	// compute refraction if it is not a case of total internal reflection
	if (kr < 1) {
		Vector3f refract_dir = Refract(dir, N, material.refractive_index).normalized();
		Vector3f refract_orig = refract_dir.dot(N) < 0.f ? (point - N * EPSILON).eval() : (point + N * EPSILON).eval();
		refract_color = CastRay(refract_orig, refract_dir, models, lights, depth + 1);
	}

	Vector3f reflect_dir = Reflect(dir, N).normalized();
	Vector3f reflect_orig = reflect_dir.dot(N) < 0.f ? (point - N * EPSILON).eval() : (point + N * EPSILON).eval(); // offset the original point to avoid occlusion by the object itself
	Vector3f reflect_color = CastRay(reflect_orig, reflect_dir, models, lights, depth + 1);

    float diffuse_light_intensity = 0, specular_light_intensity = 0;
    for (size_t i=0; i<lights.size(); i++) {
        Vector3f light_dir      = (lights[i].position - point).normalized();
        float light_distance = (lights[i].position - point).norm();

		Vector3f shadow_orig = light_dir.dot(N) < 0.f ? (point - N * EPSILON).eval() : (point + N * EPSILON).eval(); // checking if the point lies in the shadow of the lights[i]
        Vector3f shadow_pt, shadow_N;
        Material tmpmaterial;
        if (SceneIntersect(shadow_orig, light_dir, models, shadow_pt, shadow_N, tmpmaterial) && (shadow_pt-shadow_orig).norm() < light_distance)
            continue;

        diffuse_light_intensity  += lights[i].intensity * std::max(0.f, light_dir.dot(N));
        specular_light_intensity += powf(std::max(0.f, -Reflect(-light_dir, N).dot(dir)), material.specular_exponent)*lights[i].intensity;
    }

    return material.diffuse_color * diffuse_light_intensity * material.albedo[0] + Vector3f(1., 1., 1.)*specular_light_intensity * material.albedo[1] + reflect_color * material.albedo[2] * kr + refract_color * material.albedo[3] * (1 - kr);
}

void Render(const std::vector<const Model*> &models, const std::vector<Light> &lights) {
    const int   width    = 1024;
    const int   height   = 768;
    const float fov      = M_PI/3.;
    std::vector<Vector3f> framebuffer(width * height);

#pragma omp parallel for
	for (int j = 0; j < height; j++) { // actual rendering loop
        for (int i = 0; i<width; i++) {
            float dir_x =  (i + 0.5f) -  width/2.f;
            float dir_y = -(j + 0.5f) + height/2.f;    // this flips the image at the same time
            float dir_z = -height/(2.f*tan(fov/2.f));
            framebuffer[i+j*width] = CastRay(Vector3f(0.f, 0.f, 0.f), Vector3f(dir_x, dir_y, dir_z).normalized(), models, lights);
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
	RaytracingEngine::Get()->Initialize();

    Material      ivory(0.0, Vector4f(0.6f,  0.3f, 0.1f, 0.0f), Vector3f(0.4f, 0.4f, 0.3f),   50.f);
    Material      glass(1.5, Vector4f(0.0f,  0.5f, 0.1f, 0.8f), Vector3f(0.6f, 0.7f, 0.8f),  125.f);
    Material red_rubber(0.0, Vector4f(0.9f,  0.1f, 0.0f, 0.0f), Vector3f(0.3f, 0.1f, 0.1f),   10.f);
    Material     mirror(0.0, Vector4f(0.0f, 10.0f, 0.8f, 0.0f), Vector3f(1.0f, 1.0f, 1.0f), 1425.f);

    std::vector<const Model*> models;
    models.push_back(new Sphere(Vector3f(-3.f,    0.f,   -16.f), 2.f,      ivory));
	models.push_back(new Sphere(Vector3f(-1.0f, -1.5f, -12.f), 2.f,      glass));
	models.push_back(new Sphere(Vector3f( 1.5f, -0.5f, -18.f), 3.f, red_rubber));
	models.push_back(new Sphere(Vector3f(7.f, 5.f, -18.f), 4.f, mirror));

    std::vector<Light>  lights;
    lights.push_back(Light(Vector3f(-20.f, 20.f,  20.f), 1.5f));
    lights.push_back(Light(Vector3f( 30.f, 50.f, -25.f), 1.8f));
	lights.push_back(Light(Vector3f(30.f, 20.f, 30.f), 1.7f));

	Render(models, lights);

	for (auto model : models)
	{
		delete model;
	}
	models.clear();

    return 0;
}

