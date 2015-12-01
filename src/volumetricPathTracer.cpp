#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class VolumetricPathTracer: public Integrator {
public:
	VolumetricPathTracer(const PropertyList &props) {
		/* No parameters this time */
		//Right now its just path_mats
		//TODO: Vol Path Tracer, Homogeneous
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {

		float m_q = 0.05;

		Intersection its;

		// If not visible return black
		if (!scene->rayIntersect(ray, its))
			return Color3f(0.0f);

		//Get Le
		Color3f Le = Color3f(0.0f);
		const Emitter* emi2 = its.mesh->getEmitter();
		if (emi2) {
			//Its an Emitter.
			EmitterQueryRecord iRec2(ray.o, its.p, its.shFrame.n);
			Le = emi2->eval(iRec2);
		}

		// No Ld

		//Get Li
		const BSDF* bsdf = its.mesh->getBSDF();
		Vector3f toCam = -ray.d.normalized();

		BSDFQueryRecord query(its.toLocal(toCam)); //wi Camera, wo sampled ray
		query.p = its.p;
		Color3f bsdfVal = bsdf->sample(query, sampler->next2D());

		//Check if intersect with emitter
		Ray3f lightRay(its.p, its.toWorld(query.wo));

		if (sampler->next1D() > m_q)
			return Le
					+ this->Li(scene, sampler, lightRay) * bsdfVal / (1.f - m_q);
		else
			return Le;
	}

	std::string toString() const {
		return "VolumetricPathTracer[]";
	}
};

NORI_REGISTER_CLASS(VolumetricPathTracer, "volPathTracer");
NORI_NAMESPACE_END
