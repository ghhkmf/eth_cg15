#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/sampler.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class DirectMATSIntegrator: public Integrator {
public:
	DirectMATSIntegrator(const PropertyList &props) {
		/* No parameters this time */
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
		/* Find the surface that is visible in the requested direction */
		Intersection its;

		// If not visible return black
		if (!scene->rayIntersect(ray, its))
			return Color3f(0.0f);

		const std::vector<Emitter *> emis = scene->getEmitters();

		Color3f result = Color3f(0.0f);

		//Get Le
		const Emitter* emi2 = its.mesh->getEmitter();
		if (emi2) {
			EmitterQueryRecord iRec2(its.p);
			result = result + emi2->eval(iRec2);
		}

		const BSDF* bsdf = its.mesh->getBSDF();
		Vector3f toCam = -ray.d.normalized();

		BSDFQueryRecord query(its.toLocal(toCam)); //wi Camera, wo sampled ray
		query.p = its.p;
		Color3f bsdfVal = bsdf->sample(query, sampler->next2D()); //Has cos already in

		//Check if intersect with emitter
		Ray3f lightRay(its.p, its.toWorld(query.wo));

		Intersection its2;
		if (scene->rayIntersect(lightRay, its2)) {

			if (its2.mesh->getEmitter()) {
				//It intersects with light source
				const Emitter* emi = its2.mesh->getEmitter();
				EmitterQueryRecord iRec = EmitterQueryRecord(its2.p);
				Color3f Lo = emi->eval(iRec);
				result = result + Lo * bsdfVal;
			}
		}
		return result;
	}

	std::string toString() const {
		return "DirectMATSIntegrator[]";
	}
}
;

NORI_REGISTER_CLASS(DirectMATSIntegrator, "direct_mats");
NORI_NAMESPACE_END
