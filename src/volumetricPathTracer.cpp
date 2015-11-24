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
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
		/* Find the surface that is visible in the requested direction */
		Intersection its;

		if (!scene->rayIntersect(ray, its))
			return Color3f(0.0f);


		throw new NoriException("Not implemented");



		//TODO: Check if intersection is Surface, Light or Volume



		//TODO:If Volume -> Volumetric Path Tracing


		//TODO: If Surface ->


		//TODO: If Light -> Sample and return


		/*
		const std::vector<Emitter *> emis = scene->getEmitters();

		Color3f result = Color3f(0.0f);

		for (std::vector<Emitter*>::const_iterator it = emis.begin();
				it != emis.end(); ++it) {
			Emitter* emi = *it;

			EmitterQueryRecord iRec(its.p);
			Color3f Li = emi->sample(iRec,sampler->next2D());



			//Check if in Shadow
			if (!scene->rayIntersect(iRec.shadowRay)) {

				//Get BSDF
				const BSDF* bsdf = its.mesh->getBSDF();

				BSDFQueryRecord query = BSDFQueryRecord(its.toLocal(iRec.wi),
						its.toLocal(-ray.d), ESolidAngle);
				query.uv=its.uv;
				Color3f bsdfVal = bsdf->eval(query);

				// Get Cos theta
				float cos_theta_i = its.shFrame.n.dot(iRec.wi);

				if (cos_theta_i >= 0)
					result += Li * bsdfVal * cos_theta_i;
				else
					result += Li * bsdfVal * (-cos_theta_i);

			} else {
				result += Color3f(0.0f);
			}
		}
		*/
		return result;
	}

	std::string toString() const {
		return "VolumetricPathTracer[]";
	}
};

NORI_REGISTER_CLASS(VolumetricPathTracer, "volPathTracer");
NORI_NAMESPACE_END
