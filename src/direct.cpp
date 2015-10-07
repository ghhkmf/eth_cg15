#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

class DirectIntegrator: public Integrator {
public:
	DirectIntegrator(const PropertyList &props) {
		/* No parameters this time */
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
		/* Find the surface that is visible in the requested direction */
		Intersection its;

		if (!scene->rayIntersect(ray, its))
			return Color3f(0.0f);

		const std::vector<Emitter *> emis = scene->getEmitters();

		Color3f result = Color3f(0.0f);

		for (std::vector<Emitter*>::const_iterator it = emis.begin();
				it != emis.end(); ++it) {
			Emitter* emi = *it;

			EmittedValues val = emi->getEmittedValues(its.p);

			//Check if in Shadow
			if (!scene->rayIntersect(val.ray)) {

				//Get BSDF
				const BSDF* bsdf = its.mesh->getBSDF();

				BSDFQueryRecord query = BSDFQueryRecord(its.toLocal(val.wi),
						its.toLocal(-ray.d), ESolidAngle);
				query.uv=its.uv;
				Color3f bsdfVal = bsdf->eval(query);

				// Get Cos theta
				float cos_theta_i = its.shFrame.n.dot(val.wi);

				if (cos_theta_i >= 0)
					result += val.Li * bsdfVal * cos_theta_i;
				else
					result += val.Li * bsdfVal * (-cos_theta_i);

			} else {
				result += Color3f(0.0f);
			}
		}
		return result;
	}

	std::string toString() const {
		return "NormalIntegrator[]";
	}
};

NORI_REGISTER_CLASS(DirectIntegrator, "direct");
NORI_NAMESPACE_END
