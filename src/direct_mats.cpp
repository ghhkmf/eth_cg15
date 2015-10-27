#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/sampler.h>

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
			EmitterQueryRecord iRec2 = EmitterQueryRecord(its.p);
			iRec2.wi = -its.geoFrame.n;
			iRec2.n = its.geoFrame.n;
			iRec2.p = its.p + iRec2.n; // damit r = 1
			result = result + emi2->eval(iRec2);
		}

		const BSDF* bsdf = its.mesh->getBSDF();
		BSDFQueryRecord query = BSDFQueryRecord(); //TODO: Not sure about this. Does this samples wi?
		query.wo = its.toLocal(-ray.d);
		query.uv = its.uv;
		query.p = its.p;
		query.measure = ESolidAngle;
		Color3f bsdfVal = bsdf->sample(query, sampler->next2D()); //Has cos already in

		//Check if intersect with emitter
		Ray3f lightRay(query.wi);
		Intersection its2;
		if (scene->rayIntersect(lightRay,its2)) {
			if (its2.mesh->getEmitter()) {
				//It intersects with light source
				const Emitter* emi = its2.mesh->getEmitter();
				EmitterQueryRecord iRec = EmitterQueryRecord(its.p,its2.p);
				iRec.n = its2.geoFrame.n;
				Color3f Lo = emi->eval(iRec);


			//	float cos0 = its.geoFrame.n.dot(iRec.wi);
			//	if (cos0 < 0)
			//		cos0 = -cos0;
				result = result + Lo * bsdfVal;// * cos0;
			}

		}

//					BSDFQueryRecord query = BSDFQueryRecord(its.toLocal(iRec.wi),
//							its.toLocal(-ray.d), ESolidAngle);
//					query.uv = its.uv;
//					Color3f bsdfVal = bsdf->eval(query);

		/*
		 for (std::vector<Emitter*>::const_iterator it = emis.begin();
		 it != emis.end(); ++it) {

		 Emitter* emi = *it;
		 EmitterQueryRecord iRec = EmitterQueryRecord(its.p);

		 //Get BSDF
		 //For Emitter and sample
		 //Get point Le if exists
		 const Emitter* emi2 = its.mesh->getEmitter();
		 if (emi2) {
		 // Get Le
		 EmitterQueryRecord iRec2 = EmitterQueryRecord(its.p, its.p,
		 its.geoFrame.n);
		 result = result
		 + its.mesh->getEmitter()->sample(iRec2,
		 Point2f(sampler->next2D()));
		 }

		 Color3f Lo = emi->sample(iRec, sampler->next2D());

		 //Get BSDF


		 result = result + Lo + bsdfVal / emi->pdf(iRec);
		 }

		 */

		return result;
	}

	std::string toString() const {
		return "DirectEMSIntegrator[]";
	}
}
;

NORI_REGISTER_CLASS(DirectMATSIntegrator, "direct_mats");
NORI_NAMESPACE_END
