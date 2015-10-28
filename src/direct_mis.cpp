#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class DirectMISIntegrator: public Integrator {
public:
	DirectMISIntegrator(const PropertyList &props) {
		/* No parameters this time */
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
		/* Find the surface that is visible in the requested direction */
		Intersection its;

		// If not visible return black
		if (!scene->rayIntersect(ray, its))
			return Color3f(0.0f);

		const std::vector<Emitter *> emis = scene->getEmitters();

		Color3f Le = Color3f(0.0f);

		const Emitter* emi2 = its.mesh->getEmitter();
		if (emi2) {
			// Get Le
			EmitterQueryRecord iRec2(its.p);
			Le = emi2->eval(iRec2);
		}

		//Emitter part
		Color3f F_em = Color3f(0.0f);

		EmitterQueryRecord iRec = EmitterQueryRecord(its.p);

		for (std::vector<Emitter*>::const_iterator it = emis.begin();
				it != emis.end(); ++it) {

			Emitter* emi = *it;

			//Sample Emitter
			Color3f Lo_em = emi->sample(iRec, sampler->next2D()); //Div by pdf

			//Get BSDF
			const BSDF* bsdf = its.mesh->getBSDF();
			BSDFQueryRecord query = BSDFQueryRecord(its.toLocal(iRec.wi),
					its.toLocal(-ray.d), ESolidAngle);
			query.uv = its.uv;
			Color3f bsdfVal_em = bsdf->eval(query);

			//Check if something blocks the visibility
			if (!scene->rayIntersect(iRec.shadowRay)) {
				//Multiply by cos of normal of reflec. normal
				float cos_i = iRec.n.dot(-iRec.wi);
				if (cos_i < 0)
					cos_i = -cos_i;

				float cos0 = its.geoFrame.n.dot(iRec.wi);
				if (cos0 < 0)
					cos0 = -cos0;

				F_em = F_em
						+ Lo_em * bsdfVal_em * cos0 * cos_i
								/ (its.p - iRec.p).squaredNorm();
			}
		}

		// BRDF part------------------------------------
		Color3f F_mat = Color3f(0.0f);

		const BSDF* bsdf = its.mesh->getBSDF();
		Vector3f toCam = -ray.d.normalized();

		BSDFQueryRecord query(its.toLocal(toCam)); //wi Camera, wo sampled ray
		query.p = its.p;
		Color3f bsdfVal_b = bsdf->sample(query, sampler->next2D()); //Has cos already in

		float mat_pdf = bsdf->pdf(query);
		//Check if intersect with emitter
		Ray3f lightRay(its.p, its.toWorld(query.wo));

		Intersection its2;
		if (scene->rayIntersect(lightRay, its2)) {

			if (its2.mesh->getEmitter()) {
				//It intersects with light source
				const Emitter* emi = its2.mesh->getEmitter();
				EmitterQueryRecord iRec = EmitterQueryRecord(its2.p);
				Color3f Lo_b = emi->eval(iRec);

				F_mat = Lo_b * bsdfVal_b;
			}
		}

		//Same units
		float cos_i = Frame::cosTheta(query.wi);
		if (cos_i < 0)
			cos_i = -cos_i;

		float cos0 = Frame::cosTheta(query.wo);
		if (cos0 < 0)
			cos0 = -cos0;

		float ems_pdf = iRec.pdf;
		mat_pdf=mat_pdf*cos0*cos_i/(its.p - its2.p).squaredNorm();

		float w_em = ems_pdf / (ems_pdf + mat_pdf);
		float w_mat = mat_pdf / (ems_pdf + mat_pdf);

		return Le + w_em * F_em + w_mat * F_mat;
	}

	std::string toString() const {
		return "DirectMISIntegrator[]";
	}
}
;

NORI_REGISTER_CLASS(DirectMISIntegrator, "direct_mis");
NORI_NAMESPACE_END
