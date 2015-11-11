#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/sampler.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class DirectPATHMISIntegrator: public Integrator {
public:
	DirectPATHMISIntegrator(const PropertyList &props) {
		m_q = 0.4; //TODO: Chech this
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
		return this->Li(scene, sampler, ray, false);
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray,
			bool isDiffuseBounce) const {
		/* Find the surface that is visible in the requested direction */
		Intersection its;

		// If not visible return black
		if (!scene->rayIntersect(ray, its))
			return Color3f(0.0f);

		//Get Le
		Color3f Le = Color3f(0.0f);
		const Emitter* emi2 = its.mesh->getEmitter();
		if (emi2) {
			//Its an Emitter.
			EmitterQueryRecord iRec2(its.p);
			Le = emi2->eval(iRec2);
		}

		//Get Li - Infos
		const BSDF* bsdf = its.mesh->getBSDF();
		Vector3f toCam = -ray.d.normalized();

		BSDFQueryRecord query(its.toLocal(toCam)); //wi Camera, wo sampled ray
		query.p = its.p;
		Color3f bsdfVal = bsdf->sample(query, sampler->next2D());

		//Check if intersect with emitter
		Ray3f lightRay(its.p, its.toWorld(query.wo));

		// Get Ld - MULTIPLE IMPORTANCE SAMPLING

		// BSDF Part
		float mat_pdf = bsdf->pdf(query);
		Color3f F_mat(0.f);
		Intersection its2;

		if (scene->rayIntersect(lightRay, its2)) {

			if (its2.mesh->getEmitter()) {
				//It intersects with light source
				const Emitter* emi = its2.mesh->getEmitter();
				EmitterQueryRecord iRec = EmitterQueryRecord(its2.p);
				Color3f Lo_b = emi->eval(iRec);
				F_mat = Lo_b * bsdfVal;
			}
		}

		Color3f Ld(0.f);
		Color3f emiVal;
		float cos_i = 1.f;
		float cos0 = 1.f;
		Color3f emiBsdfVal;

		//Emitter Ld-------------

		const std::vector<Emitter *> emis = scene->getEmitters();
		int randomIndex = rand() % emis.size();
		Emitter* emi = emis.at(randomIndex);
		EmitterQueryRecord iRec = EmitterQueryRecord(its.p);
		emiVal = emi->sample(iRec, sampler->next2D());

		cos_i = std::abs(iRec.n.dot(-iRec.wi)); //iRec.ref = inter Point, iRec.p = point at emitter
		cos0 = std::abs(its.geoFrame.n.dot(iRec.wi));

		//Get bsdfVal
		BSDFQueryRecord emiQuery(its.toLocal(toCam), -iRec.wi, ESolidAngle);
		emiBsdfVal = bsdf->eval(emiQuery);

		Color3f F_em = emiBsdfVal * emiVal * cos_i * cos0
				/ (its.p - iRec.p).squaredNorm();

		float ems_pdf = iRec.pdf;

		//Same units--------------------

		//Same units
		cos_i = std::abs(Frame::cosTheta(query.wi));

		cos0 = std::abs(Frame::cosTheta(query.wo));

		mat_pdf = mat_pdf * cos0 * cos_i / ((its.p - its2.p).squaredNorm()*emis.size());
		float w_em = ems_pdf / (ems_pdf + mat_pdf);
		float w_mat = mat_pdf / (ems_pdf + mat_pdf);

		Ld = w_em * F_em + w_mat * F_mat;

		//----------------------------------------
		//	Russian Roulette

		Color3f Li(0.f);
		if (sampler->next1D() > m_q) {
			Li = this->Li(scene, sampler, lightRay, bsdf->isDiffuse())
					* bsdfVal;
		}

		if (isDiffuseBounce) { //== Not Specular
			return (Ld + Li) / (1 - m_q);
		} else {
			return (Le + Ld + Li) / (1 - m_q);
		}
	}

	std::string toString() const {
		return "DirectPATHMISIntegrator[]";
	}

private:
	float m_q;
}
;

NORI_REGISTER_CLASS(DirectPATHMISIntegrator, "path_mis");
NORI_NAMESPACE_END
