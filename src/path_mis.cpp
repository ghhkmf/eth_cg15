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
		m_q = 0.2;
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
		return this->Li(scene, sampler, ray, false);
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray,
			bool isDiffuseBounce) const {

		bool isNextEmitter = false;

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
			EmitterQueryRecord iRec2(ray.o, its.p, its.shFrame.n);
			Le = emi2->eval(iRec2);
		}

		//Get Li - Infos
		const BSDF* bsdf = its.mesh->getBSDF();
		Vector3f toCam = -ray.d.normalized();

		//	BSDFQueryRecord query(its.toLocal(toCam)); //wi Camera, wo sampled ray
		//	query.p = its.p;
		//	Color3f bsdfVal = bsdf->sample(query, sampler->next2D());

		//Check if intersect with emitter
		//	Ray3f lightRay(its.p, its.toWorld(query.wo));

		//Get random Emitter
		const std::vector<Emitter *> emis = scene->getEmitters();
		int randomIndex = rand() % emis.size();
		Emitter* emi = emis.at(randomIndex);

		// Get Ld - MULTIPLE IMPORTANCE SAMPLING -----------------------------------------------------------------
		Color3f Ld(0.f);
		// BSDF Part

		BSDFQueryRecord bsdfQueryWim(its.toLocal(toCam)); //wi Camera, wo sampled ray
		bsdfQueryWim.p = its.p;
		Color3f bsdfVal = bsdf->sample(bsdfQueryWim, sampler->next2D()); //Has cos already in
		Color3f bsdfVal_b = bsdf->eval(bsdfQueryWim); // No pdf

		//Check if intersect with emitter
		Ray3f lightRay(its.p, its.toWorld(bsdfQueryWim.wo));

		Color3f F_mat = Color3f(0.0f);
		float w_mat = 1;
		Intersection its_im;
		if (scene->rayIntersect(lightRay, its_im)) {
			if (its_im.mesh->getEmitter()) {
				isNextEmitter = true;
				//It intersects with light source
				const Emitter* emi = its_im.mesh->getEmitter();
				EmitterQueryRecord iRec = EmitterQueryRecord(its.p, its_im.p,
						its_im.geoFrame.n);
				Color3f Fo = emi->eval(iRec);
				F_mat = Fo * bsdfVal_b;

				EmitterQueryRecord emiRecord_im = EmitterQueryRecord(its.p,
						its_im.p, its_im.geoFrame.n);
				float cos_i_im = std::abs(its.geoFrame.n.dot(emiRecord_im.wi));
				//		float cos0_im = std::abs(
				//				its_im.geoFrame.n.dot(-emiRecord_im.wi));

				//		float dist2_im = (its.p - its_im.p).squaredNorm();

				float ems_pdf_im = emi->pdf(emiRecord_im);
				float mat_pdf_im = bsdf->pdf(bsdfQueryWim);

				w_mat = 1 / (ems_pdf_im + mat_pdf_im);
				F_mat = Fo * bsdfVal_b * cos_i_im;
			}
		}

		//Emitter Ld-----------------------------------------------------------------------------------------------------------------------------------------------

		//Sample Emitter
		EmitterQueryRecord emiRecord = EmitterQueryRecord(its.p);
		emi->sample(emiRecord, sampler->next2D());
		//Instead of iterating over all emitter
		Color3f Lo = emis.size() * emi->eval(emiRecord); //Get Lo (not divided by pdf)

		BSDFQueryRecord bsdfQueryWie = BSDFQueryRecord(its.toLocal(-ray.d),
				its.toLocal(emiRecord.wi), ESolidAngle);

		bsdfQueryWie.uv = its.uv;
		Color3f bsdfVal_em = bsdf->eval(bsdfQueryWie);

		Color3f F_em = Color3f(0.0f);
		float w_em = 1;

		//Check if something blocks the visibility
		if (!scene->rayIntersect(emiRecord.shadowRay)) {

			//	float dist2_ie = (its.p - emiRecord.p).squaredNorm();

			//	float cos0_ie = std::abs(
			//			emiRecord.n.dot(-emiRecord.wi.normalized()));
			float cos_i_ie = std::abs(
					its.geoFrame.n.dot(emiRecord.wi.normalized()));

			float ems_pdf_ie = emi->pdf(emiRecord)/emis.size();

			float mat_pdf_ie = bsdf->pdf(bsdfQueryWie);

			w_em = 1 / (ems_pdf_ie + mat_pdf_ie);

			F_em = Lo * bsdfVal_em * cos_i_ie;

		}

		//Same units------------------------------------------------------------------------------------------------------------------------------------------------------

		Ld = w_em * F_em + w_mat * F_mat;

		//----------------------------------------
		//	Russian Roulette

		Color3f Li(0.f);
		if (sampler->next1D() > m_q) {
			Li = this->Li(scene, sampler, lightRay, bsdf->isDiffuse())
								* bsdfVal;
			if (!isNextEmitter) {
				return (Le + Ld + Li) / (1 - m_q);
			} else {
				return (Le +Li) / (1 - m_q);
			}
		} else {
			return (Le + Ld) / (1 - m_q);
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
