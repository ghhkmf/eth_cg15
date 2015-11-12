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

		Color3f Le = Color3f(0.0f);
		const Emitter* emi2 = its.mesh->getEmitter();
		if (emi2) {
			// Get Le
			EmitterQueryRecord newRec(ray.o, its.p, its.shFrame.n);
			Le = emi2->eval(newRec);
		}

		//Global constants-----------------------------------------------------

		//BSDF
		const BSDF* bsdf = its.mesh->getBSDF();

		//Get random Emitter
		const std::vector<Emitter *> emis = scene->getEmitters();
		int randomIndex = rand() % emis.size();
		Emitter* emi = emis.at(randomIndex);

		//Camera Ray
		Vector3f toCam = -ray.d.normalized();

		//Emitter part - Direct Illumination ---------------------------------------------------------------------------------------------------------------------

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

			float dist2_ie = (its.p - emiRecord.p).squaredNorm();

			float cos0_ie = std::abs(
					emiRecord.n.dot(-emiRecord.wi.normalized()));
			float cos_i_ie = std::abs(
					its.geoFrame.n.dot(emiRecord.wi.normalized()));

			float ems_pdf_ie = emi->pdf(emiRecord);

			float mat_pdf_ie = bsdf->pdf(bsdfQueryWie);

			w_em = 1 / (ems_pdf_ie + mat_pdf_ie);

			F_em = Lo * bsdfVal_em * cos_i_ie;

		}

		// BRDF part---------------------------------------------------------------------------------------------------------------------------------------------

		BSDFQueryRecord bsdfQueryWim(its.toLocal(toCam)); //wi Camera, wo sampled ray
		bsdfQueryWim.p = its.p;
		bsdf->sample(bsdfQueryWim, sampler->next2D()); //Has cos already in
		Color3f bsdfVal_b = bsdf->eval(bsdfQueryWim); // No pdf

		//Check if intersect with emitter
		Ray3f lightRay(its.p, its.toWorld(bsdfQueryWim.wo));

		Color3f F_mat = Color3f(0.0f);
		float w_mat = 1;
		Intersection its_im;
		if (scene->rayIntersect(lightRay, its_im)) {
			if (its_im.mesh->getEmitter()) {
				//It intersects with light source
				const Emitter* emi = its_im.mesh->getEmitter();
				EmitterQueryRecord iRec = EmitterQueryRecord(its.p, its_im.p,
						its_im.geoFrame.n);
				Color3f Fo = emi->eval(iRec);
				F_mat = Fo * bsdfVal_b;

				EmitterQueryRecord emiRecord_im = EmitterQueryRecord(its.p,
						its_im.p, its_im.geoFrame.n);
				float cos_i_im = std::abs(its.geoFrame.n.dot(emiRecord_im.wi));
				float cos0_im = std::abs(
						its_im.geoFrame.n.dot(-emiRecord_im.wi));

				float dist2_im = (its.p - its_im.p).squaredNorm();

				float ems_pdf_im = emi->pdf(emiRecord_im);
				float mat_pdf_im = bsdf->pdf(bsdfQueryWim);

				w_mat = 1 / (ems_pdf_im + mat_pdf_im);
				F_mat = Fo * bsdfVal_b * cos_i_im;
			}
		}

		Color3f emitterPart = w_em * F_em;
		Color3f materialPart = w_mat * F_mat;

		return Le + emitterPart + materialPart;

	}

	std::string toString() const {
		return "DirectMISIntegrator[]";
	}
}
;

NORI_REGISTER_CLASS(DirectMISIntegrator, "direct_mis");
NORI_NAMESPACE_END
