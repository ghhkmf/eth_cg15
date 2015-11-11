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
			EmitterQueryRecord iRec2(its.p);
			Le = emi2->eval(iRec2);
		}

		//Global constants------------------------

		//BSDF
		const BSDF* bsdf = its.mesh->getBSDF();

		//Get random Emitter
		const std::vector<Emitter *> emis = scene->getEmitters();
		int randomIndex = rand() % emis.size();
		Emitter* emi = emis.at(randomIndex);

		//Camera Ray
		Vector3f toCam = -ray.d.normalized();

		//Emitter part - Direct Illumination --------------------------
		Color3f F_em = Color3f(0.0f);

		EmitterQueryRecord emiRecord = EmitterQueryRecord(its.p);
		Color3f Lo_em = emis.size() * emi->sample(emiRecord, sampler->next2D()); //Instead of iterating over all emitter, mulply because MC

		// W_ie = emiRecord.wi - from ref to p
		// P_ie = emiRecord.p
		// n_ie = emiRecord.n

		BSDFQueryRecord bsdfQueryWie = BSDFQueryRecord(
				its.toLocal(emiRecord.wi), its.toLocal(-ray.d), ESolidAngle);
		bsdfQueryWie.uv = its.uv;
		Color3f bsdfVal_em = bsdf->eval(bsdfQueryWie);

		//Check if something blocks the visibility
		if (!scene->rayIntersect(emiRecord.shadowRay)) {
			//Multiply by cos of normal of reflec. normal
			float cos_i = std::abs(emiRecord.n.dot(-emiRecord.wi));
			float cos0 = std::abs(its.geoFrame.n.dot(emiRecord.wi));

			F_em = Lo_em * bsdfVal_em * cos0 * cos_i
					/ (its.p - emiRecord.p).squaredNorm();
		}

		// BRDF part------------------------------------
		Color3f F_mat = Color3f(0.0f);

		BSDFQueryRecord bsdfQueryWim(its.toLocal(toCam)); //wi Camera, wo sampled ray
		bsdfQueryWim.p = its.p;
		Color3f bsdfVal_b = bsdf->sample(bsdfQueryWim, sampler->next2D()); //Has cos already in

		//W_im = bsdfQueryWim.wo
		//P_im = its_im.p
		//N_im = its_im.geoFrame.n

		//Check if intersect with emitter
		Ray3f lightRay(its.p, its.toWorld(bsdfQueryWim.wo));

		Intersection its_im;
		if (scene->rayIntersect(lightRay, its_im)) {

			if (its_im.mesh->getEmitter()) {
				//It intersects with light source
				const Emitter* emi = its_im.mesh->getEmitter();
				EmitterQueryRecord iRec = EmitterQueryRecord(its_im.p);
				Color3f Lo_b = emi->eval(iRec);
				F_mat = Lo_b * bsdfVal_b;
			}
		}

		//Same units
		EmitterQueryRecord emiRecord_im = EmitterQueryRecord(its.p, its_im.p,
						its_im.geoFrame.n);
		//w_ie part
		float cos0_ie = std::abs(emiRecord.n.dot(-emiRecord.wi.normalized()));
		float cos_i_ie = std::abs(its.geoFrame.n.dot(emiRecord.wi.normalized()));


		float cos0_im = std::abs(its.geoFrame.n.dot(emiRecord_im.wi));
		float cos_i_im = std::abs(its_im.geoFrame.n.dot(-emiRecord_im.wi));


		float ems_pdf_ie = emiRecord.pdf * cos0_ie * cos_i_ie
				/ ((its.p - emiRecord.p).squaredNorm()*emis.size());;

		float mat_pdf_ie = bsdf->pdf(
				BSDFQueryRecord(its.toLocal(toCam), emiRecord.wi, ESolidAngle));



		// w_im part
		float ems_pdf_im = emi->pdf(emiRecord_im) * cos0_im * cos_i_im
				/ ((its.p - its_im.p).squaredNorm()*emis.size());
		float mat_pdf_im = bsdf->pdf(bsdfQueryWim);





		//Total part
		float w_em = ems_pdf_ie / (ems_pdf_ie + mat_pdf_ie);
		float w_mat = mat_pdf_im / (ems_pdf_im + mat_pdf_im);

		bool isEmNan = false;
		if (w_em != w_em) {
			//	cout << "W_em is NaN" << endl;
			isEmNan = true;
		}
		bool isMatNan = false;
		if (w_mat != w_mat) {
			//	cout << "W_mat is NaN" << endl;
			isMatNan = true;
		}

		if (!isEmNan && !isMatNan)
			return Le + w_em * F_em + w_mat * F_mat;

		return Le + F_em;
	}

	std::string toString() const {
		return "DirectMISIntegrator[]";
	}
}
;

NORI_REGISTER_CLASS(DirectMISIntegrator, "direct_mis");
NORI_NAMESPACE_END
