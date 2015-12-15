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

		Color3f result;
		Color3f multiConst(1.f);

		bool run = true;

		while (run) {

			bool isNextEmitter = false;

			/* Find the surface that is visible in the requested direction */
			Intersection its;

			// If not visible return black
			if (!scene->rayIntersect(ray, its))
				return Color3f(0.0f);

			//Get Le
			Color3f Le = Color3f(0.0f);
			const Emitter* itsEmitter = its.mesh->getEmitter();
			if (itsEmitter) {
				//Its an Emitter.
				EmitterQueryRecord iRec2(ray.o, its.p, its.shFrame.n);
				Le = itsEmitter->eval(iRec2);
			}

			// Global constants
			const BSDF* currentBSDF = its.mesh->getBSDF();
			Vector3f rayToCam = -ray.d.normalized();

			//Get random Emitter
			const std::vector<Emitter *> emitterVector = scene->getEmitters();
			int randomIndex = rand() % emitterVector.size();
			Emitter* currentSampledEmitter = emitterVector.at(randomIndex);

			// Get Ld - MULTIPLE IMPORTANCE SAMPLING -----------------------------------------------------------------
			Color3f Ld(0.f);

			//Emitter Ld-----------------------------------------------------------------------------------------------------------------------------------------------

			//Sample Emitter
			EmitterQueryRecord currentSampledEmitterQuery(its.p);

			//Emitter sampling: wi , ref-> p in Emitter, sampled
			// wo is ref -> cam
			//All local

			currentSampledEmitter->sample(currentSampledEmitterQuery,
					sampler->next2D());

			//Instead of iterating over all emitter
			Color3f L_sampledEmitter = emitterVector.size()
					* currentSampledEmitter->eval(currentSampledEmitterQuery); //Get Lo (not divided by pdf)

			BSDFQueryRecord currentBSDFQueryToSampledEmitter(
					its.toLocal(rayToCam),
					its.toLocal(currentSampledEmitterQuery.wi), ESolidAngle);

			currentBSDFQueryToSampledEmitter.uv = its.uv;
			Color3f bsdfValToSampledEmitter = currentBSDF->eval(
					currentBSDFQueryToSampledEmitter);

			Color3f F_em = Color3f(0.0f);
			float w_em = 1;

			//Check if something blocks the visibility
			if (!scene->rayIntersect(currentSampledEmitterQuery.shadowRay)) {

				float cos_itsNormal_sEmitter = std::abs(
						its.shFrame.n.dot(
								currentSampledEmitterQuery.wi.normalized()));

				float ems_pdf_ie = currentSampledEmitter->pdf(
						currentSampledEmitterQuery) / emitterVector.size();

				float mat_pdf_ie = currentBSDF->pdf(
						currentBSDFQueryToSampledEmitter);

				w_em = 1 / (ems_pdf_ie + mat_pdf_ie);

				F_em = L_sampledEmitter * bsdfValToSampledEmitter
						* cos_itsNormal_sEmitter;

			} else {
				//Emitter constribution is 0, because its blocked.
			}

			// BSDF Part ----------------------------------------------------------------------------------------------

			BSDFQueryRecord currentBSDFQueryToSampledDirection(
					its.toLocal(rayToCam)); //wi Camera, wo sampled ray
			currentBSDFQueryToSampledDirection.p = its.p;

			Color3f bsdfVal_cos_pdf = currentBSDF->sample(
					currentBSDFQueryToSampledDirection, sampler->next2D()); //Has cos already in
			Color3f bsdfValToSampledDirection = currentBSDF->eval(
					currentBSDFQueryToSampledDirection); // No pdf
			//Sampled direction in wo
			float cos_itsNormal_sDirection = std::abs(
					its.shFrame.n.dot(currentBSDFQueryToSampledDirection.wo));
			float currentBSDFpdf_sDirection = currentBSDF->pdf(
					currentBSDFQueryToSampledDirection); //WARNING: PDF not in 0-1
			//Check if intersect with emitter
			Ray3f nextRay(its.p,
					its.toWorld(currentBSDFQueryToSampledDirection.wo));

			Color3f F_mat = Color3f(0.0f);
			float w_mat = 1;
			Intersection its_sampledDirection;

			if (scene->rayIntersect(nextRay, its_sampledDirection)) {
				if (its_sampledDirection.mesh->getEmitter()) {
					isNextEmitter = true;

					//It intersects with light source
					const Emitter* emi =
							its_sampledDirection.mesh->getEmitter();

					EmitterQueryRecord emitterQuerySampledDirection(its.p,
							its_sampledDirection.p,
							its_sampledDirection.geoFrame.n);
					Color3f Fo = emi->eval(emitterQuerySampledDirection);

					float ems_pdf_im = emi->pdf(emitterQuerySampledDirection);
					float mat_pdf_im = currentBSDFpdf_sDirection;

					w_mat = 1 / (ems_pdf_im + mat_pdf_im);
					F_mat = Fo * bsdfValToSampledDirection;
				} else {
					// Doesnt it doesnt hit an emitter, contribution to Ld is zero
				}
			} else {
				// Doesnt anything, ergo, it doesnt hit an emitter, contribution to Ld is zero
			}

			//Same units------------------------------------------------------------------------------------------------------------------------------------------------------

			Ld = w_em * F_em + w_mat * F_mat;
			//TODO: Problem: table_path_mis.xml

			result += (Le + Ld) * multiConst;

			//----------------------------------------
			//	Russian Roulette
			Color3f Li(0.f);
			if (sampler->next1D() > m_q) {
				Li = this->Li(scene, sampler, nextRay, currentBSDF->isDiffuse())
						* bsdfVal_cos_pdf;
				if (isNextEmitter) {
					run = false;
				} else {
					multiConst *= bsdfVal_cos_pdf / (1 - m_q);
				}
			} else {
				run = false;
			}
		}
		return result;

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
