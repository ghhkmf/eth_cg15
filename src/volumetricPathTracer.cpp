#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/medium.h>
#include <nori/sampler.h>
#include <Eigen/Dense>

NORI_NAMESPACE_BEGIN

class VolumetricPathTracer: public Integrator {
public:
	VolumetricPathTracer(const PropertyList &props) {
		m_q = 0.05;
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &_ray) const {

		bool inMedium = false;

		Color3f result;
		Color3f multiConst(1.f / (1 - m_q));
		Ray3f ray = _ray;
		Medium* medium = nullptr;

		bool run = true;

		while (run) {

			//Reset Variables
			bool isNextEmitter = false;

			Intersection its;

			// If not visible return black
			if (!scene->rayIntersect(ray, its)) {
				result += Color3f(0.0f);
				run = false;
				continue;
			}

			// Check if we are in a Medium
			medium = nullptr;
			const std::vector<Medium *> ms = scene->getMediums();
			if (ms.size() > 0) {
				for (unsigned int var = 0; var < ms.size(); ++var) {
					Medium* m = ms.at(var);
					if (m->isInside(ray.o)) {
						medium = m;
						inMedium = true;
						break;
					}
				}
			}
			if (medium) {

				float tmax = its.t;
				float t = medium->sampleFreeFlightDistance(sampler->next1D()); //Sample free flying path

				if (t < tmax) {
					//Volume interaction
					// x=itp.s
					Intersection its_x;
					its_x.p = ray.o + ray.d * t;
					its_x.t = t;

					its_x.geoFrame.n = (Normal3f) ray.d.cross(
							Vector3f(0, 0, 1));
					its_x.shFrame.n = its_x.geoFrame.n;

					MediumQueryRecord query(ray.o, its_x.p);
					float pfVal = medium->sample(query, sampler->next2D()); // this is simply sigma_t/sigma_s
					//query.wo sampled w
					//float pdf = medium->pdf(query);

					Ray3f newRay(its_x.p, query.wo);

					ray = newRay;
					//	cout << pfVal << " " << newRay.o << " " << newRay.d << endl;
					multiConst *= pfVal;
					//					return pfVal * this->Li(scene, sampler, newRay);

					//Emitter sampling?
					//
					// if (sampler->next1D() > m_q)
					//return pfVal * this->Li(scene, sampler, newRay) / (1.f - m_q);
					//else
					//return Color3f(0.f);
					//
				} else {

					//Surface interaction
					//Get Le
					Color3f Le = Color3f(0.0f);
					Ray3f shadowRay;
					const Emitter* emi2 = its.mesh->getEmitter();
					if (emi2) {
						//Its an Emitter.
						EmitterQueryRecord iRec2(ray.o, its.p, its.shFrame.n);
						Le = emi2->eval(iRec2);
					} else {

						Le = sampleRandomLight(its, scene, sampler, ray,
								shadowRay);
						//
						//Intersection itsEmi;
						//if (scene->rayIntersect(shadowRay, itsEmi)) {
						//float TrEmitter = medium->getTransmittanceValue(shadowRay.o,
						//itsEmi.p);
						//Le=Le*TrEmitter;
						//}

					}

					float Tr = medium->getTransmittanceValue(ray.o, its.p);
					float pdf1 = 1 - Tr; //PDf hitting surface

					if (Tr == 1) {
						Tr = 0;
						pdf1 = 1;
					}
					const BSDF* bsdf = its.mesh->getBSDF();
					Vector3f toCam = -ray.d.normalized();

					BSDFQueryRecord query(its.toLocal(toCam)); //wi Camera, wo sampled ray
					query.p = its.p;
					query.uv = its.uv;
					Color3f bsdfVal = bsdf->sample(query, sampler->next2D());

					//Sampled ray in query.wo

					// Next Ray
					Ray3f newRay(its.p, its.toWorld(query.wo));
					//				return bsdfVal * Tr * this->Li(scene, sampler, newRay)
					//						/ (pdf1 * pdf2);

					ray = newRay;
					result += Le * multiConst;

					if (sampler->next1D() > m_q) {
						multiConst *= bsdfVal * Tr / ((1.f - m_q) * pdf1);
						//						return Le
						//								+ bsdfVal * this->Li(scene, sampler, newRay)* Tr / ((1.f - m_q) * pdf1);
					} else {
						run = false;
						//						return Le;
					}

			}

		} else {
			// Global constants

			const BSDF* currentBSDF = its.mesh->getBSDF();
			Vector3f rayToCam = -ray.d.normalized();

			//Get random Emitter
			const std::vector<Emitter *> emitterVector = scene->getEmitters();
			int randomIndex = rand() % emitterVector.size();
			Emitter* currentSampledEmitter = emitterVector.at(randomIndex);

			//Get Le
			Color3f Le = Color3f(0.0f);
			const Emitter* itsEmitter = its.mesh->getEmitter();
			if (itsEmitter) {
				//Its an Emitter.
				EmitterQueryRecord iRec2(ray.o, its.p, its.shFrame.n);
				Le = itsEmitter->eval(iRec2);
			}

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
					* currentSampledEmitter->eval(currentSampledEmitterQuery);//Get Lo (not divided by pdf)

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

				float dist = (its.p - currentSampledEmitterQuery.p).norm();

				float ems_pdf_ie = currentSampledEmitter->pdf(
						currentSampledEmitterQuery) * cos_itsNormal_sEmitter
						/ (emitterVector.size());

				float mat_pdf_ie = currentBSDF->pdf(
						currentBSDFQueryToSampledEmitter);

				w_em = 1 / (ems_pdf_ie + mat_pdf_ie);

				F_em = L_sampledEmitter * bsdfValToSampledEmitter;

			} else {
				//Emitter constribution is 0, because its blocked.
			}

			// BSDF Part ----------------------------------------------------------------------------------------------

			BSDFQueryRecord currentBSDFQueryToSampledDirection(
					its.toLocal(rayToCam));	//wi Camera, wo sampled ray
			currentBSDFQueryToSampledDirection.p = its.p;

			Color3f bsdfVal_cos_pdf = currentBSDF->sample(
					currentBSDFQueryToSampledDirection, sampler->next2D());	//Has cos already in
			Color3f bsdfValToSampledDirection = currentBSDF->eval(
					currentBSDFQueryToSampledDirection);	// No pdf
			//Sampled direction in wo
			float cos_itsNormal_sDirection = std::abs(
					its.shFrame.n.dot(currentBSDFQueryToSampledDirection.wo));
			float currentBSDFpdf_sDirection = currentBSDF->pdf(
					currentBSDFQueryToSampledDirection);//WARNING: PDF not in 0-1
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

					float dist = (its.p - its_sampledDirection.p).norm();

					float ems_pdf_im = emi->pdf(emitterQuerySampledDirection)
							/ (dist * dist * emitterVector.size());
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

			result += (Le + Ld) * multiConst;
			ray = nextRay;

			//----------------------------------------
			//	Russian Roulette
			if (sampler->next1D() > m_q) {
				if (isNextEmitter) {
					run = false;
				} else {
					multiConst *= bsdfVal_cos_pdf / (1 - m_q);
				}
			} else {
				run = false;
			}
		}
	}
	return result;

	/*
	 if (medium) {
	 return volumetricPathTracing(scene, sampler, ray, medium);
	 } else {
	 return otherLi(scene, sampler, ray);
	 }
	 */
}
/*
 Color3f otherLi(const Scene *scene, Sampler *sampler,
 const Ray3f &ray) const {
 //This is Extended PathTracing

 Intersection its;

 // If not visible return black
 if (!scene->rayIntersect(ray, its)) {
 //This should never happen if there is a environmental Map
 return Color3f(0.3f);
 }
 //Get Le
 Color3f Le(0.0f);
 const Emitter* emi2 = its.mesh->getEmitter();
 if (emi2) {
 //Its an Emitter.
 EmitterQueryRecord iRec2(ray.o, its.p, its.shFrame.n);
 Le = emi2->eval(iRec2);
 } else {
 //Too slow if only environmental Light
 Ray3f dummy;
 //	Le = sampleRandomLight(its, scene, sampler, ray, dummy);
 }

 //Get Li
 const BSDF* bsdf = its.mesh->getBSDF();
 Vector3f toCam = -ray.d.normalized();

 BSDFQueryRecord query(its.toLocal(toCam)); //wi Camera, wo sampled ray
 query.p = its.p;
 query.uv = its.uv;
 Color3f bsdfVal = bsdf->sample(query, sampler->next2D());

 //Next path section
 Ray3f newRay(its.p, its.toWorld(query.wo));

 if (sampler->next1D() > m_q)
 return Le + this->Li(scene, sampler, newRay) * bsdfVal / (1.f - m_q);
 else
 return Le;
 }
 */

Color3f sampleRandomLight(Intersection its, const Scene* scene,
		Sampler *sampler, const Ray3f &ray, Ray3f& shadowRayOut) const {
	const std::vector<Emitter *> emis = scene->getEmitters();
	int randomIndex = rand() % emis.size();
	Emitter* emiR = emis.at(randomIndex);

	EmitterQueryRecord emiRec = EmitterQueryRecord(ray.o, its.p, its.shFrame.n);

//Instead of iterating over all emitter, mulply because MC
	Color3f Lo = emis.size() * emiR->sample(emiRec, sampler->next2D());

	BSDFQueryRecord query2(its.toLocal(emiRec.wi), its.toLocal(-ray.d),
			ESolidAngle);
	query2.uv = its.uv;
	const BSDF* bsdf1 = its.mesh->getBSDF();
	Color3f bsdfValemi = bsdf1->eval(query2);

	shadowRayOut = emiRec.shadowRay;

//Check if something blocks the visibility
	if (!scene->rayIntersect(emiRec.shadowRay)) {
		//Multiply by cos of normal of reflec. normal
		float cos0 = std::abs(its.geoFrame.n.dot(emiRec.wi));

		return Lo * bsdfValemi * cos0;
	} else {
		//IF shadowRay intersectis other object,
		//leave Ls as 0.
		return Color3f(0.f);
	}
}

Color3f volumetricPathTracing(const Scene *scene, Sampler *sampler,
		const Ray3f &ray, const Medium* medium) const {

//This is path_mats. Works for Mirror
	Intersection its;

// If not visible return black
	if (!scene->rayIntersect(ray, its)) {
		/*This should never happen if there is a environmental Map*/
		return Color3f(0.3f);
	}

	float tmax = its.t;
	float t = medium->sampleFreeFlightDistance(sampler->next1D()); //Sample free flying path

//TODO: Check fireflies in dielectric/medium -> Snell coefficient?

//TODO: Add Emitter Sampling into Rediance ->Weights?

}

std::string toString() const {
	return "VolumetricPathTracer[]";
}
private:
float m_q;
};

NORI_REGISTER_CLASS(VolumetricPathTracer, "volPathTracer");
NORI_NAMESPACE_END
