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
				////////////////////////
				// Origin is in MEDIUM, do Volumetric Path tracing
				/////////////////////////

				float tmax = its.t;
				float t = medium->sampleFreeFlightDistance(sampler->next1D()); //Sample free flying path
				//float t =tmax;
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

					// Phase function sampling
					float pfVal = medium->sample(query, sampler->next2D()); // this is simply sigma_t/sigma_s
					//query.wo sampled w
					//float pdf = medium->pdf(query);

					//Sample Emitter
					const std::vector<Emitter *> emis = scene->getEmitters();
					int randomIndex = rand() % emis.size();
					Emitter* emiR = emis.at(randomIndex);

					EmitterQueryRecord emiRec = EmitterQueryRecord(its_x.p);

					Color3f Ld = emis.size()
							* emiR->sample(emiRec, sampler->next2D());

					float cos = std::abs(emiRec.n.dot(emiRec.wi));

					//result += Ld * getTransmitance(its_x.p, emiRec.p, scene)*cos/(its_x.p-emiRec.p).squaredNorm(); //Returns 0 if blocking object

					Ray3f newRay(its_x.p, query.wo);

					ray = newRay;
					multiConst *= pfVal;

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
					}

					float Tr = medium->getTransmittanceValue(ray.o, its.p);
					float pdf1 = 1 - Tr; //PDf hitting surface

					//if (Tr == 1) {
					//	Tr = 0;
					//	pdf1 = 1;
					//}

					const BSDF* bsdf = its.mesh->getBSDF();
					Vector3f toCam = -ray.d.normalized();

					BSDFQueryRecord query(its.toLocal(toCam)); //wi Camera, wo sampled ray
					query.p = its.p;
					query.uv = its.uv;
					Color3f bsdfVal = bsdf->sample(query, sampler->next2D());
					//Sampled ray in query.wo

					// Next Ray
					Ray3f newRay(its.p, its.toWorld(query.wo));

					ray = newRay;
					result += Le * multiConst;

					if (sampler->next1D() > m_q) {
						multiConst *= bsdfVal *Tr/ ((1.f - m_q)); //*Tr
					} else {
						run = false;
					}

				}

			} else {
				////////////////////////
				// Origin of Ray is not in Medium, do MIS Path tracer
				/////////////////////////

				// Global constants

				const BSDF* currentBSDF = its.mesh->getBSDF();
				Vector3f rayToCam = -ray.d.normalized();

				//Get random Emitter
				const std::vector<Emitter *> emitterVector =
						scene->getEmitters();
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

				////////////////////////
				// EMiiter Part -------------
				//
				// Sample a random Emitter and get a direction.
				//
				// Get BSDF PDF of going through that direction
				//
				// Weight the values depending on the probabilities.
				/////////////////////////
				//Sample Emitter
				EmitterQueryRecord currentSampledEmitterQuery(its.p);

				//Emitter sampling: wi , ref-> p in Emitter, sampled
				// wo is ref -> cam
				//All local

				currentSampledEmitter->sample(currentSampledEmitterQuery,
						sampler->next2D());

				//Instead of iterating over all emitter
				Color3f L_sampledEmitter = getTransmitance(
						currentSampledEmitterQuery.ref,
						currentSampledEmitterQuery.p, scene)
						* emitterVector.size()
						* currentSampledEmitter->eval(
								currentSampledEmitterQuery);//Get Lo (not divided by pdf)

				BSDFQueryRecord currentBSDFQueryToSampledEmitter(
						its.toLocal(rayToCam),
						its.toLocal(currentSampledEmitterQuery.wi),
						ESolidAngle);

				currentBSDFQueryToSampledEmitter.uv = its.uv;
				Color3f bsdfValToSampledEmitter = currentBSDF->eval(
						currentBSDFQueryToSampledEmitter);

				Color3f F_em = Color3f(0.0f);
				float w_em = 1;

				//Check if something blocks the visibility
				if (!scene->rayIntersect(
						currentSampledEmitterQuery.shadowRay)) {

					float cos_itsNormal_sEmitter =
							std::abs(
									its.shFrame.n.dot(
											currentSampledEmitterQuery.wi.normalized()));

					float dist = (its.p - currentSampledEmitterQuery.p).norm();

					float ems_pdf_ie = currentSampledEmitter->pdf(
							currentSampledEmitterQuery)
							/ (emitterVector.size());//* cos_itsNormal_sEmitter

					float mat_pdf_ie = currentBSDF->pdf(
							currentBSDFQueryToSampledEmitter);

					w_em = 1 / (ems_pdf_ie + mat_pdf_ie);

					F_em = L_sampledEmitter * bsdfValToSampledEmitter;

				} else {
					//Emitter constribution is 0, because its blocked.
				}

				////////////////////////
				// BSDF Part -------------
				//
				// Sample the BSDF and get a direction.
				//
				// If that direction hits an emitter, get the PDF
				// to hit that point
				//
				// Weight the values depending on the probabilities.
				/////////////////////////

				BSDFQueryRecord currentBSDFQueryToSampledDirection(
						its.toLocal(rayToCam));	//wi Camera, wo sampled ray
				currentBSDFQueryToSampledDirection.p = its.p;

				Color3f bsdfVal_cos_pdf = currentBSDF->sample(
						currentBSDFQueryToSampledDirection, sampler->next2D());	//Has cos already in
				Color3f bsdfValToSampledDirection = currentBSDF->eval(
						currentBSDFQueryToSampledDirection);	// No pdf
				//Sampled direction in wo
				float cos_itsNormal_sDirection = std::abs(
						its.shFrame.n.dot(
								currentBSDFQueryToSampledDirection.wo));
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
						Color3f Fo = emi->eval(emitterQuerySampledDirection)
								* getTransmitance(
										emitterQuerySampledDirection.ref,
										emitterQuerySampledDirection.p, scene);

						float dist = (its.p - its_sampledDirection.p).norm();

						float ems_pdf_im = emi->pdf(
								emitterQuerySampledDirection)
								/ (emitterVector.size());//dist * dist *
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

		if(result.x() != result.x())
			return Color3f(0.f);

		return result;

		/*
		 if (medium) {
		 return volumetricPathTracing(scene, sampler, ray, medium);
		 } else {
		 return otherLi(scene, sampler, ray);
		 }
		 */
	}

	/**
	 * Calculates the transmitance in a Ray,
	 * Does NOT check if ray gets blocked, that depends on BSDF.
	 */
	float getTransmitance(Point3f start, Point3f endAtEmitter,
			const Scene * scene) const {

		float totalTrans = 1.f;

		const std::vector<Medium *> ms = scene->getMediums();

		Vector3f d = (start - endAtEmitter).normalized();
		Medium* medium = nullptr;

		while ((start - endAtEmitter).norm() > Epsilon) {

			medium = nullptr;
			// Check if in Medium
			for (unsigned int var = 0; var < ms.size(); ++var) {
				Medium* m = ms.at(var);
				if (m->isInside(start + d * Epsilon)) {
					medium = m;
					break;
				}
			}

			Intersection its;
			Ray3f ray(start, d, Epsilon,
					(endAtEmitter - start).norm() - Epsilon);

			if (medium) {
				//Start in Medium, find out OutgoingPoint

				if (scene->rayIntersect(ray, its)) {
					//Outgoingpoint is its.p
					totalTrans *= medium->getTransmittanceValue(start, its.p);
					start = its.p + d * Epsilon;
				} else {
					//Outgoing point not found, endAtEmitter reached
					totalTrans *= medium->getTransmittanceValue(start,
							endAtEmitter);
					start = endAtEmitter;
				}
			} else {
				//Not in medium, search next intersection
				if (scene->rayIntersect(ray, its)) {
					start = its.p;
				} else {
					//end reached, not in medium, dont do anything
					start = endAtEmitter;
				}
			}

		}
		return totalTrans;

	}
	std::string toString() const {
		return "VolumetricPathTracer[]";
	}
private:
	float m_q;
};

NORI_REGISTER_CLASS(VolumetricPathTracer, "volPathTracer");
NORI_NAMESPACE_END
