#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/Medium.h>
#include <nori/sampler.h>
#include <Eigen/Dense>

NORI_NAMESPACE_BEGIN

class VolumetricPathTracer: public Integrator {
public:
	VolumetricPathTracer(const PropertyList &props) {
		m_q = 0.05;
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {

		//Check if we are in a Medium
		// Yes -> Volumetric path tracing
		// No -> Other Integrator
		Medium* medium = nullptr;
		const std::vector<Medium *> ms = scene->getMediums();
		for (unsigned int var = 0; var < ms.size(); ++var) {
			Medium* m = ms.at(var);
			if (m->isInside(ray.o)) { //TODO: Handle border cases
				medium = m;
				break;
			}
		}
		if (medium) {
			return volumetricPathTracing(scene, sampler, ray, medium);
		} else {
			return otherLi(scene, sampler, ray);
		}

	}

	Color3f otherLi(const Scene *scene, Sampler *sampler,
			const Ray3f &ray) const {
		//This is Extended PathTracing

		Intersection its;

		// If not visible return black
		if (!scene->rayIntersect(ray, its)) {
			/*This should never happen if there is a environmental Map*/
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
		Color3f bsdfVal = bsdf->sample(query, sampler->next2D());

		//Next path section
		Ray3f newRay(its.p, its.toWorld(query.wo));

		if (sampler->next1D() > m_q)
			return Le + this->Li(scene, sampler, newRay) * bsdfVal / (1.f - m_q);
		else
			return Le;
	}

	Color3f sampleRandomLight(Intersection its, const Scene* scene,
			Sampler *sampler, const Ray3f &ray, Ray3f& shadowRayOut) const {
		const std::vector<Emitter *> emis = scene->getEmitters();
		int randomIndex = rand() % emis.size();
		Emitter* emiR = emis.at(randomIndex);

		EmitterQueryRecord emiRec = EmitterQueryRecord(ray.o, its.p,
				its.shFrame.n);

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

		if (t < tmax) {
			//Volume interaction
			// x=itp.s
			Intersection its_x;
			its_x.p = ray.o + ray.d * t;
			its_x.t = t;

			its_x.geoFrame.n = (Normal3f) ray.d.cross(Vector3f(0, 0, 1));
			its_x.shFrame.n = its_x.geoFrame.n;

			MediumQueryRecord query(ray.o, its_x.p);
			float pfVal = medium->sample(query, sampler->next2D()); // this is simply sigma_t/sigma_s
			//query.wo sampled w
			//float pdf = medium->pdf(query);

			Ray3f newRay(its_x.p, query.wo);

			//	cout << pfVal << " " << newRay.o << " " << newRay.d << endl;

			return pfVal * this->Li(scene, sampler, newRay);

			//Emitter sampling?
			/*
			 if (sampler->next1D() > m_q)
			 return pfVal * this->Li(scene, sampler, newRay) / (1.f - m_q);
			 else
			 return Color3f(0.f);
			 */
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

			/*	Le = sampleRandomLight(its, scene, sampler, ray, shadowRay);

				Intersection itsEmi;
				if (scene->rayIntersect(shadowRay, itsEmi)) {
					float TrEmitter = medium->getTransmittanceValue(shadowRay.o,
							itsEmi.p);
					Le=Le*TrEmitter;
				}
			*/
			}

			float Tr = medium->getTransmittanceValue(ray.o, its.p);
			float pdf1 = 1 - Tr; //PDf hitting surface

			if (Tr == 1){
				Tr=0;
				pdf1 = 1;
			}
			const BSDF* bsdf = its.mesh->getBSDF();
			Vector3f toCam = -ray.d.normalized();

			BSDFQueryRecord query(its.toLocal(toCam)); //wi Camera, wo sampled ray
			query.p = its.p;
			Color3f bsdfVal = bsdf->sample(query, sampler->next2D());

			//Sampled ray in query.wo

			// Next Ray
			Ray3f newRay(its.p, its.toWorld(query.wo));
//				return bsdfVal * Tr * this->Li(scene, sampler, newRay)
//						/ (pdf1 * pdf2);

			if (sampler->next1D() > m_q)
				return Le
						+ bsdfVal * this->Li(scene, sampler, newRay) * Tr
								/ ((1.f - m_q) * pdf1);
			else
				return Le;

		}
	}

	std::string toString() const {
		return "VolumetricPathTracer[]";
	}
private:
	float m_q;
};

NORI_REGISTER_CLASS(VolumetricPathTracer, "volPathTracer");
NORI_NAMESPACE_END
