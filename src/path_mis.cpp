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
		m_q = 0.9; //TODO: Chech this
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

		//Get Ld - No Ld only if emitter hit

		Color3f Ld = Color3f(0.f);
		const std::vector<Emitter *> emis = scene->getEmitters();

		int randomIndex = rand() % emis.size();

		Emitter* emi = emis.at(randomIndex);
		//Sample Emitter
		EmitterQueryRecord iRec(its.p);
		Ld = emis.size()*emi->sample(iRec, sampler->next2D()); //Instead of iterating over all emitter, mulply because MC
		float cos_i = std::abs(iRec.n.dot(-iRec.wi));
		float cos0 = std::abs(its.geoFrame.n.dot(iRec.wi));
		Ld = Ld*cos0*cos_i/ (its.p - iRec.p).squaredNorm();

		//Get Le
		Color3f Le = Color3f(0.0f);
		const Emitter* emi2 = its.mesh->getEmitter();
		if (emi2) {
			//Its an Emitter.
			EmitterQueryRecord iRec2(its.p);
			Le = emi2->eval(iRec2);
			return Le / (1 - m_q);
		}

		const BSDF* bsdf = its.mesh->getBSDF();
		Vector3f toCam = -ray.d.normalized();

		BSDFQueryRecord query(its.toLocal(toCam)); //wi Camera, wo sampled ray
		query.p = its.p;
		Color3f bsdfVal = bsdf->sample(query, sampler->next2D());

		//Check if intersect with emitter
		Ray3f lightRay(its.p, its.toWorld(query.wo));

		Color3f Li = Color3f(0.f);

		if (sampler->next1D() > m_q)
			Li = this->Li(scene, sampler, lightRay, bsdf->isDiffuse())
					* bsdfVal;

		if (isDiffuseBounce)
			return (Li + Ld) / (1 - m_q);

		return (Le + Li + Ld) / (1 - m_q);

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
