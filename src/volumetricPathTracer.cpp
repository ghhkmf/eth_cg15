#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/Medium.h>
#include <nori/sampler.h>

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
		Medium* medium=nullptr;
		const std::vector<Medium *> ms = scene->getMediums();
		for (unsigned long var = 0; var < ms.size(); ++var) {
			Medium* m = ms.at(var);
			if (m->isInside(ray.o)) {
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

	Color3f otherLi(const Scene *scene, Sampler *sampler, const Ray3f &ray) const{
		//This is path_mats. Works for Mirror
		Intersection its;

		// If not visible return black
		if (!scene->rayIntersect(ray, its)) {
			/*This should never happen if there is a environmental Map*/
			return Color3f(0.3f);
		}
		//Get Le
		Color3f Le = Color3f(0.0f);
		const Emitter* emi2 = its.mesh->getEmitter();
		if (emi2) {
			//Its an Emitter.
			EmitterQueryRecord iRec2(ray.o, its.p, its.shFrame.n);
			Le = emi2->eval(iRec2);
		}

		float tmax = its.t;

		// No Ld

		//Get Li
		const BSDF* bsdf = its.mesh->getBSDF();
		Vector3f toCam = -ray.d.normalized();

		BSDFQueryRecord query(its.toLocal(toCam)); //wi Camera, wo sampled ray
		query.p = its.p;
		Color3f bsdfVal = bsdf->sample(query, sampler->next2D());

		//Check if intersect with emitter
		Ray3f lightRay(its.p, its.toWorld(query.wo));

		if (sampler->next1D() > m_q)
			return Le
					+ this->Li(scene, sampler, lightRay) * bsdfVal / (1.f - m_q);
		else
			return Le;
	}

	Color3f volumetricPathTracing(const Scene *scene, Sampler *sampler,
			const Ray3f &ray, const Medium* medium) const{
		//TODO: DO IT
			return Color3f(0.3f);
	}

	std::string toString() const {
		return "VolumetricPathTracer[]";
	}
private:
	float m_q;
};

NORI_REGISTER_CLASS(VolumetricPathTracer, "volPathTracer");
NORI_NAMESPACE_END
