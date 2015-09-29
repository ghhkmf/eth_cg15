#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/ray.h>

NORI_NAMESPACE_BEGIN

class AverageVisibilityIntegrator : public Integrator {

public:
    AverageVisibilityIntegrator(const PropertyList &props) {
        /* No parameters this time */
    	m_rayLength = props.getFloat("length",1);
  }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
       
	

	 /* Find the surface that is visible in the requested direction */
        Intersection its;
        /* if no intersection: Set to white */
	if (!scene->rayIntersect(ray, its))
            return Color3f(1.0f);

	/* ray_direction is already in unit Hemisphere*/	
	Vector3f ray_direction = Warp::sampleUniformHemisphere(sampler, its.shFrame.n);


	/* Create ray with length: m_rayLength */	
	Ray3f second_ray(its.p,ray_direction,m_t_threshold,m_rayLength);


	Intersection its2;
	/* Traverse Ray, check if occluded*/


	if(scene->rayIntersect(second_ray,its2)){
			return Color3f(0.0f);
	}

	return Color3f(1.0f);
    }

    std::string toString() const {
        return tfm::format(
                "AverageVisibilityIntegrator[ \n"
                " Ray Length = \"%s\" \n"
                " ] ",
                m_rayLength);
    }

private:
  float m_rayLength;
  float m_t_threshold=1e-4;
};

NORI_REGISTER_CLASS(AverageVisibilityIntegrator, "av");
NORI_NAMESPACE_END
