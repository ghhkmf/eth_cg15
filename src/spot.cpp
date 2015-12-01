#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

class spotLight : public Emitter{
public:
	spotLight(const PropertyList &propList) {
		m_radiance = propList.getColor("power");
		position = propList.getPoint3("position");
		direction = propList.getVector3("direction");
		cosFalloffStart = propList.getFloat("cosFalloffStart");
		cosTotalWidth = propList.getFloat("cosTotalWidth");

	}	

	virtual std::string toString() const override {
		return "spotLight[]";
	}

	virtual Color3f samplePhoton(Ray3f &ray, const Point2f &sample1, const Point2f &sample2) const override {
		ShapeQueryRecord sRec;
		
		//sample direction from light
		Vector3f q = Warp::squareToUniformSphereCap(sample1,cosTotalWidth);
		Frame frame(direction.normalized());
		Vector3f sampleDir = frame.toWorld(q);
		ray = Ray3f(sRec.p, sampleDir);

		//compute fallOff
		float cosTheta = sampleDir.normalized().dot(direction.normalized());
		float fA = fallOff(cosTheta, cosFalloffStart, cosTotalWidth);

		ray = Ray3f(position, sampleDir);
		Color3f finals = m_radiance * fA;// / Warp::squareToUniformSphereCapPdf(Vector3f(0, 0, 1), cosTotalWidth);
		//cout << " pdf " << Warp::squareToUniformSphereCapPdf(Vector3f(0, 1, 0), cosTotalWidth);
		return finals;
	}
	
	virtual Color3f sample(EmitterQueryRecord &lRec, const Point2f &sample)const override {
		
		lRec.p = position;
		lRec.wi = (position - lRec.ref).normalized();
		lRec.n = direction.normalized();

		Ray3f shadowRay(lRec.ref, lRec.wi / lRec.wi.norm(), Epsilon, (position - lRec.ref).norm() -Epsilon);
		lRec.shadowRay = shadowRay;
		return eval(lRec) / pdf(lRec);
	}

	virtual float fallOff(float cosTh, float cosF, float cosTo) const override{
		if (cosTh < cosTo) {
			return 0;
		}
		else if (cosTh > cosF) {
			return 1;
		}
		else {
			float delta = (cosTh - cosTo) / (cosF - cosTo);
			return delta*delta*delta*delta;				
		}
	}

	virtual Color3f eval(const EmitterQueryRecord &lRec) const override {
		float cosTheta = direction.normalized().dot((lRec.ref - position).normalized());
		float fA = fallOff(cosTheta, cosFalloffStart, cosTotalWidth);
		return m_radiance * fA / (2.0f * M_PI) 
				* (1.0f - 0.5f*(cosFalloffStart + cosTotalWidth));
	}

	virtual float pdf(const EmitterQueryRecord &lRec)const override {
		float pdfVal = Warp::squareToUniformSphereCapPdf(Vector3f(0,0,1),cosTotalWidth);
		float cosThetaLight = abs(lRec.n.dot(-lRec.wi));
		return pdfVal*(lRec.p - lRec.ref).squaredNorm() / cosThetaLight;
	}

private:
	Point3f position;
	Color3f m_radiance;
	Vector3f direction;
	float cosTotalWidth;
	float cosFalloffStart;
};

NORI_REGISTER_CLASS(spotLight, "spot");
NORI_NAMESPACE_END