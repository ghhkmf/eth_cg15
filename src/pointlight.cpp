#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

class PointLight: public Emitter {
public:
	PointLight(const PropertyList &propList) {
		/* No parameters this time */
		m_position = propList.getPoint3("position", Point3f());
		m_power = propList.getColor("power", 1.f);
	}

	Color3f getRadiance(Point3f p) {
		float dist = (p - m_position).norm();
		return m_power / (4 * M_PI * dist * dist);
	}
	;

	Vector3f getWi(Point3f o) {
		return (m_position - o).normalized();
	}
	;

	Ray3f getRayToPoint(Point3f p) {
		return Ray3f(m_position, (p - m_position).normalized(), Epsilon,
				(p - m_position).norm() - Epsilon);
	}

	Color3f sample(EmitterQueryRecord &lRec, const Point2f &sample) const {
		lRec.p = m_position;
		lRec.wi = (m_position - lRec.ref).normalized();
		lRec.n=-lRec.wi;
		lRec.pdf=1.f;
		lRec.shadowRay=Ray3f(lRec.ref,lRec.wi,Epsilon,(lRec.p-lRec.ref).norm());
		return eval(lRec) / pdf(lRec);
	}
	Color3f eval(const EmitterQueryRecord &lRec) const {
		if (lRec.p == m_position) {
			return m_power / (4 * M_PI);
		} else {
			return Color3f(0.f);
		}
	}
	float pdf(const EmitterQueryRecord &lRec) const {
		return (lRec.p == m_position) ? 1.f : 0.f;
	}

	EmittedValues getEmittedValues(Point3f p) {
		EmittedValues val;
		val.p = p;
		val.Li = getRadiance(p);
		val.ray = getRayToPoint(p);
		val.wi = getWi(p);
		return val;
	}

	std::string toString() const {
		return tfm::format("PointLight[ \n"
				" Position = \"%s\" \n"
				" Power = \"%s\" \n"
				" ] ", m_position.toString(), m_power.toString());
	}

protected:
	Point3f m_position;
	Color3f m_power;
};

NORI_REGISTER_CLASS(PointLight, "point");
NORI_NAMESPACE_END
