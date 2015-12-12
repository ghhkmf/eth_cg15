/*
 This file is part of Nori, a simple educational ray tracer

 Copyright (c) 2015 by Wenzel Jakob

 Nori is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License Version 3
 as published by the Free Software Foundation.

 Nori is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <nori/bsdf.h>
#include <nori/frame.h>
NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric: public BSDF {
public:
	Dielectric(const PropertyList &propList) {
		/* Interior IOR (default: BK7 borosilicate optical glass) */
		m_intIOR = propList.getFloat("intIOR", 1.5046f);

		/* Exterior IOR (default: air) */
		m_extIOR = propList.getFloat("extIOR", 1.000277f);
	}

	virtual Color3f eval(const BSDFQueryRecord &bRec) const override {

		Vector3f reflec = getReflection(bRec);
		Vector3f refrac = getRefraction(bRec);

		bool isRefl = bRec.wo == reflec;
		bool isRefr = bRec.wo == refrac;

		float n1 = m_extIOR;
		float n2 = m_intIOR;

		float cosT = bRec.wi.z();

		if (cosT < 0) {
			n1 = m_intIOR;
			n2 = m_extIOR;
		}

		//Sample depending on Fresnel
		float snell = n1 / n2;

		if (isRefl)
			return Color3f(1.f);
		else if (isRefr)
			return snell * snell * Color3f(1.f);
		else
			return Color3f(0.0f);
	}


	virtual float pdf(const BSDFQueryRecord &bRec) const override {

		Vector3f n(0, 0, 1);
		float n1 = m_extIOR;
		float n2 = m_intIOR;

		float cosT = bRec.wi.z();

		if (cosT < 0) {
			n1 = m_intIOR;
			n2 = m_extIOR;
			n = -n;
			cosT = -cosT;
		}

		//Sample depending on Fresnel
		float F = fresnel(cosT, n1, n2);

		Vector3f refrac = getRefraction(bRec);
		Vector3f reflec = getReflection(bRec);

		float prob = 0;

		if (isEqual(bRec.wo, refrac) || isEqual(bRec.wo, reflec))
			prob = 1.f;
		else
			prob = 0.f;

		return prob;
	}

	virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const
			override {

		Vector3f n(0, 0, 1);
		float n1 = m_extIOR;
		float n2 = m_intIOR;

		float cosT = bRec.wi.z();

		if (cosT < 0) {
			n1 = m_intIOR;
			n2 = m_extIOR;
			n = -n;
			cosT = -cosT;
		}

		//Sample depending on Fresnel
		float F = fresnel(cosT, n1, n2);
		float snell = n1 / n2;

		float cons = sqrt(1.0f - (snell * snell) * (1.0f - cosT * cosT));

		// Total Reflection?
		bool isTReflec = cons > 1;

		if (sample.x() < F || isTReflec) {
			bRec.wo = getReflection(bRec);
			//return Color3f(1.f)*F;///Color3f(1.f)*F cosT;
			return Color3f(1.f); //Montecarlo deterministic

		} else {
			bRec.wo = getRefraction(bRec);
			//return snell * snell * (1 - F) * Color3f(1.f)/cos;
			return snell * snell * Color3f(1.f); //Montecarlo deterministic (1-F)/pdf=(1-F)
//			return Color3f(1.f);
		}

	}



	virtual std::string toString() const override {
		return tfm::format("Dielectric[\n"
				"  intIOR = %f,\n"
				"  extIOR = %f\n"
				"]", m_intIOR, m_extIOR);
	}
private:
	float m_intIOR, m_extIOR;


	Vector3f getOutgoingVector(const BSDFQueryRecord &bRec, bool returnRefract) const{
			Vector3f n(0, 0, 1);
			float n1 = m_extIOR;
			float n2 = m_intIOR;

			float cosT = bRec.wi.z();

			if (cosT < 0) {
				n1 = m_intIOR;
				n2 = m_extIOR;
				n = -n;
				cosT = -cosT;
			}

			//Sample depending on Fresnel
		//	float F = fresnel(cosT, n1, n2);
			float snell = n1 / n2;

			//Refraction
			Vector3f part1 = -snell * (bRec.wi - cosT * n);
			float cons = sqrt(1.0f - (snell * snell) * (1.0f - cosT * cosT));
			Vector3f part2 = n * cons;

			Vector3f refraction = part1 - part2;
			refraction = refraction.normalized();

			//Reflection
			Vector3f reflection = Vector3f(-bRec.wi.x(), -bRec.wi.y(), bRec.wi.z());

			if (returnRefract)
				return refraction;
			else
				return reflection;
		}

		Vector3f getReflection(const BSDFQueryRecord &bRec) const{
			return getOutgoingVector(bRec, false);
		}

		Vector3f getRefraction(const BSDFQueryRecord &bRec) const{
			return getOutgoingVector(bRec, true);
		}

		bool isEqual(Vector3f a, Vector3f b) const{

			bool diffA =  std::abs(a.x()-b.x())<Epsilon;
			bool diffB =  std::abs(a.y()-b.y())<Epsilon;
			bool diffC =  std::abs(a.z()-b.z())<Epsilon;
			return diffA && diffB && diffC;
		}

}
;

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END
