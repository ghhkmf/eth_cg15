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
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

/**
 * Translucency, sometimes known as semi-transparency, is a form of transparency.
 * It allows light to pass through, but unlike a transparent object, it does not
 * allow said light to form shapes. This happens because the photons scatter and
 * the image becomes less focused.
 */

class Translucent: public BSDF {
public:
	Translucent(const PropertyList &propList) { /* Interior IOR (default: BK7 borosilicate optical glass) */
		/* Interior IOR (default: BK7 borosilicate optical glass) */
		m_intIOR = propList.getFloat("intIOR", 1.5046f);

		/* Exterior IOR (default: air) */
		m_extIOR = propList.getFloat("extIOR", 1.000277f);

		m_albedo = propList.getColor("albedo");
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
			return Color3f(1.f) * m_albedo;
		else
			return snell * snell * Color3f(1.f) * m_albedo;
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

		if (isEqual(bRec.wo, reflec))
			prob = 0.5f;
		else {
			Normal3f pole = refrac;
			MatrixXf rot = Warp::getRotationMatrix(pole,
					Vector3f(0.f, 0.f, 1.f));

			prob = Warp::squareToCosineHemispherePdf(rot * bRec.wo) * 0.5f;
		}
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
			return Color3f(1.f) * m_albedo; //Montecarlo deterministic

		} else {
			bRec.wo = getRefraction(bRec);
			Point2f nSample = sample;
			nSample[0] = nSample[0] / F;

			Vector3f q = Warp::squareToCosineHemisphere(nSample);

			Vector3f pole = bRec.wo;
			MatrixXf rot = Warp::getRotationMatrix(Vector3f(0.f, 0.f, 1.f),
					pole);

			bRec.wo = rot * q;

			//return snell * snell * (1 - F) * Color3f(1.f)/cos;
			return snell * snell * Color3f(1.f) * m_albedo; //Montecarlo deterministic (1-F)/pdf=(1-F)
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
	Color3f m_albedo;

	Vector3f getOutgoingVector(const BSDFQueryRecord &bRec,
			bool returnRefract) const {
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

		refraction = -bRec.wi;

		//Reflection
		Vector3f reflection = Vector3f(-bRec.wi.x(), -bRec.wi.y(), bRec.wi.z());

		if (returnRefract)
			return refraction;
		else
			return reflection;
	}

	Vector3f getReflection(const BSDFQueryRecord &bRec) const {
		return getOutgoingVector(bRec, false);
	}

	Vector3f getRefraction(const BSDFQueryRecord &bRec) const {
		return getOutgoingVector(bRec, true);
	}

	bool isEqual(Vector3f a, Vector3f b) const {

		bool diffA = std::abs(a.x() - b.x()) < Epsilon;
		bool diffB = std::abs(a.y() - b.y()) < Epsilon;
		bool diffC = std::abs(a.z() - b.z()) < Epsilon;
		return diffA && diffB && diffC;
	}

}
;

NORI_REGISTER_CLASS(Translucent, "trans");
NORI_NAMESPACE_END
