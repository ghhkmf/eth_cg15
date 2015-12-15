/*
 9    This file is part of Nori, a simple educational ray tracer

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

class Microfacet: public BSDF {
public:
	Microfacet(const PropertyList &propList) {
		/* RMS surface roughness */
		m_alpha = propList.getFloat("alpha", 0.1f);

		/* Interior IOR (default: BK7 borosilicate optical glass) */
		m_intIOR = propList.getFloat("intIOR", 1.5046f);

		/* Exterior IOR (default: air) */
		m_extIOR = propList.getFloat("extIOR", 1.000277f);

		/* Albedo of the diffuse base material (a.k.a "kd") */
		m_kd = propList.getColor("kd", Color3f(0.5f));

		/* To ensure energy conservation, we must scale the
		 specular component by 1-kd.

		 While that is not a particularly realistic model of what
		 happens in reality, this will greatly simplify the
		 implementation. Please see the course staff if you're
		 interested in implementing a more realistic version
		 of this BRDF. */
		m_ks = 1 - m_kd.maxCoeff();
	}

	/// Evaluate the microfacet normal distribution D
	float evalBeckmann(const Normal3f &m) const {
		float temp = Frame::tanTheta(m) / m_alpha, ct = Frame::cosTheta(m),
				ct2 = ct * ct;

		return std::exp(-temp * temp) / (M_PI * m_alpha * m_alpha * ct2 * ct2);
	}

	/// Evaluate Smith's shadowing-masking function G1
	float smithBeckmannG1(const Vector3f &v, const Normal3f &m) const {
		float tanTheta = Frame::tanTheta(v);

		/* Perpendicular incidence -- no shadowing/masking */
		if (tanTheta == 0.0f)
			return 1.0f;

		/* Can't see the back side from the front and vice versa */
		if (m.dot(v) * Frame::cosTheta(v) <= 0)
			return 0.0f;

		float a = 1.0f / (m_alpha * tanTheta);
		if (a >= 1.6f)
			return 1.0f;
		float a2 = a * a;

		/* Use a fast and accurate (<0.35% rel. error) rational
		 approximation to the shadowing-masking function */
		return (3.535f * a + 2.181f * a2) / (1.0f + 2.276f * a + 2.577f * a2);
	}

	/// Evaluate the BRDF for the given pair of directions
	virtual Color3f eval(const BSDFQueryRecord &bRec) const override {

		if (bRec.measure != ESolidAngle || Frame::cosTheta(bRec.wi) <= 0
				|| Frame::cosTheta(bRec.wo) <= 0)
			return Color3f(0.0f);

		float cos_theta_i = Frame::cosTheta(bRec.wi);
		float cos_theta_o = Frame::cosTheta(bRec.wo);

		Normal3f wh = (bRec.wi + bRec.wo).normalized();

		float Dh = evalBeckmann(wh);
		float F = fresnel(wh.dot(bRec.wi), m_extIOR, m_intIOR);
		float G = smithBeckmannG1(bRec.wi, wh) * smithBeckmannG1(bRec.wo, wh);

		float cos2 = cos_theta_i * cos_theta_o;
		if (cos2 < 0) {
			cos2 = -cos2;
		}

		return m_kd / M_PI + m_ks * (Dh * F * G) / (4.f * cos2);

	}
	/// Evaluate the sampling density of \ref sample() wrt. solid angles
	virtual float pdf(const BSDFQueryRecord &bRec) const override {
		/*
		if (bRec.measure != ESolidAngle || Frame::cosTheta(bRec.wi) <= 0
				|| Frame::cosTheta(bRec.wo) <= 0) {
			return 0.0f;
		}
		*/
		Normal3f wh = (bRec.wi + bRec.wo).normalized();

		float Dh = evalBeckmann(wh);

		float cos_theta_h = Frame::cosTheta(wh);
		float cos_theta_o = Frame::cosTheta(bRec.wo);

		float J = 1 / (4 * wh.dot(bRec.wo));

		return m_ks * Dh * cos_theta_h * J + (1 - m_ks) * cos_theta_o / M_PI;
	}

	/// Sample the BRDF
	virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const
			override {

		bRec.measure = ESolidAngle;

		Point2f sample(_sample);
		if (sample(0) <= m_ks) {
			//Specular
			sample(0) = sample(0) / m_ks; // transform sample into range [0;1]
			Normal3f n = Warp::squareToBeckmann(sample, m_alpha);
			bRec.wo = 2 * n.dot(bRec.wi) * n - bRec.wi;
		} else {
			//Diffuse
			sample(0) = (sample(0) - m_ks) / (1 - m_ks); // transform sample into range [0;1]
			bRec.wo = Warp::squareToCosineHemisphere(sample);
		}
		/*
		if (Frame::cosTheta(bRec.wi) <=0 || Frame::cosTheta(bRec.wo) <= 0)
			return Color3f(0.0f);
			*/
		if (pdf(bRec) > 0)
			return eval(bRec) / pdf(bRec) * Frame::cosTheta(bRec.wo);
		else
			return Color3f(0.f);
	}

	virtual std::string toString() const override {
		return tfm::format("Microfacet[\n"
				"  alpha = %f,\n"
				"  intIOR = %f,\n"
				"  extIOR = %f,\n"
				"  kd = %s,\n"
				"  ks = %f\n"
				"]", m_alpha, m_intIOR, m_extIOR, m_kd.toString(), m_ks);
	}
private:
	float m_alpha;
	float m_intIOR, m_extIOR;
	float m_ks;
	Color3f m_kd;
};

NORI_REGISTER_CLASS(Microfacet, "microfacet");
NORI_NAMESPACE_END
