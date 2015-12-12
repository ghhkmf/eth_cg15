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

class Anisotropic: public BSDF {
public:
	Anisotropic(const PropertyList &propList) {
		ex = propList.getFloat("ex", 1.0f);
		ey = propList.getFloat("ey", 1.0f);

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
		m_ks = 1 - (m_kd.x() + m_kd.y() + m_kd.z()) / 3;
	}

	/// Evaluate the anisotropic normal distribution D
	float evalAnisotropic(const Normal3f &wh) const {
		float costhetah = Frame::cosTheta(wh);
		float d = 1.0f - costhetah*costhetah;
		if (d == 0.0f) return 0.0f;
		float e = (ex*wh.x()*wh.x() + ey*wh.y()*wh.y()) / d;

		return sqrtf((ex+2.0f)*(ey+2.0f)) * INV_TWOPI * powf(costhetah, e);
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

		float Dh = evalAnisotropic(wh); //cout << "Dh " << Dh << " ";
		float F = fresnel(wh.dot(bRec.wi), m_extIOR, m_intIOR);
		float G = smithBeckmannG1(bRec.wi, wh) * smithBeckmannG1(bRec.wo, wh);
		
		float cos2 = cos_theta_i * cos_theta_o;
		if (cos2 < 0) {
			cos2 = -cos2;
		}
		//return m_kd / M_PI + m_ks * (Dh * F ) / (4.f * wh.dot(bRec.wi));
		return m_kd / M_PI + m_ks * (Dh * F * G) / (4.f * cos2);

	}
	/// Evaluate the sampling density of \ref sample() wrt. solid angles
	virtual float pdf(const BSDFQueryRecord &bRec) const override {
		//cout << "costheta" << bRec.wo.z()<<" ";
		if (bRec.measure != ESolidAngle || Frame::cosTheta(bRec.wi) <= 0
				|| Frame::cosTheta(bRec.wo) <= 0) {
			return 0.0f;
		}
		Normal3f wh = (bRec.wi + bRec.wo).normalized();
		float costhetah = abs(wh.z());
		float ds = 1.0f - costhetah*costhetah;
		float anisoPdf = 0.0f;
		if (ds > 0.0f && bRec.wi.dot(wh) > 0.0f) {
			float e = (ex*wh.x()*wh.x() + ey*wh.y()*wh.y()) / ds;
			float d = sqrtf((ex + 1.0f) * (ey + 1.0f))*INV_TWOPI * powf(costhetah, e);
			anisoPdf = d / (4.0f*bRec.wi.dot(wh));
		}

		return anisoPdf;
	}
	virtual bool sameHemisphere(Vector3f w, Vector3f wp) const{
		return w.z()*wp.z() > 0.0f;
	}

	/// Sample the BRDF
	virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const
			override {

		bRec.measure = ESolidAngle;
		Point2f sample(_sample); 
		float u1 = sample.x();float u2 = sample.y();
		float costheta;
		float phi;
		////////////////
		//sample from first quadrant and remap to hemisphere to sample wh
		if (u1 < .25f) {
			sampleFirstQuadrant(4.0f*u1, u2, phi, costheta);
		}
		else if (u1 < 0.5f) {
			u1 = 4.0f*(.5f - u1);
			sampleFirstQuadrant(u1, u2, phi, costheta);
			phi = M_PI - phi;
		}
		else if (u1 < 0.75f) {
			u1 = 4.0f * (u1 - 0.5f);
			sampleFirstQuadrant(u1, u2, phi, costheta);
			phi += M_PI;
		}
		else {
			u1 = 4.0f * (1.0f - u1);
			sampleFirstQuadrant(u1, u2, phi, costheta);
			phi = 2.0f*M_PI - phi;
		}
		float sintheta = sqrtf(std::max(0.0f, 1.0f - costheta*costheta));
		Vector3f wh(sintheta*cosf(phi),sintheta*sinf(phi),costheta);
		//mirror wh if not same sphere as wi
		if (!sameHemisphere(bRec.wi, wh)) wh = -wh;

		///////////////
		//compute incident direction by reflecting about wh
		bRec.wo = -bRec.wi + 2.0f*bRec.wi.dot(wh)*wh;
		if (pdf(bRec) > 0)
			return eval(bRec) / pdf(bRec);// *Frame::cosTheta(bRec.wo);
		else
			return Color3f(0.f);
	}

	void sampleFirstQuadrant(float u1, float u2, float &phi, float&costheta) const {
		if (ex == ey) {
			phi = M_PI * u1 * 0.5f;
		}
		else {
			phi = atanf(sqrtf((ex + 1.0f) / (ey + 1.0f)) * tanf(M_PI * u1 * 0.5f));
		}
		float cosphi = cosf(phi);
		float sinphi = sinf(phi);
		costheta = powf(u2, 1.0f / (ex*cosphi*cosphi + ey*sinphi*sinphi + 1));
	}

	virtual std::string toString() const override {
		return tfm::format("Anisotropic[\n"
				"  alpha = %f,\n"
				"  intIOR = %f,\n"
				"  extIOR = %f,\n"
				"  kd = %s,\n"
				"  ks = %f\n"
				"]", m_alpha, m_intIOR, m_extIOR, m_kd.toString(), m_ks);
	}
private:
	float ex; float ey;
	float m_alpha;
	float m_intIOR, m_extIOR;
	float m_ks;
	Color3f m_kd;
};

NORI_REGISTER_CLASS(Anisotropic, "anisotropic");
NORI_NAMESPACE_END
