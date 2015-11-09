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
#include <math.h>

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

	virtual Color3f eval(const BSDFQueryRecord &) const override {
		/* Discrete BRDFs always evaluate to zero in Nori */
		return Color3f(0.0f);
	}

	virtual float pdf(const BSDFQueryRecord &) const override {
		/* Discrete BRDFs always evaluate to zero in Nori */
		return 0.0f;
	}

	virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const
			override {

		if (Frame::cosTheta(bRec.wi) <= 0)
			return Color3f(0.0f);

		//Sample depending on Fresnel
		float F = fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR);

		if (sample.x() <= F) {
			//Refraction
			float snell = m_extIOR/m_intIOR;
			Vector3f n(0,0,1);

			// Reflection in local coordinates
			bRec.wo = -snell*(bRec.wi-Frame::cosTheta(bRec.wi)*n)-n*(sqrt(1-snell*snell*(1-Frame::cosTheta(bRec.wi)*Frame::cosTheta(bRec.wi))));

			bRec.measure = EDiscrete;

			/* Relative index of refraction: no change */
			bRec.eta = 1.0f;
		} else {
			//Reflection
			// Reflection in local coordinates
			bRec.wo = Vector3f(-bRec.wi.x(), -bRec.wi.y(), bRec.wi.z());
			bRec.measure = EDiscrete;

			/* Relative index of reflection: no change */
			bRec.eta = 1.0f;
		}

		return Color3f(1.0f);

	}

	virtual std::string toString() const override {
		return tfm::format("Dielectric[\n"
				"  intIOR = %f,\n"
				"  extIOR = %f\n"
				"]", m_intIOR, m_extIOR);
	}
private:
	float m_intIOR, m_extIOR;
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END
