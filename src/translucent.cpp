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
		m_albedo = propList.getColor("albedo");

	}

	virtual Color3f eval(const BSDFQueryRecord &bRec) const override {
		return m_albedo;
	}

	virtual float pdf(const BSDFQueryRecord &bRec) const override {
		return Warp::squareToUniformHemispherePdf(bRec.wo);
	}

	virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const
			override {
		bRec.measure = ESolidAngle;

		//Scatter uniform inside object
		bRec.wo = -1 * Warp::squareToUniformHemisphere(sample);

		return m_albedo;
	}

	virtual std::string toString() const override {
		return tfm::format("Translucent[\n"
				"]");
	}
private:
	Color3f m_albedo;
}
;

NORI_REGISTER_CLASS(Translucent, "trans");
NORI_NAMESPACE_END
