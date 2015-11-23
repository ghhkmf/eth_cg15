/*
 This file is part of Nori, a simple educational ray tracer

 Copyright (c) 2015 by Romain Prévost

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

#include <nori/emitter.h>
#include <nori/warp.h>
#include <nori/shape.h>

NORI_NAMESPACE_BEGIN

class AreaEmitter: public Emitter {
public:
	AreaEmitter(const PropertyList &props) {
		m_radiance = props.getColor("radiance");
	}

	virtual std::string toString() const override {
		return tfm::format("AreaLight[\n"
				"  radiance = %s,\n"
				"]", m_radiance.toString());
	}

	virtual Color3f eval(const EmitterQueryRecord & lRec) const override {
		if (!m_shape)
			throw NoriException(
					"There is no shape attached to this Area light!");

		float cos_theta_i = lRec.n.dot(-lRec.wi);

		if ((cos_theta_i) >= 0) {
			return m_radiance;
		} else {
			return Color3f(0.f);
		}
	}

	virtual Color3f sample(EmitterQueryRecord & lRec,
			const Point2f & sample) const override {
		if (!m_shape)
			throw NoriException(
					"There is no shape attached to this Area light!");

		ShapeQueryRecord sRec(lRec.ref);
		m_shape->sampleSurface(sRec, sample);
		lRec.p = sRec.p;
		lRec.pdf = sRec.pdf;
		lRec.n = sRec.n;
		lRec.wi = (lRec.p - lRec.ref).normalized();

		Ray3f shadowRay(lRec.ref, lRec.wi, Epsilon,
				(lRec.p - lRec.ref).norm() - Epsilon);
		lRec.shadowRay = shadowRay;
		if (pdf(lRec) > 0)
			return eval(lRec) / pdf(lRec);
		else
			return Color3f(0.f);
	}

	virtual float pdf(const EmitterQueryRecord &lRec) const override {
		if (!m_shape)
			throw NoriException(
					"There is no shape attached to this Area light!");

		ShapeQueryRecord sRec(lRec.ref, lRec.p);
		sRec.pdf = m_shape->pdfSurface(sRec);

		//convert to solid angles
		float cos_theta_i = std::abs(lRec.n.dot(-lRec.wi));

		return sRec.pdf * (lRec.p - lRec.ref).squaredNorm() / cos_theta_i;
//return sRec.pdf
	}

	virtual Color3f samplePhoton(Ray3f &ray, const Point2f &sample1,
			const Point2f &sample2) const override {
		//ray is position and direction
		//sample 1 is for location
		//sample 2 is for direction

		//Sample position
		ShapeQueryRecord sRec(0.f);// = ShapeQueryRecord();
		m_shape->sampleSurface(sRec, sample1);
		Point3f photon_pos = sRec.p;
		Normal3f pole = sRec.n;
		//Sample direction

    	Vector3f q = Warp::squareToCosineHemisphere(sample2);
    	MatrixXf rot = Warp::getRotationMatrix(Vector3f(0.f,0.f,1.f),pole);
    	Vector3f photon_dir=(rot*q).normalized();

    	ray = Ray3f(photon_pos, photon_dir);

    	Color3f result = m_shape->getArea()*m_radiance*M_PI;
    	return result;
	}

protected:
	Color3f m_radiance;
};

NORI_REGISTER_CLASS(AreaEmitter, "area")
NORI_NAMESPACE_END
