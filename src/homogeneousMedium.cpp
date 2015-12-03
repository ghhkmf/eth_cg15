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

#include <nori/medium.h>
#include <nori/warp.h>
#include <nori/shape.h>

NORI_NAMESPACE_BEGIN

class HomogeneousMedium: public Medium {
public:
	HomogeneousMedium(const PropertyList &props) {
		m_sigma_a = props.getFloat("sigma_a");
		m_sigma_s = props.getFloat("sigma_s");
		m_sigma_t = m_sigma_a + m_sigma_s;
		m_albedo = m_sigma_s / m_sigma_t;
	}

	std::string toString() const {
		return tfm::format("HomogeneousMedium[\n"
				"  m_sigma_a = %s,\n"
				"  m_sigma_s = %s,\n"
				"]", m_sigma_a, m_sigma_s);
	}

	Color3f eval(const MediumQueryRecord & lRec) const {
		if (!m_shape)
			throw NoriException(
					"There is no shape attached to this Area light!");

		throw NoriException("To be implemented");
	}

	Color3f sample(MediumQueryRecord & lRec, const Point2f & sample) const {
		if (!m_shape)
			throw NoriException(
					"There is no shape attached to this Area light!");
		throw NoriException("To be implemented");
	}

	float pdf(const MediumQueryRecord &lRec) const {
		if (!m_shape)
			throw NoriException(
					"There is no shape attached to this Area light!");
		throw NoriException("To be implemented");
	}

	bool isInside(Point3f p) const {
		if (!m_shape)
			throw NoriException(
					"There is no shape attached to this Area light!");
		return m_shape->getBoundingBox().contains(p, 0);
	}

private:
	float phaseFunction(Point3f point, Vector3f wi, Vector3f wo) const {
		//Independent on wi, and wo. Isotropic medium
		return 1 / (4 * M_PI);
	}

protected:
	/*
	 * Absorption coefficient
	 */
	float m_sigma_a;
	/*
	 * Scattering coefficient
	 */
	float m_sigma_s;
	/*
	 * Extinction coefficient
	 */
	float m_sigma_t;
	float m_albedo;

};

NORI_REGISTER_CLASS(HomogeneousMedium, "homogen");
NORI_NAMESPACE_END
