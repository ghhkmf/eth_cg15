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

#include <nori/shape.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class Sphere : public Shape {
public:
    Sphere(const PropertyList & propList) {
        m_position = propList.getPoint3("center", Point3f());
        m_radius = propList.getFloat("radius", 1.f);

        m_bbox.expandBy(m_position - Vector3f(m_radius));
        m_bbox.expandBy(m_position + Vector3f(m_radius));
    }

    virtual BoundingBox3f getBoundingBox(uint32_t index) const override { return m_bbox; }

    virtual Point3f getCentroid(uint32_t index) const override { return m_position; }

    virtual bool rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t) const override {


    	//Analytical solution for intersection with Sphere
    	float cod =(m_position-ray.o).dot(ray.d);

    	// ray.d in opposite direction as Sphere AND not IN sphere
    	if(cod<0 && (m_position-ray.o).squaredNorm() > m_radius){
    		return false;
    	}

    	float D_2 = (m_position-ray.o).squaredNorm() - (cod*cod) ;

    	//Compare distances.
    	if (D_2 > m_radius*m_radius){
    		// No  intersection. Ray passes by
    		return false;
    	}

    	float t_cand = cod - sqrt(m_radius*m_radius-D_2);

    	//Set candidate regardless of boundaries
    	t=t_cand;

    	//Return if in boundaries
    	if(t_cand >= ray.mint && t_cand <= ray.maxt){
    		return true;
    	}

        return false;

    }

    virtual void setHitInformation(uint32_t index, const Ray3f &ray, Intersection & its) const override {


    	//Intersection Point
		its.p = ray.o+ray.d*its.t;

		Vector3f relVec = (its.p-m_position).normalized();


		// UV Coordinates SphericalCoordinates
		Point2f s_cord(
				std::atan2(relVec.y(), relVec.x()),
				std::acos(relVec.z())
		);
		its.uv[0]=0.5+(s_cord[0])/(2*M_PI); // -Pi PI - s_cord[0]
		its.uv[1]=s_cord[1]/M_PI; // 0 Pi - s_cord[1]

		// Geometric Frame
		its.geoFrame = Frame(relVec);
		its.shFrame = Frame(relVec);

    }

    virtual void sampleSurface(ShapeQueryRecord & sRec, const Point2f & sample) const override {
        Vector3f q = Warp::squareToUniformSphere(sample);
        sRec.p = m_position + m_radius * q;
        sRec.n = q;
        sRec.pdf = std::pow(1.f/m_radius,2) * Warp::squareToUniformSpherePdf(Vector3f(0));
    }
    virtual float pdfSurface(const ShapeQueryRecord & sRec) const override {
        return std::pow(1.f/m_radius,2) * Warp::squareToUniformSpherePdf(Vector3f(0));
    }


    virtual std::string toString() const override {
        return tfm::format(
                "Sphere[\n"
                "  center = %s,\n"
                "  radius = %f,\n"
                "  bsdf = %s,\n"
                "  emitter = %s\n"
                "]",
                m_position.toString(),
                m_radius,
                m_bsdf ? indent(m_bsdf->toString()) : std::string("null"),
                m_emitter ? indent(m_emitter->toString()) : std::string("null"));
    }

protected:
    Point3f m_position;
    float m_radius;
};

NORI_REGISTER_CLASS(Sphere, "sphere");
NORI_NAMESPACE_END
