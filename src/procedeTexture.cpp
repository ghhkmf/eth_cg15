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

#include <nori/object.h>
#include <nori/texture.h>
#include <fstream> // for file I/O
#include <filesystem/resolver.h>



NORI_NAMESPACE_BEGIN

template <typename T>
class ProcedTexture : public Texture<T> {
public:
	ProcedTexture(const PropertyList &props);

    virtual std::string toString() const override;

    virtual T eval(const Point2f & uv) override {
		Point2f coord = uv;
					
		double noise = ((double)rand() / (RAND_MAX)) ;
		
		double tempx = coord.x()*m_scale.x() + m_shift.x(); 
		double tempy = coord.y()*m_scale.y() + m_shift.y();

		if (abs(tempx - round(tempx))<m_delta.x() && abs(tempy - round(tempy))<m_delta.y())
			return points;

		if (abs(tempx - round(tempx))<m_delta.x() || abs(tempy - round(tempy))<m_delta.y())
			return stripes;		

		return back;
    }

protected:
    T points;
    T stripes;
	T back;

    Point2f m_delta;
	Vector2f m_shift;
	Vector2f m_scale;
	std::string filename;
};

template <>
ProcedTexture<float>::ProcedTexture(const PropertyList &props) {
    m_delta = props.getPoint2("delta", Point2f(0));
	m_shift = props.getVector2("delta", Vector2f(0));
    m_scale = props.getVector2("scale", Vector2f(1));
	points = props.getFloat("points", 0.f);
	stripes = props.getFloat("stripes", 1.f);
	back = props.getFloat("back", 1.f);
	//filename = props.getString("filename");
}

template <>
ProcedTexture<Color3f>::ProcedTexture(const PropertyList &props) {
    m_delta = props.getPoint2("delta", Point2f(0));
	m_shift = props.getVector2("shift", Vector2f(0));
    m_scale = props.getVector2("scale", Vector2f(1));
	points = props.getColor("points", Color3f(0));
	stripes = props.getColor("stripes", Color3f(1));
	back = props.getColor("back", Color3f(1));
	//filename = props.getString("filename");
}


template <>
std::string ProcedTexture<float>::toString() const {
    return tfm::format(
        "ImgTexture[\n"
                "  delta = %s,\n"
                "  scale = %s,\n"
                "  value1 = %f,\n"
                "  value2 = %f,\n"
                "]",
        m_delta.toString(),
        m_scale.toString(),
		points,
		stripes
    );
}

template <>
std::string ProcedTexture<Color3f>::toString() const {
    return tfm::format(
        "ImgTexture[\n"
                "  delta = %s,\n"
                "  scale = %s,\n"
                "  tex1 = %s,\n"
                "  tex2 = %s,\n"
                "]",
        m_delta.toString(),
        m_scale.toString(),
        points.toString(),
		stripes.toString()
    );
}

NORI_REGISTER_TEMPLATED_CLASS(ProcedTexture, float, "procedeTexture_float")
NORI_REGISTER_TEMPLATED_CLASS(ProcedTexture, Color3f, "procedeTexture_color")
NORI_NAMESPACE_END
