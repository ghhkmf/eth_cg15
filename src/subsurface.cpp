/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob, Romain Prévost

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
#include <nori/texture.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Diffuse / Lambertian BRDF model
 */
class Subsurface : public NoriObject {
public:
    Subsurface(const PropertyList &propList) : m_albedo(nullptr) {
        if(propList.has("albedo")) {
            PropertyList l;
            l.setColor("value", propList.getColor("albedo"));
            m_albedo = static_cast<Texture<Color3f> *>(NoriObjectFactory::createInstance("constant_color", l));
        }
		//scattering coeff
		sigmaS = propList.getFloat("sigmaS", 1.0f);
		//absorbtion coeff
		sigmaA = propList.getFloat("sigmaA", 1.0f);
		//refraction index
		n = propList.getFloat("n", 1.3f);
		//extinction coeff
		sigmaT = sigmaA + sigmaS;
		//scattering scale
		g = propList.getFloat("g", 0.f);
		//effective scattering coeff
		sigmaSPrime = sigmaS*(1 - g);
		//effective extinction coeff
		sigmaTPrime = sigmaSPrime + sigmaA;
		//effective transport coeff
		sigmaTr = sqrtf(3 * sigmaA*sigmaTPrime);
		//diffusion constant
		D = 1 / (3*sigmaTPrime);
		//diffuse reflection
		Fdr = -1.44 / (n*n) + .71 / n + .668 + 0.0636*n;
		A = (1 + Fdr) / (1 - Fdr);
		//real light source (negative)
		zr = 1 / sigmaTPrime;
		//positive virtual ight
		zv = zr + 4 * A*D;
		alpha = sigmaSPrime / sigmaTPrime;
    }
    virtual ~Subsurface() {
        delete m_albedo;
    }

    /// Add texture for the albedo
    virtual void addChild(NoriObject *obj) override {
        switch (obj->getClassType()) {
            case ETexture:
                if(obj->getIdName() == "albedo") {
                    if (m_albedo)
                        throw NoriException("There is already an albedo defined!");
                    m_albedo = static_cast<Texture<Color3f> *>(obj);
                }
                else {
                    throw NoriException("The name of this texture does not match any field!");
                }
                break;

            default:
                throw NoriException("Diffuse::addChild(<%s>) is not supported!",
                                    classTypeName(obj->getClassType()));
        }
    }

    virtual void activate() override {
        if(!m_albedo) {
            PropertyList l;
            l.setColor("value", Color3f(0.5f));
            m_albedo = static_cast<Texture<Color3f> *>(NoriObjectFactory::createInstance("constant_color", l));
            m_albedo->activate();
        }
    }

    /// Evaluate the BRDF model
    virtual Color3f diffusionApprox(const BSDFQueryRecord &bRec) const  {
        
		//diffuse reflectance due to dipole source
		float Rd;
		float dr; // = ||x-xr||
		float dv; // = || x − xv ||
		float first = (sigmaTr*dr + 1)*exp(-sigmaTr*dr) / (sigmaTPrime*dr*dr*dr);
		float second = zv*(sigmaTr*dv + 1)*exp(-sigmaTr*dv) / (sigmaTPrime*dv*dv*dv);
		 Rd = alpha / (4 * M_PI)*(first + second);
		return Color3f(1.0f);
    }

   
    

    /// Draw a a sample from the BRDF model
    virtual Color3f singleScatter(BSDFQueryRecord &bRec, const Point2f &sample) const  {
        

        return m_albedo->eval(bRec.uv);
    }

    bool isDiffuse() const {
        return false;
    }

	bool isSubSurface() const {
		return true;
	}

    /// Return a human-readable summary
    virtual std::string toString() const override {
        return tfm::format(
            "Diffuse[\n"
            "  albedo = %s\n"
            "]",
            m_albedo ? indent(m_albedo->toString()) : std::string("null")
        );
    }

    virtual EClassType getClassType() const override { return EBSDF; }

private:
    Texture<Color3f> * m_albedo;
	float sigmaS; float sigmaA;
	float n;
	float g;
	float sigmaSPrime;
	float sigmaTPrime;
	float sigmaT;
	float sigmaTr;
	float D;
	float Fdr;
	float A;
	float zr;
	float  zv;
	float alpha;
};

NORI_REGISTER_CLASS(Subsurface, "subsurface");
NORI_NAMESPACE_END
