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

#if !defined(__NORI_MEDIUM_H)
#define __NORI_MEDIUM_H

#include <nori/object.h>

NORI_NAMESPACE_BEGIN

struct Intersection;
/**
 * \brief Data record for conveniently querying and sampling the
 * direct illumination technique implemented by a emitter
 */
struct MediumQueryRecord {
    /// Origin point
    Point3f ref;
    //Point where Phasefunction is sampled
    Point3f p;
    /// Direction sampled
    Vector3f wi;
    // Incoming direction
    Vector3f wo;
    /// Probability of sampling wi
    float pdf;

    /// Create an unitialized query record
    MediumQueryRecord() { }

    /// Create a new query record that can be used to sample a emitter
    MediumQueryRecord(const Point3f &ref) : ref(ref) { }

    /**
     * \brief Create a query record that can be used to query the
     * sampling density after having intersected an area emitter
     */
    MediumQueryRecord(const Point3f &ref, const Point3f &p) :
        ref(ref),p(p) {
        wo = (p - ref).normalized();
    }
};

/**
 * \brief Superclass of all emitters
 */
class Medium : public NoriObject {
public:

    /**
     * \brief Sample the emitter and return the importance weight (i.e. the
     * value of the Emitter divided by the probability density
     * of the sample with respect to solid angles).
     *
     * \param lRec    An emitter query record (only ref is needed)
     * \param sample  A uniformly distributed sample on \f$[0,1]^2\f$
     *
     * \return The emitter value divided by the probability density of the sample.
     *         A zero value means that sampling failed.
     */
    virtual Color3f sample(MediumQueryRecord &lRec, const Point2f &sample) const = 0;

    /**
     * \brief Evaluate the emitter
     *
     * \param lRec
     *     A record with detailed information on the emitter query
     * \return
     *     The emitter value, evaluated for each color channel
     */
    virtual Color3f eval(const MediumQueryRecord &lRec) const = 0;

    /**
     * \brief Compute the probability of sampling \c lRec.p.
     *
     * This method provides access to the probability density that
     * is realized by the \ref sample() method.
     *
     * \param lRec
     *     A record with detailed information on the emitter query
     *
     * \return
     *     A probability/density value
     */
    virtual float pdf(const MediumQueryRecord &lRec) const = 0;

    /**
     * check if Point is in Medium
     */
    virtual bool isInside(Point3f p) const =0;

    /**
     * \brief Virtual destructor
     * */
    virtual ~Medium() {}

    /**
     * \brief Return the type of object (i.e. Mesh/Emitter/etc.)
     * provided by this instance
     * */
    virtual EClassType getClassType() const override { return EMedium; }
    /**
     * \brief Set the shape if the emitter is attached to a shape
     * */
    void setShape(Shape * shape) { m_shape = shape; }

private:
    /**
     * Evaluated the phasefunction. This is defined in the Medium class depending on the type.
     */

    virtual float phaseFunction(Point3f point,Vector3f wi, Vector3f wo) const =0;

protected:
    /// Pointer to the shape if the emitter is attached to a shape
    Shape * m_shape = nullptr;
};


NORI_NAMESPACE_END

#endif /* __NORI_MEDIUM_H */
