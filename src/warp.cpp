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

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

Vector3f Warp::sampleUniformHemisphere(Sampler *sampler, const Normal3f &pole) {
    // Naive implementation using rejection sampling
    Vector3f v;
    do {
        v.x() = 1.f - 2.f * sampler->next1D();
        v.y() = 1.f - 2.f * sampler->next1D();
        v.z() = 1.f - 2.f * sampler->next1D();
    } while (v.squaredNorm() > 1.f);

    if (v.dot(pole) < 0.f)
        v = -v;
    v /= v.norm();

    return v;
}

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
	float r = sqrtf(sample.x());
	float theta = 2*M_PI*sample.y();
	return Point2f(r*cos(theta),r*sin(theta));
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
	return (p[0]*p[0]+p[1]*p[1])<1?1/M_PI:0.0f;
}

Vector3f Warp::squareToUniformSphereCap(const Point2f &sample, float cosThetaMax) {
	Point3f cyl = squareToUniformCylinder(sample);
	float t=(cyl.z()+1.f)/2.f;
	float n_z= (cosThetaMax)*t+(1.f-t);

	float r = sqrtf(1.f-n_z*n_z);
	return Vector3f(r*cyl.x(),r*cyl.y(),n_z);
}

float Warp::squareToUniformSphereCapPdf(const Vector3f &v, float cosThetaMax) {
    float area = (1.f-cosThetaMax)*2*M_PI;
	return (v.z()>cosThetaMax)?1.f/area:0.f;
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    Point3f cyl = squareToUniformCylinder(sample);
    float r = sqrtf(1.f-cyl.z()*cyl.z());
    return Vector3f(r*cyl.x(),r*cyl.y(),cyl.z());
}

float Warp::squareToUniformSpherePdf(const Vector3f &p) {
	return 1.f/(4*M_PI);
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
	Point3f cyl = squareToUniformCylinder(sample);
	float sign = 1;
	 float r = sqrtf(1.f-cyl.z()*cyl.z());
	 if(cyl.z()<0)
		 sign=-1;
	 return Vector3f(r*cyl.x(),r*cyl.y(),sign*cyl.z());

}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
	return v.z()>0?1.f/(2*M_PI):0.0f;
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
	float theta = acos(2*sample.x()-1)/2;
	float z=cos(theta);
	float r = sqrtf(1.f-z*z);
	float phi = 2*M_PI*sample.y();
	return Vector3f(r*cos(phi),r*sin(phi),z);
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
	return v.z()>0?v.z()/M_PI:0.f;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
   float phi = 2*M_PI*sample.y();
   float theta = atan(sqrtf(-1*alpha*alpha*log(1-sample.x())));
   float z=cos(theta);
   	float r = sqrtf(1.f-z*z);
   return Vector3f(r*cos(phi),r*sin(phi),z);
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    float theta = acos(m.z());
    float D= 1/(M_PI*alpha*alpha*cos(theta)*cos(theta)*cos(theta)*cos(theta))*exp(-tan(theta)*tan(theta)/(alpha*alpha));
    return m.z()>0?D*cos(theta):0.f;
}
/**
 * Returns w_x,w_y and w_z
 */
Point3f Warp::squareToUniformCylinder(const Point2f &sample) {
	float w_z= 2*sample.x()-1;
	float phi = 2*M_PI*sample.y();
	return Point3f(cos(phi),sin(phi),w_z);
}


NORI_NAMESPACE_END
