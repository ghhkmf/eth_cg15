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

#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/scene.h>
#include <nori/photon.h>
#include <algorithm>

NORI_NAMESPACE_BEGIN

class PhotonMapper: public Integrator {
public:
	/// Photon map data structure
	typedef PointKDTree<Photon> PhotonMap;

	PhotonMapper(const PropertyList &props) {
		/* Lookup parameters */
		m_photonCount = props.getInteger("photonCount", 1000000);
		m_photonRadius = props.getFloat("photonRadius",
				0.0f /* Default: automatic */);
		m_q = 0.2;
	}

	virtual void preprocess(const Scene *scene) override {
		cout << "Gathering " << m_photonCount << " photons .. ";
		cout.flush();

		/* Create a sample generator for the preprocess step */
		Sampler *sampler =
				static_cast<Sampler *>(NoriObjectFactory::createInstance(
						"independent", PropertyList()));

		/* Allocate memory for the photon map */
		m_photonMap = std::unique_ptr<PhotonMap>(new PhotonMap());
		m_photonMap->reserve(m_photonCount);

		/* Estimate a default photon radius */
		if (m_photonRadius == 0)
			m_photonRadius = scene->getBoundingBox().getExtents().norm()
					/ 500.0f;

		/* How to add a photon?
		 * m_photonMap->push_back(Photon(
		 *	Point3f(0, 0, 0),  // Position
		 *	Vector3f(0, 0, 1), // Direction
		 *	Color3f(1, 2, 3)   // Power
		 * ));
		 */
		const std::vector<Emitter *> emis = scene->getEmitters();
		float prob_emi = 1.f / emis.size();

		m_emittedPhotonCount=0;
		for (int emittedPhotonAmount = 0; emittedPhotonAmount < m_photonCount; //TODO: Change this
				emittedPhotonAmount = m_photonMap->size()) {

			//1. Sample Photon from random Lightsource
			//Get random Emitter
			int randomIndex = rand() % emis.size();
			Emitter* emi = emis.at(randomIndex);

			Ray3f photon;
			Color3f photon_power(
					emi->samplePhoton(photon, sampler->next2D(),
							sampler->next2D()) / (prob_emi)); //*mphotoCount
			m_emittedPhotonCount++;
			tracePhoton(scene, photon, photon_power, sampler);
		}


		// Scale Power
		//m_emittedPhotonCount



		/* Build the photon map */
		m_photonMap->build();
	}

	void tracePhoton(const Scene *scene, Ray3f& photon, Color3f &photon_power,
			Sampler* sampler) { //Create loop instead of recursion

		while (true) {

			Intersection its;

			if (!scene->rayIntersect(photon, its)) {
				break;
			}
			//Store if Diffuse
			possiblyStorePhoton(its, photon, photon_power);

			//Sample BSDF
			const BSDF* bsdf = its.mesh->getBSDF();
			BSDFQueryRecord bsdfQuery(its.toLocal(-photon.d)); //wi Photon, wo sampled ray
			bsdfQuery.p = its.p;
			Color3f bsdfVal = bsdf->sample(bsdfQuery, sampler->next2D()); // bsdf*cos/pdf TODO: Check if absDot = cos

			photon = Ray3f(its.p, its.toWorld(bsdfQuery.wo));

			// Calculate new power
			Color3f photon_power = photon_power * bsdfVal;

			//Check russian roulette
			if (survivedRussianRoulette(bsdfVal, sampler->next1D())) { //bsdfVal = new/old
				//continue
				float p = 1 - std::min<float>(1, bsdfVal.maxCoeff());
				photon_power = photon_power / (1 - p);
			} else {
				//Terminate
				break;
			}
		}

	}

	bool survivedRussianRoulette(Color3f& ratio, float rand) {
		float p = 1 - std::min<float>(1, ratio.maxCoeff());

		if (rand > p) {
			return true;
		} else {
			//terminate
			return false;
		}

	}

	void possiblyStorePhoton(Intersection& its, Ray3f& photon, Color3f& power) {
		if (!its.mesh->getBSDF()->isDiffuse())
			return; //Not diffuse. Dont store photon

		m_photonMap->push_back(Photon(its.p,  // Position
				photon.d, // Direction
				power   // Power -> DIVIDE by m_emmitedPhotonCount on LUminance calculation
				));

	}

	virtual Color3f Li(const Scene *scene, Sampler *sampler,
			const Ray3f &_ray) const override {

		/* How to find photons?
		 * std::vector<uint32_t> results;
		 * m_photonMap->search(Point3f(0, 0, 0), // lookup position
		 *                     m_photonRadius,   // search radius
		 *                     results);
		 *
		 * for (uint32_t i : results) {
		 *    const Photon &photon = (*m_photonMap)[i];
		 *    cout << "Found photon!" << endl;
		 *    cout << " Position  : " << photon.getPosition().toString() << endl;
		 *    cout << " Power     : " << photon.getPower().toString() << endl;
		 *    cout << " Direction : " << photon.getDirection().toString() << endl;
		 * }
		 */

		/* Find the surface that is visible in the requested direction */
		Intersection its;

		// If not visible return black
		if (!scene->rayIntersect(_ray, its))
			return Color3f(0.0f);

		//Get Le
		Color3f Le = Color3f(0.0f);
		const Emitter* emi2 = its.mesh->getEmitter();
		if (emi2) {
			//Its an Emitter.
			EmitterQueryRecord iRec2(its.p);
			Le = emi2->eval(iRec2);
		}

		const BSDF* bsdf = its.mesh->getBSDF();
		Vector3f toRayOrigin = -_ray.d.normalized();

		if (bsdf->isDiffuse()) {
			// Query Photomap
			Color3f Li(0.f);
			std::vector < uint32_t > results;
			m_photonMap->search(its.p, // lookup position
					m_photonRadius,   // search radius
					results);
			for (uint32_t i : results) {
				const Photon &photon = (*m_photonMap)[i];

				BSDFQueryRecord query(its.toLocal(toRayOrigin),
						its.toLocal(-photon.getDirection()), ESolidAngle); //wi Camera, wo sampled ray
				query.p = its.p;
				query.uv = its.uv;

				Li = Li
						+ bsdf->eval(query) * photon.getPower()
								/ (M_PI * m_photonRadius * m_photonRadius*m_photonCount);

				//cout << "Found photon!" << endl;
				//cout << " Position  : " << photon.getPosition().toString()
				//		<< endl;
				//cout << " Power     : " << photon.getPower().toString()
				//		<< endl;
				//cout << " Direction : " << photon.getDirection().toString()
				//		<< endl;
			}

			// Terminate return
			return (Le + Li) / (1 - m_q);
		} else {

			BSDFQueryRecord query(its.toLocal(toRayOrigin)); //wi Camera, wo sampled ray
			query.p = its.p;
			Color3f bsdfVal = bsdf->sample(query, sampler->next2D());

			//Check if intersect with emitter
			Ray3f lightRay(its.p, its.toWorld(query.wo));

			if (sampler->next1D() > m_q)
				return (Le + this->Li(scene, sampler, lightRay) * bsdfVal)
						/ (1 - m_q);
			else
				return Le / (1 - m_q);
		}

	}

	virtual std::string toString() const override {
		return tfm::format("PhotonMapper[\n"
				"  photonCount = %i,\n"
				"  photonRadius = %f\n"
				"]", m_photonCount, m_photonRadius);
	}
private:
	int m_photonCount;
	int m_emittedPhotonCount;
	float m_photonRadius;
	std::unique_ptr<PhotonMap> m_photonMap;
	float m_q;
};

NORI_REGISTER_CLASS(PhotonMapper, "photonmapper");
NORI_NAMESPACE_END
