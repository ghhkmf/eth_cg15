#include <nori/emitter.h>
#include <nori/bitmap.h>
#include <Eigen/Dense>
#include <math.h>
#include <filesystem/resolver.h>
#include <fstream>
#include <nori/shape.h>
#include <iostream>


NORI_NAMESPACE_BEGIN

class EnvMapLight: public Emitter {
public:
	typedef Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> FloatMap;
	typedef Eigen::Array<float, Eigen::Dynamic, 1> FloatArray;

	EnvMapLight(const PropertyList &propList) {
		/* Load exr file */
		filesystem::path filename = getFileResolver()->resolve(
				propList.getString("filename"));

		std::ifstream is(filename.str());
		if (is.fail())
			throw NoriException("Unable to open OBJ file \"%s\"!", filename);
		m_map = Bitmap(filename.str());
		m_intensity_map = rgb2gray(m_map);

		/* Maps for PDF and CDF*/
		m_conditionalPDF_map = FloatMap(m_map.outerSize(), m_map.innerSize());
		m_conditionalCDF_map = FloatMap(m_map.outerSize(),
				m_map.innerSize() + 1);

		m_marginalPDF_map = FloatArray(m_map.outerSize());
		m_marginalCDF_map = FloatArray(m_map.outerSize() + 1);

	//	preprocessIntensityMap();
	}

	FloatMap rgb2gray(Bitmap map) {
		FloatMap res(map.outerSize(), map.innerSize());
		for (int row = 0; row < map.outerSize(); row++) {
			for (int col = 0; col < map.innerSize(); col++) {
				Color3f c = map(row, col);
				res(row, col) = sqrt(
						(c.r() * c.r() + c.g() * c.g() + c.b() * c.b()) / 3);
			}
		}
		return res;
	}

	/**
	 * Calculated the PDF and CDF of the image
	 */

	void preprocessIntensityMap() {
		int numRows = m_map.outerSize();
		int numCols = m_map.innerSize();

		FloatArray colsum(numRows);

		for (int u = 0; u < numRows; u++) { //Over rows
			precompute1D(u, numCols, m_intensity_map, m_conditionalPDF_map,
					m_conditionalCDF_map);
			for (int j = 0; j < numCols; j++) {
				colsum(u) += m_intensity_map(u, j);
			}
		}
		precompute1DArray(numCols, colsum, m_marginalPDF_map,
				m_marginalCDF_map);
	}

	void precompute1DArray(const int colTotalNum, FloatArray &values,
			FloatArray &pdf, FloatArray &cdf) {

		int I = 0;
		for (int i = 0; i < colTotalNum; i++) {
			I += values(i);
		}

		for (int i = 0; i < colTotalNum; i++) {
			pdf(i) = values(i) / I;
		}
		cdf(0) = 0;
		for (int i = 1; i < colTotalNum; i++) {
			cdf(i) = cdf(i - 1) + pdf(i - 1);
		}
		cdf(colTotalNum) = 1;
	}

	void precompute1D(const int rowNum, const int colTotalNum, FloatMap &values,
			FloatMap &pdf, FloatMap &cdf) {

		int I = 0;
		for (int i = 0; i < colTotalNum; i++) {
			I += values(rowNum, i);
		}

		for (int i = 0; i < colTotalNum; i++) {
			pdf(rowNum, i) = values(rowNum, i) / I;
		}
		cdf(0) = 0;
		for (int i = 1; i < colTotalNum; i++) {
			cdf(rowNum, i) = cdf(rowNum, i - 1) + pdf(rowNum, i - 1);
		}
		cdf(rowNum, colTotalNum) = 1;
	}

	Color3f sample(EmitterQueryRecord &lRec, const Point2f &sample) const {
		/*Uniform sampling over shape - SPHERE*/

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
	Color3f eval(const EmitterQueryRecord &lRec) const {

		if (!m_shape)
			throw NoriException(
					"There is no shape attached to this Area light!");

		//float cos_theta_i = lRec.n.dot(-lRec.wi);

		//lRec.p //Sampled point
		float r=sqrtf(m_shape->getArea()/(4*M_PI));
		float ro = acos(lRec.p.z()/r);
		float phi=atan2(lRec.p.y(),lRec.p.x());

		//Map to texture

		float col = (0.5f + phi / (2 * M_PI))*(m_map.innerSize()-1); // -Pi PI - s_cord[0]
		float row = (ro / M_PI)*(m_map.outerSize()-1); // 0 Pi - s_cord[1]

		return m_map(row,col);
	}
	float pdf(const EmitterQueryRecord &lRec) const {
		/*Uniform sampling over shape - SPHERE*/
		if (!m_shape)
			throw NoriException(
					"There is no shape attached to this Area light!");

		ShapeQueryRecord sRec(lRec.ref, lRec.p);
		sRec.pdf = m_shape->pdfSurface(sRec);

		//convert to solid angles
		float cos_theta_i = std::abs(lRec.n.dot(-lRec.wi));

		return sRec.pdf * (lRec.p - lRec.ref).squaredNorm() / cos_theta_i;
	}

	std::string toString() const {
		return tfm::format("EnvironmentalMapLight["
				" MapSize = \"%s %s\" \n"
				" IntMapSize = \"%s %s\" \n"
				" ] ", m_map.outerSize(), m_map.innerSize(),
				m_intensity_map.outerSize(), m_intensity_map.innerSize());
	}

protected:
	//Point3f m_position;
//	Color3f m_power;
	Bitmap m_map;
	FloatMap m_intensity_map;

	FloatMap m_conditionalPDF_map;
	FloatMap m_conditionalCDF_map;

	FloatArray m_marginalPDF_map;
	FloatArray m_marginalCDF_map;

};

NORI_REGISTER_CLASS(EnvMapLight, "env_map");
NORI_NAMESPACE_END
