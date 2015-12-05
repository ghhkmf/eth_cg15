#include <nori/emitter.h>
#include <nori/bitmap.h>
#include <Eigen/Dense>
#include <math.h>
#include <filesystem/resolver.h>
#include <fstream>
#include <nori/shape.h>
#include <iostream>

NORI_NAMESPACE_BEGIN

/**
 * This class is for SPHERE environmal light.
 * For cube the methods pointToIndex and indexToPint have to be changed.
 *
 * Implementation from Paper:
 * http://www.cs.virginia.edu/~gfx/courses/2007/ImageSynthesis/assignments/envsample.pdf
 */
class EnvMapLight: public Emitter {
public:
	typedef Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> FloatMap;

	EnvMapLight(const PropertyList &propList) {

		/* Load exr file */
		filesystem::path filename = getFileResolver()->resolve(
				propList.getString("filename"));
		std::ifstream is(filename.str());
		if (is.fail())
			throw NoriException("Unable to open OBJ file \"%s\"!", filename);
		m_map = Bitmap(filename.str());

		/* Maps for PDF and CDF*/
		m_conditionalPDF_map = FloatMap(m_map.rows(), m_map.cols());
		m_conditionalCDF_map = FloatMap(m_map.rows(), m_map.cols() + 1);

		m_marginalPDF_array = FloatMap(1, m_map.rows());
		m_marginalCDF_array = FloatMap(1, m_map.rows() + 1);

		// Calculate Marginal and Conditional Maps
		cout << "Calculating marginal and conditional maps" << endl;
		preprocessIntensityMap();

		// Use this to plot the maps. Check that the file
		// "filename2" exists
		//	saveFloatMap(m_conditionalCDF_map);
	}

	/*
	 * Convert the texture into a gray-Intesity map
	 */
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
	 * Prepocessing Part - Calculate the PDF and CDF of the image ------------------------------------
	 */

	void preprocessIntensityMap() {
		int numRows = m_map.rows();
		int numCols = m_map.cols();

		FloatMap m_intensity_map = rgb2gray(m_map);

		FloatMap colsum(1, numRows);
		for (int u = 0; u < numRows; u++) { //Over rows
			colsum(0, u) = precompute1D(u, numCols, m_intensity_map,
					m_conditionalPDF_map, m_conditionalCDF_map);
		}
		precompute1D(0, numRows, colsum, m_marginalPDF_array,
				m_marginalCDF_array);
	}

	float precompute1D(const int rowNum, const int colTotalNum,
			FloatMap &values, FloatMap &pdf, FloatMap &cdf) {
		float I = 0;
		for (int i = 0; i < colTotalNum; i++) {
			I += values(rowNum, i);
		}
		for (int i = 0; i < colTotalNum; i++) {
			pdf(rowNum, i) = values(rowNum, i) / I;
		}
		cdf(0, rowNum) = 0;
		for (int i = 1; i < colTotalNum; i++) {
			cdf(rowNum, i) = cdf(rowNum, i - 1) + pdf(rowNum, i - 1);
		}
		cdf(rowNum, colTotalNum) = 1;
		return I;
	}

	// Sampling Part ------------------------------------------------------------------------------

	Color3f sample(EmitterQueryRecord &lRec, const Point2f &sample) const {
		/*Sampling over the texture. Shape used for mapping from index of texture to 3D Point*/

		if (!m_shape)
			throw NoriException(
					"There is no shape attached to this Area light!");

		int rowSampledIdx = sample1D(0, m_marginalPDF_array,
				m_marginalCDF_array, (float) sample.x());
		int colSampledIdx = sample1D(rowSampledIdx, m_conditionalPDF_map,
				m_conditionalCDF_map, (float) sample.y());

		Point2f idxs(rowSampledIdx, colSampledIdx);
		indexToPoint(idxs, lRec.p);

		lRec.pdf = pdf(lRec.p);
		lRec.n = (lRec.p - m_shape->getCentroid(0)).normalized();
		lRec.wi = (lRec.p - lRec.ref).normalized();
		Ray3f shadowRay(lRec.ref, lRec.wi, Epsilon,
				(lRec.p - lRec.ref).norm() - Epsilon);
		lRec.shadowRay = shadowRay;
		if (lRec.pdf > 0)
			return eval(lRec) / lRec.pdf;
		else
			return Color3f(0.f);
	}

	/**
	 * Returns index of sampled element
	 */

	int sample1D(int rowNum, FloatMap pdf, FloatMap cdf, float sample) const {

		int index = -1;
		//Linear search for index
		for (int var = 1; var < cdf.cols(); ++var) {
			if (cdf(rowNum, var) > sample) {
				index = var - 1;
				break;
			}
		}
		if (index == -1) {
			cout << "ERROR INDEX -1. Setting it to zero" << endl;
			index = 0;
		}
		return index;
	}

	Color3f eval(const EmitterQueryRecord &lRec) const {

		if (!m_shape)
			throw NoriException(
					"There is no shape attached to this Area light!");

		Point2f indexes(0, 0);
		pointToIndex(lRec.p, indexes);


		Color3f result;
		bilinearInterpolation(indexes,result);
		return result;
	}

	/*
	 * Bilinear interpolation of pixels
	 */
	void bilinearInterpolation(Point2f idxs, Color3f& res) const{

		int x_bot = idxs.x();
		int y_bot = idxs.y();
		int x_up = (idxs.x()==(float)x_bot)?x_bot:x_bot+1;
		int y_up = (idxs.y()==(float)y_bot)?y_bot:y_bot+1;

		float x2_x = (((float) x_up)-idxs.x());
		float x_x1 = (idxs.x()-((float) x_bot));


		Color3f R1 =  (x2_x)*m_map(x_bot,y_bot) + ((x_x1))*m_map(x_up,y_bot);
		Color3f R2 = (x2_x)*m_map(x_bot,y_up) + (x_x1)*m_map(x_up,y_up);

		float y_y1=idxs.y()-((float)y_bot);
		float y2_y=((float)y_up)-idxs.y();

		Color3f r((y2_y)*R1 + (y_y1)*R2);

		res[0] =r.x();
		res[1] = r.y();
		res[2]=r.z();
	}



	float pdf(const EmitterQueryRecord &lRec) const {
		/*Uniform sampling over shape - SPHERE*/
		if (!m_shape)
			throw NoriException(
					"There is no shape attached to this Area light!");

		Point2f idxs;
		pointToIndex(lRec.p, idxs);

		return m_marginalPDF_array(0, idxs.x())
				* m_conditionalPDF_map(idxs.x(), idxs.y());
	}




	 // Helper Methods ----------------------------------------------------------------------------------

	std::string toString() const {
		return tfm::format("EnvironmentalMapLight["
				" MapSize = \"%s %s\" \n"
				" ] ", m_map.outerSize(), m_map.innerSize());
	}

	void indexToPoint(Point2f indexes, Point3f &p) const {

		float phi = 2 * M_PI * ((indexes[1] / (m_map.cols() - 1)) - 0.5f);
		float ro = indexes[0] * M_PI / (m_map.rows() - 1);

		BoundingBox3f bb = m_shape->getBoundingBox();
		float r = ((bb.getCenter() - bb.getCorner(1)).norm()) / std::sqrt(2);

		p[0] = std::cos(ro) * std::cos(phi) * r;
		p[1] = std::sin(ro) * std::cos(phi) * r;
		p[2] = std::sin(phi) * r;
	}

	void pointToIndex(Point3f point, Point2f &idx) const {
		//float r = (m_shape->getCentroid(0) - point).norm();
		BoundingBox3f bb = m_shape->getBoundingBox();
		float r = ((bb.getCenter() - bb.getCorner(1)).norm()) / std::sqrt(2);
		float ro = acos(point.z() / r);
		float phi = atan2(point.y(), point.x());

		float col = (0.5f + phi / (2 * M_PI)) * (m_map.cols() - 1); // -PI PI -> 0 cols-1
		float row = (ro / M_PI) * (m_map.rows() - 1); // 0 PI -> 0 rows-1

		idx[0] = row;
		idx[1] = col;
	}

	void saveFloatMap(FloatMap map) {
		int numCols = map.cols();
		int numRows = map.rows();

		Bitmap m(Vector2i(numCols, numRows)); //Invert to match size

		for (int r = 0; r < numRows; r++) {
			for (int c = 0; c < numCols; c++) {
				m(r, c) = Color3f(map(r, c) * 200);
			}
		}

		filesystem::path filename2 = getFileResolver()->resolve("test.exr"); // For controlling

		m.save(filename2.str());
	}

protected:
	//Point3f m_position;
//	Color3f m_power;
	Bitmap m_map;

	FloatMap m_conditionalPDF_map;
	FloatMap m_conditionalCDF_map;

	FloatMap m_marginalPDF_array; //#rows=1 //element at (0,X) is marginal(Row x in Map)
	FloatMap m_marginalCDF_array; // #rows=1

};

NORI_REGISTER_CLASS(EnvMapLight, "env_map");
NORI_NAMESPACE_END
