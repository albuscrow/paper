/********************************************************
 * 读取、解析obj文件，将其中信息存放到相应的数据结构中
 ********************************************************/

#ifndef __OBJDATA_H__
#define __OBJDATA_H__

#include <string>
#include <vector>

namespace objdata
{
	class Vector3
	{
		double component0_, component1_, component2_;
	public:
		inline Vector3(double component0, double component1, double component2);
		double value0() const { return component0_; }
		double value1() const { return component1_; }
		double value2() const { return component2_; }
		void value0(double value) { component0_ = value; }
		void value1(double value) { component1_ = value; }
		void value2(double value) { component2_ = value; }
		Vector3 &operator+=(const Vector3 &v);
		friend inline const Vector3 operator+(const Vector3 &v0, const Vector3 &v1);
		friend inline const Vector3 operator-(const Vector3 &v0, const Vector3 &v1);
		friend inline const Vector3 operator*(const Vector3 &v, double k);
		friend inline const Vector3 operator*(double k, const Vector3 &v);
		friend inline double operator*(const Vector3 &v0, const Vector3 &v1);
		friend const Vector3 cross(const Vector3 &v0, const Vector3 &v1);
		friend inline std::ostream &operator<<(std::ostream &, const Vector3 &v);
		double inline norm() const;
		void normalize();
	};

	class TextureCoord : public Vector3
	{
	public:
		inline TextureCoord(double u = 0.0, double v = 0.0, double w = 0.0);
		inline TextureCoord(const Vector3 &v);
		double u() const { return value0(); }
		double v() const { return value1(); }
		double w() const { return value2(); }
		void u(double value) { value0(value); }
		void v(double value) { value1(value); }
		void w(double value) { value2(value); }
	};

	class NormalCoord : public Vector3
	{
	public:
		inline NormalCoord(double i = 0.0, double j = 0.0, double k = 0.0);
		inline NormalCoord(const Vector3 &v);
		double i() const { return value0(); }
		double j() const { return value1(); }
		double k() const { return value2(); }
		void i(double value) { value0(value); }
		void j(double value) { value1(value); }
		void k(double value) { value2(value); }
	};

	class VertexCoord : public Vector3
	{
	public:
		inline VertexCoord(double x = 0.0, double y = 0.0, double z = 0.0);
		inline VertexCoord(const Vector3 &v);
		inline VertexCoord(const TextureCoord &t);
		inline VertexCoord(const NormalCoord &n);
		double x() const { return value0(); }
		double y() const { return value1(); }
		double z() const { return value2(); }
		void x(double value) { value0(value); }
		void y(double value) { value1(value); }
		void z(double value) { value2(value); }
	};

	/*
	 * 面片使用纹理、法向的情况
	 */
	enum FaceCase
	{
		V_0_0,			// 只有顶点
		V_0_N,			// 有顶点和法向
		V_T_N			// 有顶点、纹理和法向
	};

	/*
	 * 面片信息
	 */
	struct Face
	{
		std::vector<int> vertexCoordIndex;		// 顶点坐标的索引
		std::vector<int> textureCoordIndex;		// 纹理坐标的索引
		std::vector<int> normalCoordIndex;		// 法向坐标的索引
		FaceCase m_eFaceCase;					// 该面片使用纹理、法向的情况
		int m_nMtlIdx;							// 使用材质的编号
		Face();
		Face(const Face& face);
	};

	/*
	 * 解析obj文件时记录每个mtl的适用范围
	 */
	struct MtlRange
	{
		std::string m_sMtlName;		// mtl的名字
		int m_nBegin, m_nEnd;		// mtl开始、结束的面片编号
		MtlRange(const std::string &mtlName, const int begin);
		MtlRange(const MtlRange &mtlRange);
	};

	/*
	 * 解析mtl文件时记录所有mtl的信息
	 */
	struct MtlTex
	{
		std::string m_sMtlTexName;				// mtl的名字
		double m_fKa[3], m_fKd[3], m_fKs[3];	// 材质属性
		std::string m_sTexFileName;				// 纹理文件名
		unsigned int m_nTexBindingIdx;			// 纹理绑定时的编号
		MtlTex(const std::string mtlTexName);
		MtlTex(const MtlTex &mtlTex);
	};

	/*
	 * 存放从obj文件中读取的数据
	 */
	class ObjData
	{
		void parseVertex(std::istringstream &iss);
		void parseTexture(std::istringstream &iss);
		void parseNormal(std::istringstream &iss);
		void parseFace(std::istringstream &iss, const std::string &str);
		std::string parseMtlLib(std::istringstream &iss) const;
		void parseUseMtl(std::istringstream &iss, std::vector<MtlRange> &mtlRangeList);

		void parseNewMtl(std::istringstream &iss);
		void parseTexFile(std::istringstream &iss, const std::string &filePath);
		void parseKa(std::istringstream &iss);
		void parseKd(std::istringstream &iss);
		void parseKs(std::istringstream &iss);
	public:
		std::vector<VertexCoord> vertexCoordList;		// 存放顶点坐标
		std::vector<TextureCoord> textureCoordList;		// 存放纹理坐标
		std::vector<NormalCoord> normalCoordList;		// 存放法向坐标
		std::vector<Face> faceList;						// 存放面片信息
		std::vector<MtlTex> mtlTexList;

		void readObj(const char *filePath);		// 读取obj文件数据并存放
	};
}

#include "obj_data.inl"

#endif
