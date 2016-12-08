#include <cmath>
#include <iostream>

namespace objdata
{
	Vector3::Vector3(double component0, double component1, double component2)
	{
		component0_ = component0;
		component1_ = component1;
		component2_ = component2;
	}

	const Vector3 operator+(const Vector3 &v0, const Vector3 &v1)
	{
		return Vector3(v0.component0_ + v1.component0_,
					   v0.component1_ + v1.component1_,
					   v0.component2_ + v1.component2_);
	}

	const Vector3 operator-(const Vector3 &v0, const Vector3 &v1)
	{
		return Vector3(v0.component0_ - v1.component0_,
					   v0.component1_ - v1.component1_,
					   v0.component2_ - v1.component2_);
	}

	double operator*(const Vector3 &v0, const Vector3 &v1)
	{
		return v0.component0_ * v1.component0_ + 
			   v0.component1_ * v1.component1_ + 
			   v0.component2_ * v1.component2_;
	}

	const Vector3 operator*(const Vector3 &v, double k)
	{
		return Vector3(k * v.component0_, k * v.component1_, k * v.component2_);
	}

	const Vector3 operator*(double k, const Vector3 &v)
	{
		return Vector3(k * v.component0_, k * v.component1_, k * v.component2_);
	}

	std::ostream &operator<<(std::ostream &out, const Vector3 &v)
	{
		// 这个函数的修改一定要慎重，因为它涉及到.edit文件中的控制顶点的读写。
		// 如果加上了括号和逗号，在读取edit文件时就一定要考虑到它们

		//return (out << "(" << v.component0_ << ", " << v.component1_ << ", " << v.component2_ << ")");
		return (out << v.component0_ << " " << v.component1_ << " " << v.component2_);
	}

	double inline Vector3::norm() const
	{
		return sqrt(component0_ * component0_ +
					component1_ * component1_ +
					component2_ * component2_);
	}

	/***************************************************************************/

	TextureCoord::TextureCoord(double u, double v, double w)
	:Vector3(u, v, w) {}

	TextureCoord::TextureCoord(const Vector3 &v)
	:Vector3(v) {}

	/***************************************************************************/

	NormalCoord::NormalCoord(double i, double j, double k)
	:Vector3(i, j, k) {}

	NormalCoord::NormalCoord(const Vector3 &v)
	:Vector3(v) {}

	/***************************************************************************/

	VertexCoord::VertexCoord(double x, double y, double z)
	:Vector3(x, y, z) {}

	VertexCoord::VertexCoord(const Vector3 &v)
	:Vector3(v) {}

	VertexCoord::VertexCoord(const TextureCoord &t)
	:Vector3(t.u(), t.v(), t.w()) {}

	VertexCoord::VertexCoord(const NormalCoord &n)
	:Vector3(n.i(), n.j(), n.k()) {}
}
