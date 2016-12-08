#ifndef __MATRIXSTACK_H__
#define __MATRIXSTACK_H__

#include <stack>
#include <iostream>

namespace matrix_stack
{
	class Vector3
	{
		float m_fVector[3];
	public:
		Vector3(float x = 0.0, float y = 0.0, float z = 0.0);
		float &operator[](int i) { return m_fVector[i]; }
		const float &operator[](int i) const { return m_fVector[i]; }
		friend const Vector3 operator-(const Vector3 &v0, const Vector3 &v1);
		friend const Vector3 operator*(const Vector3 &v0, const Vector3 &v1);
		void normalize();
	};

/************************************************************************/

	class Vector4
	{
		float m_fVector[4];
	public:
		Vector4(float val = 0.0);
		float &operator[](int i) { return m_fVector[i]; }
		const float &operator[](int i) const { return m_fVector[i]; }
	};

/************************************************************************/

	class Matrix4x4
	{
		Vector4 m_fMatrix[4];
	public:
		Matrix4x4(float val = 0.0f);
		Matrix4x4(float *mat);
		//Matrix4x4(float val, float val2);
		Vector4 &operator[](int i) { return m_fMatrix[i]; }
		const Vector4 &operator[](int i) const { return m_fMatrix[i]; }
		friend const Matrix4x4 operator*(const Matrix4x4 &m0, const Matrix4x4 &m1);
		//void toPtr3x3(float *m3x3) const;
		void inverseRotateMatrix(float *m3x3) const;
		//friend std::ostream &operator<<(const std::ostream &, const Matrix4x4 &);
	};

/************************************************************************/

	class MatrixStack
	{
	protected:
		static const double PI;
		static const double ZERO;
		std::stack<Matrix4x4> m_matrixStack;
		mutable float m_topMatrixArray[16];
		friend std::ostream& operator<<(std::ostream &, const MatrixStack &);
	public:
		MatrixStack() { m_matrixStack.push(Matrix4x4(1.0f)); }

		void loadMatrix(const Matrix4x4 &matrix);
		void multMatrix(const Matrix4x4 &matrix);
		void loadIdentity() { loadMatrix(Matrix4x4(1.0f)); }

		void pushMatrix() { m_matrixStack.push(m_matrixStack.top()); }
		void popMatrix() { m_matrixStack.pop(); }

		const float *top() const;
		void top(Matrix4x4 &result) const { result = m_matrixStack.top(); }
	};

/************************************************************************/

	class ModelViewMatrixStack : public MatrixStack
	{
	public:
		void translate(float x, float y, float z);
		void rotate(float angle, float x, float y, float z);
		void scale(float x, float y, float z);
		void lookAt(float eyeX, float eyeY, float eyeZ,
					float centerX, float centerY, float centerZ,
					float upX, float upY, float upZ);
		const float *normalMatrix();
	};

/************************************************************************/

	class ProjectionMatrixStack : public MatrixStack
	{
	public:
		void frustum(float left, float right, float bottom, float top,
					 float near, float far);
		void perspective(float fovy, float aspect, float near, float far);
		void ortho(float left, float right, float bottom, float top,
				   float near, float far);
		void ortho2D(float left, float right, float bottom, float top)
		{ ortho(left, right, bottom, top, -1.0f, 1.0f); }
	};
}

#endif
