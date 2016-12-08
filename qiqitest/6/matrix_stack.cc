#include <cmath>
#include <iostream>
#include <iomanip>
#include "matrix_stack.h"

namespace matrix_stack
{
	Vector3::Vector3(float x, float y, float z)
	{
		m_fVector[0] = x;
		m_fVector[1] = y;
		m_fVector[2] = z;
	}

	const Vector3 operator-(const Vector3 &v0, const Vector3 &v1)
	{
		return Vector3(v0[0] - v1[0], v0[1] - v1[1], v0[2] - v1[2]);
	}

	const Vector3 operator*(const Vector3 &v0, const Vector3 &v1)
	{
		float x = v0[1] * v1[2] - v0[2] * v1[1];
		float y = v0[2] * v1[0] - v0[0] * v1[2];
		float z = v0[0] * v1[1] - v0[1] * v1[0];

		return Vector3(x, y, z);
	}

	void Vector3::normalize()
	{
		float norm = sqrt(m_fVector[0] * m_fVector[0] +
						  m_fVector[1] * m_fVector[1] +
						  m_fVector[2] * m_fVector[2]);
		m_fVector[0] /= norm;
		m_fVector[1] /= norm;
		m_fVector[2] /= norm;
	}

	/************************************************************************/

	Vector4::Vector4(float val)
	{
		m_fVector[0] = m_fVector[1] = m_fVector[2] = m_fVector[3] = val;
	}

	/************************************************************************/

	Matrix4x4::Matrix4x4(float val)
	{
		m_fMatrix[0][0] = m_fMatrix[1][1] = m_fMatrix[2][2] = m_fMatrix[3][3] = val;
	}

	Matrix4x4::Matrix4x4(float *mat)
	{
		for (int i = 0; i < 16; ++i)
		{
			m_fMatrix[i / 4][i % 4] = mat[i];
		}
	}

	//Matrix4x4::Matrix4x4(float val, float val2)
	//{
		//m_fMatrix[0][0] = -0.3494;
		//m_fMatrix[0][1] = -0.5322;
		//m_fMatrix[0][2] = 0.7712;
		//m_fMatrix[0][3] = 0;

		//m_fMatrix[1][0] = -0.4789;
		//m_fMatrix[1][1] = 0.8088;
		//m_fMatrix[1][2] = 0.3413;
		//m_fMatrix[1][3] = 0;

		//m_fMatrix[2][0] = -0.8054;
		//m_fMatrix[2][1] = -0.2501;
		//m_fMatrix[2][2] = -0.5374;
		//m_fMatrix[2][3] = 0;

		//m_fMatrix[3][0] = 0;
		//m_fMatrix[3][1] = 0;
		//m_fMatrix[3][2] = 0;
		//m_fMatrix[3][3] = 1;
	//}

	const Matrix4x4 operator*(const Matrix4x4 &m0, const Matrix4x4 &m1)
	{
		Matrix4x4 result;
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				for (int k = 0; k < 4; ++k)
					result[i][j] += m0[k][j] * m1[i][k];
		return result;
	}

	/*void Matrix4x4::toPtr3x3(float *m3x3) const
	{
		for (int col = 0; col < 3; ++col)
			for (int row = 0; row < 3; ++row)
				m3x3[col * 3 + row] = m_fMatrix[col][row];
	}*/

	void Matrix4x4::inverseRotateMatrix(float *m3x3) const
	{
		float m00 = m_fMatrix[1][1] * m_fMatrix[2][2] - m_fMatrix[2][1] * m_fMatrix[1][2];
		float m01 = m_fMatrix[2][1] * m_fMatrix[0][2] - m_fMatrix[0][1] * m_fMatrix[2][2];
		float m02 = m_fMatrix[0][1] * m_fMatrix[1][2] - m_fMatrix[1][1] * m_fMatrix[0][2];

		float m10 = m_fMatrix[2][0] * m_fMatrix[1][2] - m_fMatrix[1][0] * m_fMatrix[2][2];
		float m11 = m_fMatrix[0][0] * m_fMatrix[2][2] - m_fMatrix[2][0] * m_fMatrix[0][2];
		float m12 = m_fMatrix[1][0] * m_fMatrix[0][2] - m_fMatrix[0][0] * m_fMatrix[1][2];

		float m20 = m_fMatrix[1][0] * m_fMatrix[2][1] - m_fMatrix[2][0] * m_fMatrix[1][1];
		float m21 = m_fMatrix[0][1] * m_fMatrix[2][0] - m_fMatrix[0][0] * m_fMatrix[2][1];
		float m22 = m_fMatrix[0][0] * m_fMatrix[1][1] - m_fMatrix[1][0] * m_fMatrix[0][1];

		m3x3[0] = m00;		m3x3[3] = m10;		m3x3[6] = m20;
		m3x3[1] = m01;		m3x3[4] = m11;		m3x3[7] = m21;
		m3x3[2] = m02;		m3x3[5] = m12;		m3x3[8] = m22;
	}

	/*std::ostream &operator<<(ostream &out, const Matrix4x4 &m)
	{
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				out << m[i][j] << " ";
			}
			out << endl;
		}
		return out;
	}*/

	/************************************************************************/

	const double MatrixStack::PI = 3.1415926535;
	const double MatrixStack::ZERO = 0.0000001;

	void MatrixStack::loadMatrix(const Matrix4x4 &matrix)
	{
		m_matrixStack.pop();
		m_matrixStack.push(matrix);
	}

	void MatrixStack::multMatrix(const Matrix4x4 &matrix)
	{
		Matrix4x4 topMatrix = m_matrixStack.top();
		m_matrixStack.pop();
		m_matrixStack.push(topMatrix * matrix);
	}

	const float *MatrixStack::top() const
	{
		Matrix4x4 topMatrix = m_matrixStack.top();
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				m_topMatrixArray[i * 4 + j] = topMatrix[i][j];
		return m_topMatrixArray;
	}

	std::ostream &operator<<(std::ostream &os, const MatrixStack &matrixStack)
	{
		os << "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" << std::endl;
		os.precision(4);
		MatrixStack tempStack = matrixStack;
		for (unsigned int idx = 0; idx < tempStack.m_matrixStack.size(); tempStack.popMatrix(), ++idx)
		{
			os << idx << ":" << std::endl;
			const float *topMatrix = matrixStack.top();
			for (int i = 0; i < 4; ++i)
			{
				for (int j = 0; j < 4; ++j)
				{
					os << std::setw(8) << std::fixed << topMatrix[j * 4 + i];
				}
				os << std::endl;
			}
			os << std::endl;
		}
		os.precision(6);
		os << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;

		return os;
	}

	/************************************************************************/

	void ModelViewMatrixStack::translate(float x, float y, float z)
	{
		Matrix4x4 translateMatrix(1.0f);
		translateMatrix[3][0] = x;
		translateMatrix[3][1] = y;
		translateMatrix[3][2] = z;

		multMatrix(translateMatrix);
	}

	void ModelViewMatrixStack::rotate(float angle, float xUnnorm, float yUnnorm, float zUnnorm)
	{
		if (sqrt(angle) < ZERO)
			return;

		float a = angle * PI / 180.0f;
		float c = cos(a);
		float s = sin(a);

		float norm = sqrt(xUnnorm * xUnnorm + yUnnorm * yUnnorm + zUnnorm * zUnnorm);
		if (sqrt(norm) < ZERO)
			return;

		float x = xUnnorm / norm;
		float y = yUnnorm / norm;
		float z = zUnnorm / norm;

		float tempX = (1.0f - c) * x;
		float tempY = (1.0f - c) * y;
		float tempZ = (1.0f - c) * z;

		Matrix4x4 rotateMatrix(1.0f);
		rotateMatrix[0][0] = tempX * x + c;
		rotateMatrix[0][1] = tempX * y + s * z;
		rotateMatrix[0][2] = tempX * z - s * y;

		rotateMatrix[1][0] = tempY * x - s * z;
		rotateMatrix[1][1] = tempY * y + c;
		rotateMatrix[1][2] = tempY * z + s * x;

		rotateMatrix[2][0] = tempZ * x + s * y;
		rotateMatrix[2][1] = tempZ * y - s * x;
		rotateMatrix[2][2] = tempZ * z + c;

		multMatrix(rotateMatrix);
	}

	void ModelViewMatrixStack::scale(float x, float y, float z)
	{
		Matrix4x4 scaleMatrix(1.0f);
		scaleMatrix[0][0] = x;
		scaleMatrix[1][1] = y;
		scaleMatrix[2][2] = z;

		multMatrix(scaleMatrix);
	}

	void ModelViewMatrixStack::lookAt(float eyeX, float eyeY, float eyeZ,
									  float centerX, float centerY, float centerZ,
									  float upX, float upY, float upZ)
	{
		Vector3 eye(eyeX, eyeY, eyeZ);
		Vector3 center(centerX, centerY, centerZ);
		Vector3 f(center - eye);
		f.normalize();

		Vector3 up(upX, upY, upZ);
		up.normalize();

		Vector3 s(f * up);
		Vector3 u(s * f);

		Matrix4x4 lookAtMatrix(1.0f);
		for (int i = 0; i < 3; ++i)
		{
			lookAtMatrix[i][0] = s[i];
			lookAtMatrix[i][1] = u[i];
			lookAtMatrix[i][2] = -f[i];
		}

		multMatrix(lookAtMatrix);
		translate(-eyeX, -eyeY, -eyeZ);
	}

	const float *ModelViewMatrixStack::normalMatrix()
	{
		//Matrix4x4 m = m_matrixStack.back();
		Matrix4x4 m = m_matrixStack.top();
		float determinant = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
						  - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
						  + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

		m_topMatrixArray[0 * 3 + 0] =  (m[1][1] * m[2][2] - m[2][1] * m[1][2]) / determinant;
		m_topMatrixArray[0 * 3 + 1] = -(m[1][0] * m[2][2] - m[2][0] * m[1][2]) / determinant;
		m_topMatrixArray[0 * 3 + 2] =  (m[1][0] * m[2][1] - m[2][0] * m[1][1]) / determinant;
		m_topMatrixArray[1 * 3 + 0] = -(m[0][1] * m[2][2] - m[2][1] * m[0][2]) / determinant;
		m_topMatrixArray[1 * 3 + 1] =  (m[0][0] * m[2][2] - m[2][0] * m[0][2]) / determinant;
		m_topMatrixArray[1 * 3 + 2] = -(m[0][0] * m[2][1] - m[2][0] * m[0][1]) / determinant;
		m_topMatrixArray[2 * 3 + 0] =  (m[0][1] * m[1][2] - m[1][1] * m[0][2]) / determinant;
		m_topMatrixArray[2 * 3 + 1] = -(m[0][0] * m[1][2] - m[1][0] * m[0][2]) / determinant;
		m_topMatrixArray[2 * 3 + 2] =  (m[0][0] * m[1][1] - m[1][0] * m[0][1]) / determinant;

		return m_topMatrixArray;
	}

	/************************************************************************/

	void ProjectionMatrixStack::frustum(float left, float right, float bottom, float top,
										float near, float far)
	{
		Matrix4x4 frustumMatrix;
		frustumMatrix[0][0] = (2.0f * near) / (right - left);
		frustumMatrix[1][1] = (2.0f * near) / (top - bottom);
		frustumMatrix[2][0] = (right + left) / (right - left);
		frustumMatrix[2][1] = (top + bottom) / (top - bottom);
		frustumMatrix[2][2] = (near + far) / (near - far);
		frustumMatrix[2][3] = -1.0f;
		frustumMatrix[3][2] = (2.0f * near * far) / (near - far);

		multMatrix(frustumMatrix);
	}

	void ProjectionMatrixStack::perspective(float fovy, float aspect, float near, float far)
	{
		float f = 1.0f / tan(fovy * PI / 360);

		Matrix4x4 perspectiveMatrix;
		perspectiveMatrix[0][0] = f / aspect;
		perspectiveMatrix[1][1] = f;
		perspectiveMatrix[2][2] = (near + far) / (near - far);
		perspectiveMatrix[2][3] = -1.0f;
		perspectiveMatrix[3][2] = (2.0f * near * far) / (near - far);

		multMatrix(perspectiveMatrix);
	}

	void ProjectionMatrixStack::ortho(float left, float right, float bottom, float top,
									  float near, float far)
	{
		Matrix4x4 orthoMatrix(1.0f);
		orthoMatrix[0][0] = 2.0f / (right - left);
		orthoMatrix[1][1] = 2.0f / (top - bottom);
		orthoMatrix[2][2] = 2.0f / (near - far);
		orthoMatrix[3][0] = (left + right) / (left - right);
		orthoMatrix[3][1] = (bottom + top) / (bottom - top);
		orthoMatrix[3][2] = (near + far) / (near - far);

		multMatrix(orthoMatrix);
	}
}
