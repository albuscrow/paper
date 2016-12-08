#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cuda.h>
#include <cuda_gl_interop.h>
#include <cublas_v2.h>
#include "common_data.h"

using namespace std;

int totalMemD = 0;
int permanentMemD = 0;			// 与模型无关的内存使用量（只需要程序退出时释放）
int modelMemD = 0;				// 仅和模型相关的内存使用量（重新载入模型时释放）
int degreeMemD = 0;				// 仅和B样条体次数相关的内存使用量（重新设定B样条体、重新载入模型时释放）
int tessMemD = 0;				// 和细分程度相关的内存使用量（重新设定采样点数、重新设定B样条体、重新载入模型时释放）
int viewMemD = 0;				// 显示函数申请的显存量

//ofstream fout("cuda.txt");

void callCudaThreadSynchronize()
{
	cudaThreadSynchronize();
}

/* B 样条体求值所需的矩阵 */
extern float matrix_b_spline_f[185];
static __device__ float matrix_b_spline_d[185];

/* 根据阶数、控制顶点数、左端节点的编号返回相应的 B 样条矩阵（用于 B 样条体求值） */
template <typename T>
__host__ __device__ T *matrixCase(T *matrix_b_spline, int order, int ctrlPointNum, int leftIdx)
{
	if (order == 1)
		return matrix_b_spline;					// MB1
	else if (order == 2)
		return matrix_b_spline + 1;				// MB2
	else if (order == 3)
	{
		if (ctrlPointNum == 3)
			return matrix_b_spline + 5;			// MB30
		else
		{
			if (leftIdx == 2)
				return matrix_b_spline + 14;	// MB31
			else if (leftIdx == ctrlPointNum - 1)
				return matrix_b_spline + 23;	// MB32
			else
				return matrix_b_spline + 32;	// MB33
		}
	}
	else
	{
		if (ctrlPointNum == 4)
			return matrix_b_spline + 41;		// MB40
		else if (ctrlPointNum == 5)
		{
			if (leftIdx == 3)
				return matrix_b_spline + 57;	// MB41
			else
				return matrix_b_spline + 73;	// MB42
		}
		else if (ctrlPointNum == 6)
		{
			if (leftIdx == 3)
				return matrix_b_spline + 89;	// MB43
			else if (leftIdx == 4)
				return matrix_b_spline + 105;	// MB44
			else
				return matrix_b_spline + 121;	// MB45
		}
		else
		{
			if (leftIdx == 3)
				return matrix_b_spline + 89;	// MB43
			else if (leftIdx == 4)
				return matrix_b_spline + 137;	// MB46
			else if (leftIdx == ctrlPointNum - 2)
				return matrix_b_spline + 153;	// MB47
			else if (leftIdx == ctrlPointNum - 1)
				return matrix_b_spline + 121;	// MB45
			else
				return matrix_b_spline + 169;	// MB48
		}
	}
}

// 便于CPU端调用的一个代理函数
double *matrixCaseHost(double *matrix_b_spline, int order, int ctrlPointNum, int leftIdx)
{
	return matrixCase(matrix_b_spline, order, ctrlPointNum, leftIdx);
}

static __device__ float3 ctrlPointD[15][15][15];	// 原始控制顶点，目前只用于求truth或者FFD结果
static __device__ float knotListD[3 * 20];			// 节点序列

/*
 * 使用矩阵乘法求 B 样条体的值
 * 仅用于 FFD 算法
 */
__device__ float3 BSplineVolumeValueMatrixD(float u, float v, float w,
											int leftUIdx, int leftVIdx, int leftWIdx,
											int orderU, int orderV, int orderW,
											int ctrlPointNumU, int ctrlPointNumV, int ctrlPointNumW)
{
	float3 result;
	float3 tempCtrlPoint1[4];
	float3 tempCtrlPoint2[4][4];

	float *M, temp[4], mul1[4];

	float tempKnot = knotListD[leftUIdx];
	u = (u - tempKnot) / (knotListD[leftUIdx + 1] - tempKnot);
	tempKnot = knotListD[20 + leftVIdx];
	v = (v - tempKnot) / (knotListD[20 + leftVIdx + 1] - tempKnot);
	tempKnot = knotListD[40 + leftWIdx];
	w = (w - tempKnot) / (knotListD[40 + leftWIdx + 1] - tempKnot);

	// 由三维控制顶点算出二维临时控制顶点
	temp[0] = 1.0f;
	temp[1] = w;
	temp[2] = w * w;
	temp[3] = temp[2] * w;

	M = matrixCase(matrix_b_spline_d, orderW, ctrlPointNumW, leftWIdx);

	for (int i = 0; i < orderW; ++i)
	{
		mul1[i] = 0.0f;
		for (int j = 0; j < orderW; ++j)
		{
			mul1[i] += temp[j] * M[j * orderW + i];
		}
	}
	for (int i = 0; i < orderU; ++i)
	{
		for (int j = 0; j < orderV; ++j)
		{
			tempCtrlPoint2[i][j].x = 0.0f;
			tempCtrlPoint2[i][j].y = 0.0f;
			tempCtrlPoint2[i][j].z = 0.0f;
			for (int k = 0; k < orderW; ++k)
			{
				float3 cp = ctrlPointD[leftUIdx - i][leftVIdx - j][leftWIdx - k];
				tempCtrlPoint2[i][j].x += cp.x * mul1[orderW - 1 - k];
				tempCtrlPoint2[i][j].y += cp.y * mul1[orderW - 1 - k];
				tempCtrlPoint2[i][j].z += cp.z * mul1[orderW - 1 - k];
			}
		}
	}

	// 由二维临时控制顶点算出一维临时控制顶点
	temp[1] = v;
	temp[2] = v * v;
	temp[3] = temp[2] * v;

	M = matrixCase(matrix_b_spline_d, orderV, ctrlPointNumV, leftVIdx);

	for (int i = 0; i < orderV; ++i)
	{
		mul1[i] = 0.0;
		for (int j = 0; j < orderV; ++j)
		{
			mul1[i] += temp[j] * M[j * orderV + i];
		}
	}
	for (int i = 0; i < orderU; ++i)
	{
		tempCtrlPoint1[i].x = 0.0f;
		tempCtrlPoint1[i].y = 0.0f;
		tempCtrlPoint1[i].z = 0.0f;
		for (int j = 0; j < orderV; ++j)
		{
			tempCtrlPoint1[i].x += tempCtrlPoint2[i][j].x * mul1[orderV - 1 - j];
			tempCtrlPoint1[i].y += tempCtrlPoint2[i][j].y * mul1[orderV - 1 - j];
			tempCtrlPoint1[i].z += tempCtrlPoint2[i][j].z * mul1[orderV - 1 - j];
		}
	}

	// 由一维临时控制顶点算出结果
	temp[1] = u;
	temp[2] = u * u;
	temp[3] = temp[2] * u;

	M = matrixCase(matrix_b_spline_d, orderU, ctrlPointNumU, leftUIdx);

	for (int i = 0; i < orderU; ++i)
	{
		mul1[i] = 0.0;
		for (int j = 0; j < orderU; ++j)
		{
			mul1[i] += temp[j] * M[j * orderU + i];
		}
	}
	result.x = 0.0f;
	result.y = 0.0f;
	result.z = 0.0f;
	for (int i = 0; i < orderU; ++i)
	{
		result.x += tempCtrlPoint1[i].x * mul1[orderU - 1 - i];
		result.y += tempCtrlPoint1[i].y * mul1[orderU - 1 - i];
		result.z += tempCtrlPoint1[i].z * mul1[orderU - 1 - i];
	}
	return result;
}

/*
 * kernel，计算三个方向参数分别为 u, v, w 的点的 B 样条体值
 * 仅用于 FFD 算法
 */
__global__ void fromParamToCoordOnePoint(float3 *vertexCoordListD, float3 *vertexParamListD,
										 int vertexCount, int orderU, int orderV, int orderW,
										 int ctrlPointNumU, int ctrlPointNumV, int ctrlPointNumW,
										 int knotIntervalCountU, int knotIntervalCountV, int knotIntervalCountW)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx >= vertexCount)
		return;

	float3 tempVertexParam = vertexParamListD[idx];
	float u = tempVertexParam.x;
	float v = tempVertexParam.y;
	float w = tempVertexParam.z;

	// 预先将其值设为最大，将末端点归入最后一段
	int leftUIdx, leftVIdx, leftWIdx;
	leftUIdx = orderU - 1 + knotIntervalCountU - 1;
	leftVIdx = orderV - 1 + knotIntervalCountV - 1;
	leftWIdx = orderW - 1 + knotIntervalCountW - 1;

	// 沿 U 方向查找当前点所在的节点区间
	for (int i = orderU - 1; i <= orderU - 1 + knotIntervalCountU - 1; ++i)
	{
		if (u >= knotListD[i] && u < knotListD[i + 1])
		{
			leftUIdx = i;
			break;
		}
	}
	// 沿 V 方向查找当前点所在的节点区间
	for (int j = orderV - 1; j <= orderV - 1 + knotIntervalCountV - 1; ++j)
	{
		if (v >= knotListD[20 + j] && v < knotListD[20 + j + 1])
		{
			leftVIdx = j;
			break;
		}
	}
	// 沿 W 方向查找当前点所在的节点区间
	for (int k = orderW - 1; k <= orderW - 1 + knotIntervalCountW - 1; ++k)
	{
		if (w >= knotListD[40 + k] && w < knotListD[40 + k + 1])
		{
			leftWIdx = k;
			break;
		}
	}
	vertexCoordListD[idx] = BSplineVolumeValueMatrixD(u, v, w, leftUIdx, leftVIdx, leftWIdx,
													  orderU, orderV, orderW,
													  ctrlPointNumU, ctrlPointNumV, ctrlPointNumW);
}

float3 *vertexParamListD = 0;					// 模型顶点参数序列
float3 *vertexCoordListD = 0;					// 模型顶点坐标序列

//float3 *vertexParamListD_teapot = 0;					// 模型顶点参数序列
//float3 *normalParamListD_teapot = 0;					// 模型顶点参数序列
//float3 *vertexCoordListD_teapot = 0;					// 模型顶点参数序列
int vertexCount_teapot;

int order[3], ctrlPointNum[3], knotIntervalCount[3], knotCount[3];		// 三个方向的阶数、控制顶点数、节点区间数、节点数
float knotList[3][20];														// 三个方向的节点向量
float3 ctrlPoint[15][15][15];												// B样条体的控制顶点

/*
 * 根据所有顶点的参数，计算出相应的 B 样条体值
 * 仅用于 FFD 算法
 */
void fromParamToCoordD(CommonData *commonData)
{
	int vertexCount = commonData->vertexCount();
	int threadCount = commonData->ffdThreadCount();
	fromParamToCoordOnePoint<<<vertexCount / threadCount + 1, threadCount>>>(
													vertexCoordListD, vertexParamListD,
													vertexCount, order[U], order[V], order[W],
													ctrlPointNum[U], ctrlPointNum[V], ctrlPointNum[W],
													knotIntervalCount[U], knotIntervalCount[V], knotIntervalCount[W]);
	float3 *vertexCoordList = new float3[vertexCount];
	cudaMemcpy(vertexCoordList, vertexCoordListD, sizeof(float3) * vertexCount, cudaMemcpyDeviceToHost);
	for (int i = 0; i < vertexCount; ++i)
		commonData->setVertexCoord(i, vertexCoordList[i].x, vertexCoordList[i].y, vertexCoordList[i].z);
	delete []vertexCoordList;
}

/*------------------------------------------------------- 上面是FFD算法部分  ---------------------------------------------------------*/
/*------------------------------------------------------- 下面是AFFD算法部分 ---------------------------------------------------------*/

/* 把数字a转换成一个逗号分节的string */
string longNumber(int a)
{
	string result;
	do
	{
		ostringstream oss;

		int remainder = a % 1000;
		if (a >= 1000)
		{
			if (remainder < 10)
				oss << "00" << remainder;
			else if (remainder >= 10 && remainder < 100)
				oss << "0" << remainder;
			else 
				oss << remainder;
		}
		else
			oss << remainder;

		if (result.size() == 0)
			result = oss.str();
		else
			result = oss.str() + "," + result;

		a /= 1000;
	}while(a > 0);

	return result;
}

/* 打印显存使用量 */
void printMemD(const char *file, const char *function, int line, int memSize, string info)
{
	/* 只取文件名部分，路径舍弃 */
	string fileName(file);
	int lastSlashPos = fileName.rfind('/');
	fileName = fileName.substr(lastSlashPos + 1, fileName.size());

/*#define PRINT_MEM*/
#ifdef PRINT_MEM
	/*作废totalMemD += memSize;*/
	cout << info << "\n"
		 << "\t文件" << fileName << "，函数" << function << ", 第" << line << "行，申请显存" << longNumber(memSize) << "字节, "
		 << "目前累计使用显存" << longNumber(permanentMemD + modelMemD + degreeMemD + tessMemD + viewMemD) << "字节\n"
		 << "\t其中permanent = " << longNumber(permanentMemD) << ", model = " << longNumber(modelMemD)
		 << ", degreeMemD = " << longNumber(degreeMemD) << ", tessMemD = " << longNumber(tessMemD)
		 << ", view = " << longNumber(viewMemD) << endl;
#endif
}

void printCudaError(const char *file, const char *function, int line)
{
	/* 只取文件名部分，路径舍弃 */
	string fileName(file);
	int lastSlashPos = fileName.rfind('/');
	fileName = fileName.substr(lastSlashPos + 1, fileName.size());

	cudaError_t cymError = cudaGetLastError();
	if (cymError)
		cout << fileName << "第" << line << "行, 错误代码" << cymError << ": " << cudaGetErrorString(cymError) << endl;
}

__host__ __device__ inline const float3 operator+(const float3 &a, const float3 &b)
{
	return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}

__host__ __device__ inline const float3 operator-(const float3 &a, const float3 &b)
{
	return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

__host__ __device__ inline const float3 operator-(const float3 &a)
{
	return make_float3(-a.x, -a.y, -a.z);
}

__host__ __device__ inline const float3 operator*(float a, const float3 &b)
{
	return make_float3(a * b.x, a * b.y, a * b.z);
}

__host__ __device__ inline const float3 operator*(const float3 &a, float b)
{
	return make_float3(a.x * b, a.y * b, a.z * b);
}

__host__ __device__ inline const float3 operator/(const float3 &a, float b)
{
	return make_float3(a.x / b, a.y / b, a.z / b);
}

__host__ __device__ inline float operator*(const float3 &a, const float3 &b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

__device__ float3 cross(const float3 &a, const float3 &b)
{
	return make_float3(a.y * b.z - a.z * b.y,
					   a.z * b.x - a.x * b.z,
					   a.x * b.y - a.y * b.x);
}

__host__ __device__ inline void operator*=(float3 &a, float b)
{
	a.x *= b;
	a.y *= b;
	a.z *= b;
}

__host__ __device__ inline void operator/=(float3 &a, float b)
{
	a.x /= b;
	a.y /= b;
	a.z /= b;
}

__host__ __device__ inline void operator+=(float3 &a, const float3 &b)
{
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
}

__host__ __device__ inline void operator-=(float3 &a, const float3 &b)
{
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
}

__device__ inline float length(const float3 &v)
{
	return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

__device__ inline void normalize(float3 &v)
{
	float length_inverse = 1.0 / length(v);
	v *= length_inverse;
}

cublasHandle_t cublas_handle = 0;

/* 将 B 样条矩阵载入显存 */
void loadMatrixBSplineD()
{
	cudaMemcpyToSymbol(matrix_b_spline_d, matrix_b_spline_f, sizeof(float) * 185);

	cublasCreate(&cublas_handle);
}

static __device__ float3 newCtrlPointD[15][15][15][4][4][4];	// 使用皮本上的优化算法之后新生成的控制顶点，每个节点盒都有一个4x4x4的控制顶点

/* 
 * 计算采样点值的优化算法，事先对每个节点盒分别计算B样条体的控制顶点乘以Mu, Mv, Mw的结果并存储
 * 本函数就是进行这个计算，算法具体思路可以看皮本
 */
__global__ void calcNewCtrlPointD(int order_u, int order_v, int order_w,
								  int ctrlPointNum_u, int ctrlPointNum_v, int ctrlPointNum_w)
{
	int ii = blockIdx.x;
	int jj = blockIdx.y;
	int kk = blockIdx.z;

	int leftUIdx = ii + order_u - 1;
	int leftVIdx = jj + order_v - 1;
	int leftWIdx = kk + order_w - 1;
	float *Mu = matrixCase(matrix_b_spline_d, order_u, ctrlPointNum_u, leftUIdx);
	float *Mv = matrixCase(matrix_b_spline_d, order_v, ctrlPointNum_v, leftVIdx);
	float *Mw = matrixCase(matrix_b_spline_d, order_w, ctrlPointNum_w, leftWIdx);

	// 第一个矩阵乘法
	int base_i = leftUIdx - order_u + 1;
	int base_j = leftVIdx - order_v + 1;
	int base_k = leftWIdx - order_w + 1;

	for (int k = 0; k < order_w; ++k)
		for (int i = 0; i < order_u; ++i)
			for (int j = 0; j < order_v; ++j)
			{
				newCtrlPointD[ii][jj][kk][i][j][k] = make_float3(0.0, 0.0, 0.0);
				for (int l = 0; l < order_u; ++l)
				{
					float3 cp = ctrlPointD[base_i + l][base_j + j][base_k + k];
					newCtrlPointD[ii][jj][kk][i][j][k] += Mu[i * order_u + l] * cp;
				}
			}

	// 第二个矩阵乘法
	float3 box[4][4][4];
	for (int i = 0; i < order_u; ++i)
		for (int j = 0; j < order_v; ++j)
			for (int k = 0; k < order_w; ++k)
			{
				box[i][j][k] = make_float3(0.0, 0.0, 0.0);
				for (int l = 0; l < order_v; ++l)
				{
					float3 cp = newCtrlPointD[ii][jj][kk][i][l][k];
					box[i][j][k] += Mv[j * order_v + l] * cp;
				}
			}

	// 第三个矩阵乘法
	for (int j = 0; j < order_v; ++j)
		for (int k = 0; k < order_w; ++k)
			for (int i = 0; i < order_u; ++i)
			{
				newCtrlPointD[ii][jj][kk][i][j][k] = make_float3(0.0, 0.0, 0.0);
				for (int l = 0; l < order_w; ++l)
				{
					float3 cp = box[i][j][l];
							newCtrlPointD[ii][jj][kk][i][j][k] += Mw[k * order_w + l] * cp;
						}
					}
}

/*
 * 将B样条体控制顶点拷贝到显存
 * 另外，将皮本上新算法中的控制顶点拷贝到显存
 */
void copyCtrlPointD(CommonData *commonData)
{
	for (int i = 0; i < ctrlPointNum[U]; ++i)
	{
		for (int j = 0; j < ctrlPointNum[V]; ++j)
		{
			for (int k = 0; k < ctrlPointNum[W]; ++k)
			{
				ctrlPoint[i][j][k].x = (float)commonData->getCtrlPoint(i, j, k).x();
				ctrlPoint[i][j][k].y = (float)commonData->getCtrlPoint(i, j, k).y();
				ctrlPoint[i][j][k].z = (float)commonData->getCtrlPoint(i, j, k).z();
			}
		}
	}
	cudaMemcpyToSymbol(ctrlPointD, &ctrlPoint[0][0][0], sizeof(float3) * 15 * 15 * 15);
	printCudaError(__FILE__, __FUNCTION__, __LINE__);

	calcNewCtrlPointD<<<dim3(knotIntervalCount[U], knotIntervalCount[V], knotIntervalCount[W]), 1>>>
		(order[U], order[V], order[W], ctrlPointNum[U], ctrlPointNum[V], ctrlPointNum[W]);
	printCudaError(__FILE__, __FUNCTION__, __LINE__);
}

extern float teapot_ctrl_point[3 * 4 * 4 * 36];
static __device__ float3 teapot_ctrl_pointD[4 * 4 * 36];

int u_seg, v_seg;

/* 预计算，将内存中的数据拷贝到相应的显存空间中 */
void preCalcD(CommonData *commonData)
{
	u_seg = commonData->u_seg();
	v_seg = commonData->v_seg();

	float3 temp[16 * 36];
	for (int i = 0; i < 16 * 36; ++i)
	{
		temp[i].x = teapot_ctrl_point[i * 3];
		temp[i].y = teapot_ctrl_point[i * 3 + 1];
		temp[i].z = teapot_ctrl_point[i * 3 + 2];
	}
	cudaMemcpyToSymbol(teapot_ctrl_pointD, temp, sizeof(float3) * 16 * 36);

	for (int i = 0; i < 3; ++i)
	{
		order[i] = commonData->order(i);
		ctrlPointNum[i] = commonData->ctrlPointCount(i);
		knotIntervalCount[i] = commonData->knotIntervalCount(i);
		knotCount[i] = order[i] + ctrlPointNum[i];
	}
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < knotCount[i]; ++j)
			knotList[i][j] = (float)commonData->getKnot(i, j);
	cudaMemcpyToSymbol(knotListD, &knotList[0][0], sizeof(float) * 3 * 20);
	printCudaError(__FILE__, __FUNCTION__, __LINE__);

	int vertexCount = commonData->vertexCount();
	float3 *vertexParamListAlloc = new float3[vertexCount];
	for (int i = 0; i < vertexCount; ++i)
	{
		vertexParamListAlloc[i].x = (float)commonData->vertexParam(i).u();
		vertexParamListAlloc[i].y = (float)commonData->vertexParam(i).v();
		vertexParamListAlloc[i].z = (float)commonData->vertexParam(i).w();
	}
	modelMemD += sizeof(float3) * vertexCount;
	printMemD(__FILE__, __FUNCTION__, __LINE__, sizeof(float3) * vertexCount, "@原始模型上所有顶点的参数，仅用于FFD");

	cudaMalloc((void**)&vertexParamListD, sizeof(float3) * vertexCount);
	printCudaError(__FILE__, __FUNCTION__, __LINE__);

	cudaMemcpy(vertexParamListD, vertexParamListAlloc, sizeof(float3) * vertexCount, cudaMemcpyHostToDevice);
	printCudaError(__FILE__, __FUNCTION__, __LINE__);

	delete []vertexParamListAlloc;
	vertexParamListAlloc = 0;

	modelMemD += sizeof(float3) * vertexCount;
	printMemD(__FILE__, __FUNCTION__, __LINE__, sizeof(float3) * vertexCount, "@原始模型上所有顶点的坐标，仅用于FFD");
	cudaMalloc((void**)&vertexCoordListD, sizeof(float3) * vertexCount);
	printCudaError(__FILE__, __FUNCTION__, __LINE__);

	copyCtrlPointD(commonData);
	printCudaError(__FILE__, __FUNCTION__, __LINE__);

	// teapot的顶点数据拷贝到显存
	//vertexCount_teapot = commonData->vertexCount_teapot();
	//vertexParamListAlloc = new float3[vertexCount_teapot];
	//for (int i = 0; i < vertexCount_teapot; ++i)
	//{
		//vertexParamListAlloc[i].x = (float)commonData->vertexParam_teapot(i).x();
		//vertexParamListAlloc[i].y = (float)commonData->vertexParam_teapot(i).y();
		//vertexParamListAlloc[i].z = (float)commonData->vertexParam_teapot(i).z();
	//}
	//modelMemD += sizeof(float3) * vertexCount_teapot;
	//cudaMalloc((void**)&vertexParamListD_teapot, sizeof(float3) * vertexCount_teapot);
	//printCudaError(__FILE__, __FUNCTION__, __LINE__);
	//cudaMemcpy(vertexParamListD_teapot, vertexParamListAlloc, sizeof(float3) * vertexCount_teapot, cudaMemcpyHostToDevice);

	//for (int i = 0; i < vertexCount_teapot; ++i)
	//{
		//vertexParamListAlloc[i].x = (float)commonData->normalParam_teapot(i).i();
		//vertexParamListAlloc[i].y = (float)commonData->normalParam_teapot(i).j();
		//vertexParamListAlloc[i].z = (float)commonData->normalParam_teapot(i).k();
	//}
	//cudaMalloc((void**)&normalParamListD_teapot, sizeof(float3) * vertexCount_teapot);
	//printCudaError(__FILE__, __FUNCTION__, __LINE__);
	//cudaMemcpy(normalParamListD_teapot, vertexParamListAlloc, sizeof(float3) * vertexCount_teapot, cudaMemcpyHostToDevice);

	//delete []vertexParamListAlloc;
}

int *matrixFittingIdxD;
float *matrixFittingD;

void loadTriangleMatrixD()
{
	extern int matrixFittingIdx[100];
	cudaMalloc((void**)&matrixFittingIdxD, sizeof(int) * 100);
	permanentMemD += sizeof(int) * 100;
	printMemD(__FILE__, __FUNCTION__, __LINE__, sizeof(int) * 100, "@拟合矩阵的索引矩阵");
	cudaMemcpy(matrixFittingIdxD, matrixFittingIdx, sizeof(int) * 100, cudaMemcpyHostToDevice);
	printCudaError(__FILE__, __FUNCTION__, __LINE__);

	extern float matrixFitting[39417];
	cudaMalloc((void**)&matrixFittingD, sizeof(float) * 39417);
	permanentMemD += sizeof(float) * 39417;
	printMemD(__FILE__, __FUNCTION__, __LINE__, sizeof(float) * 39417, "@拟合矩阵");
	cudaMemcpy(matrixFittingD, matrixFitting, sizeof(float) * 39417, cudaMemcpyHostToDevice);
	printCudaError(__FILE__, __FUNCTION__, __LINE__);
}

struct TriangleD
{
	float3 v[3], n[3], n_adj_origin[3], n_adj[3];
#ifdef LINE
	float3 bary_origin[3];
#endif
	int nc[3];		// nc0, nc1, nc2分别代表v2v0, v0v1, v1v2边的法向数量
	float3 vt[3];
};

TriangleD *triangleListD;
float *sampleValueD, *triangleCtrlPointD;
float3 *sampleValueD_PN;
float *triangleCtrlPointD_PN, *triangleNormalCtrlPointD_PN;
int *triangle_adjacent_tableD;
int degree, degree_lower, triangleCtrlPointNum, triangleCtrlPointNum_lower, triangleNum, constrait_point_num;

int blockSizeStep0 = 128, activeThreadNumStep0, blockNumStep0;
int blockSizeStep1 = 128, activeThreadNumStep1, blockNumStep1;
int blockSizeAdjNormal = 128, activeThreadNumAdjNormal, blockNumAdjNormal;
int blockSizeStep0_PN = 128, blockNumStep0_PN;
#ifdef TRUTH
float *B_1D_truth, *sampleValueD_truth;
int activeThreadNumStep0_truth, blockNumStep0_truth;
#endif

int matrixStartIdxFitting;

__host__ __device__ inline int index2c(int i, int j, int stride)
{
	return j * stride + i;
}

#ifdef DRAW_TRIANGULAR_CTRL_POINTS
float *triangular_ctrl_points;
#endif

void loadTriangleListD(const vector<Triangle> &triangleList, int *triangle_adjacent_table, int deg)
{
	triangleNum = triangleList.size();
	degree = deg;
	/*degree_lower = deg;*/
	degree_lower = 3;
	triangleCtrlPointNum = (degree + 1) * (degree + 2) / 2;
	triangleCtrlPointNum_lower = (degree_lower + 1) * (degree_lower + 2) / 2;
	constrait_point_num = 3 * degree_lower;

	TriangleD *tempTriangleList = new TriangleD[triangleNum];
	for (vector<Triangle>::size_type i = 0; i < triangleNum; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			tempTriangleList[i].v[j].x = triangleList[i].v[j].x();
			tempTriangleList[i].v[j].y = triangleList[i].v[j].y();
			tempTriangleList[i].v[j].z = triangleList[i].v[j].z();

			tempTriangleList[i].n[j].x = triangleList[i].n[j].i();
			tempTriangleList[i].n[j].y = triangleList[i].n[j].j();
			tempTriangleList[i].n[j].z = triangleList[i].n[j].k();

			tempTriangleList[i].n_adj_origin[j].x = triangleList[i].n_adj[j].i();
			tempTriangleList[i].n_adj_origin[j].y = triangleList[i].n_adj[j].j();
			tempTriangleList[i].n_adj_origin[j].z = triangleList[i].n_adj[j].k();

			tempTriangleList[i].n_adj[j].x = triangleList[i].n_adj[j].i();
			tempTriangleList[i].n_adj[j].y = triangleList[i].n_adj[j].j();
			tempTriangleList[i].n_adj[j].z = triangleList[i].n_adj[j].k();

#ifdef LINE
			tempTriangleList[i].bary_origin[j].x = triangleList[i].bary_origin[j].x();
			tempTriangleList[i].bary_origin[j].y = triangleList[i].bary_origin[j].y();
			tempTriangleList[i].bary_origin[j].z = triangleList[i].bary_origin[j].z();
#endif

			tempTriangleList[i].nc[j] = triangleList[i].n_count[j];

			tempTriangleList[i].vt[j].x = triangleList[i].vt[j].u();
			tempTriangleList[i].vt[j].y = triangleList[i].vt[j].v();
			tempTriangleList[i].vt[j].z = triangleList[i].vt[j].w();
		}
	}
	cudaMalloc((void**)&triangleListD, sizeof(TriangleD) * triangleNum);
	degreeMemD += sizeof(TriangleD) * triangleNum;
	printMemD(__FILE__, __FUNCTION__, __LINE__, sizeof(TriangleD) * triangleNum, "@原始模型上所有三角形信息");

	cudaMemcpy(triangleListD, tempTriangleList, sizeof(TriangleD) * triangleNum, cudaMemcpyHostToDevice);

	delete []tempTriangleList;

	cudaMalloc(&sampleValueD, sizeof(float) * (triangleCtrlPointNum + constrait_point_num) * triangleNum * 6);
	degreeMemD += sizeof(float) * (triangleCtrlPointNum + constrait_point_num) * triangleNum * 6;
	printMemD(__FILE__, __FUNCTION__, __LINE__, sizeof(float) * (triangleCtrlPointNum + constrait_point_num) * triangleNum * 6,
			  "@为了求Bezier曲面片的控制顶点，需要在其上进行采样，结果放在这里。即第二个矩阵乘法用到的矩阵T");

	cudaMalloc(&sampleValueD_PN, sizeof(float3) * triangleNum * 3 * 2);
	cudaMalloc(&triangleCtrlPointD_PN, sizeof(float) * (1 + 2 + 3 + 4) * triangleNum * 3);
	cudaMalloc(&triangleNormalCtrlPointD_PN, sizeof(float) * (1 + 2 + 3) * triangleNum * 3);

	cudaMalloc(&triangleCtrlPointD, sizeof(float) * triangleCtrlPointNum_lower * triangleNum * 6);

	cudaMalloc(&triangle_adjacent_tableD, sizeof(int) * triangleNum * 3);
	cudaMemcpy(triangle_adjacent_tableD, triangle_adjacent_table, sizeof(int) * triangleNum * 3, cudaMemcpyHostToDevice);

#ifdef TRUTH
	cudaMalloc(&sampleValueD_truth, sizeof(float) * triangleCtrlPointNum * triangleNum * 3);
	degreeMemD += sizeof(float) * triangleCtrlPointNum * triangleNum * 3;
	printMemD(__FILE__, __FUNCTION__, __LINE__, sizeof(float) * triangleCtrlPointNum * triangleNum * 3,
			  "@为了求精确Bezier曲面片的控制顶点，需要在其上进行采样，结果放在这里。即第二个矩阵乘法用到的矩阵T");

	activeThreadNumStep0_truth = triangleCtrlPointNum * triangleNum;
	blockNumStep0_truth = ceil(static_cast<double>(activeThreadNumStep0_truth) / blockSizeStep0);
#endif

	activeThreadNumStep0 = triangleCtrlPointNum * triangleNum;
	blockNumStep0 = ceil(static_cast<double>(activeThreadNumStep0) / blockSizeStep0);

	activeThreadNumStep1 = constrait_point_num * triangleNum;
	blockNumStep1 = ceil(static_cast<double>(activeThreadNumStep1) / blockSizeStep1);

	activeThreadNumAdjNormal = triangleNum * 3;
	blockNumAdjNormal = ceil(static_cast<double>(activeThreadNumAdjNormal) / blockSizeAdjNormal);

	blockNumStep0_PN = ceil(static_cast<double>(3 * triangleNum) / blockSizeStep0_PN);

#ifdef TRUTH
	extern float matrixTriangle[9][55*55];
	float *temp = new float[triangleCtrlPointNum * triangleCtrlPointNum];
	for (int i = 0; i < triangleCtrlPointNum; ++i)
	{
		for (int j = 0; j < triangleCtrlPointNum; ++j)
		{
			temp[index2c(i, j, triangleCtrlPointNum)] = matrixTriangle[degree - 1][i * triangleCtrlPointNum + j];
		}
	}
	cudaMalloc(&B_1D_truth, sizeof(float) * triangleCtrlPointNum * triangleCtrlPointNum);
	degreeMemD += sizeof(float) * triangleCtrlPointNum * triangleCtrlPointNum;
	printMemD(__FILE__, __FUNCTION__, __LINE__, sizeof(float) * triangleCtrlPointNum * triangleCtrlPointNum, "@第一个矩阵乘法用到的矩阵(B-1)T存放在这里");
	cudaMemcpy(B_1D_truth, temp, sizeof(float) * triangleCtrlPointNum * triangleCtrlPointNum, cudaMemcpyHostToDevice);
	delete temp;
#endif
	/***************************************************************************/

	extern int matrixFittingIdx[100];
	matrixStartIdxFitting = matrixFittingIdx[degree * 10 + degree_lower];

	cout << "triangleNum = " << triangleNum << endl;
	cout << "degree = " << degree << ", degree_lower = " << degree_lower << ", constrait_point_num = " << constrait_point_num << endl;
	cout << "triangleCtrlPointNum = " << triangleCtrlPointNum << ", triangleCtrlPointNum_lower = " << triangleCtrlPointNum_lower << endl;
	cout << "activeThreadNumStep1 = " << activeThreadNumStep1 << ", blockNumStep1 = " << blockNumStep1 << endl;

#ifdef DRAW_TRIANGULAR_CTRL_POINTS
	triangular_ctrl_points = new float[3 * triangleCtrlPointNum_lower * triangleNum];
#endif
}

double power(double a, int n)
{
	if (n <= 0)
		return 1.0;
	double result = a;
	for (int i = 1; i < n; ++i)
		result *= a;
	return result;
}

int factorial(int n)
{
	int result = 1;
	for (int i = 1; i <= n; ++i)
		result *= i;
	return result;
}

float B(double u, double v, double w, int n, int3 c)
{
	return factorial(n) / factorial(c.x) / factorial(c.y) / factorial(c.z) * power(u, c.x) * power(v, c.y) * power(w, c.z);
}

float *BqD, *BqD_PN, *BBD, *RD;
int *my_to_truth_tableD;
float3 *parameter3D, *parameterND;
#ifdef TRUTH
float *BqD_truth, *BBD_truth, *RD_truth;
#endif
int segmentPerEdge, samplePointPerTriangle;
int blockSizeCopy = 256, activeThreadNumCopy, blockNumCopy;
int *my_to_truth_table;

void generateUVW(int samplePointPerEdge)
{
	segmentPerEdge = samplePointPerEdge - 1;
	samplePointPerTriangle = (samplePointPerEdge + 1) * samplePointPerEdge / 2;

	activeThreadNumCopy = samplePointPerTriangle * triangleNum;
	blockNumCopy = ceil(static_cast<double>(activeThreadNumCopy) / blockSizeCopy);

	double *a = new double[samplePointPerTriangle * 3];
	int idx = 0;
	for (int i = segmentPerEdge; i >= 0; --i)
	{
		for (int j = segmentPerEdge - i; j >= 0; --j)
		{
			int k = segmentPerEdge - i - j;
			a[idx++] = (double)i / segmentPerEdge;
			a[idx++] = (double)j / segmentPerEdge;
			a[idx++] = (double)k / segmentPerEdge;
		}
	}

	float *b = new float[samplePointPerTriangle * triangleCtrlPointNum_lower];
	for (int row = 0; row < samplePointPerTriangle; ++row)
	{
		int idx = 0;
		for (int i = degree_lower; i >= 0; --i)
		{
			for (int j = degree_lower - i; j >= 0; --j)
			{
				int k = degree_lower - i - j;
				double u = a[row * 3 + 0];
				double v = a[row * 3 + 1];
				double w = a[row * 3 + 2];
				b[index2c(row, idx, samplePointPerTriangle)] = B(u, v, w, degree_lower, make_int3(i, j, k));
				++idx;
			}
		}
	}

	float *b_PN = new float[samplePointPerTriangle * 6];
	for (int row = 0; row < samplePointPerTriangle; ++row)
	{
		int idx = 0;
		for (int i = 2; i >= 0; --i)
		{
			for (int j = 2 - i; j >= 0; --j)
			{
				int k = 2 - i - j;
				double u = a[row * 3 + 0];
				double v = a[row * 3 + 1];
				double w = a[row * 3 + 2];
				b_PN[index2c(row, idx, samplePointPerTriangle)] = B(u, v, w, 2, make_int3(i, j, k));
				++idx;
			}
		}
	}

	cudaMalloc(&BqD, sizeof(float) * samplePointPerTriangle * triangleCtrlPointNum_lower);
	tessMemD += sizeof(float) * samplePointPerTriangle * triangleCtrlPointNum_lower;
	printMemD(__FILE__, __FUNCTION__, __LINE__, sizeof(float) * samplePointPerTriangle * triangleCtrlPointNum_lower, "@第一个矩阵乘法用到的矩阵Bq存放在这里");
	cudaMemcpy(BqD, b, sizeof(float) * samplePointPerTriangle * triangleCtrlPointNum_lower, cudaMemcpyHostToDevice);

	cudaMalloc(&BqD_PN, sizeof(float) * samplePointPerTriangle * 6);
	tessMemD += sizeof(float) * samplePointPerTriangle * 6;
	printMemD(__FILE__, __FUNCTION__, __LINE__, sizeof(float) * samplePointPerTriangle * 6, "@第一个矩阵乘法用到的矩阵Bq存放在这里");
	cudaMemcpy(BqD_PN, b_PN, sizeof(float) * samplePointPerTriangle * 6, cudaMemcpyHostToDevice);

	/***********************************************************************************************************************************/
	cudaMalloc(&BBD, sizeof(float) * samplePointPerTriangle * (triangleCtrlPointNum + constrait_point_num));
	tessMemD += sizeof(float) * samplePointPerTriangle * (triangleCtrlPointNum + constrait_point_num);
	printMemD(__FILE__, __FUNCTION__, __LINE__, sizeof(float) * samplePointPerTriangle * (triangleCtrlPointNum + constrait_point_num), "@第二个矩阵乘法用到的矩阵BB存放在这里");

	cudaMalloc(&RD, sizeof(float) * samplePointPerTriangle * triangleNum * 6);
	tessMemD += sizeof(float) * samplePointPerTriangle * triangleNum * 6;
	printMemD(__FILE__, __FUNCTION__, __LINE__, sizeof(float) * samplePointPerTriangle * triangleNum * 6, "@第二个矩阵乘法的结果RD存放在这里");

	cudaMalloc(&parameter3D, sizeof(float3) * samplePointPerTriangle * triangleNum);
	cudaMalloc(&parameterND, sizeof(float3) * samplePointPerTriangle * triangleNum);

	delete []a;
	delete []b;

	cudaMalloc(&my_to_truth_tableD, sizeof(int) * samplePointPerTriangle * triangleNum);
	my_to_truth_table = new int[samplePointPerTriangle * triangleNum];
	fill(my_to_truth_table, my_to_truth_table + samplePointPerTriangle * triangleNum, 0);
}

#ifdef TRUTH
void generateUVW_truth(int samplePointPerEdge)
{
	double *a = new double[samplePointPerTriangle * 3];
	int idx = 0;
	for (int i = segmentPerEdge; i >= 0; --i)
	{
		for (int j = segmentPerEdge - i; j >= 0; --j)
		{
			int k = segmentPerEdge - i - j;
			a[idx++] = (double)i / segmentPerEdge;
			a[idx++] = (double)j / segmentPerEdge;
			a[idx++] = (double)k / segmentPerEdge;
		}
	}

	float *b = new float[samplePointPerTriangle * triangleCtrlPointNum * 3];
	for (int row = 0; row < samplePointPerTriangle; ++row)
	{
		int idx = 0;
		for (int i = degree; i >= 0; --i)
		{
			for (int j = degree - i; j >= 0; --j)
			{
				int k = degree - i - j;
				double u = a[row * 3 + 0];
				double v = a[row * 3 + 1];
				double w = a[row * 3 + 2];
				b[index2c(row, idx, samplePointPerTriangle * 3)] = B(u, v, w, degree, make_int3(i, j, k));
				//b[row * triangleCtrlPointNum + idx] = B(u, v, w, degree, make_int3(i, j, k));
				++idx;
			}
		}
	}

	/***********************************************************************************************************************************/

	for (int row = 0; row < samplePointPerTriangle; ++row)
	{
		int idx = 0;
		for (int i = degree; i >= 0; --i)
		{
			for (int j = degree - i; j >= 0; --j)
			{
				int k = degree - i - j;
				double u = a[row * 3 + 0];
				double v = a[row * 3 + 1];
				double w = a[row * 3 + 2];
				b[index2c(row + samplePointPerTriangle, idx, samplePointPerTriangle * 3)] = factorial(degree) / (factorial(i) * factorial(j) * factorial(k)) *
												 (i * power(u, i - 1) * power(v, j) * power(w, k) - k * power(u, i) * power(v, j) * power(w, k - 1));
				++idx;
			}
		}
	}

	/***********************************************************************************************************************************/

	for (int row = 0; row < samplePointPerTriangle; ++row)
	//for (int row = samplePointPerTriangle * 2; row < samplePointPerTriangle * 3; ++row)
	{
		int idx = 0;
		for (int i = degree; i >= 0; --i)
		{
			for (int j = degree - i; j >= 0; --j)
			{
				int k = degree - i - j;
				double u = a[row * 3 + 0];
				double v = a[row * 3 + 1];
				double w = a[row * 3 + 2];
				b[index2c(row + samplePointPerTriangle * 2, idx, samplePointPerTriangle * 3)] = factorial(degree) / (factorial(i) * factorial(j) * factorial(k)) *
												 (j * power(u, i) * power(v, j - 1) * power(w, k) - k * power(u, i) * power(v, j) * power(w, k - 1));
				++idx;
			}
		}
	}
	cudaMalloc(&BqD_truth, sizeof(float) * samplePointPerTriangle * triangleCtrlPointNum * 3);
	tessMemD += sizeof(float) * samplePointPerTriangle * triangleCtrlPointNum * 3;
	printMemD(__FILE__, __FUNCTION__, __LINE__, sizeof(float) * samplePointPerTriangle * triangleCtrlPointNum * 3, "@第一个矩阵乘法用到的矩阵Bq存放在这里");
	cudaMemcpy(BqD_truth, b, sizeof(float) * samplePointPerTriangle * triangleCtrlPointNum * 3, cudaMemcpyHostToDevice);

	/***********************************************************************************************************************************/
	cudaMalloc(&BBD_truth, sizeof(float) * samplePointPerTriangle * triangleCtrlPointNum * 3);
	tessMemD += sizeof(float) * samplePointPerTriangle * triangleCtrlPointNum * 3;
	printMemD(__FILE__, __FUNCTION__, __LINE__, sizeof(float) * samplePointPerTriangle * triangleCtrlPointNum * 3, "@第二个矩阵乘法用到的矩阵BB存放在这里");

	cudaMalloc(&RD_truth, sizeof(float) * samplePointPerTriangle * 3 * triangleNum * 3);
	tessMemD += sizeof(float) * samplePointPerTriangle * 3 * triangleNum * 3;
	printMemD(__FILE__, __FUNCTION__, __LINE__, sizeof(float) * samplePointPerTriangle * 3 * triangleNum * 3, "@第二个矩阵乘法的结果RD存放在这里");

	delete []a;
	delete []b;
}
#endif

/* 
 * 使用矩阵乘法求 B 样条体的值，和上面一个类似函数的区别在于不负责 u、v、w 重新参数化的工作，
 * 而且也不负责求合适的 B 样条矩阵，这两项工作均需调用函数之前完成，参数列表得到简化
 * 目前仅用于求truth
 */
__device__ float3 BSplineVolumeValueMatrixD2(float *Mu, float *Mv, float *Mw,
											 float u, float v, float w, float *shared_array,
											 int leftUIdx, int leftVIdx, int leftWIdx,
											 int orderU, int orderV, int orderW)
{
#define NB		// NB表示使用比较好的算法，如果不define NB，则使用最原始的算法，逻辑也相对清晰
#ifdef NB
	float *mul1 = (float *)shared_array;
	float *mul2 = (float *)&mul1[blockDim.x * 4];
	float *temp = (float *)&mul2[blockDim.x * 4];

	// 由三维控制顶点算出二维临时控制顶点
	temp[3 * threadIdx.x + 0] = w;
	temp[3 * threadIdx.x + 1] = w * w;
	temp[3 * threadIdx.x + 2] = w * w * w;

	for (int i = 0; i < orderW; ++i)
	{
		mul1[4 * threadIdx.x + i] = Mw[i];
		for (int j = 1; j < orderW; ++j)
			mul1[4 * threadIdx.x + i] += temp[3 * threadIdx.x + j - 1] * Mw[j * orderW + i];
	}

	// 由二维临时控制顶点算出一维临时控制顶点
	temp[3 * threadIdx.x + 0] = v;
	temp[3 * threadIdx.x + 1] = v * v;
	temp[3 * threadIdx.x + 2] = v * v * v;

	for (int i = 0; i < orderV; ++i)
	{
		mul2[4 * threadIdx.x + i] = Mv[i];
		for (int j = 1; j < orderV; ++j)
			mul2[4 * threadIdx.x + i] += temp[3 * threadIdx.x + j - 1] * Mv[j * orderV + i];
	}

	float3 tempCtrlPoint2[4];
	float3 tempCtrlPoint1[4];
	for (int i = 0; i < orderU; ++i)
	{
		for (int j = 0; j < orderV; ++j)
		{
			tempCtrlPoint2[j] = make_float3(0.0f, 0.0f, 0.0f);
			for (int k = 0; k < orderW; ++k)
			{
				float3 cp = ctrlPointD[leftUIdx - i][leftVIdx - j][leftWIdx - k];
				tempCtrlPoint2[j] += cp * mul1[4 * threadIdx.x + orderW - 1 - k];
			}
		}
		tempCtrlPoint1[i] = make_float3(0.0f, 0.0f, 0.0f);
		for (int j = 0; j < orderV; ++j)
			tempCtrlPoint1[i] += tempCtrlPoint2[j] * mul2[4 * threadIdx.x + orderV - 1 - j];
	}

	// 由一维临时控制顶点算出结果
	temp[3 * threadIdx.x + 0] = u;
	temp[3 * threadIdx.x + 1] = u * u;
	temp[3 * threadIdx.x + 2] = u * u * u;

	for (int i = 0; i < orderU; ++i)
	{
		mul1[4 * threadIdx.x + i] = Mu[i];
		for (int j = 1; j < orderU; ++j)
			mul1[4 * threadIdx.x + i] += temp[3 * threadIdx.x + j - 1] * Mu[j * orderU + i];
	}
	float3 result = make_float3(0.0f, 0.0f, 0.0f);
	for (int i = 0; i < orderU; ++i)
		result += tempCtrlPoint1[i] * mul1[4 * threadIdx.x + orderU - 1 - i];

	return result;

	/*-------------------------------------------------------------------------------------------------*/

#else

	// 第一个矩阵乘法
	int base_i = leftUIdx - orderU + 1;
	int base_j = leftVIdx - orderV + 1;
	int base_k = leftWIdx - orderW + 1;
	float3 box[4][4][4], temp;
	for (int k = 0; k < orderW; ++k)
		for (int i = 0; i < orderU; ++i)
			for (int j = 0; j < orderV; ++j)
			{
				temp = make_float3(0.0, 0.0, 0.0);
				for (int l = 0; l < orderU; ++l)
				{
					float3 cp = ctrlPointD[base_i + l][base_j + j][base_k + k];
					temp += Mu[i * orderU + l] * cp;
				}
				box[i][j][k] = temp;
			}

	// 第二个矩阵乘法
	float3 box1[4][4][4];
	for (int i = 0; i < orderU; ++i)
		for (int j = 0; j < orderV; ++j)
			for (int k = 0; k < orderW; ++k)
			{
				temp = make_float3(0.0, 0.0, 0.0);
				for (int l = 0; l < orderV; ++l)
				{
					float3 cp = box[i][l][k];
					temp += Mv[j * orderV + l] * cp;
				}
				box1[i][j][k] = temp;
			}

	// 第三个矩阵乘法
	for (int j = 0; j < orderV; ++j)
		for (int k = 0; k < orderW; ++k)
			for (int i = 0; i < orderU; ++i)
			{
				temp = make_float3(0.0, 0.0, 0.0);
				for (int l = 0; l < orderW; ++l)
				{
					float3 cp = box1[i][j][l];
					temp += Mw[k * orderW + l] * cp;
				}
				box[i][j][k] = temp;
			}

	// 由三维控制顶点算出二维临时控制顶点
	float t[4];
	t[0] = 1.0f;
	t[1] = u;
	t[2] = u * u;
	t[3] = t[2] * u;

	float3 cp2D[4][4];
	for (int j = 0; j < orderV; ++j)
		for (int k = 0; k < orderW; ++k)
		{
			cp2D[j][k] = make_float3(0.0f, 0.0f, 0.0f);
			for (int i = 0; i < orderU; ++i)
			{
				cp2D[j][k] += t[i] * box[i][j][k];
			}
		}

	// 由二维临时控制顶点算出一维临时控制顶点
	t[1] = v;
	t[2] = v * v;
	t[3] = t[2] * v;

	float3 cp1D[4];
	for (int k = 0; k < orderW; ++k)
	{
		cp1D[k] = make_float3(0.0f, 0.0f, 0.0f);
		for (int j = 0; j < orderV; ++j)
			cp1D[k] += t[j] * cp2D[j][k];
	}

	// 由一维临时控制顶点算出结果
	t[1] = w;
	t[2] = w * w;
	t[3] = t[2] * w;

	temp = make_float3(0.0f, 0.0f, 0.0f);
	for (int k = 0; k < orderW; ++k)
		temp += t[k] * cp1D[k];

	return temp;
#endif
}

/* 新的合并算法 */
__device__ void BSplineVolumeValueMatrixD_combine(float u, float v, float w, float *shared_array,
											 int i_idx, int j_idx, int k_idx,
											 int orderU, int orderV, int orderW,
											 float3 &f, float3 &fu, float3 &fv)
{
	int base2 = 2 * threadIdx.x;
	int base3 = 3 * threadIdx.x;
	float *tu = &shared_array[base3];
	float *tu_ = &shared_array[blockDim.x * 3 + base2];
	float *tv = &shared_array[blockDim.x * 5 + base3];
	float *tv_ = &shared_array[blockDim.x * 8 + base2];
	float *tw = &shared_array[blockDim.x * 10 + base3];

	tu[0] = u; tu[1] = u * u, tu[2] = u * tu[1];
	tu_[0] = 2 * u; tu_[1] = 3 * tu[1];

	tv[0] = v; tv[1] = v * v; tv[2] = v * tv[1];
	tv_[0] = 2 * v; tv_[1] = 3 * tv[1];

	tw[0] = w; tw[1] = w * w; tw[2] = w * tw[1];

	/************* 将i = 0 提到前面，减少tu和tu_数组的大小 ****************/
	/******** orderU至少是2，所以这里可以将i = 0的情况提到for之外 *********/
	float3 cp2D[4];
	for (int j = 0; j < orderV; ++j)
	{
		cp2D[j] = newCtrlPointD[i_idx][j_idx][k_idx][0][j][0];
		for (int k = 1; k < orderW; ++k)
			cp2D[j] += tw[k - 1] * newCtrlPointD[i_idx][j_idx][k_idx][0][j][k];
	}
	// orderV至少是2，所以这里可以将tv[0] * cp2D[1]提到for之外
	float3 cp1D = cp2D[0] + tv[0] * cp2D[1], cp1Dv = cp2D[1];
	for (int j = 2; j < orderV; ++j)
	{
		cp1D += tv[j - 1] * cp2D[j];
		cp1Dv += tv_[j - 2] * cp2D[j];
	}
	f = cp1D;
	fv = cp1Dv;

	/*************** 将i = 1 提到前面，减少tu_数组的大小 ******************/
	/******** orderU至少是2，所以这里可以将i = 1的情况提到for之外 *********/
	for (int j = 0; j < orderV; ++j)
	{
		cp2D[j] = newCtrlPointD[i_idx][j_idx][k_idx][1][j][0];
		for (int k = 1; k < orderW; ++k)
			cp2D[j] += tw[k - 1] * newCtrlPointD[i_idx][j_idx][k_idx][1][j][k];
	}
	// orderV至少是2，所以这里可以将tv[0] * cp2D[1]提到for之外
	cp1D = cp2D[0] + tv[0] * cp2D[1];
	cp1Dv = cp2D[1];
	for (int j = 2; j < orderV; ++j)
	{
		cp1D += tv[j - 1] * cp2D[j];
		cp1Dv += tv_[j - 2] * cp2D[j];
	}
	f += tu[0] * cp1D;
	fu = cp1D;
	fv += tu[0] * cp1Dv;

	/*********************************************************************/
	for (int i = 2; i < orderU; ++i)
	{
		for (int j = 0; j < orderV; ++j)
		{
			cp2D[j] = newCtrlPointD[i_idx][j_idx][k_idx][i][j][0];
			for (int k = 1; k < orderW; ++k)
				cp2D[j] += tw[k - 1] * newCtrlPointD[i_idx][j_idx][k_idx][i][j][k];
		}
		// orderV至少是2，所以这里可以将tv[0] * cp2D[1]提到for之外
		cp1D = cp2D[0] + tv[0] * cp2D[1];
		cp1Dv = cp2D[1];
		for (int j = 2; j < orderV; ++j)
		{
			cp1D += tv[j - 1] * cp2D[j];
			cp1Dv += tv_[j - 2] * cp2D[j];
		}
		f += tu[i - 1] * cp1D;
		fu += tu_[i - 2] * cp1D;
		fv += tu[i - 1] * cp1Dv;
	}
}

/* 求采样点在u方向的偏导，由优化之后的采样点求值算法改造而来 */
__device__ float3 BSplineVolumeValueMatrixDu(float u, float v, float w, float *shared_array,
											 int i_idx, int j_idx, int k_idx,
											 int orderU, int orderV, int orderW)
{
	float *tu = (float *)shared_array;
	float *tv = (float *)&tu[blockDim.x * 2];
	float *tw = (float *)&tv[blockDim.x * 3];
	int base2 = 2 * threadIdx.x;
	int base3 = 3 * threadIdx.x;

	tu[base2] = 2 * u; tu[base2 + 1] = 3 * u * u;
	tv[base3] = v; tv[base3 + 1] = v * v; tv[base3 + 2] = v * v * v;
	tw[base3] = w; tw[base3 + 1] = w * w, tw[base3 + 2] = w * w * w;

	// 一步完成三维控制顶点->二维临时控制顶点->一维临时控制顶点->结果
	float3 cp2D[4], cp1D, result;
	for (int j = 0; j < orderV; ++j)
	{
		cp2D[j] = newCtrlPointD[i_idx][j_idx][k_idx][1][j][0];
		for (int k = 1; k < orderU; ++k)
			cp2D[j] += tw[base3 + k - 1] * newCtrlPointD[i_idx][j_idx][k_idx][1][j][k];
	}
	cp1D = cp2D[0];
	for (int j = 1; j < orderV; ++j)
		cp1D += tv[base3 + j - 1] * cp2D[j];
	result = cp1D;

	// 为了把tu从[3]缩成[2]，将i=0的情况提到了前面
	for (int i = 2; i < orderU; ++i)
	{
		for (int j = 0; j < orderV; ++j)
		{
			cp2D[j] = newCtrlPointD[i_idx][j_idx][k_idx][i][j][0];
			for (int k = 1; k < orderW; ++k)
				cp2D[j] += tw[base3 + k - 1] * newCtrlPointD[i_idx][j_idx][k_idx][i][j][k];
		}
		cp1D = cp2D[0];
		for (int j = 1; j < orderV; ++j)
			cp1D += tv[base3 + j - 1] * cp2D[j];

		result += tu[base2 + i - 2] * cp1D;
	}
	return result;

	/*-----------------------------------------------------*/

	/*float tu[4], tv[4], tw[4];*/
	/*tu[0] = 0; tu[1] = 1; tu[2] = 2 * u; tu[3] = 3 * u * u;*/
	/*tv[0] = 1; tv[1] = v; tv[2] = v * v; tv[3] = v * v * v;*/
	/*tw[0] = 1; tw[1] = w; tw[2] = w * w, tw[3] = w * w * w;*/

	/*// 一步完成三维控制顶点->二维临时控制顶点->一维临时控制顶点->结果*/
	/*float3 cp2D[4], cp1D, result = make_float3(0.0f, 0.0f, 0.0f);*/
	/*for (int k = 0; k < orderW; ++k)*/
	/*{*/
		/*for (int j = 0; j < orderV; ++j)*/
		/*{*/
			/*cp2D[j] = make_float3(0.0f, 0.0f, 0.0f);*/
			/*for (int i = 0; i < orderU; ++i)*/
				/*cp2D[j] += tu[i] * newCtrlPointD[i_idx][j_idx][k_idx][i][j][k];*/
		/*}*/
		/*cp1D = make_float3(0.0f, 0.0f, 0.0f);*/
		/*for (int j = 0; j < orderV; ++j)*/
			/*cp1D += tv[j] * cp2D[j];*/

		/*result += tw[k] * cp1D;*/
	/*}*/
	/*return result;*/
}

/* 求采样点在v方向的偏导，由优化之后的采样点求值算法改造而来 */
__device__ float3 BSplineVolumeValueMatrixDv(float u, float v, float w, float *shared_array,
											 int i_idx, int j_idx, int k_idx,
											 int orderU, int orderV, int orderW)
{
	float *tu = (float *)shared_array;
	float *tv = (float *)&tu[blockDim.x * 3];
	float *tw = (float *)&tv[blockDim.x * 2];
	int base2 = 2 * threadIdx.x;
	int base3 = 3 * threadIdx.x;
	tu[base3] = u; tu[base3 + 1] = u * u; tu[base3 + 2] = u * u * u;
	tv[base2] = 2 * v; tv[base2 + 1] = 3 * v * v;
	tw[base3] = w; tw[base3 + 1] = w * w, tw[base3 + 2] = w * w * w;

	// 一步完成三维控制顶点->二维临时控制顶点->一维临时控制顶点->结果
	float3 cp2D[4], cp1D, result;
	for (int k = 0; k < orderW; ++k)
	{
		cp2D[k] = newCtrlPointD[i_idx][j_idx][k_idx][0][1][k];
		for (int i = 1; i < orderU; ++i)
			cp2D[k] += tu[base3 + i - 1] * newCtrlPointD[i_idx][j_idx][k_idx][i][1][k];
	}
	cp1D = cp2D[0];
	for (int k = 1; k < orderW; ++k)
		cp1D += tw[base3 + k - 1] * cp2D[k];
	result = cp1D;

	// 为了把tv从[3]缩成[2]，将j=0的情况提到了前面
	for (int j = 2; j < orderV; ++j)
	{
		for (int k = 0; k < orderW; ++k)
		{
			cp2D[k] = newCtrlPointD[i_idx][j_idx][k_idx][0][j][k];
			for (int i = 1; i < orderU; ++i)
				cp2D[k] += tu[base3 + i - 1] * newCtrlPointD[i_idx][j_idx][k_idx][i][j][k];
		}
		cp1D = cp2D[0];
		for (int k = 1; k < orderW; ++k)
			cp1D += tw[base3 + k - 1] * cp2D[k];

		result += tv[base2 + j - 2] * cp1D;
	}
	return result;

	/*-----------------------------------------------------*/

	/*float tu[4], tv[4], tw[4];*/
	/*tu[0] = 1; tu[1] = u; tu[2] = u * u; tu[3] = u * u * u;*/
	/*tv[0] = 0; tv[1] = 1; tv[2] = 2 * v; tv[3] = 3 * v * v;*/
	/*tw[0] = 1; tw[1] = w; tw[2] = w * w, tw[3] = w * w * w;*/

	/*// 一步完成三维控制顶点->二维临时控制顶点->一维临时控制顶点->结果*/
	/*float3 cp2D[4], cp1D, result = make_float3(0.0f, 0.0f, 0.0f);*/
	/*for (int k = 0; k < orderW; ++k)*/
	/*{*/
		/*for (int j = 0; j < orderV; ++j)*/
		/*{*/
			/*cp2D[j] = make_float3(0.0f, 0.0f, 0.0f);*/
			/*for (int i = 0; i < orderU; ++i)*/
				/*cp2D[j] += tu[i] * newCtrlPointD[i_idx][j_idx][k_idx][i][j][k];*/
		/*}*/
		/*cp1D = make_float3(0.0f, 0.0f, 0.0f);*/
		/*for (int j = 0; j < orderV; ++j)*/
			/*cp1D += tv[j] * cp2D[j];*/

		/*result += tw[k] * cp1D;*/
	/*}*/
	/*return result;*/
}

/* 求采样点在w方向的偏导，由优化之后的采样点求值算法改造而来 */
__device__ float3 BSplineVolumeValueMatrixDw(float u, float v, float w, float *shared_array,
											 int i_idx, int j_idx, int k_idx,
											 int orderU, int orderV, int orderW)
{
	int base2 = 2 * threadIdx.x;
	int base3 = 3 * threadIdx.x;
	float *tu = &shared_array[base3];
	float *tv = &shared_array[blockDim.x * 3 + base3];
	float *tw = &shared_array[blockDim.x * 6 + base2];
	tu[0] = u; tu[1] = u * u; tu[2] = u * tu[1];
	tv[0] = v; tv[1] = v * v; tv[2] = v * tv[1];
	tw[0] = 2 * w; tw[1] = 3 * w * w;

	float3 cp2D[4], cp1D, result;
	for (int i = 0; i < orderU; ++i)
	{
		cp2D[i] = newCtrlPointD[i_idx][j_idx][k_idx][i][0][1];
		for (int j = 1; j < orderV; ++j)
			cp2D[i] += tv[j - 1] * newCtrlPointD[i_idx][j_idx][k_idx][i][j][1];
	}
	cp1D = cp2D[0];
	for (int i = 1; i < orderU; ++i)
		cp1D += tu[i - 1] * cp2D[i];
	result = cp1D;

	// 为了把tw从[3]缩成[2]，将k = 1的情况提到了前面
	for (int k = 2; k < orderW; ++k)
	{
		for (int i = 0; i < orderU; ++i)
		{
			cp2D[i] = newCtrlPointD[i_idx][j_idx][k_idx][i][0][k];
			for (int j = 1; j < orderV; ++j)
				cp2D[i] += tv[j - 1] * newCtrlPointD[i_idx][j_idx][k_idx][i][j][k];
		}
		cp1D = cp2D[0];
		for (int i = 1; i < orderU; ++i)
			cp1D += tu[i - 1] * cp2D[i];

		result += tw[k - 2] * cp1D;
	}
	return result;

	/*-----------------------------------------------------*/

	//float tu[4], tv[4], tw[4];
	//tu[0] = 1; tu[1] = u; tu[2] = u * u; tu[3] = u * u * u;
	//tv[0] = 1; tv[1] = v; tv[2] = v * v, tv[3] = v * v * v;
	//tw[0] = 0; tw[1] = 1; tw[2] = 2 * w; tw[3] = 3 * w * w;

	//float3 cp2D[4], cp1D, result = make_float3(0.0f, 0.0f, 0.0f);
	//for (int i = 0; i < orderU; ++i)
	//{
		//for (int j = 0; j < orderV; ++j)
		//{
			//cp2D[j] = make_float3(0.0f, 0.0f, 0.0f);
			//for (int k = 0; k < orderW; ++k)
				//cp2D[j] += tw[k] * newCtrlPointD[i_idx][j_idx][k_idx][i][j][k];
		//}
		//cp1D = make_float3(0.0f, 0.0f, 0.0f);
		//for (int j = 0; j < orderV; ++j)
			//cp1D += tv[j] * cp2D[j];

		//result += tu[i] * cp1D;
	//}
	//return result;

	/*-----------------------------------------------------*/

	// 最原始的算法
	/*float tu[4], tv[4], tw[4];*/
	/*tu[0] = 1; tu[1] = u; tu[2] = u * u; tu[3] = u * u * u;*/
	/*tv[0] = 1; tv[1] = v; tv[2] = v * v, tv[3] = v * v * v;*/
	/*tw[0] = 0; tw[1] = 1; tw[2] = 2 * w; tw[3] = 3 * w * w;*/

	/*float3 cp2D[4], cp1D, result = make_float3(0.0f, 0.0f, 0.0f);*/
	/*for (int k = 0; k < orderW; ++k)*/
	/*{*/
		/*for (int j = 0; j < orderV; ++j)*/
		/*{*/
			/*cp2D[j] = make_float3(0.0f, 0.0f, 0.0f);*/
			/*for (int i = 0; i < orderU; ++i)*/
				/*cp2D[j] += tu[i] * newCtrlPointD[i_idx][j_idx][k_idx][i][j][k];*/
		/*}*/
		/*cp1D = make_float3(0.0f, 0.0f, 0.0f);*/
		/*for (int j = 0; j < orderV; ++j)*/
			/*cp1D += tv[j] * cp2D[j];*/

		/*result += tw[k] * cp1D;*/
	/*}*/
	/*return result;*/
}

__global__ void calcCtrlPoint_PN(TriangleD *triangleListD, int *triangle_adjacent_tableD, float3 *sampleValueD_PN, float *triangleCtrlPointD_PN, float *triangleNormalCtrlPointD_PN, int f, int m_)
{
	int triangleIdx = blockDim.x * blockIdx.x + threadIdx.x;
	if (triangleIdx >= f)
		return;

	int adj_face_idx[3];
	adj_face_idx[0] = triangle_adjacent_tableD[triangleIdx * 3];
	adj_face_idx[1] = triangle_adjacent_tableD[triangleIdx * 3 + 1];
	adj_face_idx[2] = triangle_adjacent_tableD[triangleIdx * 3 + 2];

	//int adj_edge_idx[3] = { -1, -1, -1 };		// 实际上应该使用这一句，以此判断有没有相邻三角形
	int adj_edge_idx[3] = { 0, 0, 0 };			// 但是对于某些模型使用上一句会出现内存越界，所以暂且使用这一句，权宜之计
	//bool handle[3] = { false, false, false };	// 在该边有一个以上法向时，是否处理这条边。当这条边有相邻面片才需要处理
	int adj_corner_ctrlpoint_idx[3][2] = { { 0, 2 }, { 1, 0 }, { 2, 1 } };		// 相邻三角形0, 1, 2号边上的控制顶点编号（仅有角点，没有边点）
	for (int i = 0; i < 3; ++i)
		if (adj_face_idx[i] >= 0)
		{
			adj_edge_idx[i] = adj_face_idx[i] & 0x3;
			adj_face_idx[i] = adj_face_idx[i] >> 2;
			//handle[i] = true;
		}
	//printf("edge_id = (%d, %d, %d), face_id = (%d, %d, %d)\n", adj_edge_idx[0], adj_edge_idx[1], adj_edge_idx[2],
																//adj_face_idx[0], adj_face_idx[1], adj_face_idx[2]);

	int n_count[3];
	n_count[0] = triangleListD[triangleIdx].nc[0];
	n_count[1] = triangleListD[triangleIdx].nc[1];
	n_count[2] = triangleListD[triangleIdx].nc[2];

	float *p_x  = &triangleCtrlPointD_PN[m_ * triangleIdx];
	float *p_y  = &triangleCtrlPointD_PN[m_ * (f + triangleIdx)];
	float *p_z  = &triangleCtrlPointD_PN[m_ * (f * 2 + triangleIdx)];

	float3 v0 = sampleValueD_PN[triangleIdx * 3];
	float3 v1 = sampleValueD_PN[triangleIdx * 3 + 1];
	float3 v2 = sampleValueD_PN[triangleIdx * 3 + 2];
	float3 n0 = sampleValueD_PN[(f + triangleIdx) * 3];
	float3 n1 = sampleValueD_PN[(f + triangleIdx) * 3 + 1];
	float3 n2 = sampleValueD_PN[(f + triangleIdx) * 3 + 2];
	normalize(n0);
	normalize(n1);
	normalize(n2);

	/*********************** 计算几何控制顶点 **********************/
	p_x[0] = v0.x; p_y[0] = v0.y; p_z[0] = v0.z; // 控制顶点0
	p_x[6] = v1.x; p_y[6] = v1.y; p_z[6] = v1.z; // 控制顶点6
	p_x[9] = v2.x; p_y[9] = v2.y; p_z[9] = v2.z; // 控制顶点9

	float3 e = make_float3(0.0f, 0.0f, 0.0f);

	float3 v01 = v1 - v0;
	float3 result;
	if (n_count[1] < 2)		// 该条边只有一个法向
	{
		result = (v0 * 2 + v1 - n0 * (v01 * n0)) / 3;
	}
	else
	{
		float3 n_oppo = triangleListD[adj_face_idx[1]].n_adj[adj_corner_ctrlpoint_idx[adj_edge_idx[1]][0]];
		float3 n_ave = cross(n0, n_oppo);
		normalize(n_ave);
		result = v0 + v01 * n_ave / 3 * n_ave;
	}
	e += result;
	p_x[1] = result.x; p_y[1] = result.y; p_z[1] = result.z; // 控制顶点1

	float3 v02 = v2 - v0;
	if (n_count[0] < 2)		// 该条边只有一个法向
	{
		result = (v0 * 2 + v2 - n0 * (v02 * n0)) / 3;
	}
	else
	{
		float3 n_oppo = triangleListD[adj_face_idx[0]].n_adj[adj_corner_ctrlpoint_idx[adj_edge_idx[0]][1]];
		float3 n_ave = cross(n0, n_oppo);
		normalize(n_ave);
		result = v0 + v02 * n_ave / 3 * n_ave;
	}
	e += result;
	p_x[2] = result.x; p_y[2] = result.y; p_z[2] = result.z; // 控制顶点2

	float3 v10 = v0 - v1;
	if (n_count[1] < 2)		// 该条边只有一个法向
	{
		result = (v1 * 2 + v0 - n1 * (v10 * n1)) / 3;
	}
	else
	{
		float3 n_oppo = triangleListD[adj_face_idx[1]].n_adj[adj_corner_ctrlpoint_idx[adj_edge_idx[1]][1]];
		float3 n_ave = cross(n1, n_oppo);
		normalize(n_ave);
		result = v1 + v10 * n_ave / 3 * n_ave;
	}
	e += result;
	p_x[3] = result.x; p_y[3] = result.y; p_z[3] = result.z; // 控制顶点3

	float3 v12 = v2 - v1;
	if (n_count[2] < 2)		// 该条边只有一个法向
	{
		result = (v1 * 2 + v2 - n1 * (v12 * n1)) / 3;
	}
	else
	{
		float3 n_oppo = triangleListD[adj_face_idx[2]].n_adj[adj_corner_ctrlpoint_idx[adj_edge_idx[2]][0]];
		float3 n_ave = cross(n1, n_oppo);
		normalize(n_ave);
		result = v1 + v12 * n_ave / 3 * n_ave;
	}
	e += result;
	p_x[7] = result.x; p_y[7] = result.y; p_z[7] = result.z; // 控制顶点7

	float3 v20 = v0 - v2;
	if (n_count[0] < 2)		// 该条边只有一个法向
	{
		result = (v2 * 2 + v0 - n2 * (v20 * n2)) / 3;
	}
	else
	{
		float3 n_oppo = triangleListD[adj_face_idx[0]].n_adj[adj_corner_ctrlpoint_idx[adj_edge_idx[0]][0]];
		float3 n_ave = cross(n2, n_oppo);
		normalize(n_ave);
		result = v2 + v20 * n_ave / 3 * n_ave;
	}
	e += result;
	p_x[5] = result.x; p_y[5] = result.y; p_z[5] = result.z; // 控制顶点5

	float3 v21 = v1 - v2;
	if (n_count[2] < 2)		// 该条边只有一个法向
	{
		result = (v2 * 2 + v1 - n2 * (v21 * n2)) / 3;
	}
	else
	{
		float3 n_oppo = triangleListD[adj_face_idx[2]].n_adj[adj_corner_ctrlpoint_idx[adj_edge_idx[2]][1]];
		float3 n_ave = cross(n2, n_oppo);
		normalize(n_ave);
		result = v2 + v21 * n_ave / 3 * n_ave;
	}
	e += result;
	p_x[8] = result.x; p_y[8] = result.y; p_z[8] = result.z; // 控制顶点8

	e /= 6;
	float3 v_total = (v0 + v1 + v2) / 3;
	result = e + (e - v_total) / 2;
	p_x[4] = result.x; p_y[4] = result.y; p_z[4] = result.z; // 控制顶点4


	/*********************** 计算法向控制顶点 **********************/
	p_x  = &triangleNormalCtrlPointD_PN[6 * triangleIdx];
	p_y  = &triangleNormalCtrlPointD_PN[6 * (f + triangleIdx)];
	p_z  = &triangleNormalCtrlPointD_PN[6 * (f * 2 + triangleIdx)];

	p_x[0] = n0.x; p_y[0] = n0.y; p_z[0] = n0.z; // 控制顶点0
	p_x[3] = n1.x; p_y[3] = n1.y; p_z[3] = n1.z; // 控制顶点3
	p_x[5] = n2.x; p_y[5] = n2.y; p_z[5] = n2.z; // 控制顶点5

	float value01 = 2 * v01 * (n0 + n1) / (v01 * v01);
	result = n0 + n1 - value01 * v01;
	normalize(result);
	p_x[1] = result.x; p_y[1] = result.y; p_z[1] = result.z; // 控制顶点1

	float value12 = 2 * v12 * (n1 + n2) / (v12 * v12);
	result = n1 + n2 - value12 * v12;
	normalize(result);
	p_x[4] = result.x; p_y[4] = result.y; p_z[4] = result.z; // 控制顶点4

	float value20 = 2 * v20 * (n2 + n0) / (v20 * v20);
	result = n2 + n0 - value20 * v20;
	normalize(result);
	p_x[2] = result.x; p_y[2] = result.y; p_z[2] = result.z; // 控制顶点2
}

__global__ void calcSampleValueThread_PN(TriangleD *triangleListD, float3 *sampleValueD_PN,
									  int f, int n, int orderU, int orderV, int orderW,
									  int ctrlPointNumU, int ctrlPointNumV, int ctrlPointNumW)
{
	int globalIdx = blockDim.x * blockIdx.x + threadIdx.x;
	if (globalIdx >= 3 * f)
		return;
	int triangleIdx = globalIdx / 3;
	TriangleD &triangle = triangleListD[triangleIdx];

	int localIdx = globalIdx % 3;
	float3 vertex = triangle.v[localIdx];
	float u = vertex.x;
	float v = vertex.y;
	float w = vertex.z;

	//float tempFloorFloat = (sqrtf((float)localIdx * 8 + 9) - 3) * 0.5;
	//int floor = rintf(tempFloorFloat);
	//if ((floor * 2 + 3) * (floor * 2 + 3) != localIdx * 8 + 9)
		//floor = ceilf(tempFloorFloat);
	//int room = localIdx - (floor + 1) * floor * 0.5;
	//float3 barycentric_coord;
	//barycentric_coord.x = (float)(n - floor) / n;
	//barycentric_coord.y = (float)(floor - room) / n;
	//barycentric_coord.z = 1.0f - barycentric_coord.x - barycentric_coord.y;

	//float3 v0 = triangle.v[0];
	//float3 v1 = triangle.v[1];
	//float3 v2 = triangle.v[2];

	//// u, v, w 表示经过重心坐标插值之后的采样点的x, y, z分量
	//float u = v0.x * barycentric_coord.x + v1.x * barycentric_coord.y + v2.x * barycentric_coord.z;
	//float v = v0.y * barycentric_coord.x + v1.y * barycentric_coord.y + v2.y * barycentric_coord.z;
	//float w = v0.z * barycentric_coord.x + v1.z * barycentric_coord.y + v2.z * barycentric_coord.z;

	// u, v, w方向节点区间数量
	int knot_interval_count_u = orderU + ctrlPointNumU - (orderU - 1) * 2 - 1;
	int knot_interval_count_v = orderV + ctrlPointNumV - (orderV - 1) * 2 - 1;
	int knot_interval_count_w = orderW + ctrlPointNumW - (orderW - 1) * 2 - 1;

	// 预先将其值设为最大，将末端点归入最后一段 
	int left_idx_u = orderU - 1 + knot_interval_count_u - 1;
	int left_idx_v = orderV - 1 + knot_interval_count_v - 1;
	int left_idx_w = orderW - 1 + knot_interval_count_w - 1;

	// 沿 U 方向查找当前点所在的节点区间 
	for (int ii = orderU - 1; ii <= orderU - 1 + knot_interval_count_u - 1; ++ii)
	{
		if (u >= knotListD[ii] && u < knotListD[ii + 1])
		{
			left_idx_u = ii;
			break;
		}
	}
	// 沿 V 方向查找当前点所在的节点区间 
	for (int jj = orderV - 1; jj <= orderV - 1 + knot_interval_count_v - 1; ++jj)
	{
		if (v >= knotListD[20 + jj] && v < knotListD[20 + jj + 1])
		{
			left_idx_v = jj;
			break;
		}
	}
	// 沿 W 方向查找当前点所在的节点区间 
	for (int kk = orderW - 1; kk <= orderW - 1 + knot_interval_count_w - 1; ++kk)
	{
		if (w >= knotListD[40 + kk] && w < knotListD[40 + kk + 1])
		{
			left_idx_w = kk;
			break;
		}
	}

	float tmpKnot = knotListD[left_idx_u];
	float tmpKnot1 = knotListD[left_idx_u + 1];
	float x_stride = tmpKnot1 - tmpKnot;
	u = (u - tmpKnot) / x_stride;

	tmpKnot = knotListD[20 + left_idx_v];
	tmpKnot1 = knotListD[20 + left_idx_v + 1];
	float y_stride = tmpKnot1 - tmpKnot;
	v = (v - tmpKnot) / y_stride;

	tmpKnot = knotListD[40 + left_idx_w];
	tmpKnot1 = knotListD[40 + left_idx_w + 1];
	float z_stride = tmpKnot1 - tmpKnot;
	w = (w - tmpKnot) / z_stride;

	extern __shared__ float shared_array[];
	// 算出该线程负责的采样点的 B 样条体值
	// fu 表示J_bar矩阵第一列三个元素：偏F_bar_x偏u、偏F_bar_y偏u、偏F_bar_z偏u
	// fv 表示J_bar矩阵第二列三个元素：偏F_bar_x偏v、偏F_bar_y偏v、偏F_bar_z偏v
	float3 result, fu, fv;
	BSplineVolumeValueMatrixD_combine(u, v, w, shared_array,
	left_idx_u - (orderU - 1), left_idx_v - (orderV - 1), left_idx_w - (orderW - 1),
											   orderU, orderV, orderW,
											   result, fu, fv);
	__syncthreads();

	//sampleValueD[index2c(localIdx, triangleIdx		  , 3)] = result.x;
	//sampleValueD[index2c(localIdx, triangleIdx + f	  , 3)] = result.y;
	//sampleValueD[index2c(localIdx, triangleIdx + f * 2, 3)] = result.z;
	sampleValueD_PN[3 * triangleIdx + localIdx].x = result.x;
	sampleValueD_PN[3 * triangleIdx + localIdx].y = result.y;
	sampleValueD_PN[3 * triangleIdx + localIdx].z = result.z;
	//printf("%d: result = (%f, %f, %f)\n", globalIdx, result.x, result.y, result.z);

	//printf("%d: result = (%f, %f, %f)\n", threadIdx.x, result.x, result.y, result.z);

	///////////////////////////////////////////////////////////////////////////////

	// fw 表示J_bar矩阵第三列三个元素：偏F_bar_x偏w、偏F_bar_y偏w、偏F_bar_z偏w
	float3 fw = BSplineVolumeValueMatrixDw(u, v, w, shared_array,
	left_idx_u - (orderU - 1), left_idx_v - (orderV - 1), left_idx_w - (orderW - 1),
										   orderU, orderV, orderW);
	//__syncthreads();

	//v0 = triangle.n[0];
	//v1 = triangle.n[1];
	//v2 = triangle.n[2];

	//// u, v, w 表示经过重心坐标插值之后的法向的x, y, z分量
	//u = v0.x * barycentric_coord.x + v1.x * barycentric_coord.y + v2.x * barycentric_coord.z;
	//v = v0.y * barycentric_coord.x + v1.y * barycentric_coord.y + v2.y * barycentric_coord.z;
	//w = v0.z * barycentric_coord.x + v1.z * barycentric_coord.y + v2.z * barycentric_coord.z;

	vertex = triangle.n[localIdx];
	u = vertex.x;
	v = vertex.y;
	w = vertex.z;

	float3 *sampleNormalD_PN = sampleValueD_PN + 3 * f;
	// J_bar_star_T_[012]表示J_bar的伴随矩阵的转置(即J_bar*T)的第一行三个元素
	float J_bar_star_T_0 = fv.y * fw.z - fw.y * fv.z;
	float J_bar_star_T_1 = fw.y * fu.z - fu.y * fw.z;
	float J_bar_star_T_2 = fu.y * fv.z - fv.y * fu.z;
	//sampleNormalD[index2c(localIdx, triangleIdx, 3)] = u * J_bar_star_T_0 * x_stride + v * J_bar_star_T_1 * y_stride + w * J_bar_star_T_2 * z_stride;
	sampleNormalD_PN[3 * triangleIdx + localIdx].x = u * J_bar_star_T_0 * x_stride + v * J_bar_star_T_1 * y_stride + w * J_bar_star_T_2 * z_stride;

	// J_bar_star_T_[012]表示J_bar的伴随矩阵的转置(即J_bar*T)的第二行三个元素
	J_bar_star_T_0 = fv.z * fw.x - fw.z * fv.x;
	J_bar_star_T_1 = fw.z * fu.x - fu.z * fw.x;
	J_bar_star_T_2 = fu.z * fv.x - fv.z * fu.x;
	sampleNormalD_PN[3 * triangleIdx + localIdx].y = u * J_bar_star_T_0 * x_stride + v * J_bar_star_T_1 * y_stride + w * J_bar_star_T_2 * z_stride;

	// J_bar_star_T_[012]表示J_bar的伴随矩阵的转置(即J_bar*T)的第三行三个元素
	J_bar_star_T_0 = fv.x * fw.y - fw.x * fv.y;
	J_bar_star_T_1 = fw.x * fu.y - fu.x * fw.y;
	J_bar_star_T_2 = fu.x * fv.y - fv.x * fu.y;
	sampleNormalD_PN[3 * triangleIdx + localIdx].z = u * J_bar_star_T_0 * x_stride + v * J_bar_star_T_1 * y_stride + w * J_bar_star_T_2 * z_stride;
}

__global__ void calcSampleValueThread(TriangleD *triangleListD, float *sampleValueD,
									  int activeThreadNum, int m, int f, int c, int n,
									  int orderU, int orderV, int orderW,
									  int ctrlPointNumU, int ctrlPointNumV, int ctrlPointNumW)
{
	int globalIdx = blockDim.x * blockIdx.x + threadIdx.x;
	if (globalIdx >= activeThreadNum)
		return;
	int triangleIdx = globalIdx / m;
	TriangleD &triangle = triangleListD[triangleIdx];

	int localIdx = globalIdx % m;

	float tempFloorFloat = (sqrtf((float)localIdx * 8 + 9) - 3) * 0.5;
	int floor = rintf(tempFloorFloat);
	if ((floor * 2 + 3) * (floor * 2 + 3) != localIdx * 8 + 9)
		floor = ceilf(tempFloorFloat);
	int room = localIdx - ((floor + 1) * floor >> 1);
	float3 barycentric_coord;
	barycentric_coord.x = (float)(n - floor) / n;
	barycentric_coord.y = (float)(floor - room) / n;
	barycentric_coord.z = 1.0f - barycentric_coord.x - barycentric_coord.y;

	float3 v0 = triangle.v[0];
	float3 v1 = triangle.v[1];
	float3 v2 = triangle.v[2];

	// u, v, w 表示经过重心坐标插值之后的采样点的x, y, z分量
	float u = v0.x * barycentric_coord.x + v1.x * barycentric_coord.y + v2.x * barycentric_coord.z;
	float v = v0.y * barycentric_coord.x + v1.y * barycentric_coord.y + v2.y * barycentric_coord.z;
	float w = v0.z * barycentric_coord.x + v1.z * barycentric_coord.y + v2.z * barycentric_coord.z;

	// u, v, w方向节点区间数量
	int knot_interval_count_u = orderU + ctrlPointNumU - (orderU - 1) * 2 - 1;
	int knot_interval_count_v = orderV + ctrlPointNumV - (orderV - 1) * 2 - 1;
	int knot_interval_count_w = orderW + ctrlPointNumW - (orderW - 1) * 2 - 1;

	// 预先将其值设为最大，将末端点归入最后一段 
	int left_idx_u = orderU - 1 + knot_interval_count_u - 1;
	int left_idx_v = orderV - 1 + knot_interval_count_v - 1;
	int left_idx_w = orderW - 1 + knot_interval_count_w - 1;

	// 沿 U 方向查找当前点所在的节点区间 
	for (int ii = orderU - 1; ii <= orderU - 1 + knot_interval_count_u - 1; ++ii)
	{
		if (u >= knotListD[ii] && u < knotListD[ii + 1])
		{
			left_idx_u = ii;
			break;
		}
	}
	// 沿 V 方向查找当前点所在的节点区间 
	for (int jj = orderV - 1; jj <= orderV - 1 + knot_interval_count_v - 1; ++jj)
	{
		if (v >= knotListD[20 + jj] && v < knotListD[20 + jj + 1])
		{
			left_idx_v = jj;
			break;
		}
	}
	// 沿 W 方向查找当前点所在的节点区间 
	for (int kk = orderW - 1; kk <= orderW - 1 + knot_interval_count_w - 1; ++kk)
	{
		if (w >= knotListD[40 + kk] && w < knotListD[40 + kk + 1])
		{
			left_idx_w = kk;
			break;
		}
	}

	float tmpKnot = knotListD[left_idx_u];
	float tmpKnot1 = knotListD[left_idx_u + 1];
	float x_stride = tmpKnot1 - tmpKnot;
	u = (u - tmpKnot) / x_stride;

	tmpKnot = knotListD[20 + left_idx_v];
	tmpKnot1 = knotListD[20 + left_idx_v + 1];
	float y_stride = tmpKnot1 - tmpKnot;
	v = (v - tmpKnot) / y_stride;

	tmpKnot = knotListD[40 + left_idx_w];
	tmpKnot1 = knotListD[40 + left_idx_w + 1];
	float z_stride = tmpKnot1 - tmpKnot;
	w = (w - tmpKnot) / z_stride;

	extern __shared__ float shared_array[];
	// 算出该线程负责的采样点的 B 样条体值
	// fu 表示J_bar矩阵第一列三个元素：偏F_bar_x偏u、偏F_bar_y偏u、偏F_bar_z偏u
	// fv 表示J_bar矩阵第二列三个元素：偏F_bar_x偏v、偏F_bar_y偏v、偏F_bar_z偏v
	float3 result, fu, fv;
	BSplineVolumeValueMatrixD_combine(u, v, w, shared_array,
	left_idx_u - (orderU - 1), left_idx_v - (orderV - 1), left_idx_w - (orderW - 1),
											   orderU, orderV, orderW,
											   result, fu, fv);
	__syncthreads();

	sampleValueD[index2c(localIdx, triangleIdx		  , m + c)] = result.x;
	sampleValueD[index2c(localIdx, triangleIdx + f	  , m + c)] = result.y;
	sampleValueD[index2c(localIdx, triangleIdx + f * 2, m + c)] = result.z;

	///////////////////////////////////////////////////////////////////////////////

	// fw 表示J_bar矩阵第三列三个元素：偏F_bar_x偏w、偏F_bar_y偏w、偏F_bar_z偏w
	float3 fw = BSplineVolumeValueMatrixDw(u, v, w, shared_array,
	left_idx_u - (orderU - 1), left_idx_v - (orderV - 1), left_idx_w - (orderW - 1),
										   orderU, orderV, orderW);
	//__syncthreads();

	v0 = triangle.n[0];
	v1 = triangle.n[1];
	v2 = triangle.n[2];

	// u, v, w 表示经过重心坐标插值之后的法向的x, y, z分量
	u = v0.x * barycentric_coord.x + v1.x * barycentric_coord.y + v2.x * barycentric_coord.z;
	v = v0.y * barycentric_coord.x + v1.y * barycentric_coord.y + v2.y * barycentric_coord.z;
	w = v0.z * barycentric_coord.x + v1.z * barycentric_coord.y + v2.z * barycentric_coord.z;

	float *sampleNormalD = sampleValueD + 3 * f * (m + c);
	// J_bar_star_T_[012]表示J_bar的伴随矩阵的转置(即J_bar*T)的第一行三个元素
	float J_bar_star_T_0 = fv.y * fw.z - fw.y * fv.z;
	float J_bar_star_T_1 = fw.y * fu.z - fu.y * fw.z;
	float J_bar_star_T_2 = fu.y * fv.z - fv.y * fu.z;
	sampleNormalD[index2c(localIdx, triangleIdx, m + c)] = u * J_bar_star_T_0 * x_stride + v * J_bar_star_T_1 * y_stride + w * J_bar_star_T_2 * z_stride;

	// J_bar_star_T_[012]表示J_bar的伴随矩阵的转置(即J_bar*T)的第二行三个元素
	J_bar_star_T_0 = fv.z * fw.x - fw.z * fv.x;
	J_bar_star_T_1 = fw.z * fu.x - fu.z * fw.x;
	J_bar_star_T_2 = fu.z * fv.x - fv.z * fu.x;
	sampleNormalD[index2c(localIdx, triangleIdx + f, m + c)] = u * J_bar_star_T_0 * x_stride + v * J_bar_star_T_1 * y_stride + w * J_bar_star_T_2 * z_stride;

	// J_bar_star_T_[012]表示J_bar的伴随矩阵的转置(即J_bar*T)的第三行三个元素
	J_bar_star_T_0 = fv.x * fw.y - fw.x * fv.y;
	J_bar_star_T_1 = fw.x * fu.y - fu.x * fw.y;
	J_bar_star_T_2 = fu.x * fv.y - fv.x * fu.y;
	sampleNormalD[index2c(localIdx, triangleIdx + f * 2, m + c)] = u * J_bar_star_T_0 * x_stride + v * J_bar_star_T_1 * y_stride + w * J_bar_star_T_2 * z_stride;
}

__global__ void calcConstraitSampleValueThread(TriangleD *triangleListD, float *sampleValueD,
											   int activeThreadNum, int m, int f, int c, int n_,
											   int orderU, int orderV, int orderW,
											   int ctrlPointNumU, int ctrlPointNumV, int ctrlPointNumW)
{
	int globalIdx = blockDim.x * blockIdx.x + threadIdx.x;
	if (globalIdx >= activeThreadNum)
		return;
	int triangleIdx = globalIdx / c;
	TriangleD &triangle = triangleListD[triangleIdx];

	int localIdx = globalIdx % c;

	int floor = -1, room = -1;
	if (localIdx < 2 * n_ - 1)
	{
		floor = (localIdx + 1) / 2;
		if (localIdx % 2 == 1)
			room = 0;
		else
			room = floor;
	}
	else
	{
		floor = n_;
		room = localIdx - (2 * n_ - 1);
	}
	float3 barycentric_coord;
	barycentric_coord.x = (float)(n_ - floor) / n_;
	barycentric_coord.y = (float)(floor - room) / n_;
	barycentric_coord.z = 1.0f - barycentric_coord.x - barycentric_coord.y;

	float3 v0 = triangle.v[0];
	float3 v1 = triangle.v[1];
	float3 v2 = triangle.v[2];

	// u, v, w 表示经过重心坐标插值之后的采样点的x, y, z分量
	float u = v0.x * barycentric_coord.x + v1.x * barycentric_coord.y + v2.x * barycentric_coord.z;
	float v = v0.y * barycentric_coord.x + v1.y * barycentric_coord.y + v2.y * barycentric_coord.z;
	float w = v0.z * barycentric_coord.x + v1.z * barycentric_coord.y + v2.z * barycentric_coord.z;

	// u, v, w方向节点区间数量
	int knot_interval_count_u = orderU + ctrlPointNumU - (orderU - 1) * 2 - 1;
	int knot_interval_count_v = orderV + ctrlPointNumV - (orderV - 1) * 2 - 1;
	int knot_interval_count_w = orderW + ctrlPointNumW - (orderW - 1) * 2 - 1;

	// 预先将其值设为最大，将末端点归入最后一段 
	int left_idx_u = orderU - 1 + knot_interval_count_u - 1;
	int left_idx_v = orderV - 1 + knot_interval_count_v - 1;
	int left_idx_w = orderW - 1 + knot_interval_count_w - 1;

	// 沿 U 方向查找当前点所在的节点区间 
	for (int ii = orderU - 1; ii <= orderU - 1 + knot_interval_count_u - 1; ++ii)
	{
		if (u >= knotListD[ii] && u < knotListD[ii + 1])
		{
			left_idx_u = ii;
			break;
		}
	}
	// 沿 V 方向查找当前点所在的节点区间 
	for (int jj = orderV - 1; jj <= orderV - 1 + knot_interval_count_v - 1; ++jj)
	{
		if (v >= knotListD[20 + jj] && v < knotListD[20 + jj + 1])
		{
			left_idx_v = jj;
			break;
		}
	}
	// 沿 W 方向查找当前点所在的节点区间 
	for (int kk = orderW - 1; kk <= orderW - 1 + knot_interval_count_w - 1; ++kk)
	{
		if (w >= knotListD[40 + kk] && w < knotListD[40 + kk + 1])
		{
			left_idx_w = kk;
			break;
		}
	}

	float tmpKnot = knotListD[left_idx_u];
	float tmpKnot1 = knotListD[left_idx_u + 1];
	float x_stride = tmpKnot1 - tmpKnot;
	u = (u - tmpKnot) / x_stride;

	tmpKnot = knotListD[20 + left_idx_v];
	tmpKnot1 = knotListD[20 + left_idx_v + 1];
	float y_stride = tmpKnot1 - tmpKnot;
	v = (v - tmpKnot) / y_stride;

	tmpKnot = knotListD[40 + left_idx_w];
	tmpKnot1 = knotListD[40 + left_idx_w + 1];
	float z_stride = tmpKnot1 - tmpKnot;
	w = (w - tmpKnot) / z_stride;

	extern __shared__ float shared_array[];
	// 算出该线程负责的采样点的 B 样条体值f
	// fu 表示J_bar矩阵第一列三个元素：偏F_bar_x偏u、偏F_bar_y偏u、偏F_bar_z偏u
	// fv 表示J_bar矩阵第二列三个元素：偏F_bar_x偏v、偏F_bar_y偏v、偏F_bar_z偏v
	float3 result, fu, fv;
	BSplineVolumeValueMatrixD_combine(u, v, w, shared_array,
	left_idx_u - (orderU - 1), left_idx_v - (orderV - 1), left_idx_w - (orderW - 1),
											   orderU, orderV, orderW,
											   result, fu, fv);
	__syncthreads();

	sampleValueD[index2c(localIdx + m, triangleIdx		  , m + c)] = result.x;
	sampleValueD[index2c(localIdx + m, triangleIdx + f	  , m + c)] = result.y;
	sampleValueD[index2c(localIdx + m, triangleIdx + f * 2, m + c)] = result.z;

	////////////////////////////////////////////////////////////////////////////

	// fw 表示J_bar矩阵第三列三个元素：偏F_bar_x偏w、偏F_bar_y偏w、偏F_bar_z偏w
	float3 fw = BSplineVolumeValueMatrixDw(u, v, w, shared_array,
	left_idx_u - (orderU - 1), left_idx_v - (orderV - 1), left_idx_w - (orderW - 1),
										   orderU, orderV, orderW);
	//__syncthreads();

	v0 = triangle.n[0];
	v1 = triangle.n[1];
	v2 = triangle.n[2];

	// u, v, w 表示经过重心坐标插值之后的法向的x, y, z分量
	u = v0.x * barycentric_coord.x + v1.x * barycentric_coord.y + v2.x * barycentric_coord.z;
	v = v0.y * barycentric_coord.x + v1.y * barycentric_coord.y + v2.y * barycentric_coord.z;
	w = v0.z * barycentric_coord.x + v1.z * barycentric_coord.y + v2.z * barycentric_coord.z;

	float *sampleNormalD = sampleValueD + 3 * f * (m + c);
	// J_bar_star_T_[012]表示J_bar的伴随矩阵的转置(即J_bar*T)的第一行三个元素
	float J_bar_star_T_0 = fv.y * fw.z - fw.y * fv.z;
	float J_bar_star_T_1 = fw.y * fu.z - fu.y * fw.z;
	float J_bar_star_T_2 = fu.y * fv.z - fv.y * fu.z;
	sampleNormalD[index2c(localIdx + m, triangleIdx, m + c)] = u * J_bar_star_T_0 * x_stride + v * J_bar_star_T_1 * y_stride + w * J_bar_star_T_2 * z_stride;

	// J_bar_star_T_[012]表示J_bar的伴随矩阵的转置(即J_bar*T)的第二行三个元素
	J_bar_star_T_0 = fv.z * fw.x - fw.z * fv.x;
	J_bar_star_T_1 = fw.z * fu.x - fu.z * fw.x;
	J_bar_star_T_2 = fu.z * fv.x - fv.z * fu.x;
	sampleNormalD[index2c(localIdx + m, triangleIdx + f, m + c)] = u * J_bar_star_T_0 * x_stride + v * J_bar_star_T_1 * y_stride + w * J_bar_star_T_2 * z_stride;

	// J_bar_star_T_[012]表示J_bar的伴随矩阵的转置(即J_bar*T)的第三行三个元素
	J_bar_star_T_0 = fv.x * fw.y - fw.x * fv.y;
	J_bar_star_T_1 = fw.x * fu.y - fu.x * fw.y;
	J_bar_star_T_2 = fu.x * fv.y - fv.x * fu.y;
	sampleNormalD[index2c(localIdx + m, triangleIdx + f * 2, m + c)] = u * J_bar_star_T_0 * x_stride + v * J_bar_star_T_1 * y_stride + w * J_bar_star_T_2 * z_stride;
}

__global__ void calcAdjustNormal(TriangleD *triangleListD, int f,
								 int orderU, int orderV, int orderW,
								 int ctrlPointNumU, int ctrlPointNumV, int ctrlPointNumW)
{
	int globalIdx = blockDim.x * blockIdx.x + threadIdx.x;
	int triangleIdx = globalIdx / 3;
	if (triangleIdx >= f)
		return;
	int i = globalIdx % 3;

	float3 vertex = triangleListD[triangleIdx].v[i];
	float u = vertex.x;
	float v = vertex.y;
	float w = vertex.z;

	// u, v, w方向节点区间数量
	int knot_interval_count_u = orderU + ctrlPointNumU - (orderU - 1) * 2 - 1;
	int knot_interval_count_v = orderV + ctrlPointNumV - (orderV - 1) * 2 - 1;
	int knot_interval_count_w = orderW + ctrlPointNumW - (orderW - 1) * 2 - 1;

	// 预先将其值设为最大，将末端点归入最后一段 
	int left_idx_u = orderU - 1 + knot_interval_count_u - 1;
	int left_idx_v = orderV - 1 + knot_interval_count_v - 1;
	int left_idx_w = orderW - 1 + knot_interval_count_w - 1;

	// 沿 U 方向查找当前点所在的节点区间 
	for (int ii = orderU - 1; ii <= orderU - 1 + knot_interval_count_u - 1; ++ii)
	{
		if (u >= knotListD[ii] && u < knotListD[ii + 1])
		{
			left_idx_u = ii;
			break;
		}
	}
	// 沿 V 方向查找当前点所在的节点区间 
	for (int jj = orderV - 1; jj <= orderV - 1 + knot_interval_count_v - 1; ++jj)
	{
		if (v >= knotListD[20 + jj] && v < knotListD[20 + jj + 1])
		{
			left_idx_v = jj;
			break;
		}
	}
	// 沿 W 方向查找当前点所在的节点区间 
	for (int kk = orderW - 1; kk <= orderW - 1 + knot_interval_count_w - 1; ++kk)
	{
		if (w >= knotListD[40 + kk] && w < knotListD[40 + kk + 1])
		{
			left_idx_w = kk;
			break;
		}
	}
	float tmpKnot = knotListD[left_idx_u];
	float tmpKnot1 = knotListD[left_idx_u + 1];
	float x_stride = tmpKnot1 - tmpKnot;
	u = (u - tmpKnot) / x_stride;

	tmpKnot = knotListD[20 + left_idx_v];
	tmpKnot1 = knotListD[20 + left_idx_v + 1];
	float y_stride = tmpKnot1 - tmpKnot;
	v = (v - tmpKnot) / y_stride;

	tmpKnot = knotListD[40 + left_idx_w];
	tmpKnot1 = knotListD[40 + left_idx_w + 1];
	float z_stride = tmpKnot1 - tmpKnot;
	w = (w - tmpKnot) / z_stride;

	extern __shared__ float shared_array[];

	// fu 表示J_bar矩阵第一列三个元素：偏F_bar_x偏u、偏F_bar_y偏u、偏F_bar_z偏u
	float3 fu = BSplineVolumeValueMatrixDu(u, v, w, shared_array,
	left_idx_u - (orderU - 1), left_idx_v - (orderV - 1), left_idx_w - (orderW - 1),
										   orderU, orderV, orderW);
	__syncthreads();

	// fv 表示J_bar矩阵第二列三个元素：偏F_bar_x偏v、偏F_bar_y偏v、偏F_bar_z偏v
	float3 fv = BSplineVolumeValueMatrixDv(u, v, w, shared_array,
	left_idx_u - (orderU - 1), left_idx_v - (orderV - 1), left_idx_w - (orderW - 1),
										   orderU, orderV, orderW);
	__syncthreads();

	// fw 表示J_bar矩阵第三列三个元素：偏F_bar_x偏w、偏F_bar_y偏w、偏F_bar_z偏w
	float3 fw = BSplineVolumeValueMatrixDw(u, v, w, shared_array,
	left_idx_u - (orderU - 1), left_idx_v - (orderV - 1), left_idx_w - (orderW - 1),
										   orderU, orderV, orderW);
	__syncthreads();

	vertex = triangleListD[triangleIdx].n_adj_origin[i];
	u = vertex.x;
	v = vertex.y;
	w = vertex.z;

	// J_bar_star_T_[012]表示J_bar的伴随矩阵的转置(即J_bar*T)的第一行三个元素
	float J_bar_star_T_0 = fv.y * fw.z - fw.y * fv.z;
	float J_bar_star_T_1 = fw.y * fu.z - fu.y * fw.z;
	float J_bar_star_T_2 = fu.y * fv.z - fv.y * fu.z;
	triangleListD[triangleIdx].n_adj[i].x = u * J_bar_star_T_0 * x_stride + v * J_bar_star_T_1 * y_stride + w * J_bar_star_T_2 * z_stride;

	// J_bar_star_T_[012]表示J_bar的伴随矩阵的转置(即J_bar*T)的第二行三个元素
	J_bar_star_T_0 = fv.z * fw.x - fw.z * fv.x;
	J_bar_star_T_1 = fw.z * fu.x - fu.z * fw.x;
	J_bar_star_T_2 = fu.z * fv.x - fv.z * fu.x;
	triangleListD[triangleIdx].n_adj[i].y = u * J_bar_star_T_0 * x_stride + v * J_bar_star_T_1 * y_stride + w * J_bar_star_T_2 * z_stride;

	// J_bar_star_T_[012]表示J_bar的伴随矩阵的转置(即J_bar*T)的第三行三个元素
	J_bar_star_T_0 = fv.x * fw.y - fw.x * fv.y;
	J_bar_star_T_1 = fw.x * fu.y - fu.x * fw.y;
	J_bar_star_T_2 = fu.x * fv.y - fv.x * fu.y;
	triangleListD[triangleIdx].n_adj[i].z = u * J_bar_star_T_0 * x_stride + v * J_bar_star_T_1 * y_stride + w * J_bar_star_T_2 * z_stride;
}

void calcSampleValue(AlgorithmType algo_type)
{
	if (algo_type == CYM)
	{
		// 计算采样点的值和法向
		//calcSampleValueThread<<<blockNumStep0, blockSizeStep0, sizeof(float) * blockSizeStep0 * 9>>>
		calcSampleValueThread<<<blockNumStep0, blockSizeStep0, sizeof(float) * blockSizeStep0 * 13>>>
											(triangleListD, sampleValueD,
											 activeThreadNumStep0, triangleCtrlPointNum, triangleNum, constrait_point_num,
											 degree, order[U], order[V], order[W],
											 ctrlPointNum[U], ctrlPointNum[V], ctrlPointNum[W]);

		// 计算约束点的值和法向
		//calcConstraitSampleValueThread<<<blockNumStep1, blockSizeStep1, sizeof(float) * blockSizeStep1 * 9>>>
		calcConstraitSampleValueThread<<<blockNumStep1, blockSizeStep1, sizeof(float) * blockSizeStep1 * 13>>>
											(triangleListD, sampleValueD,
											 activeThreadNumStep1, triangleCtrlPointNum, triangleNum, constrait_point_num,
											 degree_lower, order[U], order[V], order[W],
											 ctrlPointNum[U], ctrlPointNum[V], ctrlPointNum[W]);
	}
	else
	{
		calcSampleValueThread_PN<<<blockNumStep0_PN, blockSizeStep0_PN, sizeof(float) * blockSizeStep1 * 13>>>
										(triangleListD, sampleValueD_PN,
										  triangleNum, degree, order[U], order[V], order[W],
										  ctrlPointNum[U], ctrlPointNum[V], ctrlPointNum[W]);
	}
	//cudaError_t error = cudaGetLastError();
	//if (error != 0)
	//{
		//cout << "第0步出错\t";
		//printf("CUDA error: %s\n", cudaGetErrorString(error));
	//}

	//error = cudaGetLastError();
	//if (error != 0)
	//{
		//cout << "第一步出错\t";
		//printf("CUDA error: %s\n", cudaGetErrorString(error));
	//}
	//float3 *test = new float3[triangleNum * 3];
	//for (int i = 0; i < triangleNum * 3; ++i)
		//test[i] = make_float3(1.0f, 2.0f, 3.0f);
	//cudaMemcpy(test, sampleValueD_PN, sizeof(float3) * triangleNum * 3, cudaMemcpyDeviceToHost);
	//for (int i = 0; i < triangleNum; ++i)
	//{
		//cout << test[i * 3].x << ", " << test[i * 3].y << ", " << test[i * 3].z << endl;
		//cout << test[i * 3 + 1].x << ", " << test[i * 3 + 1].y << ", " << test[i * 3 + 1].z << endl;
		//cout << test[i * 3 + 2].x << ", " << test[i * 3 + 2].y << ", " << test[i * 3 + 2].z << endl;
		//cout << "==============" << endl;
	//}
	//delete []test;

	//float *test = new float[(triangleCtrlPointNum + constrait_point_num) * triangleNum * 6];
	//cudaMemcpy(test, sampleValueD, sizeof(float) * (triangleCtrlPointNum + constrait_point_num) * triangleNum * 6, cudaMemcpyDeviceToHost);
	//float *n = test + (triangleCtrlPointNum + constrait_point_num) * triangleNum * 3;
	//for (int i = 0; i < triangleNum; ++i)
	//{
		//cout << "i = " << i << endl;
		//for (int j = 0; j < constrait_point_num; ++j)
		//{
			//cout << "\t" << j << " " << test[i * (triangleCtrlPointNum + constrait_point_num) + triangleCtrlPointNum + j] << " "
				 //<< test[(i + triangleNum) * (triangleCtrlPointNum + constrait_point_num) + triangleCtrlPointNum + j] << " "
				 //<< test[(i + triangleNum * 2) * (triangleCtrlPointNum + constrait_point_num) + triangleCtrlPointNum + j] << "\t";
			//cout << n[i * (triangleCtrlPointNum + constrait_point_num) + triangleCtrlPointNum + j] << " "
				 //<< n[(i + triangleNum) * (triangleCtrlPointNum + constrait_point_num) + triangleCtrlPointNum + j] << " "
				 //<< n[(i + triangleNum * 2) * (triangleCtrlPointNum + constrait_point_num) + triangleCtrlPointNum + j] << endl;
		//}
	//}
	//delete []test;
}

#ifdef TRUTH
__global__ void calcSampleValueThread_truth(TriangleD *triangleListD, float *sampleValueD_truth,
									  int activeThreadNum, int m, int f, int n,
									  int orderU, int orderV, int orderW,
									  int ctrlPointNumU, int ctrlPointNumV, int ctrlPointNumW)
{
	int globalIdx = blockDim.x * blockIdx.x + threadIdx.x;
	if (globalIdx >= activeThreadNum)
		return;

	int triangleIdx = globalIdx / m;

	int localIdx = globalIdx % m;

	float tempFloorFloat = (sqrtf((float)localIdx * 8 + 9) - 3) / 2;
	int floor = rintf(tempFloorFloat);
	if ((floor * 2 + 3) * (floor * 2 + 3) != localIdx * 8 + 9)
		floor = ceilf(tempFloorFloat);
	int room = localIdx - (floor + 1) * floor / 2;
	float3 barycentric_coord;
	barycentric_coord.x = (float)(n - floor) / n;
	barycentric_coord.y = (float)(floor - room) / n;
	barycentric_coord.z = 1.0f - barycentric_coord.x - barycentric_coord.y;

	TriangleD &triangle = triangleListD[triangleIdx];
	float3 v0 = triangle.v[0];
	float3 v1 = triangle.v[1];
	float3 v2 = triangle.v[2];

	float u = v0.x * barycentric_coord.x + v1.x * barycentric_coord.y + v2.x * barycentric_coord.z;
	float v = v0.y * barycentric_coord.x + v1.y * barycentric_coord.y + v2.y * barycentric_coord.z;
	float w = v0.z * barycentric_coord.x + v1.z * barycentric_coord.y + v2.z * barycentric_coord.z;

	int i = (u - knotListD[0]) / (knotListD[orderU] - knotListD[0]);
	int j = (v - knotListD[20 + 0]) / (knotListD[20 + orderV] - knotListD[20 + 0]);
	int k = (w - knotListD[40 + 0]) / (knotListD[40 + orderW] - knotListD[40 + 0]);
	if (i >= ctrlPointNumU + orderU - 2 * (orderU - 1) - 1)
		--i;
	if (j >= ctrlPointNumV + orderV - 2 * (orderV - 1) - 1)
		--j;
	if (k >= ctrlPointNumW + orderW - 2 * (orderW - 1) - 1)
		--k;

	/* 确定此 block 需要的 u、v、w 三个方向的 B 样条矩阵 */
	float *Mu = matrixCase(matrix_b_spline_d, orderU, ctrlPointNumU, i + orderU - 1);
	float *Mv = matrixCase(matrix_b_spline_d, orderV, ctrlPointNumV, j + orderV - 1);
	float *Mw = matrixCase(matrix_b_spline_d, orderW, ctrlPointNumW, k + orderW - 1);

	float tmpKnot = knotListD[i + orderU - 1];
	float tmpKnot1 = knotListD[i + orderU];
	u = (u - tmpKnot) / (tmpKnot1 - tmpKnot);

	tmpKnot = knotListD[20 + j + orderV - 1];
	tmpKnot1 = knotListD[20 + j + orderV];
	v = (v - tmpKnot) / (tmpKnot1 - tmpKnot);

	tmpKnot = knotListD[40 + k + orderW - 1];
	tmpKnot1 = knotListD[40 + k + orderW];
	w = (w - tmpKnot) / (tmpKnot1 - tmpKnot);

	extern __shared__ float shared_array[];
	/* 算出该线程负责的采样点的 B 样条体值 */
	float3 result = BSplineVolumeValueMatrixD2(Mu, Mv, Mw,
											   u, v, w, shared_array,
											   i + orderU - 1, j + orderV - 1, k + orderW - 1,
											   orderU, orderV, orderW);

	sampleValueD_truth[index2c(localIdx, triangleIdx, m)] = result.x;
	sampleValueD_truth[index2c(localIdx, triangleIdx + f, m)] = result.y;
	sampleValueD_truth[index2c(localIdx, triangleIdx + f * 2, m)] = result.z;
}

void calcSampleValue_truth()
{
	calcSampleValueThread_truth<<<blockNumStep0_truth, blockSizeStep0, sizeof(float) * blockSizeStep0 * 11>>>
										(triangleListD, sampleValueD_truth,
										 activeThreadNumStep0_truth, triangleCtrlPointNum, triangleNum,
										 degree, order[U], order[V], order[W],
										 ctrlPointNum[U], ctrlPointNum[V], ctrlPointNum[W]);
}
#endif

/************************************************************************************************************/
#define NEW_MOVE			// 定义这个表示使用最终的move函数，否则使用最原始的PN算法
#ifdef NEW_MOVE
// 调整拟合出来的控制顶点
__global__ void move(TriangleD *triangleListD, float *triangleCtrlPointD, int *triangle_adjacent_tableD,
					 int m_, int f, float center_factor, bool use_pn)
{
	int triangleIdx = blockIdx.x * blockDim.x + threadIdx.x;
	if (triangleIdx >= f)
		return;

	int adj_face_idx[3];
	adj_face_idx[0] = triangle_adjacent_tableD[triangleIdx * 3];
	adj_face_idx[1] = triangle_adjacent_tableD[triangleIdx * 3 + 1];
	adj_face_idx[2] = triangle_adjacent_tableD[triangleIdx * 3 + 2];

	//int adj_edge_idx[3] = { -1, -1, -1 };		// 实际上应该使用这一句，以此判断有没有相邻三角形
	int adj_edge_idx[3] = { 0, 0, 0 };			// 但是对于某些模型使用上一句会出现内存越界，所以暂且使用这一句，权宜之计
	//bool handle[3] = { false, false, false };	// 在该边有一个以上法向时，是否处理这条边。当这条边有相邻面片才需要处理
	for (int i = 0; i < 3; ++i)
		if (adj_face_idx[i] >= 0)
		{
			adj_edge_idx[i] = adj_face_idx[i] & 0x3;
			adj_face_idx[i] = adj_face_idx[i] >> 2;
			//handle[i] = true;
		}

	int n_count[3];
	n_count[0] = triangleListD[triangleIdx].nc[0];
	n_count[1] = triangleListD[triangleIdx].nc[1];
	n_count[2] = triangleListD[triangleIdx].nc[2];
	//printf("ncount = (%d, %d, %d): triangleIdx = %d\n", n_count[0], n_count[1], n_count[2], triangleIdx);

	float *p_x  = &triangleCtrlPointD[m_ * triangleIdx];
	float *p_y  = &triangleCtrlPointD[m_ * (f + triangleIdx)];
	float *p_z  = &triangleCtrlPointD[m_ * (f * 2 + triangleIdx)];

	int edge_ctrlpoint_idx[6] = { 5, 2,			1, 3,		7, 8 };				// 依次处理的边控制顶点的序号(0,1号属于0号边， 2,3号属于1号边， 4,5号属于2号边)
	int corner_ctrlpoint_idx[6] = { 9, 0,		0, 6,		6, 9 };				// 上面每个边控制顶点对应的角控制顶点序号(0,1号属于0号边， 2,3号属于1号边， 4,5号属于2号边)
	int oppo_corner_ctrlpoint_idx[6] = { 0, 9,	6, 0,		9, 6 };				// 上面每个角控制顶点所在边的另一个角控制顶点编号(0,1号属于0号边， 2,3号属于1号边， 4,5号属于2号边)
	//int adj_corner_ctrlpoint_idx[3][2] = { { 0, 9 }, { 6, 0 }, { 9, 6 } };		// 相邻三角形0, 1, 2号边上的控制顶点编号（仅有角点，没有边点）
	//int adj_edge_ctrlpoint_idx[3][2] = { { 2, 5 }, { 3, 1 }, { 8, 7 } };		// 相邻三角形0, 1, 2号边上的控制顶点编号（仅有边点，没有角点）
	int adjust_normal_idx[6] = { 2, 0,		0, 1,		1, 2 };
	int adj_corner_ctrlpoint_idx[3][2] = { { 0, 2 }, { 1, 0 }, { 2, 1 } };		// 相邻三角形0, 1, 2号边上的控制顶点编号（仅有角点，没有边点）

	//const float ZERO = 10e-6;
	float3 delta = make_float3(0.0f, 0.0f, 0.0f), sum = make_float3(0.0f, 0.0f, 0.0f);
	// 六个边点，按0 1 2号边的顺序处理，即5, 2, 1, 3, 7, 8号控制顶点
	//printf("for开始, triangleIdx = %d\n", triangleIdx);
	for (int i = 0; i < 6; ++i)
	{
		float3 v_ctrlpoint_corner = make_float3(*(p_x + corner_ctrlpoint_idx[i]), *(p_y + corner_ctrlpoint_idx[i]), *(p_z + corner_ctrlpoint_idx[i]));
		float3 v_ctrlpoint_corner_oppo = make_float3(*(p_x + oppo_corner_ctrlpoint_idx[i]), *(p_y + oppo_corner_ctrlpoint_idx[i]), *(p_z + oppo_corner_ctrlpoint_idx[i]));
		float3 v01 = v_ctrlpoint_corner_oppo - v_ctrlpoint_corner;
		float3 v_mid = 0.5 * (v_ctrlpoint_corner + v_ctrlpoint_corner_oppo);
		//float3 n_ctrlpoint_corner = make_float3(*(pn_x + corner_ctrlpoint_idx[i]), *(pn_y + corner_ctrlpoint_idx[i]), *(pn_z + corner_ctrlpoint_idx[i]));
		float3 n_ctrlpoint_corner = triangleListD[triangleIdx].n_adj[adjust_normal_idx[i]];
		normalize(n_ctrlpoint_corner);

		// p 是要处理的边控制顶点
		float3 p = make_float3(*(p_x + edge_ctrlpoint_idx[i]), *(p_y + edge_ctrlpoint_idx[i]), *(p_z + edge_ctrlpoint_idx[i]));
		if (n_count[i / 2] < 2)		// 该条边只有一个法向
		{
			//if (adj_face_idx[i / 2] >= 0)		// 只有当这条边的另一侧有面片时才会处理
			//{
				float3 result = p - ((p - v_ctrlpoint_corner) * n_ctrlpoint_corner) * n_ctrlpoint_corner;

#ifdef RE_LENGTH
				float len0 = length(result);
				float3 result_vector = result - v_ctrlpoint_corner;
				float l_origin = length(p - v_ctrlpoint_corner);
				float l_current = length(result_vector);
				result_vector *= l_origin / l_current;
				result = v_ctrlpoint_corner + result_vector;
				float len1 = length(result);
				printf("delta_leng_1_normal = %f\n", len1 - len0);
#endif

				delta += (result - p);

				*(p_x + edge_ctrlpoint_idx[i]) = result.x;
				*(p_y + edge_ctrlpoint_idx[i]) = result.y;
				*(p_z + edge_ctrlpoint_idx[i]) = result.z;

				sum += result;
			//}
			//printf("only one : result_%d = (%f, %f, %f)\n", edge_ctrlpoint_idx[i], result.x, result.y, result.z);
		}
		//else if (handle[i / 2])		// 该条边有一个以上法向，且需要处理
		else							// 该条边有一个以上法向
		{
			float3 n1 = triangleListD[adj_face_idx[i / 2]].n_adj[adj_corner_ctrlpoint_idx[adj_edge_idx[i / 2]][i % 2]];
			//printf("else开始, triangleIdx = %d, adj_face = %d, cp = %d, n1 = (%f, %f, %f)\n", triangleIdx, adj_face_idx[i / 2], edge_ctrlpoint_idx[i], n1.x, n1.y, n1.z);
			normalize(n1);
			//if (use_pn)
			//{
				float3 n_ave = cross(n_ctrlpoint_corner, n1);
				normalize(n_ave);
				//printf("t = %d, n_cross = %f, %f, %f\n", triangleIdx, n_ave.x, n_ave.y, n_ave.z);
				//float3 result = v_ctrlpoint_corner + v01 * n_ave * 0.333333 * n_ave;				// 原始的pn尖锐边算法，将1/3点往法向上投影，效果不佳
				float3 result = v_ctrlpoint_corner + ((p - v_ctrlpoint_corner) * n_ave) * n_ave;	// 由我的算法改良而来，将差的控制顶点往法向上投影，效果很好
#ifdef RE_LENGTH
				float len0 = length(result);
				float3 result_vector = result - v_ctrlpoint_corner;
				float l_origin = length(p - v_ctrlpoint_corner);
				float l_current = length(result_vector);
				result_vector *= l_origin / l_current;
				result = v_ctrlpoint_corner + result_vector;
				float len1 = length(result);
				printf("delta_leng_pn = %f\n", len1 - len0);
#endif
				delta += (result - p);
				*(p_x + edge_ctrlpoint_idx[i]) = result.x;
				*(p_y + edge_ctrlpoint_idx[i]) = result.y;
				*(p_z + edge_ctrlpoint_idx[i]) = result.z;
				sum += result;
				//printf("2 : result_%d = (%f, %f, %f)\n", edge_ctrlpoint_idx[i], result.x, result.y, result.z);
			//}
			//else
			//{
				//float t0 = 1.2345f, t1 = 2.3456f;
				//float3 center0, center1;
				//bool t0_exist = false, t1_exist = false;
				//if (fabs(n_ctrlpoint_corner * v01) > ZERO)
				//{
					//t0 = (v_mid - v_ctrlpoint_corner) * v01 / (n_ctrlpoint_corner * v01);
					//center0 = v_ctrlpoint_corner + t0 * n_ctrlpoint_corner;
					//t0_exist = true;
					////if (triangleIdx == 10 && i == 5)
					////{
						////printf("n0 = (%f, %f, %f), t0 = %f, center0 = (%f, %f, %f)\n",
								////n_ctrlpoint_corner.x, n_ctrlpoint_corner.y, n_ctrlpoint_corner.z, t0, center0.x, center0.y, center0.z);
					////}
				//}
				//if (fabs(n1 * v01) > ZERO)
				//{
					//t1 = (v_mid - v_ctrlpoint_corner) * v01 / (n1 * v01);
					//center1 = v_ctrlpoint_corner + t1 * n1;
					//t1_exist = true;
					////if (triangleIdx == 10 && i == 5)
					////{
						////printf("n1 = (%f, %f, %f), t1 = %f, center1 = (%f, %f, %f)\n",
								////n1.x, n1.y, n1.z, t1, center1.x, center1.y, center1.z);
					////}
				//}

				////printf("t0 = %f, t1 = %f, triangleIdx = %d, cp = %d\n", t0, t1, triangleIdx, edge_ctrlpoint_idx[i]);
				//float3 center_mid;
				//if (t0_exist && t1_exist)	// 当前三角形和相邻三角形都不精确
				//{
					//float3 center_delta = center0 - center1;
					//float t = (v_ctrlpoint_corner - center0) * center_delta / (center_delta * center_delta);
					//center_mid = center0 + t * center_delta;


					//float3 rad0 = v_ctrlpoint_corner - center0;
					//float r0 = sqrt(rad0.x * rad0.x + rad0.y * rad0.y + rad0.z * rad0.z);
					//float3 rad1 = v_ctrlpoint_corner - center1;
					//float r1 = sqrt(rad1.x * rad1.x + rad1.y * rad1.y + rad1.z * rad1.z);
					////printf("都不精确, 三角形=%d, cp=%d, t = %f, r0 = %f, r1 = %f\n", triangleIdx, edge_ctrlpoint_idx[i], t, r0, r1);
				//}
				//else if (t0_exist)			// 当前三角形不精确，相邻三角形精确
				//{
					//float t = (v_ctrlpoint_corner - center0) * n1 / (n1 * n1);
					//center_mid = center0 + t * n1;
					////printf("当前三角形不精确，相邻三角形精确, 三角形=%d, cp=%d, n1 = (%f, %f, %f), t = %f\n", triangleIdx, edge_ctrlpoint_idx[i], n1.x, n1.y, n1.z, t);
				//}
				//else if (t1_exist)			// 当前三角形精确，相邻三角形不精确
				//{
					//float t = (v_ctrlpoint_corner - center1) * n_ctrlpoint_corner / (n_ctrlpoint_corner * n_ctrlpoint_corner);
					//center_mid = center1 + t * n_ctrlpoint_corner;
					////printf("当前三角形精确，相邻三角形不精确, 三角形=%d, cp=%d, t = %f\n", triangleIdx, edge_ctrlpoint_idx[i], t);
				//}
				//else						// 当前三角形和相邻三角形都精确
				//{
					////printf("两个都精确, 三角形=%d, cp=%d\n", triangleIdx, edge_ctrlpoint_idx[i]);
					//continue;
				//}
				//float3 n_ave = v_ctrlpoint_corner - center_mid;
				//normalize(n_ave);
				//float3 result = p - ((p - v_ctrlpoint_corner) * n_ave) * n_ave;
				////printf("t = %d, n_ave = %f, %f, %f\tp = %f, %f, %f\t, result=%f, %f, %f\n", triangleIdx, n_ave.x, n_ave.y, n_ave.z, p.x, p.y, p.z, result.x, result.y, result.z);
//#ifdef RE_LENGTH
				//float len0 = length(result);
				//float3 result_vector = result - v_ctrlpoint_corner;
				//float l_origin = length(p - v_ctrlpoint_corner);
				//float l_current = length(result_vector);
				//result_vector *= l_origin / l_current;
				//result = v_ctrlpoint_corner + result_vector;
				//float len1 = length(result);
				//printf("delta_leng_my = %f\n", len1 - len0);
//#endif
				//delta += (result - p);
				//*(p_x + edge_ctrlpoint_idx[i]) = result.x;
				//*(p_y + edge_ctrlpoint_idx[i]) = result.y;
				//*(p_z + edge_ctrlpoint_idx[i]) = result.z;
				//sum += result;


				//float3 n_pn = cross(n_ctrlpoint_corner, n1);
				//normalize(n_pn);
				//float3 result_pn = v_ctrlpoint_corner + ((p - v_ctrlpoint_corner) * n_pn) * n_pn;

				//float3 del = result_pn - result;
				//float dot = n_ave * n_pn;
				////printf("del = %f, %f, %f\t\tdot = %f\n", del.x, del.y, del.z, dot);



				////printf("2 : result_%d = (%f, %f, %f)\n", edge_ctrlpoint_idx[i], result.x, result.y, result.z);
			//}
		}
	}

	// 中间控制顶点，即4号控制顶点
#ifdef LESS_THAN_2
	if (n_count[0] < 2 && n_count[1] < 2 && n_count[2] < 2)
#endif
	{
		float3 p = make_float3(*(p_x + 4), *(p_y + 4), *(p_z + 4));

		/******** 平均顶点位置，PN-Triangle方法 *********/
		//sum *= 1.0 / 6;
		//float3 result = sum + (sum - p) * 0.5;

		/******** 平均delta *********/
		delta *= center_factor / 6;
		//delta *= 1.5 / 6;
		float3 result = p + delta;

		/******** 写结果 *********/
		*(p_x + 4) = result.x;
		*(p_y + 4) = result.y;
		*(p_z + 4) = result.z;
		//printf("result_4 = (%f, %f, %f)\n", result.x, result.y, result.z);
	}
}
#else
__global__ void	move(TriangleD *triangleListD, float *triangleCtrlPointD, float *triangleNormalCtrlPointD, int m_, int f)
{
	int triangleIdx = blockIdx.x * blockDim.x + threadIdx.x;
	if (triangleIdx >= f)
		return;

	float *p_x = &triangleCtrlPointD[m_ * triangleIdx];
	float *p_y = &triangleCtrlPointD[m_ * (f + triangleIdx)];
	float *p_z = &triangleCtrlPointD[m_ * (f * 2 + triangleIdx)];
	float *pn_x = &triangleNormalCtrlPointD[m_ * triangleIdx];
	float *pn_y = &triangleNormalCtrlPointD[m_ * (f + triangleIdx)];
	float *pn_z = &triangleNormalCtrlPointD[m_ * (f * 2 + triangleIdx)];

	/******* 点1 *******/
/*#define MOVE1*/ 	// MOVE1被定义表示永远使用原始三角形的信息进行调整，理论上是错误的，只有初始情况下正确，仅供调试时用
#ifdef MOVE1
	float3 v = triangleListD[triangleIdx].v[0];
	float3 n = triangleListD[triangleIdx].n[0];
#else
	float3 v = make_float3(*p_x, *p_y, *p_z);
	float3 n = make_float3(*pn_x, *pn_y, *pn_z);
	float length = sqrt(n.x * n.x + n.y * n.y + n.z * n.z);
	normalize(n);
#endif
	float3 p = make_float3(*(p_x + 1), *(p_y + 1), *(p_z + 1));
	float3 result = p - ((p - v) * n) * n;
	float3 delta = result - p;
	/*if (threadIdx.x == 0)*/
	/*{*/
		/*printf("triangleIdx = %d\n", triangleIdx);*/
		/*printf("待投点 = (%f, %f, %f), 法向 = (%f, %f, %f), 角点 = (%f, %f, %f),\n结果 = (%f, %f, %f), 差值 = (%f, %f, %f)\n",*/
				/*p.x, p.y, p.z,			n.x, n.y, n.z,		v.x, v.y, v.z,		result.x, result.y, result.z,		delta.x, delta.y, delta.z);*/
	/*}*/
	*(p_x + 1) = result.x;
	*(p_y + 1) = result.y;
	*(p_z + 1) = result.z;
	float3 sum = result;

	// 点2
	p = make_float3(*(p_x + 2), *(p_y + 2), *(p_z + 2));
	result = p - ((p - v) * n) * n;
	*(p_x + 2) = result.x;
	*(p_y + 2) = result.y;
	*(p_z + 2) = result.z;
	sum += result;

	/******* 点3 *******/
#ifdef MOVE1
	v = triangleListD[triangleIdx].v[1];
	n = triangleListD[triangleIdx].n[1];
#else
	v = make_float3(*(p_x + 6), *(p_y + 6), *(p_z + 6));
	n = make_float3(*(pn_x + 6), *(pn_y + 6), *(pn_z + 6));
	normalize(n);
#endif
	p = make_float3(*(p_x + 3), *(p_y + 3), *(p_z + 3));
	result = p - ((p - v) * n) * n;
	*(p_x + 3) = result.x;
	*(p_y + 3) = result.y;
	*(p_z + 3) = result.z;
	sum += result;

	// 点7
	p = make_float3(*(p_x + 7), *(p_y + 7), *(p_z + 7));
	result = p - ((p - v) * n) * n;
	*(p_x + 7) = result.x;
	*(p_y + 7) = result.y;
	*(p_z + 7) = result.z;
	sum += result;

	/******* 点8 *******/
#ifdef MOVE1
	v = triangleListD[triangleIdx].v[2];
	n = triangleListD[triangleIdx].n[2];
#else
	v = make_float3(*(p_x + 9), *(p_y + 9), *(p_z + 9));
	n = make_float3(*(pn_x + 9), *(pn_y + 9), *(pn_z + 9));
	normalize(n);
#endif
	p = make_float3(*(p_x + 8), *(p_y + 8), *(p_z + 8));
	result = p - ((p - v) * n) * n;
	*(p_x + 8) = result.x;
	*(p_y + 8) = result.y;
	*(p_z + 8) = result.z;
	sum += result;

	// 点5
	p = make_float3(*(p_x + 5), *(p_y + 5), *(p_z + 5));
	result = p - ((p - v) * n) * n;
	*(p_x + 5) = result.x;
	*(p_y + 5) = result.y;
	*(p_z + 5) = result.z;
	sum += result;

	/******* 点4 *******/
	p = make_float3(*(p_x + 4), *(p_y + 4), *(p_z + 4));
	sum *= 1.0 / 6;
	result = sum + (sum - p) * 0.5;
	*(p_x + 4) = result.x;
	*(p_y + 4) = result.y;
	*(p_z + 4) = result.z;
}
#endif

float center_factor = 1.5f;

void calcTriangleCtrlPoint(bool adjust_silhouette, bool use_pn, AlgorithmType algo_type)
{
	if (algo_type == CYM)
	{
		float alpha = 1.0f, beta = 0.0f;
		/* 计算面片和法向的控制顶点*/
		cublasStatus_t stat = cublasSgemm(cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N,
										  triangleCtrlPointNum_lower, triangleNum * 6, triangleCtrlPointNum + constrait_point_num,
										  &alpha,
										  matrixFittingD + matrixStartIdxFitting, triangleCtrlPointNum_lower,
										  sampleValueD, triangleCtrlPointNum + constrait_point_num,
										  &beta,
										  triangleCtrlPointD, triangleCtrlPointNum_lower);
		if (stat != CUBLAS_STATUS_SUCCESS)
		{
			cout << "triangleCtrlPointD fail!!!!!!!!!!!!!\tstat = " << stat << endl;
			printCudaError(__FILE__, __FUNCTION__, __LINE__);
			return;
		}
	}

	// 计算每个三角片的三个用于调整控制顶点的法向
	calcAdjustNormal<<<blockNumAdjNormal, blockSizeAdjNormal, sizeof(float) * blockSizeAdjNormal * 8>>>
										(triangleListD, triangleNum,
										 order[U], order[V], order[W],
										 ctrlPointNum[U], ctrlPointNum[V], ctrlPointNum[W]);

	if (algo_type == CYM)
	{
		// 调整侧影轮廓线
		if (adjust_silhouette)
		{
			const int move_block_size = 256;
			int move_block_num = ceil(static_cast<double>(triangleNum) / move_block_size);
			//cout << "move 开始" << endl;
#ifdef NEW_MOVE
			move<<<move_block_num, move_block_size>>>(triangleListD, triangleCtrlPointD, triangle_adjacent_tableD,
					triangleCtrlPointNum_lower, triangleNum, center_factor, use_pn);
#else
			move<<<move_block_num, move_block_size>>>(triangleListD, triangleCtrlPointD, triangleCtrlPointD + 3 * triangleNum * triangleCtrlPointNum_lower, triangleCtrlPointNum_lower, triangleNum);
#endif

//cudaThreadSynchronize();

			//move<<<move_block_num, move_block_size>>>(triangleListD, triangleCtrlPointD, triangleCtrlPointD + 3 * triangleNum * triangleCtrlPointNum_lower, triangle_adjacent_tableD,
					//triangleCtrlPointNum_lower, triangleNum, center_factor);
#ifndef MORPH
			cout << "center_factor = " << center_factor << endl;
#endif
			printCudaError(__FILE__, __FUNCTION__, __LINE__);
		}
#ifdef DRAW_TRIANGULAR_CTRL_POINTS
		// 将计算好的控制顶点传回内存，仅用于调试，在最终结果上显示控制顶点，测效率时需删除
		cudaMemcpy(triangular_ctrl_points, triangleCtrlPointD, sizeof(float) * 3 * triangleNum * triangleCtrlPointNum_lower, cudaMemcpyDeviceToHost);
#endif
	}
	else
	{
		int blockNum = ceil(static_cast<double>(triangleNum) / 128);
		calcCtrlPoint_PN<<<blockNum, 128>>>(triangleListD, triangle_adjacent_tableD, sampleValueD_PN, triangleCtrlPointD_PN, triangleNormalCtrlPointD_PN, triangleNum, triangleCtrlPointNum_lower);
#ifdef DRAW_TRIANGULAR_CTRL_POINTS
		// 将计算好的控制顶点传回内存，仅用于调试，在最终结果上显示控制顶点，测效率时需删除
		cudaMemcpy(triangular_ctrl_points, triangleCtrlPointD_PN, sizeof(float) * 3 * triangleNum * triangleCtrlPointNum_lower, cudaMemcpyDeviceToHost);
#endif
	}


#ifndef MORPH
	cout << "triangleNum = " << triangleNum << endl;
#endif
	//float *test = new float[6 * triangleNum * 3];
	//cudaMemcpy(test, triangleNormalCtrlPointD_PN, sizeof(float) * 6 * triangleNum * 3, cudaMemcpyDeviceToHost);
	//for (int i = 0; i < triangleNum; ++i)
	//{
		//for (int j = 0; j < 6; ++j)
		//{
			////cout << i * 10 + j << ", " << (i + triangleNum) * 10 + j << ", " << (i + triangleNum * 2) * 10 + j << endl;
			//cout
				//<< test[i * 6 + j] << ", "
				//<< test[(i + triangleNum) * 6 + j] << ", "
				//<< test[(i + triangleNum * 2) * 6 + j] << endl;
		//}
		//cout << "================" << endl;
	//}

	//float *test = new float[triangleCtrlPointNum_lower * triangleNum * 6];
	//cudaMemcpy(test, triangleCtrlPointD, sizeof(float) * triangleCtrlPointNum_lower * triangleNum * 6, cudaMemcpyDeviceToHost);
	//float *v = test, *n = test + triangleCtrlPointNum_lower * triangleNum * 3;
	//for (int i = 0; i < triangleNum; ++i)
	////for (int i = 24; i < 25; ++i)
	//{
		//cout << "i = " << i << endl;
		//for (int j = 0; j < triangleCtrlPointNum_lower; ++j)
		//{
			////if (j != 0 && j != 6) continue;
			////cout << i * 10 + j << ", " << (i + triangleNum) * 10 + j << ", " << (i + triangleNum * 2) * 10 + j << endl;
			//cout << "\t" << j << ": " << v[i * triangleCtrlPointNum_lower + j] << ", "
				 //<< v[(i + triangleNum) * triangleCtrlPointNum_lower + j] << ", "
				 //<< v[(i + triangleNum * 2) * triangleCtrlPointNum_lower + j] << "\t";
			//double x = n[i * triangleCtrlPointNum_lower + j];
			//double y = n[(i + triangleNum) * triangleCtrlPointNum_lower + j];
			//double z = n[(i + triangleNum * 2) * triangleCtrlPointNum_lower + j];
			//double length = sqrt(x * x + y * y + z * z);
			//cout << "\t" << x / length << ", " << y / length << ", " << z / length << endl;
		//}
		//cout << "================" << endl;
	//}
	//delete []test;
}

#ifdef TRUTH
void matrixMul1_truth()
{
	float alpha = 1.0f, beta = 0.0f;
	cublasStatus_t stat = cublasSgemm(cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N,
									  samplePointPerTriangle * 3, triangleCtrlPointNum, triangleCtrlPointNum, 
									  &alpha,
									  BqD_truth, samplePointPerTriangle * 3,
									  B_1D_truth, triangleCtrlPointNum,
									  &beta,
									  BBD_truth, samplePointPerTriangle * 3);
	if (stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CtrlPoint_truth fail!!!!!!!!!!!!!\tstat = " << stat << endl;
		printCudaError(__FILE__, __FUNCTION__, __LINE__);
		return;
	}
}
#endif

/************************************************************************************************************/

__global__ void copy(float *RD, int u_seg, int v_seg,
					 int activeThreadNumCopy, bool firstLoad, float maxX, float maxY, float maxZ,
					 TriangleD *triangleListD, int segmentPerEdge, int f, int q,
					 float *normalPtrVBO, float *texCoordPtrVBO, float *texCoord3DPtrVBO, float *vertexPtrVBO)
{
	int globalIdx = blockDim.x * blockIdx.x + threadIdx.x;
	if (globalIdx >= activeThreadNumCopy)
		return;

	int triangleIdx = globalIdx / q;
	int localIdx = globalIdx % q;

	vertexPtrVBO[globalIdx * 3 + 0] = RD[triangleIdx * q + localIdx];
	vertexPtrVBO[globalIdx * 3 + 1] = RD[(triangleIdx + f) * q + localIdx];
	vertexPtrVBO[globalIdx * 3 + 2] = RD[(triangleIdx + f * 2) * q + localIdx];

	float *ND = RD + 3 * f * q;
	normalPtrVBO[globalIdx * 3 + 0] = ND[triangleIdx * + q + localIdx];
	normalPtrVBO[globalIdx * 3 + 1] = ND[(triangleIdx + f) * + q + localIdx];
	normalPtrVBO[globalIdx * 3 + 2] = ND[(triangleIdx + f * 2) * + q + localIdx];

	if (firstLoad)
	{
		// 计算纹理坐标
		float tempFloorFloat = (sqrtf((float)(localIdx) * 8 + 9) - 3) / 2;
		int floor = rintf(tempFloorFloat);
		if ((floor * 2 + 3) * (floor * 2 + 3) != localIdx * 8 + 9)
			floor = ceilf(tempFloorFloat);
		int room = localIdx - (floor + 1) * floor / 2;

		float3 barycentric_coord;
		float3 vt0 = triangleListD[triangleIdx].vt[0];
		float3 vt1 = triangleListD[triangleIdx].vt[1];
		float3 vt2 = triangleListD[triangleIdx].vt[2];
		float u, v, w;

		barycentric_coord.x = (float)(segmentPerEdge - floor) / segmentPerEdge;
		barycentric_coord.y = (float)(floor - room) / segmentPerEdge;
		barycentric_coord.z = 1.0f - barycentric_coord.x - barycentric_coord.y;
		u = vt0.x * barycentric_coord.x + vt1.x * barycentric_coord.y + vt2.x * barycentric_coord.z;
		v = vt0.y * barycentric_coord.x + vt1.y * barycentric_coord.y + vt2.y * barycentric_coord.z;
		//w = vt0.z * barycentric_coord.x + vt1.z * barycentric_coord.y + vt2.z * barycentric_coord.z;
		w = vt0.z;

		// 这是对前8个原始三角形面片的特殊处理，因为它们的第一个顶点映射到多个纹理坐标
		if ((triangleIdx < (u_seg * v_seg * 2 - v_seg) * 8) &&		// 首先保证属于前8个原始面片
			(triangleIdx / v_seg % (u_seg * 2 - 1) == 0))			// 子面片属于原始面片产生的最开始的v_seg个面片
		{
			if (floor == 0)
				u = v = 0;
			else
			{
				float v0 = vt2.y;
				float v1 = vt1.y;
				float k = (v1 - v0) / u;
				float b = v0 * (1.0 - k);
				v = k * v + b;
			}
		}

		// 存储二维纹理坐标
		texCoordPtrVBO[globalIdx * 3 + 0] = u;
		texCoordPtrVBO[globalIdx * 3 + 1] = v;
		texCoordPtrVBO[globalIdx * 3 + 2] = w;

		// 存储三维纹理坐标
		//float minMax = maxX;
		//if (minMax > maxY)
			//minMax = maxY;
		//if (minMax > maxZ)
			//minMax = maxZ;
		////texCoord3DPtrVBO[globalIdx * 3 + 0] = vertexPtrVBO[globalIdx * 3 + 0] / maxX;
		////texCoord3DPtrVBO[globalIdx * 3 + 1] = vertexPtrVBO[globalIdx * 3 + 1] / maxY;
		////texCoord3DPtrVBO[globalIdx * 3 + 2] = vertexPtrVBO[globalIdx * 3 + 2] / maxZ;
		//texCoord3DPtrVBO[globalIdx * 3 + 0] = vertexPtrVBO[globalIdx * 3 + 0] / minMax;
		//texCoord3DPtrVBO[globalIdx * 3 + 1] = vertexPtrVBO[globalIdx * 3 + 1] / minMax;
		//texCoord3DPtrVBO[globalIdx * 3 + 2] = vertexPtrVBO[globalIdx * 3 + 2] / minMax;
	}
}

#ifdef LINE
__global__ void make_bary(TriangleD *triangleListD, float *baryPtrVBO, float *oriBaryPtrVBO, int n, int q)
{
	/***************************** 生成切割后每个顶点的重心坐标 ******************************/
	int localIdx = threadIdx.x;
	float tempFloorFloat = (sqrtf((float)localIdx * 8 + 9) - 3) / 2;
	int floor = rintf(tempFloorFloat);
	if ((floor * 2 + 3) * (floor * 2 + 3) != localIdx * 8 + 9)
		floor = ceilf(tempFloorFloat);
	int room = localIdx - (floor + 1) * floor / 2;
	float3 barycentric_coord;
	barycentric_coord.x = (float)(n - floor) / n;
	barycentric_coord.y = (float)(floor - room) / n;
	barycentric_coord.z = 1.0f - barycentric_coord.x - barycentric_coord.y;

	int globalIdx = blockDim.x * blockIdx.x + threadIdx.x;
	baryPtrVBO[globalIdx * 3 + 0] = (float)(n - floor) / n;
	baryPtrVBO[globalIdx * 3 + 1] = barycentric_coord.y;
	baryPtrVBO[globalIdx * 3 + 2] = barycentric_coord.z;

	/***************************** 生成切割前每个顶点的重心坐标 ******************************/
	int triangleIdx = blockIdx.x;

	// 切割后三角形三个顶点在原始三角形中的重心坐标
	float3 bary_origin0 = triangleListD[triangleIdx].bary_origin[0];
	float3 bary_origin1 = triangleListD[triangleIdx].bary_origin[1];
	float3 bary_origin2 = triangleListD[triangleIdx].bary_origin[2];

	// 当前点在原始三角形上的重心坐标
	float3 bary_origin = bary_origin0 * barycentric_coord.x + bary_origin1 * barycentric_coord.y + bary_origin2 * barycentric_coord.z;

	// 存储目前处理的采样点在原始三角片上的重心坐标
	oriBaryPtrVBO[globalIdx * 3 + 0] = bary_origin.x;
	oriBaryPtrVBO[globalIdx * 3 + 1] = bary_origin.y;
	oriBaryPtrVBO[globalIdx * 3 + 2] = bary_origin.z;
}
#endif

#ifdef TRUTH
__global__ void copy_truth(float *RD_truth,
					 int activeThreadNumCopy, bool firstLoad,
					 TriangleD *triangleListD, int segmentPerEdge, int f, int q,
					 float *normalPtrVBO_truth, float *vertexPtrVBO_truth)
{
	int globalIdx = blockDim.x * blockIdx.x + threadIdx.x;
	if (globalIdx >= activeThreadNumCopy)
		return;

	int triangleIdx = globalIdx / q;
	int localIdx = globalIdx % q;

	vertexPtrVBO_truth[globalIdx * 3 + 0] = RD_truth[triangleIdx * q * 3 + localIdx];
	vertexPtrVBO_truth[globalIdx * 3 + 1] = RD_truth[(triangleIdx + f) * q * 3 + localIdx];
	vertexPtrVBO_truth[globalIdx * 3 + 2] = RD_truth[(triangleIdx + f * 2) * q * 3 + localIdx];

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////

	float ux = RD_truth[triangleIdx * q * 3 + q + localIdx];
	float uy = RD_truth[(triangleIdx + f) * q * 3 + q + localIdx];
	float uz = RD_truth[(triangleIdx + f * 2) * q * 3 + q + localIdx];

	float vx = RD_truth[triangleIdx * q * 3 + q * 2 + localIdx];
	float vy = RD_truth[(triangleIdx + f) * q * 3 + q * 2 + localIdx];
	float vz = RD_truth[(triangleIdx + f * 2) * q * 3 + q * 2 + localIdx];

	float nx = uy * vz - uz * vy;
	float ny = uz * vx - ux * vz;
	float nz = ux * vy - uy * vx;
	float l = sqrtf(nx * nx + ny * ny + nz * nz);
	nx /= l;
	ny /= l;
	nz /= l;

	normalPtrVBO_truth[globalIdx * 3 + 0] = nx;
	normalPtrVBO_truth[globalIdx * 3 + 1] = ny;
	normalPtrVBO_truth[globalIdx * 3 + 2] = nz;
}
#endif

bool registered = false;
GLuint normalVBO = 0, texCoordVBO = 0, texCoord3DVBO = 0, vertexVBO = 0;
#ifdef LINE
GLuint baryVBO = 0, oriBaryVBO = 0;
#endif
float *normalPtrVBO;							// 读写缓冲区对象所用的指针
float *texCoordPtrVBO;							// 读写缓冲区对象所用的指针
float *texCoord3DPtrVBO;						// 读写缓冲区对象所用的指针
float *vertexPtrVBO;							// 读写缓冲区对象所用的指针
#ifdef LINE
float *baryPtrVBO, *oriBaryPtrVBO;							// 读写缓冲区对象所用的指针
#endif

struct cudaGraphicsResource *normalVBO_CUDA, *texCoordVBO_CUDA, *texCoord3DVBO_CUDA, *vertexVBO_CUDA;
#ifdef LINE
struct cudaGraphicsResource *baryVBO_CUDA, *oriBaryVBO_CUDA;
#endif


void tessellateD(bool firstLoad, float maxX, float maxY, float maxZ, AlgorithmType algo_type)
{
	float alpha = 1.0f, beta = 0.0f;
	// 计算三角化点的坐标和法向
	cublasStatus_t stat;
	if (algo_type == CYM)
	{
		stat = cublasSgemm(cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N,
										  samplePointPerTriangle, triangleNum * 6, triangleCtrlPointNum_lower,
										  &alpha,
										  BqD, samplePointPerTriangle,
										  triangleCtrlPointD, triangleCtrlPointNum_lower,
										  &beta,
										  RD, samplePointPerTriangle);
	}
	//else if (algo_type == PN_CUTTING)
	else
	{
		stat = cublasSgemm(cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N,
										  samplePointPerTriangle, triangleNum * 3, triangleCtrlPointNum_lower,
										  &alpha,
										  BqD, samplePointPerTriangle,
										  triangleCtrlPointD_PN, triangleCtrlPointNum_lower,
										  &beta,
										  RD, samplePointPerTriangle);
		stat = cublasSgemm(cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N,
										  samplePointPerTriangle, triangleNum * 3, 6,
										  &alpha,
										  BqD_PN, samplePointPerTriangle,
										  triangleNormalCtrlPointD_PN, 6,
										  &beta,
										  RD + samplePointPerTriangle * triangleNum * 3, samplePointPerTriangle);
	}
	if (stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "RD fail!!!!!!!!!!!!!\tstat = " << stat << "\t\t";
		printCudaError(__FILE__, __FUNCTION__, __LINE__);
		cout << endl;
		return;
	}

	cudaGraphicsMapResources(1, &normalVBO_CUDA, 0);
	cudaGraphicsMapResources(1, &texCoordVBO_CUDA, 0);
	cudaGraphicsMapResources(1, &texCoord3DVBO_CUDA, 0);
	cudaGraphicsMapResources(1, &vertexVBO_CUDA, 0);
#ifdef LINE
	cudaGraphicsMapResources(1, &baryVBO_CUDA, 0);
	cudaGraphicsMapResources(1, &oriBaryVBO_CUDA, 0);
#endif
	//size_t size2 = sizeof(float) * samplePointPerTriangle * triangleNum * 2;
	size_t size3 = sizeof(float) * samplePointPerTriangle * triangleNum * 3;
	cudaGraphicsResourceGetMappedPointer((void**)&normalPtrVBO, &size3, normalVBO_CUDA);
	cudaGraphicsResourceGetMappedPointer((void**)&texCoordPtrVBO, &size3, texCoordVBO_CUDA);
	cudaGraphicsResourceGetMappedPointer((void**)&texCoord3DPtrVBO, &size3, texCoord3DVBO_CUDA);
	cudaGraphicsResourceGetMappedPointer((void**)&vertexPtrVBO, &size3, vertexVBO_CUDA);
#ifdef LINE
	cudaGraphicsResourceGetMappedPointer((void**)&baryPtrVBO, &size3, baryVBO_CUDA);
	cudaGraphicsResourceGetMappedPointer((void**)&oriBaryPtrVBO, &size3, oriBaryVBO_CUDA);
#endif

	copy<<<blockNumCopy, blockSizeCopy>>>(RD, u_seg, v_seg,
										  activeThreadNumCopy, firstLoad, maxX, maxY, maxZ, triangleListD, segmentPerEdge, triangleNum, samplePointPerTriangle,
										  normalPtrVBO, texCoordPtrVBO, texCoord3DPtrVBO, vertexPtrVBO);

#ifdef LINE
	make_bary<<<triangleNum, samplePointPerTriangle>>>(triangleListD, baryPtrVBO, oriBaryPtrVBO, segmentPerEdge, samplePointPerTriangle);
#endif

	cudaGraphicsUnmapResources(1, &normalVBO_CUDA, 0);
	cudaGraphicsUnmapResources(1, &texCoordVBO_CUDA, 0);
	cudaGraphicsUnmapResources(1, &texCoord3DVBO_CUDA, 0);
	cudaGraphicsUnmapResources(1, &vertexVBO_CUDA, 0);
#ifdef LINE
	cudaGraphicsUnmapResources(1, &baryVBO_CUDA, 0);
	cudaGraphicsUnmapResources(1, &oriBaryVBO_CUDA, 0);
#endif
}

//#ifdef TRUTH
GLuint normalVBO_truth = 0, vertexVBO_truth = 0;
float *normalPtrVBO_truth;						// 读写缓冲区对象所用的指针
float *vertexPtrVBO_truth;						// 读写缓冲区对象所用的指针

struct cudaGraphicsResource* normalVBO_CUDA_truth;
struct cudaGraphicsResource* vertexVBO_CUDA_truth;

//double vertex_error_ave_max = 0.0, vertex_error_max_max = 0.0;
//double normal_error_ave_max = 0.0, normal_error_max_max = 0.0;

//int triangleCoord(int floor, int room)
//{
	//return (1 + floor) * floor / 2 + room;
//}

//__global__ void my_to_truth(int f, int q, int point_per_real_face_u, int point_per_real_face_v,
		//float *myV, float *realV, float *parameterD, int *my_to_truth_tableD, int *belongs_to_originD)
//{
	//int globalIdx = blockIdx.x * blockDim.x + threadIdx.x;
	//int triangleIdx = globalIdx / q;
	//if (triangleIdx >= f)
		//return;
	//realV += belongs_to_originD[triangleIdx] * point_per_real_face_u * point_per_real_face_v * 3;

	//float u = parameterD[globalIdx * 3];
	//float v = parameterD[globalIdx * 3 + 1];
	//int u_pre = (int)(u * (point_per_real_face_u - 1));
	//int v_pre = (int)(v * (point_per_real_face_v - 1));
	//if (u_pre < 0)
		//u_pre = 0;
	//if (v_pre < 0)
		//v_pre = 0;
	//int u_next = u_pre + 1;
	//int v_next = v_pre + 1;
	//if (u_next > point_per_real_face_u - 1)
		//u_next = point_per_real_face_u - 1;
	//if (v_next > point_per_real_face_v - 1)
		//v_next = point_per_real_face_v - 1;
	////if (blockIdx.x == 0)
	////{
		////printf("u, v = %f, %f\n", u, v);
	////}

	//float dist_min = 999999, idx_min = -1;
	//// 编号u，v的点
		//float dx = myV[globalIdx * 3] - realV[(u_pre * point_per_real_face_v + v_pre) * 3];
		//float dy = myV[globalIdx * 3 + 1] - realV[(u_pre * point_per_real_face_v + v_pre) * 3 + 1];
		//float dz = myV[globalIdx * 3 + 2] - realV[(u_pre * point_per_real_face_v + v_pre) * 3 + 2];
		//float dist = sqrt(dx * dx + dy * dy + dz * dz);
		//if (dist < dist_min)
		//{
			//dist_min = dist;
			//idx_min = u_pre * point_per_real_face_v + v_pre;
		//}

	//// 编号u+1，v的点
		//dx = myV[globalIdx * 3] - realV[(u_next * point_per_real_face_v + v_pre) * 3];
		//dy = myV[globalIdx * 3 + 1] - realV[(u_next * point_per_real_face_v + v_pre) * 3 + 1];
		//dz = myV[globalIdx * 3 + 2] - realV[(u_next * point_per_real_face_v + v_pre) * 3 + 2];
		//dist = sqrt(dx * dx + dy * dy + dz * dz);
		//if (dist < dist_min)
		//{
			//dist_min = dist;
			//idx_min = u_next * point_per_real_face_v + v_pre;
		//}

	//// 编号u，v+1的点
		//dx = myV[globalIdx * 3] - realV[(u_pre * point_per_real_face_v + v_next) * 3];
		//dy = myV[globalIdx * 3 + 1] - realV[(u_pre * point_per_real_face_v + v_next) * 3 + 1];
		//dz = myV[globalIdx * 3 + 2] - realV[(u_pre * point_per_real_face_v + v_next) * 3 + 2];
		//dist = sqrt(dx * dx + dy * dy + dz * dz);
		//if (dist < dist_min)
		//{
			//dist_min = dist;
			//idx_min = u_pre * point_per_real_face_v + v_next;
		//}

	//// 编号u+1，v+1的点
		//dx = myV[globalIdx * 3] - realV[(u_next * point_per_real_face_v + v_next) * 3];
		//dy = myV[globalIdx * 3 + 1] - realV[(u_next * point_per_real_face_v + v_next) * 3 + 1];
		//dz = myV[globalIdx * 3 + 2] - realV[(u_next * point_per_real_face_v + v_next) * 3 + 2];
		//dist = sqrt(dx * dx + dy * dy + dz * dz);
		//if (dist < dist_min)
		//{
			//dist_min = dist;
			//idx_min = u_next * point_per_real_face_v + v_next;
		//}
	//my_to_truth_tableD[globalIdx] = idx_min + belongs_to_originD[triangleIdx] * point_per_real_face_u * point_per_real_face_v;
//}

__global__ void deformTeapot(float3 *vertexParamListD_teapot, float3 *normalParamListD_teapot,
		int activeThreadNum,
		float *normalPtrVBO_truth, float *vertexPtrVBO_truth,
		int orderU, int orderV, int orderW,
		int ctrlPointNumU, int ctrlPointNumV, int ctrlPointNumW)
{
	// u, v, w 表示经过重心坐标插值之后的采样点的x, y, z分量
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= activeThreadNum)
		return;
	float u = vertexParamListD_teapot[i].x;
	float v = vertexParamListD_teapot[i].y;
	float w = vertexParamListD_teapot[i].z;

	// u, v, w方向节点区间数量
	int knot_interval_count_u = orderU + ctrlPointNumU - (orderU - 1) * 2 - 1;
	int knot_interval_count_v = orderV + ctrlPointNumV - (orderV - 1) * 2 - 1;
	int knot_interval_count_w = orderW + ctrlPointNumW - (orderW - 1) * 2 - 1;

	// 预先将其值设为最大，将末端点归入最后一段 
	int left_idx_u = orderU - 1 + knot_interval_count_u - 1;
	int left_idx_v = orderV - 1 + knot_interval_count_v - 1;
	int left_idx_w = orderW - 1 + knot_interval_count_w - 1;

	// 沿 U 方向查找当前点所在的节点区间 
	for (int ii = orderU - 1; ii <= orderU - 1 + knot_interval_count_u - 1; ++ii)
	{
		if (u >= knotListD[ii] && u < knotListD[ii + 1])
		{
			left_idx_u = ii;
			break;
		}
	}
	// 沿 V 方向查找当前点所在的节点区间 
	for (int jj = orderV - 1; jj <= orderV - 1 + knot_interval_count_v - 1; ++jj)
	{
		if (v >= knotListD[20 + jj] && v < knotListD[20 + jj + 1])
		{
			left_idx_v = jj;
			break;
		}
	}
	// 沿 W 方向查找当前点所在的节点区间 
	for (int kk = orderW - 1; kk <= orderW - 1 + knot_interval_count_w - 1; ++kk)
	{
		if (w >= knotListD[40 + kk] && w < knotListD[40 + kk + 1])
		{
			left_idx_w = kk;
			break;
		}
	}

	float tmpKnot = knotListD[left_idx_u];
	float tmpKnot1 = knotListD[left_idx_u + 1];
	float x_stride = tmpKnot1 - tmpKnot;
	u = (u - tmpKnot) / x_stride;

	tmpKnot = knotListD[20 + left_idx_v];
	tmpKnot1 = knotListD[20 + left_idx_v + 1];
	float y_stride = tmpKnot1 - tmpKnot;
	v = (v - tmpKnot) / y_stride;

	tmpKnot = knotListD[40 + left_idx_w];
	tmpKnot1 = knotListD[40 + left_idx_w + 1];
	float z_stride = tmpKnot1 - tmpKnot;
	w = (w - tmpKnot) / z_stride;

	extern __shared__ float shared_array[];
	// 算出该线程负责的采样点的 B 样条体值
	// fu 表示J_bar矩阵第一列三个元素：偏F_bar_x偏u、偏F_bar_y偏u、偏F_bar_z偏u
	// fv 表示J_bar矩阵第二列三个元素：偏F_bar_x偏v、偏F_bar_y偏v、偏F_bar_z偏v
	float3 result, fu, fv;
	BSplineVolumeValueMatrixD_combine(u, v, w, shared_array,
	left_idx_u - (orderU - 1), left_idx_v - (orderV - 1), left_idx_w - (orderW - 1),
											   orderU, orderV, orderW,
											   result, fu, fv);
	__syncthreads();

	//sampleValueD[index2c(localIdx, triangleIdx		  , m + c)] = result.x;
	//sampleValueD[index2c(localIdx, triangleIdx + f	  , m + c)] = result.y;
	//sampleValueD[index2c(localIdx, triangleIdx + f * 2, m + c)] = result.z;
	vertexPtrVBO_truth[i * 3 + 0] = result.x;
	vertexPtrVBO_truth[i * 3 + 1] = result.y;
	vertexPtrVBO_truth[i * 3 + 2] = result.z;

	//vertexPtrVBO_truth[0] = -1.0;
	//vertexPtrVBO_truth[1] = 0.0;
	//vertexPtrVBO_truth[2] = 1.0;

	//vertexPtrVBO_truth[3] = 0.0;
	//vertexPtrVBO_truth[4] = -1.0;
	//vertexPtrVBO_truth[5] = 1.0;

	//printf("%f, %f, %f = %f, %f, %f\n", u, v, w, vertexPtrVBO_truth[i * 3], vertexPtrVBO_truth[i * 3 + 1], vertexPtrVBO_truth[i * 3 + 2]);

	///////////////////////////////////////////////////////////////////////////////

	// fw 表示J_bar矩阵第三列三个元素：偏F_bar_x偏w、偏F_bar_y偏w、偏F_bar_z偏w
	float3 fw = BSplineVolumeValueMatrixDw(u, v, w, shared_array,
	left_idx_u - (orderU - 1), left_idx_v - (orderV - 1), left_idx_w - (orderW - 1),
										   orderU, orderV, orderW);
	//__syncthreads();

	u = normalParamListD_teapot[i].x;
	v = normalParamListD_teapot[i].y;
	w = normalParamListD_teapot[i].z;
	//printf("%f, %f, %f\n", u, v, w);

	//float *sampleNormalD = sampleValueD + 3 * f * (m + c);
	// J_bar_star_T_[012]表示J_bar的伴随矩阵的转置(即J_bar*T)的第一行三个元素
	float J_bar_star_T_0 = fv.y * fw.z - fw.y * fv.z;
	float J_bar_star_T_1 = fw.y * fu.z - fu.y * fw.z;
	float J_bar_star_T_2 = fu.y * fv.z - fv.y * fu.z;
	normalPtrVBO_truth[i * 3 + 0] = u * J_bar_star_T_0 * x_stride + v * J_bar_star_T_1 * y_stride + w * J_bar_star_T_2 * z_stride;
	//sampleNormalD[index2c(localIdx, triangleIdx, m + c)] = u * J_bar_star_T_0 * x_stride + v * J_bar_star_T_1 * y_stride + w * J_bar_star_T_2 * z_stride;

	// J_bar_star_T_[012]表示J_bar的伴随矩阵的转置(即J_bar*T)的第二行三个元素
	J_bar_star_T_0 = fv.z * fw.x - fw.z * fv.x;
	J_bar_star_T_1 = fw.z * fu.x - fu.z * fw.x;
	J_bar_star_T_2 = fu.z * fv.x - fv.z * fu.x;
	normalPtrVBO_truth[i * 3 + 1] = u * J_bar_star_T_0 * x_stride + v * J_bar_star_T_1 * y_stride + w * J_bar_star_T_2 * z_stride;
	//sampleNormalD[index2c(localIdx, triangleIdx + f, m + c)] = u * J_bar_star_T_0 * x_stride + v * J_bar_star_T_1 * y_stride + w * J_bar_star_T_2 * z_stride;

	// J_bar_star_T_[012]表示J_bar的伴随矩阵的转置(即J_bar*T)的第三行三个元素
	J_bar_star_T_0 = fv.x * fw.y - fw.x * fv.y;
	J_bar_star_T_1 = fw.x * fu.y - fu.x * fw.y;
	J_bar_star_T_2 = fu.x * fv.y - fv.x * fu.y;
	normalPtrVBO_truth[i * 3 + 2] = u * J_bar_star_T_0 * x_stride + v * J_bar_star_T_1 * y_stride + w * J_bar_star_T_2 * z_stride;
	//sampleNormalD[index2c(localIdx, triangleIdx + f * 2, m + c)] = u * J_bar_star_T_0 * x_stride + v * J_bar_star_T_1 * y_stride + w * J_bar_star_T_2 * z_stride;

	//printf("%f, %f, %f = %f, %f, %f\n", u, v, w, normalPtrVBO_truth[i * 3], normalPtrVBO_truth[i * 3 + 1], normalPtrVBO_truth[i * 3 + 2]);
}





using namespace objdata;

double color_map_vertex(const VertexCoord &v0, const VertexCoord v1, double range)
{
	return (v0 - v1).norm() / range;
}

double color_map_normal(const NormalCoord &n0, const NormalCoord &n1, double range)
{
	double result = 2 * asin((n0 - n1).norm() * 0.5) / range;
	if (result < 0)
		result = 0;
	else if (result > 1)
		result = 1;
	return result;
	//return 2 * asin((n0 - n1).norm() * 0.5) / range;
}

__device__ float power(float a, int x)
{
	float result = 1;
	for (int i = 0; i < x; ++i)
		result *= a;
	return result;
}

__device__ float B(float t, int n, int i)
{
	int factorial[] = {1, 1, 2, 6};
	double factorial_1[] = {1, 1, 0.5, 0.166666666666666666};

	return factorial[n] * factorial_1[i] * factorial_1[n - i]
			* power(t, i) * power(1 - t, n - i);
}

__global__ void calc3Dparameter(float *parameterD, float3 *parameter3D, float3 *parameterND, int total)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx >= total)
		return;

	float u = parameterD[idx * 3 + 0];
	float v = parameterD[idx * 3 + 1];
	int surface_id = parameterD[idx * 3 + 2];

	float3 *cp = teapot_ctrl_pointD + surface_id * 16;

	float3 result = make_float3(0, 0, 0);
	for (int i = 0; i < 4; ++i)
	{
		float3 temp_ctrlpoint = make_float3(0, 0, 0);
		for (int j = 0; j < 4; ++j)
		{
			temp_ctrlpoint += B(v, 3, j) * cp[i * 4 + j];
		}
		result += B(u, 3, i) * temp_ctrlpoint;
	}
	parameter3D[idx] = result;
	//printf("3d = %f, %f, %f\n", result.x, result.y, result.z);
	//if (parameter3D[idx].x > 1.0)
		//parameter3D[idx].x = 1.0;
	//else if (parameter3D[idx].x < -1.0)
		//parameter3D[idx].x = -1.0;

	//if (parameter3D[idx].y > 1.0)
		//parameter3D[idx].y = 1.0;
	//else if (parameter3D[idx].y < -1.0)
		//parameter3D[idx].y = -1.0;

	//if (parameter3D[idx].z > 1.0)
		//parameter3D[idx].z = 1.0;
	//else if (parameter3D[idx].z < -1.0)
		//parameter3D[idx].z = -1.0;
	//if (blockIdx.x == 0)
		//printf("%d : (%f, %f, %f) = (%f, %f, %f)\n", idx, u, v, surface_id, result.x, result.y, result.z);

	if (surface_id < 4 && u < 0.0001)				// 处理壶盖顶部的三角形面片
	{
		parameterND[idx] = make_float3(0, 0, 1);
		//printf("xiaoyu\n");
	}
	else if (surface_id < 8 && u < 0.0001)			// 处理壶身底部的三角形面片
	{
		parameterND[idx] = make_float3(0, 0, -1);
		//printf("xiaoyu2\n");
	}
	else
	{
		// 计算(u, v)点u方向导矢的值
		float3 result_u = make_float3(0, 0, 0);
		for (int j = 0; j < 4; ++j)
		{
			float3 temp_ctrlpoint = make_float3(0, 0, 0);
			for (int i = 0; i < 3; ++i)
			{
				temp_ctrlpoint += B(u, 2, i) * (cp[i * 4 + j] - cp[(i + 1) * 4 + j]);
			}
			result_u += B(v, 3, j) * temp_ctrlpoint;
		}

		// 计算(u, v)点v方向导矢的值
		float3 result_v = make_float3(0, 0, 0);
		for (int i = 0; i < 4; ++i)
		{
			float3 temp_ctrlpoint = make_float3(0, 0, 0);
			for (int j = 0; j < 3; ++j)
			{
				temp_ctrlpoint += B(v, 2, j) * (cp[i * 4 + j] - cp[i * 4 + j + 1]);
			}
			result_v += B(u, 3, i) * temp_ctrlpoint;
		}
		//Point3 normal = cross(result_u, result_v);
		float3 normal = cross(result_v, result_u);
		normalize(normal);
		parameterND[idx] = normal;
		//if (blockIdx.x == 0)
		//printf("(%f, %f, %f) = (%f, %f, %f)\n", u, v, surface_id, normal.x, normal.y, normal.z);

		//if (surface_id < 4 && u < 0.0001)
		//{
			////parameterND[idx] = make_float3(0, 0, 1);
			//printf("xiaoyu, ND = %f, %f, %f\n", parameterND[idx].x, parameterND[idx].y, parameterND[idx].z);
		//}
		//else if (surface_id < 8 && u < 0.0001)
		//{
			////parameterND[idx] = make_float3(0, 0, -1);
			//printf("xiaoyu2, ND = %f, %f, %f\n", parameterND[idx].x, parameterND[idx].y, parameterND[idx].z);
		//}
	}
}

float *texture_coord;

void tessellateD_truth(bool adjust_silhouette, bool firstLoad)
{
	cudaError_t cymError;
	//cymError = cudaMemcpy(my_to_truth_table, my_to_truth_tableD, sizeof(int) * samplePointPerTriangle * triangleNum, cudaMemcpyDeviceToHost);
	//if (cymError)
		//cout << __FILE__ << "第" << __LINE__ << "行, 错误代码" << cymError << ": " << cudaGetErrorString(cymError) << endl;
	//size_t size2 = sizeof(float) * samplePointPerTriangle * triangleNum * 2;
	size_t size3 = sizeof(float) * samplePointPerTriangle * triangleNum * 3;
	//size_t size3_truth = sizeof(float) * samplePointPerTriangle * triangleNum * 3;
	cymError = cudaGetLastError();
	if (cymError)
		cout << __FILE__ << "第" << __LINE__ << "行, 错误代码" << cymError << ": " << cudaGetErrorString(cymError) << endl;
	cudaGraphicsMapResources(1, &normalVBO_CUDA, 0);
	cudaGraphicsMapResources(1, &vertexVBO_CUDA, 0);
	cudaGraphicsMapResources(1, &texCoordVBO_CUDA, 0);
	cudaGraphicsMapResources(1, &texCoord3DVBO_CUDA, 0);
	cudaGraphicsMapResources(1, &normalVBO_CUDA_truth, 0);
	cudaGraphicsMapResources(1, &vertexVBO_CUDA_truth, 0);
	cudaGraphicsResourceGetMappedPointer((void**)&normalPtrVBO, &size3, normalVBO_CUDA);
	cudaGraphicsResourceGetMappedPointer((void**)&vertexPtrVBO, &size3, vertexVBO_CUDA);
	cudaGraphicsResourceGetMappedPointer((void**)&texCoordPtrVBO, &size3, texCoordVBO_CUDA);
	cudaGraphicsResourceGetMappedPointer((void**)&texCoord3DPtrVBO, &size3, texCoord3DVBO_CUDA);
	cudaGraphicsResourceGetMappedPointer((void**)&normalPtrVBO_truth, &size3, normalVBO_CUDA_truth);
	cudaGraphicsResourceGetMappedPointer((void**)&vertexPtrVBO_truth, &size3, vertexVBO_CUDA_truth);

	// 计算参数
	if (firstLoad)
	{
		int block_size = 128;
		int block_num = ceil(static_cast<double>(triangleNum * samplePointPerTriangle) / block_size);
		calc3Dparameter<<<block_num, block_size>>>(texCoordPtrVBO, parameter3D, parameterND, samplePointPerTriangle * triangleNum);
		cymError = cudaGetLastError();
		if (cymError)
			cout << __FILE__ << "第" << __LINE__ << "行, 错误代码" << cymError << ": " << cudaGetErrorString(cymError) << endl;

		texture_coord = new float[size3 / sizeof(float)];
		//cout << "第一次，找对应" << endl;
		//cout << "block_num = " << block_num << endl;
		//int size = belongs_to_origin.size();
		//int *belongs_to_originD;
		//cudaMalloc((void**)&belongs_to_originD, sizeof(int) * size);
		////cout << "belongsize = " << size << endl;
		//cymError = cudaMemcpy(belongs_to_originD, &belongs_to_origin[0], sizeof(int) * size, cudaMemcpyHostToDevice);
		//if (cymError)
			//cout << __FILE__ << "第" << __LINE__ << "行, 错误代码" << cymError << ": " << cudaGetErrorString(cymError) << endl;
		//block_size = 128;
		//block_num = ceil(static_cast<double>(triangleNum * samplePointPerTriangle) / block_size);
		//cout << "my_to_truth.blockNum = " << block_num << endl;
		//my_to_truth<<<block_num, block_size>>>(triangleNum, samplePointPerTriangle,
				//(u_seg + 1), (v_seg + 1), vertexPtrVBO, vertexPtrVBO_truth, parameterD,
				//my_to_truth_tableD, belongs_to_originD);
		//cymError = cudaGetLastError();
		//if (cymError)
			//cout << __FILE__ << "第" << __LINE__ << "行, 错误代码" << cymError << ": " << cudaGetErrorString(cymError) << endl;

		//cymError = cudaMemcpy(my_to_truth_table, my_to_truth_tableD, sizeof(int) * samplePointPerTriangle * triangleNum, cudaMemcpyDeviceToHost);
		//if (cymError)
			//cout << __FILE__ << "第" << __LINE__ << "行, 错误代码" << cymError << ": " << cudaGetErrorString(cymError) << endl;


		//cout << "triangleNum = " << triangleNum << ", samplePointPerTriangle = "
			 //<< samplePointPerTriangle << endl;
		//for (int i = 0; i < samplePointPerTriangle * triangleNum; ++i)
		//{
			//cout << my_to_truth_table[i] << " ";
			//if (i % 20 == 19)
				//cout << endl;
		//}

		//int *tttt = new int[size];
		//cudaMemcpy(tttt, belongs_to_originD, sizeof(int) * size, cudaMemcpyDeviceToHost);
		//for (int i = 0; i < size; ++i)
		//{
			//cout << tttt[i] << " ";
			//if (i % 20 == 19)
				//cout << endl;
		//}
		//delete []tttt;
	}
	//cout << "计算参数完成" << endl;

	// 变形基准茶壶
	int block_size = 128;
	int block_num = ceil(static_cast<double>(triangleNum * samplePointPerTriangle) / block_size);
	deformTeapot<<<block_num, block_size, sizeof(float) * block_size * 13>>>
					(parameter3D, parameterND,
					triangleNum * samplePointPerTriangle,
					normalPtrVBO_truth, vertexPtrVBO_truth,
					order[U], order[V], order[W],
					ctrlPointNum[U], ctrlPointNum[V], ctrlPointNum[W]);
	cymError = cudaGetLastError();
	if (cymError)
		cout << __FILE__ << "第" << __LINE__ << "行, 错误代码" << cymError << ": " << cudaGetErrorString(cymError) << endl;
	//cout << "基准茶壶变形完成" << endl;

	//cudaMemcpy(vertexPtrVBO_truth, vertexPtrVBO, sizeof(float) * 3 * samplePointPerTriangle * 8, cudaMemcpyDeviceToDevice);

	/*------------------------ 测量误差 ----------------------------*/
	/* 顶点误差 */
	float *result = new float[size3 / sizeof(float)];
	float *result_truth = new float[size3 / sizeof(float)];

	cymError = cudaMemcpy(result, vertexPtrVBO, size3, cudaMemcpyDeviceToHost);
	if (cymError)
		cout << __FILE__ << "第" << __LINE__ << "行, 错误代码" << cymError << ": " << cudaGetErrorString(cymError) << endl;
	cymError = cudaMemcpy(result_truth, vertexPtrVBO_truth, size3, cudaMemcpyDeviceToHost);
	if (cymError)
		cout << __FILE__ << "第" << __LINE__ << "行, 错误代码" << cymError << ": " << cudaGetErrorString(cymError) << endl;

	double vertex_error_ave_max = 0.0, vertex_error_max_max = 0.0;
	double normal_error_ave_max = 0.0, normal_error_max_max = 0.0;
	double error_ave = 0.0, error_max = 0.0;
	//cout << "准备进入for循环计算每个点的误差, size3 = " << size3 << ", size3 = " << size3 << endl;
	//double maxZ = -100000, minZ = 100000, minX, maxX, minY, maxY;
	for (int i = 0; i < samplePointPerTriangle * triangleNum; ++i)
	{
		double x0 = result[i * 3];
		double y0 = result[i * 3 + 1];
		double z0 = result[i * 3 + 2];
		//if (z0 > maxZ)
		//{
			//maxX = x0;
			//maxY = y0;
			//maxZ = z0;
		//}
		//if (z0 < minZ)
		//{
			//minX = x0;
			//minY = y0;
			//minZ = z0;
		//}
		//int real_idx = my_to_truth_table[i];
		//double x1 = result_truth[real_idx * 3];
		//double y1 = result_truth[real_idx * 3 + 1];
		//double z1 = result_truth[real_idx * 3 + 2];
		double x1 = result_truth[i * 3];
		double y1 = result_truth[i * 3 + 1];
		double z1 = result_truth[i * 3 + 2];
		//cout << "顶点 " << x0 << ", " << y0 << ", " << z0 <<  "====" << x1 << ", " << y1 << ", " << z1 << endl;
		double error = sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1) + (z0 - z1) * (z0 - z1));
		//cout << "error = " << error << endl;
		error_ave += error;
		if (error_max < error)
			error_max = error;

		float vertex_diff = color_map_vertex(VertexCoord(x0, y0, z0), VertexCoord(x1, y1, z1), 0.04);
		texture_coord[i * 3] = vertex_diff;
		texture_coord[i * 3 + 1] = 0.5;
	}
	cout << "for循环完成" << endl;
	//cout << "maxZ = " << maxZ << ", minZ = " << minZ << endl;
	//cout << "min = " << minX << ", " << minY << ", " << minZ << endl;
	//cout << "max = " << maxX << ", " << maxY << ", " << maxZ << endl;
	cudaMemcpy(texCoordPtrVBO, texture_coord, size3, cudaMemcpyHostToDevice);
	/*cout << "eeeeee samplePonitPerTriangle = " << samplePointPerTriangle << endl;*/
	/*cout << "eeeeee triangleNum = " << triangleNum << endl;*/
	/*cout << "eeeeee samplePonitPerTriangle * triangleNum = " << samplePointPerTriangle * triangleNum << endl;*/

	/*cout << "eeeeee error = " << error_ave / (samplePointPerTriangle * triangleNum)*/
		 /*<< ", error_max = " << error_max << endl;*/

	if (error_ave > vertex_error_ave_max)
		vertex_error_ave_max = error_ave;
	if (error_max > vertex_error_max_max)
		vertex_error_max_max = error_max;
	if (adjust_silhouette)
		cout << "调整过，误差大" << endl;
	else
		cout << "未调整，误差小" << endl;
	cout << "eeeeee 平均顶点误差 = " << vertex_error_ave_max / (samplePointPerTriangle * triangleNum) << ", 最大顶点误差 = " << vertex_error_max_max << endl;
	cout << "error_ave = " << error_ave / (samplePointPerTriangle * triangleNum) << endl;

	//[> 体积误差 <]
	//double volume = 0.0;
	//for (int f = 0; f < triangleNum; ++f)
	//{
		//for (int i = 0; i < segmentPerEdge; ++i)
		//{
			//for (int j = 0; j <= i; ++j)
			//{
				//// smooth FFD算法结果
				//double v0x = result[samplePointPerTriangle * 3 * f + triangleCoord(i, j) * 3 + 0];
				//double v0y = result[samplePointPerTriangle * 3 * f + triangleCoord(i, j) * 3 + 1];
				//double v0z = result[samplePointPerTriangle * 3 * f + triangleCoord(i, j) * 3 + 2];
				//double v1x = result[samplePointPerTriangle * 3 * f + triangleCoord(i + 1, j) * 3 + 0];
				//double v1y = result[samplePointPerTriangle * 3 * f + triangleCoord(i + 1, j) * 3 + 1];
				//double v1z = result[samplePointPerTriangle * 3 * f + triangleCoord(i + 1, j) * 3 + 2];
				//double v2x = result[samplePointPerTriangle * 3 * f + triangleCoord(i + 1, j + 1) * 3 + 0];
				//double v2y = result[samplePointPerTriangle * 3 * f + triangleCoord(i + 1, j + 1) * 3 + 1];
				//double v2z = result[samplePointPerTriangle * 3 * f + triangleCoord(i + 1, j + 1) * 3 + 2];
				//volume += (v0z + v1z + v2z) * ((v1x - v0x) * (v2y - v0y) - (v2x - v0x) * (v1y - v0y));
				//if (i < segmentPerEdge - 1)
				//{
					//double v0x = result[samplePointPerTriangle * 3 * f + triangleCoord(i + 1, j) * 3 + 0];
					//double v0y = result[samplePointPerTriangle * 3 * f + triangleCoord(i + 1, j) * 3 + 1];
					//double v0z = result[samplePointPerTriangle * 3 * f + triangleCoord(i + 1, j) * 3 + 2];
					//double v1x = result[samplePointPerTriangle * 3 * f + triangleCoord(i + 2, j + 1) * 3 + 0];
					//double v1y = result[samplePointPerTriangle * 3 * f + triangleCoord(i + 2, j + 1) * 3 + 1];
					//double v1z = result[samplePointPerTriangle * 3 * f + triangleCoord(i + 2, j + 1) * 3 + 2];
					//double v2x = result[samplePointPerTriangle * 3 * f + triangleCoord(i + 1, j + 1) * 3 + 0];
					//double v2y = result[samplePointPerTriangle * 3 * f + triangleCoord(i + 1, j + 1) * 3 + 1];
					//double v2z = result[samplePointPerTriangle * 3 * f + triangleCoord(i + 1, j + 1) * 3 + 2];
					//volume += (v0z + v1z + v2z) * ((v1x - v0x) * (v2y - v0y) - (v2x - v0x) * (v1y - v0y));
				//}
			//}
		//}
	//}
	//volume /= 6;
	//double volume_truth = 0.0;
	//for (vector<int>::size_type i = 0; i < teapotFaceList.size() / 3; ++i)
	//{
		//int id0 = teapotFaceList[i * 3];
		//int id1 = teapotFaceList[i * 3 + 1];
		//int id2 = teapotFaceList[i * 3 + 2];
		//double v0x = result_truth[id0 * 3];
		//double v0y = result_truth[id0 * 3 + 1];
		//double v0z = result_truth[id0 * 3 + 2];
		//double v1x = result_truth[id1 * 3];
		//double v1y = result_truth[id1 * 3 + 1];
		//double v1z = result_truth[id1 * 3 + 2];
		//double v2x = result_truth[id2 * 3];
		//double v2y = result_truth[id2 * 3 + 1];
		//double v2z = result_truth[id2 * 3 + 2];
		//volume_truth += (v0z + v1z + v2z) * ((v1x - v0x) * (v2y - v0y) - (v2x - v0x) * (v1y - v0y));
	//}
	//volume_truth /= 6;
	//cout << "eeeeee 近似体积 = " << volume << ", 真实体积 = " << volume_truth << endl;
	//cout << "eeeeee 体积误差 = " << volume - volume_truth << ", 误差率 = " << fabs(volume - volume_truth) / volume_truth << endl;

	/* 法向误差 */
	cymError = cudaMemcpy(result, normalPtrVBO, size3, cudaMemcpyDeviceToHost);
	if (cymError)
		cout << __FILE__ << "第" << __LINE__ << "行, 错误代码" << cymError << ": " << cudaGetErrorString(cymError) << endl;
	cymError = cudaMemcpy(result_truth, normalPtrVBO_truth, size3, cudaMemcpyDeviceToHost);
	if (cymError)
		cout << __FILE__ << "第" << __LINE__ << "行, 错误代码" << cymError << ": " << cudaGetErrorString(cymError) << endl;

	error_ave = 0.0, error_max = 0.0;
	float x0_max, y0_max, z0_max, x1_max, y1_max, z1_max;
	const float PI = 3.14159265358979;
	for (int i = 0; i < samplePointPerTriangle * triangleNum; ++i)
	{
		double x0 = result[i * 3];
		double y0 = result[i * 3 + 1];
		double z0 = result[i * 3 + 2];
		//int real_idx = my_to_truth_table[i];
		//double x1 = result_truth[real_idx * 3];
		//double y1 = result_truth[real_idx * 3 + 1];
		//double z1 = result_truth[real_idx * 3 + 2];
		double x1 = result_truth[i * 3];
		double y1 = result_truth[i * 3 + 1];
		double z1 = result_truth[i * 3 + 2];
		//printf("(%f, %f, %f) = (%f, %f, %f)\n", x0, y0, z0, x1, y1, z1);
		double length = sqrt(x0 * x0 + y0 * y0 + z0 * z0);
		if (fabs(length) < 0.00001)
			cout << "长度异常：" << length << endl;
		x0 /= length; y0 /= length; z0 /= length;
		length = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
		if (fabs(length) < 0.00001)
			cout << "长度异常2：" << length << endl;
		x1 /= length; y1 /= length; z1 /= length;
		//cout << "ori = " << x0 << ", " << y0 << ", " << z0 << "\t"
			 //<< "deform = " << x1 << ", " << y1 << ", " << z1 << endl;
		double error = sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1) + (z0 - z1) * (z0 - z1));
		//cout << "before error * 0.5 = " << error * 0.5 << "\t";
		error = 2 * asin(error * 0.5);
		//error = 1 * asin(error / 1);
		error_ave += error;
		//cout << "after error = " << error << endl;
		if (error_max < error)
		{
			error_max = error;
			x0_max = x0;
			y0_max = y0;
			z0_max = z0;
			x1_max = x1;
			y1_max = y1;
			z1_max = z1;
		}
		//double normal_diff = color_map_normal(VertexCoord(x0, y0, z0), VertexCoord(x1, y1, z1), PI / 3);
		double normal_diff = color_map_normal(VertexCoord(x0, y0, z0), VertexCoord(x1, y1, z1), PI / 30);
		texture_coord[i * 3] = normal_diff;
		texture_coord[i * 3 + 1] = 0.5;
	}
	cymError = cudaGetLastError();
	if (cymError)
		cout << __FILE__ << "第" << __LINE__ << "行, 错误代码" << cymError << ": " << cudaGetErrorString(cymError) << endl;
	cudaMemcpy(texCoord3DPtrVBO, texture_coord, size3, cudaMemcpyHostToDevice);
	cymError = cudaGetLastError();
	if (cymError)
		cout << __FILE__ << "第" << __LINE__ << "行, 错误代码" << cymError << ": " << cudaGetErrorString(cymError) << endl;
	if (error_ave > normal_error_ave_max)
		normal_error_ave_max = error_ave;
	if (error_max > normal_error_max_max)
		normal_error_max_max = error_max;
	cout << "max0 = " << x0_max << ", " << y0_max << ", " << z0_max;
	cout << "\tmax1 = " << x1_max << ", " << y1_max << ", " << z1_max << endl;
	cout << "总法向误差 = " << normal_error_ave_max << endl;
	cout << "eeeeee 平均法向误差（角度） = " << normal_error_ave_max / (samplePointPerTriangle * triangleNum) * 180 / PI
		 << ", 最大法向误差（角度） = " << normal_error_max_max * 180 / PI << endl << endl;
	cout << "error_ave = " << error_ave / (samplePointPerTriangle * triangleNum) * 180 / PI << endl;

	cudaGraphicsUnmapResources(1, &normalVBO_CUDA, 0);
	cudaGraphicsUnmapResources(1, &vertexVBO_CUDA, 0);
	cudaGraphicsUnmapResources(1, &texCoordVBO_CUDA, 0);
	cudaGraphicsUnmapResources(1, &texCoord3DVBO_CUDA, 0);
	cudaGraphicsUnmapResources(1, &normalVBO_CUDA_truth, 0);
	cudaGraphicsUnmapResources(1, &vertexVBO_CUDA_truth, 0);
	cymError = cudaGetLastError();
	if (cymError)
		cout << __FILE__ << "第" << __LINE__ << "行, 错误代码" << cymError << ": " << cudaGetErrorString(cymError) << endl;

	delete []result;
	delete []result_truth;
	/*---------------------- 测量误差完成 --------------------------*/
}
//#endif

/************************************************************************************************************/

void setGLDevice()
{
	cudaGLSetGLDevice(0);
}

/* 使用缓冲区对象进行 cuda 和 OpenGL 协同工作之前，需要进行一些初始化 */
void regGLBuffer()
{
	printCudaError(__FILE__, __FUNCTION__, __LINE__);
	if (registered)
	{
		cudaGraphicsUnregisterResource(normalVBO_CUDA);
		cudaGraphicsUnregisterResource(texCoordVBO_CUDA);
		cudaGraphicsUnregisterResource(texCoord3DVBO_CUDA);
		cudaGraphicsUnregisterResource(vertexVBO_CUDA);
#ifdef LINE
		cudaGraphicsUnregisterResource(baryVBO_CUDA);
		cudaGraphicsUnregisterResource(oriBaryVBO_CUDA);
#endif

//#ifdef TRUTH
		cudaGraphicsUnregisterResource(normalVBO_CUDA_truth);
		cudaGraphicsUnregisterResource(vertexVBO_CUDA_truth);
//#endif
		registered = false;
	}
	printCudaError(__FILE__, __FUNCTION__, __LINE__);
	cudaGraphicsGLRegisterBuffer(&normalVBO_CUDA, normalVBO, cudaGraphicsMapFlagsWriteDiscard);
	printCudaError(__FILE__, __FUNCTION__, __LINE__);
	cudaGraphicsGLRegisterBuffer(&texCoordVBO_CUDA, texCoordVBO, cudaGraphicsMapFlagsWriteDiscard);
	cudaGraphicsGLRegisterBuffer(&texCoord3DVBO_CUDA, texCoord3DVBO, cudaGraphicsMapFlagsWriteDiscard);
	cudaGraphicsGLRegisterBuffer(&vertexVBO_CUDA, vertexVBO, cudaGraphicsMapFlagsWriteDiscard);
#ifdef LINE
	cudaGraphicsGLRegisterBuffer(&baryVBO_CUDA, baryVBO, cudaGraphicsMapFlagsWriteDiscard);
	cudaGraphicsGLRegisterBuffer(&oriBaryVBO_CUDA, oriBaryVBO, cudaGraphicsMapFlagsWriteDiscard);
#endif
//#ifdef TRUTH
	cudaGraphicsGLRegisterBuffer(&normalVBO_CUDA_truth, normalVBO_truth, cudaGraphicsMapFlagsWriteDiscard);
	cudaGraphicsGLRegisterBuffer(&vertexVBO_CUDA_truth, vertexVBO_truth, cudaGraphicsMapFlagsWriteDiscard);
//#endif
	registered = true;
	printCudaError(__FILE__, __FUNCTION__, __LINE__);
}

/************************************************************************************************************/

void cudaFreeNonZero(void **ptr)
{
	if (*ptr)
	{
		cudaFree(*ptr);
		*ptr = 0;
	}
}

void freeTessMemD()
{
	cudaFreeNonZero((void**)&BqD);
	cudaFreeNonZero((void**)&BqD_PN);
	cudaFreeNonZero((void**)&RD);
	cudaFreeNonZero((void**)&my_to_truth_tableD);
	cudaFreeNonZero((void**)&parameter3D);
	cudaFreeNonZero((void**)&parameterND);
	delete []my_to_truth_table;
#ifdef TRUTH
	cudaFreeNonZero((void**)&BqD_truth);
	cudaFreeNonZero((void**)&BBD_truth);
	cudaFreeNonZero((void**)&RD_truth);
#endif
}

void freeModelMemD()
{
	cudaFreeNonZero((void**)&vertexParamListD);
	cudaFreeNonZero((void**)&vertexCoordListD);
	//cudaFreeNonZero((void**)&vertexParamListD_teapot);
	//cudaFreeNonZero((void**)&vertexCoordListD_teapot);

	cudaFreeNonZero((void**)&triangleListD);
	cudaFreeNonZero((void**)&sampleValueD);
	cudaFreeNonZero((void**)&sampleValueD_PN);
	cudaFreeNonZero((void**)&triangleCtrlPointD);
	cudaFreeNonZero((void**)&triangleCtrlPointD_PN);
	cudaFreeNonZero((void**)&triangleNormalCtrlPointD_PN);
	cudaFreeNonZero((void**)&triangle_adjacent_tableD);
#ifdef TRUTH
	cudaFreeNonZero((void**)&sampleValueD_truth);
	cudaFreeNonZero((void**)&B_1D_truth);
#endif

	degreeMemD = 0;
	modelMemD = 0;

#ifdef DRAW_TRIANGULAR_CTRL_POINTS
	delete []triangular_ctrl_points;
#endif

	freeTessMemD();
}

void freeMemD()
{
	if (registered)
	{
		cudaGraphicsUnregisterResource(normalVBO_CUDA);
		cudaGraphicsUnregisterResource(texCoordVBO_CUDA);
		cudaGraphicsUnregisterResource(texCoord3DVBO_CUDA);
		cudaGraphicsUnregisterResource(vertexVBO_CUDA);
#ifdef LINE
		cudaGraphicsUnregisterResource(baryVBO_CUDA);
		cudaGraphicsUnregisterResource(oriBaryVBO_CUDA);
#endif
//#ifdef TRUTH
		cudaGraphicsUnregisterResource(normalVBO_CUDA_truth);
		cudaGraphicsUnregisterResource(vertexVBO_CUDA_truth);
//#endif
		registered = false;
	}
	if (cublas_handle)
	{
		cublasDestroy(cublas_handle);
	}
	cudaFreeNonZero((void**)&matrixFittingIdxD);
	cudaFreeNonZero((void**)&matrixFittingD);
	permanentMemD = 0;
	freeModelMemD();
}
