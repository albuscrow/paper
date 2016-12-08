#ifndef __COMMONDATA_H__
#define __COMMONDATA_H__

#include <fstream>
#include <cstring>
#include <list>
#include <lapackpp/gmd.h>
#include "obj_data.h"
#include "matrix_stack.h"

#define NORMALIZE_TO_1	// 模型读进来之后单位化到[-1,1]x[-1,1]x[-1,1]区间
//#define MORPH			// 表示进行变形动画，因此不进行一些时间的打印
/* 
 * 注意：TRUTH模式和LINE模式不能同时使用
 * 如果同时define，只能得到其中之一的效果
 */
//#define TRUTH
#define LINE

//#define RE_LENGTH		// 对边控制顶点的调整要不要保持长度

//#define DRAW_TRIANGULAR_CTRL_POINTS			// 绘制三角形Bezier曲面片控制定点，仅调试用

#define IS_BORDER_INFO		// 有时需要根据一条边是不是边界进行调整

/*---------------------------------------------------------*/

/*
 * 参数方向
 */
enum Direction
{
	U = 0,
	V = 1,
	W = 2
};

/*---------------------------------------------------------*/

/*
 * 坐标轴方向
 */
enum Coord
{
	X = 0,
	Y = 1,
	Z = 2
};

/*---------------------------------------------------------*/

/*
 * 使用的变形算法
 */
enum Algorithm
{
	FFD,
	AFFD
};

/*---------------------------------------------------------*/

/*
 * 控制顶点的编辑方向（可能有多个方向）
 */
enum EditDirection
{
	NONE_AXIS = 0,
	X_AXIS = 1,
	Y_AXIS = 2,
	Z_AXIS = 4
};

/*---------------------------------------------------------*/

enum AlgorithmType
{
	CYM,
	PN_CUTTING,
	PN_NO_CUTTING
};

/*---------------------------------------------------------*/

/*
 * 控制顶点类
 */
class CtrlPoint : public objdata::VertexCoord		// 控制顶点类
{
protected:
	bool m_bSelected;					// 是否被选中
public:
	CtrlPoint(double x = 0.0, double y = 0.0, double z = 0.0, bool selected = false) : VertexCoord(x, y, z)
	{m_bSelected = selected;}
	bool selected() const {return m_bSelected;}
	void setSelect(bool selected) {m_bSelected = selected;}
};

/*---------------------------------------------------------*/

/*
 * 直接编辑点类
 */
class DirectPoint : public CtrlPoint		// 控制顶点类
{
	double m_fU, m_fV, m_fW;
public:
	DirectPoint(const VertexCoord &p);
	DirectPoint(double u, double v, double w, double x, double y, double z);
	double u() const {return m_fU;}
	double v() const {return m_fV;}
	double w() const {return m_fW;}
	void setU(double u) {m_fU = u;}
	void setV(double v) {m_fV = v;}
	void setW(double w) {m_fW = w;}
	friend std::ostream& operator<<(std::ostream &, const DirectPoint &);
};

/*---------------------------------------------------------*/

template <int degree>
class BezierTriangle
{
	objdata::VertexCoord ctrl_points[(degree + 2) * (degree + 1) / 2];
public:
	objdata::VertexCoord &operator[](int i) { return ctrl_points[i]; }
	objdata::VertexCoord calcValue(const objdata::VertexCoord &baryCoord);
};

/*---------------------------------------------------------*/

struct EdgeWithNormal
{
	int v_max, n_count, n_min[2], n_max[2], face_id[2];
#ifdef IS_BORDER_INFO
	bool border;
#endif
	EdgeWithNormal();
	EdgeWithNormal(int vmax, int nmin, int nmax, int faceid);
	void setXthNormal(int edge_id, int nmin, int nmax, int faceid);
};

std::ostream &operator<<(std::ostream &out, const EdgeWithNormal &edge);

/*---------------------------------------------------------*/

/*
 * 裁剪后的多边形类
 * 每个对象代表一个多边形
 */
class TrimmedPolygon
{
public:
	enum OnFace {space = 0, x_min = 1, y_min = 2, z_min = 4, x_max = 8,
				 y_max = 16, z_max = 32};
	TrimmedPolygon(bool useTexture, int mtlIdx, const objdata::VertexCoord *v_origin, int origin_face_idx);
	objdata::VertexCoord &operator[](int idx) {return vertexList[idx];}
	objdata::NormalCoord &getNormal(int idx) {return normalList[idx];}
	objdata::NormalCoord &getNormalAdj(int idx) {return normal_adj_list[idx];}
	objdata::TextureCoord &getTexture(int idx) {return textureList[idx];}
#ifdef LINE
	objdata::VertexCoord &getBary(int idx) {return bary_origin_list[idx];}
#endif
	int normalCount(int id) { return normal_count_list[id]; }
#ifdef IS_BORDER_INFO
	bool border(int id) { return is_border_list[id]; }
#endif
	bool adjust(int id) { return adjust_list[id]; }
	int size() const {return vertexList.size();}
	int texSize() const {return textureList.size();}
	void push_back(const objdata::VertexCoord &v) {vertexList.push_back(v);}
	void push_back_n(const objdata::NormalCoord &n) {normalList.push_back(n);}
	void push_back_n_adj(const objdata::NormalCoord &n) {normal_adj_list.push_back(n);}
#ifdef LINE
	void push_back_bary(const objdata::VertexCoord &v) { bary_origin_list.push_back(v); }
#endif
	void push_back_texture(const objdata::TextureCoord &vt) {textureList.push_back(vt);}
	//void clear()
	//{
		//vertexList.clear();
		//normalList.clear();
		//textureList.clear();
		//normal_count_list.clear();
		//adjust_list.clear();
	//#ifdef LINE
		//bary_origin_list.clear();
	//#endif
	//}
	bool useTexture() const {return m_bUseTexture;}
	int mtlIdx() const {return m_nMtlIdx;}
	objdata::VertexCoord v_origin[3];
	int origin_face_idx_;
	void push_back_n_count(int nc) {normal_count_list.push_back(nc);}
#ifdef IS_BORDER_INFO
	void push_back_border(int b) { is_border_list.push_back(b); }
#endif
	void push_back_adjust(bool ad) { adjust_list.push_back(ad); }
	void setNormalAdj(int idx, const objdata::NormalCoord &n) {normal_adj_list[idx] = n;}
	void normalizeAllNormals();
private:
	std::vector<objdata::VertexCoord> vertexList;		// 顶点列表
	std::vector<objdata::NormalCoord> normalList;		// 法向列表
	std::vector<objdata::NormalCoord> normal_adj_list;	// 此法向用于调整控制顶点
	std::vector<objdata::TextureCoord> textureList;

	// normal_count_list[i]表示vertexList[i-1]到vertexList[i]那条边的法向数量
	std::vector<int> normal_count_list;
#ifdef IS_BORDER_INFO
	// is_border_list[i]表示vertexList[i-1]到vertexList[i]那条边是否是边界
	std::vector<bool> is_border_list;
#endif
	// adjust_list[i]表示vertexList[i]是否进行第一阶段调整
	std::vector<bool> adjust_list;
#ifdef LINE
	std::vector<objdata::VertexCoord> bary_origin_list;
#endif
	//std::vector<int> on_face_list;
	bool m_bUseTexture;
	int m_nMtlIdx;
};

/*---------------------------------------------------------*/

class Triangle
{
public:
	objdata::VertexCoord v[3], v_origin[3];
	objdata::NormalCoord n[3];
	objdata::NormalCoord n_adj[3];
	objdata::TextureCoord vt[3];
#ifdef LINE
	objdata::VertexCoord bary_origin[3];
#endif
	int n_count[3];	// nc0, nc1, nc2分别代表v2v0, v0v1, v1v2边的法向数量
	int origin_face_idx_;
	Triangle(const objdata::VertexCoord &v0, const objdata::VertexCoord &v1, const objdata::VertexCoord &v2,
			 const objdata::NormalCoord &n0, const objdata::NormalCoord &n1, const objdata::NormalCoord &n2,
			 const objdata::NormalCoord &n_adj0, const objdata::NormalCoord &n_adj1, const objdata::NormalCoord &n_adj2,
#ifdef LINE
			 const objdata::VertexCoord &bary_origin0, const objdata::VertexCoord &bary_origin1, const objdata::VertexCoord &bary_origin2,
#endif
			 const objdata::VertexCoord *v_origin, int origin_face_idx,
			 int nc0, int nc1, int nc2,
			 const objdata::TextureCoord &vt0 = 0.0f, const objdata::TextureCoord &vt1 = 0.0f, const objdata::TextureCoord &vt2 = 0.0f);
};

typedef objdata::TextureCoord VertexParam;

/*---------------------------------------------------------*/

/*
 * 数据类
 * 存储各个类共用的数据
 * 完成核心算法
 */
class CommonData
{
	static const int MAXVERTEX  = 100000;		// 预估的最大顶点数
	static const int MAXFACE = 10000;			// 预估的最大面片数
	static const int MAXCTRLPOINT = 15;			// B 样条体各方向最大控制顶点数
	static const int MAXORDER = 4;				// B 样条体各方向最大阶数
	static const int MINORDER = 2;				// B 样条体各方向最小阶数
	static const int MAXKNOTINTERVAL = MAXCTRLPOINT - MINORDER + 1;		// 最大节点区间数
	static const int MAXKNOT = MAXCTRLPOINT + MAXORDER;					// 最大节点数
	static const int MAXDIRECTPOINT = 20;		// 直接编辑点的最大数量
	static const int MAXANIMATIONSTEP = 1000;	// 关键帧最大数量
	static const double ZERO;					// 一个很小的正数
	static const double EXPAND;					// 防止调整侧影轮廓线后某些点跑出侧影轮廓线，所以进行一些扩展
	static const double PI;
	//static const double TAN_FACTOR;

	objdata::ObjData *objData;
	double m_fMin[3], m_fMax[3];				// X, Y, Z 三个方向坐标的最大最小值
	double m_fLength[3];						// X, Y, Z 三个方向跨度的一半
	double m_fTotalLength;						// 上面三者平方和的开方（简称模型半径）
	int m_nOrder[3];							// B 样条体三个方向的阶数
	int m_nCtrlPointCount[3];					// B 样条体三个方向控制顶点的数量
	int m_nKnotIntervalCount[3];				// B 样条体三个方向节点区间的数量
	//float m_fTrimmedBezierKnotPoint[10][20];	// 各种次数 Bézier 曲线的节点向量
	int m_nSamplePointCount;					// 三角化时每个参数方向上的采样点数量

	std::vector<VertexParam> vertexParamList;	// 调整好之后的顶点参数坐标
	std::vector<objdata::VertexCoord> vertexCoordList;	// 调整好之后的顶点几何坐标
	std::vector<objdata::Face> faceList;					// 调整好后的面片信息
	std::vector<BezierTriangle<3> > cubicBezierTriangleList;	// faceList中每个片面对应的三次三角Bezier曲面片，用PN Triangle方法生成。这个数据结构仅用于计算对应点，无需传给CudaCalc.cu
	std::vector<BezierTriangle<2> > quadraticBezierTriangleList;	// faceList中每个片面的法向对应的二次三角Bezier曲面片，用PN Triangle方法生成。这个数据结构仅用于计算对应点，无需传给CudaCalc.cu

	double knotList[3][MAXKNOT];				// B 样条体三个方向的节点向量
	/* B 样条体的控制顶点 */
	CtrlPoint ctrlPoint[MAXCTRLPOINT][MAXCTRLPOINT][MAXCTRLPOINT];
	//CtrlPoint sourceCtrlPoint[MAXCTRLPOINT][MAXCTRLPOINT][MAXCTRLPOINT];
	//CtrlPoint targetCtrlPoint[MAXCTRLPOINT][MAXCTRLPOINT][MAXCTRLPOINT];
	CtrlPoint animationList[MAXANIMATIONSTEP][MAXCTRLPOINT][MAXCTRLPOINT][MAXCTRLPOINT];
	int animationTimeList[MAXANIMATIONSTEP];
	int animationStep;
	std::vector<DirectPoint> directPointVector;
	/* 每一节点盒当中裁剪后的多边形 */
	std::vector<TrimmedPolygon> trimmedPolygonList[MAXKNOTINTERVAL][MAXKNOTINTERVAL][MAXKNOTINTERVAL];
	std::vector<Triangle> triangleList;

	Algorithm m_eAlgorithm;						// 使用的变形算法
	int m_nFFDThreadCount;						// 使用 GPU 加速 FFD 算法时线程块大小
	bool m_bGPU;								// 是否使用 GPU 加速
	bool m_bLoaded;								// 是否已载入模型
	bool m_bFirstLoad;
	bool m_bLoadingEdit;
	bool adjust_silhouette_, adjust_split_points_, use_pn_;
	AlgorithmType algorithm_type_;
	int m_nTessPointCount;		// 模型三角化后，生成的三角化点的数量
	std::vector<int *> faceMtlList;	// 模型经过三角化后，所有的三角形顶点统一编号，而所有的面片按照材质的不同放入不同的列表
	std::vector<int> faceMtlCountList;	// 上面各个列表的大小
	std::list<EdgeWithNormal> *edge_table_;
	EdgeWithNormal edge_illegal_;
	std::list<int> *face_adjacent_table_;
	std::list<int> *face_children_table_;
	int *triangle_adjacent_table_;
	int findNormalCount(int v0, int v1);
#ifdef IS_BORDER_INFO
	int findIsBorder(int v0, int v1);
#endif
	const EdgeWithNormal &findEdge(int v0, int v1);
	objdata::VertexCoord BSplineVolumeValueMatrix(double u, double v, double w, int leftUIdx, int leftVIdx, int leftWIdx);
	void calcKnotCtrlPoint();
	void directLinearSolve(const LaGenMatDouble &D);
	void generateEdge();
	void PN();
	objdata::VertexCoord calcBaryCoord(const objdata::VertexCoord &v0, const objdata::VertexCoord &v1, const objdata::VertexCoord &v2, const objdata::VertexCoord &v);
	//objdata::VertexCoord bezierTriangleValue2(const objdata::VertexCoord &baryCoord, int n, const objdata::VertexCoord *ctrl_points);
	//void bezierTriangleValue(const objdata::VertexCoord &baryCoord, int n, const objdata::VertexCoord *ctrl_points, objdata::VertexCoord &value, objdata::NormalCoord &normal);
	void clipPolygon();
	int clipPolygonAgainstPlane(TrimmedPolygon &polygon, double a, double b, double c, double d, bool isPositive);
	void triangulatePolygon();
	void triangulatePolygon_PN_NO_CUTTING();
	/*
	 * 简单三角剖分的坐标转换，将整个三角形看成一座楼房
	 * floor 代表层数，最上面是第0层，向下依次是第1层、第2层……
	 * room 代表房间号，每层最左面是第0间，往右依次是第1间、第2间……
	 * triangleCoord 返回的是每个三角化点的一维编号，编号从0开始，
	 * 从上到下，从左到右，遍历所有三角化点
	 */
	int triangleCoord(int floor, int room) const
	{return (1 + floor) * floor / 2 + room;}
	double BSplineBase(int i, int k, double t, const Direction direction);
	int triangularAdjacent(int idx0, int idx1);
public:
	std::vector<matrix_stack::Matrix4x4> cubemap_matrix_list_;
	int morphStep, sourceIdx, targetIdx;
	CommonData();
	~CommonData();
	void setAlgorithmType();
	AlgorithmType getAlgorithmType() { return algorithm_type_; }
	void increaseCenterCtrlPointFactor();
	void decreaseCenterCtrlPointFactor();
	void setUsePN();
	bool usePN() { return use_pn_; }
	void testRotate(bool clockwise);
	void loadObj(const char *fileName);
	/* 返回模型半径 */
	double length() const {return m_fTotalLength;}
	double minX() const {return m_fMin[0];}
	double minY() const {return m_fMin[1];}
	double minZ() const {return m_fMin[2];}
	double maxX() const {return m_fMax[0];}
	double maxY() const {return m_fMax[1];}
	double maxZ() const {return m_fMax[2];}
	/* 返回模型的顶点数量 */
	unsigned int vertexCount() const {return vertexParamList.size();}
	/* 返回模型面片数量 */
	unsigned int faceCount() const {return faceList.size();}
	/* 返回第 i 个面片的顶点数量 */
	unsigned int vertexCoordIdxCount(unsigned int i) const
	{return faceList[i].vertexCoordIndex.size();}
	/* 返回第 i 个面片第 j 个顶点的编号 */
	int vertexCoordIdx(unsigned int i, unsigned int j) const
	{return faceList[i].vertexCoordIndex[j];}
	int textureVertexIdx(unsigned int i, unsigned int j) const
	{return faceList[i].textureCoordIndex[j];}
	/* 返回编号为 idx 的顶点 */
	objdata::VertexCoord vertexCoord(int idx) const
	{return vertexCoordList[idx];}
	objdata::TextureCoord textureVertex(int idx) const
	{return objData->textureCoordList[idx];}
	/* 将编号为 idx 的顶点的坐标设置为 (x, y, z) */
	void setVertexCoord(int idx, float x, float y, float z)
	{
		vertexCoordList[idx].x(x);
		vertexCoordList[idx].y(y);
		vertexCoordList[idx].z(z);
	}
	//{vertexCoordList[idx] = objdata::VertexCoord((double)x, (double)y, (double)z);}

	/* 返回编号为 idx 的顶点的参数坐标 */
	VertexParam vertexParam(int idx) const {return vertexParamList[idx];}

	/* 返回 B 样条体 direction 方向（可能为 U, V, W）的阶数 */
	int order(int direction) const {return m_nOrder[direction];}
	/* 设置 B 样条体 direction 方向（可能为 U, V, W）的阶数 */
	void setOrder(int order, int direction) {m_nOrder[direction] = order;}
	/* 返回 B 样条体 direction 方向（可能为U, V, W）的控制顶点数 */
	int ctrlPointCount(int direction) const
	{return m_nCtrlPointCount[direction];}
	/* 设置 B 样条体 direction 方向（可能为U, V, W）的控制顶点数 */
	void setCtrlPointCount(int ctrlPointCount, int direction)
	{m_nCtrlPointCount[direction] = ctrlPointCount;}
	void cancelAllSelection();
	/* 返回相应控制顶点，i, j, k 分别为 U, V, W 方向的序号 */
	CtrlPoint getCtrlPoint(int i, int j, int k) const
	{return ctrlPoint[i][j][k];}
	/* 返回相应节点，i 表示（U, V, W）方向，j 表示序号*/
	double getKnot(int i, int j) const {return knotList[i][j];}
	void setCtrlPointSelected(int i, int j, int k);
	objdata::VertexCoord fromParamToCoord(double u, double v, double w);
	//void fromParamToCoordCopy();
	void setAllDirectPoint(std::vector<DirectPoint> &result)
	{directPointVector = result;}
	void setAllCtrlPoint(std::vector<CtrlPoint> &result);
	void setTargetCtrlPoint(std::vector<std::vector<CtrlPoint> > &, std::vector<int> &, std::vector<matrix_stack::Matrix4x4> &);
	void morph();
	void ctrlPointTranslate(double delta, int editDirection);
	void ctrlPointRotate(double theta, int editDirection);
	void ctrlPointScale(double factor, int editDirection);
	void resetFirstLoad() {m_bFirstLoad = true;}
	void setLoadingEdit(bool loading) {m_bLoadingEdit = loading;}
	bool loadingEdit() const {return m_bLoadingEdit;}
	void directPointTranslate(double delta, int editDirection);
	void directPointRotate(double theta, int editDirection);
	void directPointScale(double factor, int editDirection);
	void calcFinalResult();
	/* 返回三角化时每个参数方向上的采样点数量 */
	int samplePointCount() const {return m_nSamplePointCount;}
	void setSamplePointCount(int count);
	/* 返回 GPU 加速 FFD 算法时线程块大小 */
	int ffdThreadCount() const {return m_nFFDThreadCount;}

	/* 返回当前使用的变形算法（FFD 或 AFFD） */
	Algorithm algorithm() const {return m_eAlgorithm;}
	void setAlgorithm(const Algorithm algo);
	/* 是否使用 GPU 计算 */
	bool isGPU() const {return m_bGPU;};
	/* 设置计算工具（CPU 或者 GPU） */
	void setCalcDevice(bool isGPU) {m_bGPU = isGPU;};

	void preCalc(bool reset_ctrl_point = true);
	void execute();
	int objMtlTexListSize() const {return objData->mtlTexList.size();}
	//int *faceMtl(const int i) {return faceMtlList[i];}
	int *faceMtl(int i) {return faceMtlList[i];}
	int faceMtlListSize() const {return faceMtlList.size();}
	int faceMtlCount(int i) const {return faceMtlCountList[i];}
	objdata::MtlTex getMtl(int i) const;
	/* 返回d方向（可能为U, V, W）的节点区间数 */
	int knotIntervalCount(int d) const {return m_nKnotIntervalCount[d];}
	/* 返回次数为 degree 的 Bézier 曲线的节点向量 */
	//float *trimmedBezierKnotPoint(int degree)
	//{return &m_fTrimmedBezierKnotPoint[degree][0];}

	int tessPointCount() const {return m_nTessPointCount;}

	void newSamplePointTesslate();
	void callTesslateD();

	int directPointVectorSize() const {return directPointVector.size();}
	const DirectPoint &directPoint(int idx) const
	{return directPointVector[idx];}
	void directPointVectorPush(DirectPoint &p);
	//{directPointVector.push_back(p);}
	//bool directPointVectorEmpty() const {return directPointVector.empty();}
	void setDirectPointSelected(int idx);
	int directPointSelectedCount();
	void saveEdit(std::ofstream &fout);
	void testMove(bool out);
	void testMoveTorus(bool out);

	bool adjustSilhouette() const { return adjust_silhouette_; }
	bool adjustSplitPoints() const { return adjust_split_points_; }
	void setAdjustSilhouette(bool adjust) {adjust_silhouette_ = adjust;}
	void setAdjustSplitPoints(bool adjust) {adjust_split_points_ = adjust;}

	//void moveShip(int target, int id0, int id1, int xyz);
};

#endif
