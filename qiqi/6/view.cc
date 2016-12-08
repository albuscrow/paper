#include "widget.h"
#include "view.h"
#include <vector>
#include <fstream>
#include <cmath>
#include <lapackpp/gmd.h>
#include <lapackpp/blas3pp.h>
#include <lapackpp/laslv.h>

#include <iostream>
using std::cout;
using std::endl;

using namespace objdata;

#ifndef CALLBACK
#define CALLBACK
#endif

using std::string;
using std::vector;

extern GLuint normalVBO;				// 结果法向量的缓冲区对象标识符
extern GLuint texCoordVBO;				// 结果纹理坐标的缓冲区对象标识符
extern GLuint texCoord3DVBO;			// 结果纹理三维坐标的缓冲区对象标识符
extern GLuint vertexVBO;				// 结果顶点坐标的缓冲区对象标识符
#ifdef LINE
extern GLuint baryVBO;					// 结果顶点在切割后的三角形中的坐标的缓冲区对象标识符
extern GLuint oriBaryVBO;				// 结果顶点在切割前的三角形中的坐标的缓冲区对象标识符
#endif

//#ifdef TRUTH
extern GLuint normalVBO_truth;			// 结果法向量的缓冲区对象标识符
extern GLuint vertexVBO_truth;			// 结果顶点坐标的缓冲区对象标识符
//#endif



void setGLDevice();

const int View::EYE_Z_FACTOR = 2;
const double View::ZERO = 0.000001;
const double View::PI = 3.14159265358979;
const double View::RIGHTFAC = 0.8;
const double View::BOTTOMFAC = 0.8;

void drawSphere(double r, double ox, double oy, double oz)
{
	const double PI = 3.14159265358979;
	vector<VertexCoord> v;
	vector<NormalCoord> vn;
	v.push_back(VertexCoord(0, 0, 0));
	vn.push_back(NormalCoord(0, 1, 0));
	// 顶端和底端两个点
	//fout << "v 0 " << r * scale_y << " 0\nvn 0 1 0\n"
		 //<< "v 0 " << -r * scale_y << " 0\nvn 0 -1 0" << endl;
	v.push_back(VertexCoord(ox, r + oy, oz));
	vn.push_back(NormalCoord(0, 1, 0));
	v.push_back(VertexCoord(ox, -r + oy, oz));
	vn.push_back(NormalCoord(0, -1, 0));

	//int m = 4, n = 4;
	//int m = 8, n = 8;
	//int m = 12, n = 12;
	int m = 30, n = 30;
	//int m = 40, n = 40;

	// 球上除顶端和底端外所有的顶点
	for (int i = 1; i < m; ++i)
	{
		double phi = PI * static_cast<double>(i) / m;
		for (int j = 0; j < n; ++j)
		{
			double theta = 2 * PI * static_cast<double>(j) / n;
			double x = r * sin(phi) * cos(theta);
			double y = r * cos(phi);
			double z = r * sin(phi) * sin(theta);
			//fout << "v " << x * scale_x << " " << y * scale_y << " " << z * scale_z << endl;
			//fout << "vn " << x * scale_x_1 << " " << y * scale_y_1 << " " << z * scale_z_1 << endl;
			v.push_back(VertexCoord(x + ox, y + oy, z + oz));
			NormalCoord n(x, y, z);
			n.normalize();
			vn.push_back(n);
		}
	}

	// 球顶端的一圈面片
	for (int j = 0; j < n; ++j)
	{
		int next_idx = j + 1 < n ? j + 1 : 0;
		//fout << "f 1//1 " << next_idx + 3 << "//" << next_idx + 3 << " "
			 //<< j + 3 << "//" << j + 3 << " " << endl;
		glBegin(GL_TRIANGLES);

		glNormal3f(vn[1].i(), vn[1].j(), vn[1].k());
		//cout << "vn1.x = " << vn[1].i() << ", vn1.y = " << vn[1].j()
			 //<< ", vn1.z = " << vn[1].k() << endl;
		glVertex3f(v[1].x(), v[1].y(), v[1].z());

		glNormal3f(vn[next_idx + 3].i(), vn[next_idx + 3].j(), vn[next_idx + 3].k());
		glVertex3f(v[next_idx + 3].x(), v[next_idx + 3].y(), v[next_idx + 3].z());

		glNormal3f(vn[j + 3].i(), vn[j + 3].j(), vn[j + 3].k());
		glVertex3f(v[j + 3].x(), v[j + 3].y(), v[j + 3].z());

		glEnd();
	}

	// 球底端的一圈面片
	int base = 2 + n * (m - 2);
	for (int j = 0; j < n; ++j)
	{
		int next_idx = j + 1 < n ? j + 1 : 0;
		//fout << "f 2//2 " << base + j + 1 << "//" << base + j + 1 << " "
			 //<< base + next_idx + 1 << "//" << base + next_idx + 1 << endl;
		glBegin(GL_TRIANGLES);

		glNormal3f(vn[2].i(), vn[2].j(), vn[2].k());
		glVertex3f(v[2].x(), v[2].y(), v[2].z());

		glNormal3f(vn[base + j + 1].i(), vn[base + j + 1].j(), vn[base + j + 1].k());
		glVertex3f(v[base + j + 1].x(), v[base + j + 1].y(), v[base + j + 1].z());

		glNormal3f(vn[base + next_idx + 1].i(), vn[base + next_idx + 1].j(), vn[base + next_idx + 1].k());
		glVertex3f(v[base + next_idx + 1].x(), v[base + next_idx + 1].y(), v[base + next_idx + 1].z());

		glEnd();
	}

	// 球中间的部分
	for (int i = 1; i <= m - 2; ++i)
	//for (int i = 1; i < 5; ++i)
	{
		int base_upper = 2 + n * (i - 1);
		int base_lower = 2 + n * i;
		for (int j = 0; j < n; ++j)
		{
			int j_plus_1 = j + 1 < n ? j + 1 : 0;
			int upperleft = base_upper + j_plus_1 + 1;
			int upperright = base_upper + j + 1;
			int lowerleft = base_lower + j_plus_1 + 1;
			int lowerright = base_lower + j + 1;
			//fout << "f " << upperright << "//" << upperright << " "
						 //<< upperleft << "//" << upperleft << " "
						 //<< lowerright << "//" << lowerright << "\n"
				 //<< "f " << lowerright << "//" << lowerright << " "
						 //<< upperleft << "//" << upperleft << " "
						 //<< lowerleft << "//" << lowerleft << endl;
			glBegin(GL_TRIANGLES);
			glNormal3f(vn[upperright].i(), vn[upperright].j(), vn[upperright].k());
			glVertex3f(v[upperright].x(), v[upperright].y(), v[upperright].z());

			glNormal3f(vn[upperleft].i(), vn[upperleft].j(), vn[upperleft].k());
			glVertex3f(v[upperleft].x(), v[upperleft].y(), v[upperleft].z());

			glNormal3f(vn[lowerright].i(), vn[lowerright].j(), vn[lowerright].k());
			glVertex3f(v[lowerright].x(), v[lowerright].y(), v[lowerright].z());
			glEnd();

			glBegin(GL_TRIANGLES);
			glNormal3f(vn[lowerright].i(), vn[lowerright].j(), vn[lowerright].k());
			glVertex3f(v[lowerright].x(), v[lowerright].y(), v[lowerright].z());

			glNormal3f(vn[upperleft].i(), vn[upperleft].j(), vn[upperleft].k());
			glVertex3f(v[upperleft].x(), v[upperleft].y(), v[upperleft].z());

			glNormal3f(vn[lowerleft].i(), vn[lowerleft].j(), vn[lowerleft].k());
			glVertex3f(v[lowerleft].x(), v[lowerleft].y(), v[lowerleft].z());
			glEnd();
		}
	}
}

View::View(QWidget *parent, CommonData *commonData) : QGLWidget(parent)
{
	//fout.open("View.txt");
	m_pParent = reinterpret_cast<Widget *>(parent);

	setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
	setMouseTracking(true);		// 即使不按键也会响应鼠标的移动，如果影响动画效率，可以关掉，但此时就无法进行编辑了

	m_bLocalViewer = false;
	//m_fDist = 0.0f;

	m_bLoaded = false;								// 初始化为“没有载入obj文件”
	m_bIsFlat = false;								// 默认绘制模式为 GL_SMOOTH
	//m_bIsFlat = true;								// 测使用
	m_bIsOrtho = true;								// 默认使用平行投影
	//m_bIsOrtho = false;								// 默认使用平行投影
	m_bCtrlPressed = false;							// 默认没有按下Ctrl
	m_bShiftPressed = false;						// 默认没有按下Shift
	m_eMouseButton = NO_BUTTON;						// 默认没有鼠标键按下
	m_nCtrlSelectStartX = m_nCtrlSelectStartY = 0;	// 选取框开始坐标
	m_nCtrlSelectEndX = m_nCtrlSelectEndY = 0;		// 选取框结束坐标
	m_nDisplayMode = DISPLAY_MODEL;					// 默认显示模型
	//m_nDisplayMode = 0;								// 测试
	//m_nDisplayMode |= DISPLAY_CTRL_POINTS;			// 默认显示控制网
	//m_nDisplayMode |= DISPLAY_KNOT;				// 默认显示节点
	//m_nDisplayMode |= DISPLAY_DIRECT_POINT;			// 默认显示直接编辑点
	//m_nDisplayMode |= DISPLAY_EDIT_ICON;				// 默认显示编辑图标
	//m_nDisplayMode |= DISPLAY_TRIANGULAR_CTRL_POINTS;				// 默认显示编辑图标
	m_nTexCase = 2;
	m_fScale = 1.0;
	m_fDeltaX = m_fDeltaY = 0.0;
	m_bUseCubeMap = false;
	m_bUseEnvMap = false;
//#ifdef TRUTH
	m_eErrorDisplayMode = MESH;
	m_bDisplayTruth = false;
	error_texture_ = true;
	colormap_name = 0;
//#endif
#ifdef LINE
	m_bLine = false;
#endif
	tex2DNameList = 0;
	tex3DName = 0;
	texCubeName = 0;
	m_pCommonData = commonData;
	m_nEditDirection = NONE_AXIS;
	//m_eEditMode = TRANSLATE;
	m_eEditMode = ROTATE;
	m_nSelectedOriginalFace = 0;

	axisConvert123To124[0] = 0;
	axisConvert123To124[1] = 1;
	axisConvert123To124[2] = 2;
	axisConvert123To124[3] = 4;

	light_position[0] = 0.0f;
	light_position[1] = 0.0f;
	light_position[2] = 1.0f;
	light_position[3] = 0.0f;

	setFocusPolicy(Qt::StrongFocus);				// 对键盘事件的响应方式

	setFormat(QGLFormat(QGL::DoubleBuffer | QGL::DepthBuffer));

	/* 不重置 */
	for (int i = 0; i < 3; ++i)
	{
		rotateAxis[i] = 1.0;
		rotateAxisCubeMap[i] = 1.0;
	}
	lastMatrix = matrix_stack::Matrix4x4(1.0f);
	//lastMatrix = matrix_stack::Matrix4x4(1.0f, 2.0f);

	lastMatrixCubeMap = matrix_stack::Matrix4x4(1.0f);

	queryBackBuffer = 0;
	queryFrontBuffer = 1;
#ifdef DRAW_TRIANGULAR_CTRL_POINTS
	//draw_which_ = 22;
	draw_which_ = 23;
	//draw_which_ = 0;
	//draw_which_ = 63;
#endif
}

View::~View()
{
	//fout.close();			// 测试用

	delTextures();
	delBuffers();
	//delGlobalIdxVBOList();
}

void View::genTextures()
{
//#ifdef TRUTH
	glGenTextures(1, &colormap_name);
	cout << "colormap_name_index = "<< colormap_name << endl;
	loadColormap();
//#endif

	/*-----------------------------------------------------------------*/

	tex2DNameList = new GLuint[m_pCommonData->objMtlTexListSize() + 1];
	glGenTextures(m_pCommonData->objMtlTexListSize() + 1, tex2DNameList);

	//cout << "size = " << m_pCommonData->objMtlTexListSize() << endl;
	//for (int i = 0; i < m_pCommonData->objMtlTexListSize() + 1; ++i)
		//cout << "i = " << i << ", texcym = " << tex2DNameList[i] << endl;
}

/* 对模型材质、颜色不同的各部分各生成一个VBO */
void View::genBuffers()
{
	glGenBuffers(1, &normalVBO);
	glGenBuffers(1, &texCoordVBO);
	glGenBuffers(1, &texCoord3DVBO);
	glGenBuffers(1, &vertexVBO);
#ifdef LINE
	glGenBuffers(1, &baryVBO);
	glGenBuffers(1, &oriBaryVBO);
#endif

//#ifdef TRUTH
	glGenBuffers(1, &normalVBO_truth);
	glGenBuffers(1, &vertexVBO_truth);
//#endif

	for (int i = 0; i < m_pCommonData->objMtlTexListSize() + 1; ++i)
	{
		globalIdxVBOList.push_back(0);
		glGenBuffers(1, &globalIdxVBOList[globalIdxVBOList.size() - 1]);
	}

	glGenBuffers(1, &teapotIdxVBO);
}

void printMemD(const char *file, const char *function, int line, int memSize, string info);
extern int viewMemD;

void View::delTextures()
{
//#ifdef TRUTH
	if (colormap_name)
		glDeleteTextures(1, &colormap_name);
//#endif
	if(tex2DNameList)
	{
		glDeleteTextures(globalIdxVBOList.size(), tex2DNameList);
		delete[] tex2DNameList;
		tex2DNameList = 0;
	}
}

void View::delBuffers()
{
	glBindBuffer(GL_ARRAY_BUFFER, normalVBO);
	glDeleteBuffers(1, &normalVBO);
	glBindBuffer(GL_ARRAY_BUFFER, texCoordVBO);
	glDeleteBuffers(1, &texCoordVBO);
	glBindBuffer(GL_ARRAY_BUFFER, texCoord3DVBO);
	glDeleteBuffers(1, &texCoord3DVBO);
	glBindBuffer(GL_ARRAY_BUFFER, vertexVBO);
	glDeleteBuffers(1, &vertexVBO);
#ifdef LINE
	glBindBuffer(GL_ARRAY_BUFFER, baryVBO);
	glDeleteBuffers(1, &baryVBO);
	glBindBuffer(GL_ARRAY_BUFFER, oriBaryVBO);
	glDeleteBuffers(1, &oriBaryVBO);
#endif

//#ifdef TRUTH
	glBindBuffer(GL_ARRAY_BUFFER, normalVBO_truth);
	glDeleteBuffers(1, &normalVBO_truth);
	glBindBuffer(GL_ARRAY_BUFFER, vertexVBO_truth);
	glDeleteBuffers(1, &vertexVBO_truth);
//#endif

	for (vector<GLuint>::size_type i = 0; i < globalIdxVBOList.size(); ++i)
	{
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, globalIdxVBOList[i]);
		glDeleteBuffers(1, &globalIdxVBOList[i]);
	}
	globalIdxVBOList.clear();

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, teapotIdxVBO);
	glDeleteBuffers(1, &teapotIdxVBO);

	viewMemD = 0;
}

void View::setBufferData()
{
	int tessPointCount = m_pCommonData->tessPointCount();
	glBindBuffer(GL_ARRAY_BUFFER, normalVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * tessPointCount * 3, 0, GL_DYNAMIC_DRAW);
	//cout << "@法向VBO" << endl << "\t";
	viewMemD += sizeof(float) * tessPointCount * 3;
	printMemD(__FILE__, __FUNCTION__, __LINE__, sizeof(float) * tessPointCount * 3, "@法向VBO");

	glBindBuffer(GL_ARRAY_BUFFER, texCoordVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * tessPointCount * 2, 0, GL_DYNAMIC_DRAW);
	//cout << "@二维纹理坐标VBO" << endl << "\t";
	viewMemD += sizeof(float) * tessPointCount * 2;
	printMemD(__FILE__, __FUNCTION__, __LINE__, sizeof(float) * tessPointCount * 2, "@二维纹理坐标VBO");

	glBindBuffer(GL_ARRAY_BUFFER, texCoord3DVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * tessPointCount * 2, 0, GL_DYNAMIC_DRAW);
	//cout << "@三维纹理坐标VBO" << endl << "\t";
	viewMemD += sizeof(float) * tessPointCount * 2;
	printMemD(__FILE__, __FUNCTION__, __LINE__, sizeof(float) * tessPointCount * 2, "@三维纹理坐标VBO");

	glBindBuffer(GL_ARRAY_BUFFER, vertexVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * tessPointCount * 3, 0, GL_DYNAMIC_DRAW);
	//cout << "@顶点坐标VBO" << endl << "\t";
	viewMemD += sizeof(float) * tessPointCount * 3;
	printMemD(__FILE__, __FUNCTION__, __LINE__, sizeof(float) * tessPointCount * 3, "@顶点坐标VBO");

	//cout << "nVBO = " << normalVBO << ", texVBO = " << texCoordVBO
		 //<< ", vVBO = " << vertexVBO << endl;
#ifdef LINE
	glBindBuffer(GL_ARRAY_BUFFER, baryVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * tessPointCount * 3, 0, GL_DYNAMIC_DRAW);
	//cout << "@重心坐标VBO" << endl << "\t";
	viewMemD += sizeof(float) * tessPointCount * 3;
	printMemD(__FILE__, __FUNCTION__, __LINE__, sizeof(float) * tessPointCount * 3, "@重心坐标VBO");

	glBindBuffer(GL_ARRAY_BUFFER, oriBaryVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * tessPointCount * 3, 0, GL_DYNAMIC_DRAW);
	//cout << "@重心坐标VBO" << endl << "\t";
	viewMemD += sizeof(float) * tessPointCount * 3;
	printMemD(__FILE__, __FUNCTION__, __LINE__, sizeof(float) * tessPointCount * 3, "@原始重心坐标VBO");
#endif


//#ifdef TRUTH
	int vertexCount_teapot = m_pCommonData->vertexCount_teapot();
	glBindBuffer(GL_ARRAY_BUFFER, normalVBO_truth);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * vertexCount_teapot * 3, 0, GL_DYNAMIC_DRAW);
	//cout << "@法向VBO_truth" << endl << "\t";
	viewMemD += sizeof(float) * vertexCount_teapot * 3;
	printMemD(__FILE__, __FUNCTION__, __LINE__, sizeof(float) * vertexCount_teapot * 3, "@法向VBO_truth");

	glBindBuffer(GL_ARRAY_BUFFER, vertexVBO_truth);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * vertexCount_teapot * 3, 0, GL_DYNAMIC_DRAW);
	//cout << "@顶点坐标VBO_truth" << endl << "\t";
	viewMemD += sizeof(float) * vertexCount_teapot * 3;
	printMemD(__FILE__, __FUNCTION__, __LINE__, sizeof(float) * vertexCount_teapot * 3, "@顶点坐标VBO_truth");
//#endif

	for (int i = 0; i < m_pCommonData->faceMtlListSize(); ++i)
	{
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, globalIdxVBOList[i]);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(int) * m_pCommonData->faceMtlCount(i) * 3,
					 m_pCommonData->faceMtl(i), GL_STATIC_DRAW);
		//cout << "@索引VBO_" << i << endl << "\t";
		viewMemD += sizeof(float) * m_pCommonData->faceMtlCount(i) * 3;
		printMemD(__FILE__, __FUNCTION__, __LINE__, sizeof(int) * m_pCommonData->faceMtlCount(i) * 3, "@索引VBO_");
	}

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, teapotIdxVBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(int) * m_pCommonData->teapotFaceList.size(),
					 &m_pCommonData->teapotFaceList[0], GL_STATIC_DRAW);
}

/* 载入 OBJ 模型*/
void View::loadObj(const QString &texName)
{
	//cout << "load length = " << m_pCommonData->length() << endl;
	for (int i = 0; i < 3; ++i)
	{
		lastPos[i] = 0.0;
		curPos[i] = 0.0;
		//不重置rotateAxis[i] = 1.0;
		lastPosCubeMap[i] = 0.0;
		curPosCubeMap[i] = 0.0;
		//不重置rotateAxisCubeMap[i] = 1.0;
	}
	m_fTheta = 0.0;
	m_fThetaCubeMap = 0.0;
	//m_fScale = 1.0;
	//m_fDeltaX = m_fDeltaY = 0.0;
	//不重置lastMatrix = Matrix4x4(1.0f);
	//不重置lastMatrixCubeMap = Matrix4x4(1.0f);

	m_bLoaded = true;
	setProjection();

	if (texName.toUtf8().constData() == string("自带2D纹理"))
	{
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		load2DTexture();
	}

	//float min[3], delta_inverse[3];
	//min[0] = m_pCommonData->minX();
	//min[1] = m_pCommonData->minY();
	//min[2] = m_pCommonData->minZ();
	//delta_inverse[0] = 1.0f / (m_pCommonData->maxX() - min[0]);
	//delta_inverse[1] = 1.0f / (m_pCommonData->maxY() - min[1]);
	//delta_inverse[2] = 1.0f / (m_pCommonData->maxZ() - min[2]);
	//glUseProgram(prog);
	//glUniform3fv(min_vertex_id_, 1, min);
	//glUniform3fv(delta_vertex_inverse_id_, 1, delta_inverse);

	//glUseProgram(0);
}

void View::load2DTexture()
{
	for (int i = 0; i < m_pCommonData->objMtlTexListSize(); ++i)
	{
		MtlTex mtlTex = m_pCommonData->getMtl(i);
		if (mtlTex.m_sTexFileName != "")
		{
			//glActiveTexture(GL_TEXTURE0 + mtlTex.m_nTexBindingIdx);
//#define QT_BIND
#ifdef QT_BIND
			int aa = bindTexture(QPixmap(QString(tr(mtlTex.m_sTexFileName.c_str()))));
			cout << "QT_BIND  = " << aa << endl;
#else
			glBindTexture(GL_TEXTURE_2D, tex2DNameList[i]);
			//cout << "size2 = " << m_pCommonData->objMtlTexListSize() << endl;
			//cout << "i = " << i << ", bind = " << tex2DNameList[i] << endl;
			//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
			//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
			//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
			QImage image(QString(tr(mtlTex.m_sTexFileName.c_str())));
			cout << mtlTex.m_sTexFileName.c_str() << endl;
			//unsigned char *a = image.bits();
			/*
			 * QRgb类型是一个typedef，等于一个unsigned int
			 * 如果使用pixel(int, int)函数取出QRgb值，则它的值是0xAARRGGBB
			 * 但是如果按照unsigned char一个一个取值，则上述一个QRgb对象对应
			 * 的四个值依次是：0xBB, 0xGG, 0xRR, 0xAA。
			 * glTexImage2D最后一个数组要求按照0xRR, 0xGG, 0xBB, 0xAA存放，
			 * 所以需要转换，0xAA和0xGG的位置不变，0xRR和0xBB的位置对调
			 */
			//image = image.convertToFormat(QImage::Format_ARGB32);
			/*
			 * QImage 类存储图像数据是按照图像在看图软件中显示的那样，
			 * 从上到下、从左到右一行一行存储的。
			 * 而OpenGL要求 glTexImage2D 最后一个参数指向的数据是从下到上、
			 * 从左到右存储，所以需要使用
			 * QImage QImage::mirrored(bool horizontal = false, bool vertical = true) const
			 * 函数
			 */
			image = image.mirrored();
			//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, image.width(), image.height(),
						 //0, GL_BGRA, GL_UNSIGNED_BYTE, image.bits());
			gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, image.width(), image.height(),
							  GL_BGRA, GL_UNSIGNED_BYTE, image.bits());
			//glGenerateMipmap(GL_TEXTURE_2D);
#endif
		}
	}
}

//#ifdef TRUTH
void View::loadColormap()
{
	QString colormap_path("/home/cym/program/OBJ/colormap.png");
	glBindTexture(GL_TEXTURE_2D, colormap_name);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	QImage image(colormap_path);
	//image.convertToFormat(QImage::Format_ARGB32);
	image = image.mirrored();
	gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, image.width(), image.height(),
					  GL_BGRA, GL_UNSIGNED_BYTE, image.bits());
}
//#endif

struct VolumeHeader
{
	char magic[4];
	int  version;
	char texName[256];
	bool wrap;
	int  volSize;
	int  numChannels;
	int  bytesPerChannel;
};

void View::load3DTexture(const string &textureName)
{
	string texFileName = string("/home/cym/program/textures/3D_textures/") + textureName + ".vol";
	std::ifstream fin(texFileName.c_str());
	char buf[4096];
	fin.read(buf, 4096);

	VolumeHeader * header = (VolumeHeader *)buf;
	if ((header->magic[0] != 'V') || (header->magic[1] != 'O') ||
		(header->magic[2] != 'L') || (header->magic[3] != 'U'))
		std::cerr << "bad header: invalid magic" << endl;
	if (header->version != 4)
		std::cerr << "bad header: version != 4" << endl;
	if (header->bytesPerChannel != 1)
		std::cerr << "bad header: only byte textures supported" << endl;

	int volBytes = header->volSize * header->volSize * header->volSize * header->numChannels;
	//cout << "volSize = " << header->volSize << endl;
	//cout << "numChannels = " << header->numChannels << endl;
	//cout << "volBytes = " << volBytes << endl;
	unsigned char *data = new unsigned char[volBytes];
	fin.read((char *)data, volBytes);

	glDeleteTextures(1, &tex3DName);
	glGenTextures(1, &tex3DName);
	cout << "tex3DName = " << tex3DName << endl;
	glBindTexture(GL_TEXTURE_3D, tex3DName);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexImage3D(GL_TEXTURE_3D, 0, GL_RGB, header->volSize, header->volSize,
				 header->volSize, 0, GL_RGB, GL_UNSIGNED_BYTE, data);

	fin.close();
	fin.clear();

	delete []data;
}

void View::loadCubeMap(const QString &textureName)
{
	glDeleteTextures(1, &texCubeName);
	glGenTextures(1, &texCubeName);
	cout << "texCubeName = " << texCubeName << endl;
	glBindTexture(GL_TEXTURE_CUBE_MAP, texCubeName);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	QString path = tr("/home/cym/program/textures/CubeMaps/") + textureName + "/";
	QString cubeMapFilename[] = {path + "posx.jpg", path + "negx.jpg",
								 path + "posy.jpg", path + "negy.jpg",
								 path + "posz.jpg", path + "negz.jpg"};
	for (int i = 0; i < 6; ++i)
	{
		QImage image(cubeMapFilename[i]);
		//image.convertToFormat(QImage::Format_ARGB32);
		glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_RGBA, image.width(), image.height(),
					 0, GL_BGRA, GL_UNSIGNED_BYTE, image.bits());
	}
}

/* 进行投影变换 */
void View::setProjection()
{
	double r_origin = m_pCommonData->length();
	double r_scale = r_origin / m_fScale;
	if (m_bLoaded)				// 模型已载入
	{
		m_projectionMatrixStack.loadIdentity();
		/* 保证整个模型都能显示在窗口内 */
		if (width() >= height())
		{
			GLfloat x = GLfloat(width()) / height();
			if (m_bIsOrtho)
			{
				m_projectionMatrixStack.ortho(-x * r_scale, x * r_scale,
											  -r_scale, r_scale,
											  r_origin * (EYE_Z_FACTOR - 2), r_origin * (EYE_Z_FACTOR + 1));
				//m_projectionMatrixStack.ortho(-x * r_scale, x * r_scale,
											  //-r_scale, r_scale,
											  //r_origin * (EYE_Z_FACTOR - 1), r_origin * (EYE_Z_FACTOR + 1));
				m_fRight = x * r_scale;
				m_fBottom = -r_scale;
			}
			else
			{
				// 下式括号中的内容通过(length / m_fScale) / (length * EYE_Z_FACTOR)简化而来
				// length是物体的半径，对一个物体来说是常量
				// 分子表示物体中心和视锥上边缘的距离；分母表示视点和物体的距离
				double fovy = 2 * atan(1.0 / (m_fScale * EYE_Z_FACTOR)) * 180 / PI;
				cout << "fovy = " << fovy << endl;
				m_projectionMatrixStack.perspective(fovy, x, r_origin, 8);
				m_fRight = x * r_scale;
				m_fBottom = -r_scale;
			}
		}
		else
		{
			GLfloat y = GLfloat(height()) / width();
			if (m_bIsOrtho)
			{
				m_projectionMatrixStack.ortho(-r_scale, r_scale,
											  -y * r_scale, y * r_scale,
											  r_origin * (EYE_Z_FACTOR - 1), r_origin * (EYE_Z_FACTOR + 1));
				//m_projectionMatrixStack.ortho(-r_scale, r_scale,
											  //-y * r_scale, y * r_scale,
											  //r_origin * (EYE_Z_FACTOR - 1), r_origin * (EYE_Z_FACTOR + 1));
				m_fRight = r_scale;
				m_fBottom = -y * r_scale;
			}
			else
			{
				// 下式括号中的内容通过(length * y / m_fScale) / (length * EYE_Z_FACTOR)简化而来
				// length是物体的半径，对一个物体来说是常量
				// 分子表示物体中心和视锥上边缘的距离；分母表示视点和物体的距离
				double fovy = 2 * atan(y / (m_fScale * EYE_Z_FACTOR)) * 180 / PI;
				m_projectionMatrixStack.perspective(fovy, 1.0 / y, r_origin, 8);
				m_fRight = r_scale;
				m_fBottom = -y * r_scale;
			}
		}
		glMatrixMode(GL_PROJECTION);						// deprecated
		glLoadMatrixf(m_projectionMatrixStack.top());		// deprecated
	}
}

/* 完成OpenGL的初始化工作 */
void View::initializeGL()
{
	setGLDevice();

	//glGenBuffers(1, &normalVBO);
	//glGenBuffers(1, &texCoordVBO);
	//glGenBuffers(1, &vertexVBO);

	GLfloat light_ambient[] = { 0.5, 0.5, 0.5, 1.0 };		// 光源的环境光分量
	GLfloat light_diffuse[] = { 0.8, 0.8, 0.8, 1.0 };		// 光源的漫射分量
	GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };		// 光源的镜面分量

	/* 设置光源0的环境分量、漫射分量、镜面分量和位置 */
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	//GLfloat lmodel_ambient[] = {1.0, 1.0, 1.0, 1.0};
	//glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);

	//glClearColor(0.796875, 0.90625, 0.80859375, 0.0);		// 清除屏幕的颜色（淡绿色）
	//glClearColor(0.0, 0.0, 0.0, 0.0);						// 清除屏幕的颜色（黑色）
	//glClearColor(1.0, 0.0, 0.0, 0.0);						// 清除屏幕的颜色（红色）
	glClearColor(1.0, 1.0, 1.0, 0.0);						// 清除屏幕的颜色（白色）

	//glClearColor(0.91015625, 0.90234375, 0.8984375, 0.0);		// 最常用，清除屏幕的颜色（灰色，kde）

	//glClearColor(0.9490196, 0.945098, 0.941176, 0.0);			// 清除屏幕的颜色（灰色，gnome）
	glEnable(GL_CULL_FACE);
	//glDisable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);								// 打开深度测试

	glEnable(GL_LIGHTING);									// 打开光照
	glEnable(GL_LIGHT0);									// 打开光源0
	//glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

	string vertexShader, fragmentShader;
//#ifdef LINE
	//readShaderSource("affd_line", vertexShader, fragmentShader);
//#else
//#ifdef TRUTH
	//readShaderSource("affd_truth", vertexShader, fragmentShader);
//#else
	readShaderSource("affd", vertexShader, fragmentShader);
//#endif
//#endif

	bool success = installShaders(vertexShader.c_str(), fragmentShader.c_str());
	if (!success)
	{
		cout << "installShader未成功！" << endl;
		return;
	}

	glGenQueries(1, &queryID[queryBackBuffer]);
	glGenQueries(1, &queryID[queryFrontBuffer]);
	// dummy query to prevent OpenGL errors from popping out
	glQueryCounter(queryID[queryFrontBuffer], GL_TIMESTAMP);
}

void View::swapQueryBuffers()
{
	if (queryBackBuffer)
	{
		queryBackBuffer = 0;
		queryFrontBuffer = 1;
	}
	else
	{
		queryBackBuffer = 1;
		queryFrontBuffer = 0;
	}
}

/* 当窗口大小发生变化时需要调用的函数 */
void View::resizeGL(int width, int height)
{
	glViewport(0, 0, (GLsizei)width, (GLsizei)height);		// 将视口设置为和窗口同样大小
	setProjection();							// 进行投影变换
}

/* 完成绘制工作的函数 */
void View::paintGL()										// 完成opengl绘制的函数
{
	if (m_bIsFlat)											// 绘制模式（gl_flat 和 gl_smooth）
		glShadeModel(GL_FLAT);
	else
		glShadeModel(GL_SMOOTH);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);		// 清空颜色缓存和深度缓存

	glBeginQuery(GL_TIME_ELAPSED, queryID[queryBackBuffer]);

	if (m_bLoaded)											// 如果obj文件已经载入
	{
		//m_calcTime.start();

		/* 将负责旋转的矩阵保存，不影响模型显示效果 */
		m_modelViewMatrixStack.loadIdentity();

		m_modelViewMatrixStack.rotate(m_fTheta, rotateAxis[0], rotateAxis[1], rotateAxis[2]);	// 计算新的旋转矩阵，即：m = e·r = r
		m_modelViewMatrixStack.multMatrix(lastMatrix);		// 左乘前一次的矩阵，即：m = r·l
		m_modelViewMatrixStack.top(lastMatrix);				// 此次处理结果保存在lastmatrix中，即：l = m

		//cout << m_modelViewMatrixStack << endl;

		/* 进行实际的旋转、平移和缩放 */
		m_modelViewMatrixStack.loadIdentity();
		double length = m_pCommonData->length();
		// 视点设置在Z轴正半轴
		m_modelViewMatrixStack.lookAt(0.0, 0.0, length * EYE_Z_FACTOR, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);	

		double r = 0.0;
		if (width() > height())
			r = 2.0 * length / (m_fScale * height());
		else
			r = 2.0 * length / (m_fScale * width());

		//if (m_eMouseButton == LEFT_BUTTON)								// 鼠标左键被按下
		//{
		//	if (m_nEditDirection == NONE_AXIS)							// 不对控制定点进行编辑
		//	{
				if (m_bShiftPressed)								// Shift键被按下，此时要旋转模型和cubemap
				{
					m_modelViewMatrixStack.multMatrix(lastMatrix);			// 旋转
					m_modelViewMatrixStack.translate(m_fDeltaX * r, m_fDeltaY * r, 0.0);
					light_position[0] = m_fDeltaX * r;
					light_position[1] = m_fDeltaY * r;
				}
				else
				{
					//m_modelViewMatrixStack.translate(m_fDeltaX * r, m_fDeltaY * r, 0.0);
					//m_modelViewMatrixStack.multMatrix(lastMatrix);			// 旋转
					m_modelViewMatrixStack.translate(m_fDeltaX * r, m_fDeltaY * r, 0.0);
					m_modelViewMatrixStack.multMatrix(lastMatrix);			// 旋转
					light_position[0] = m_fDeltaX * r;
					light_position[1] = m_fDeltaY * r;
				}
		//	}
		//}

		glMatrixMode(GL_MODELVIEW);								// deprecated
		glLoadMatrixf(m_modelViewMatrixStack.top());			// deprecated

		m_fTheta = 0.0;

		/* 算出垂直于屏幕指向外的向量 */
		m_fVector001[0] = lastMatrix[1][0] * lastMatrix[2][1] - lastMatrix[1][1] * lastMatrix[2][0];
		m_fVector001[1] = lastMatrix[0][1] * lastMatrix[2][0] - lastMatrix[0][0] * lastMatrix[2][1];
		m_fVector001[2] = lastMatrix[0][0] * lastMatrix[1][1] - lastMatrix[0][1] * lastMatrix[1][0];

		/* 开始模型的绘制 */
		if ((m_nDisplayMode & DISPLAY_MODEL) != 0)			// 需要绘制模型
		{
			if (m_pCommonData->algorithm() == FFD)			// 基本FFD算法
			{
				glUseProgram(0);
				//GLfloat mat_diffuse[] = {0.745, 0.478, 0.196, 1.0};

				float amb = 0.1;
				GLfloat mat_ambient[] = {amb, amb, amb, 1.0};

				float dif = 1.0;
				GLfloat mat_diffuse[] = {dif, dif, dif, 1.0};

				float spe = 0.9;
				GLfloat mat_specular[] = {spe, spe, spe, 1.0};

				GLfloat mat_shiness[] = { 40 };

				glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
				glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
				glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
				glMaterialfv(GL_FRONT, GL_SHININESS, mat_shiness);
				//glEnable(GL_TEXTURE_2D);
				unsigned int faceCount = m_pCommonData->faceCount();
				for (unsigned int i = 0; i < faceCount; ++i)
				{
					VertexCoord normal;
					unsigned int totalVertex = m_pCommonData->vertexCoordIdxCount(i);
					glBegin(GL_POLYGON);
						for (unsigned int j = 0; j < totalVertex; ++j)
						{
							int index = m_pCommonData->vertexCoordIdx(i, j);
							VertexCoord v0 = m_pCommonData->vertexCoord(index);
							//int texIdx = m_pCommonData->textureVertexIdx(i, j);
							//TextureCoord vt0 = m_pCommonData->textureVertex(texIdx);

							if (j == 0)						// 每个面片的法向量都在访问第一个顶点时计算
							{
								int index1 = m_pCommonData->vertexCoordIdx(i, 1);
								int index2 = m_pCommonData->vertexCoordIdx(i, 2);

								VertexCoord v1 = m_pCommonData->vertexCoord(index1);
								VertexCoord v2 = m_pCommonData->vertexCoord(index2);
								VertexCoord v10 = v0 - v1;
								VertexCoord v12 = v2 - v1;

								//normal = v12 * v10;
								normal = cross(v12, v10);
								normal.normalize();
							}
							glNormal3d((GLdouble)normal.x(), (GLdouble)normal.y(), (GLdouble)normal.z());
							//glTexCoord2f((GLdouble)vt0.u(), (GLdouble)vt0.v());
							glVertex3d((GLdouble)v0.x(), (GLdouble)v0.y(), (GLdouble)v0.z());
						}
					glEnd();
				}
				glDisable(GL_TEXTURE_2D);
			}
			else											// AFFD算法
			{
				if (m_pCommonData->isGPU())					// 使用GPU计算
				{
					glUseProgram(prog);
					glUniformMatrix4fv(MVMatrix_id_, 1, false, m_modelViewMatrixStack.top());
					glUniformMatrix4fv(PMatrix_id_, 1, false, m_projectionMatrixStack.top());
					glUniformMatrix3fv(NormalMatrix_id_, 1, false, m_modelViewMatrixStack.normalMatrix());
					if (m_bDisplayTruth)
					{
						glBindBuffer(GL_ARRAY_BUFFER, vertexVBO_truth);
						const GLuint index_mPosition = glGetAttribLocation(prog, "MCvertex");
						glVertexAttribPointer(index_mPosition, 3, GL_FLOAT, GL_FALSE, 0, 0);
						glEnableVertexAttribArray(index_mPosition);

						glBindBuffer(GL_ARRAY_BUFFER, normalVBO_truth);
						const GLuint index_mNormal = glGetAttribLocation(prog, "MCnormal");
						glVertexAttribPointer(index_mNormal, 3, GL_FLOAT, GL_FALSE, 0, 0);
						glEnableVertexAttribArray(index_mNormal);

						glUniform1i(local_viewer_id_, m_bLocalViewer);
						//glUniform4f(light_position_id_, 0, 0, m_fDist, 1.0);
						glUniform4f(light_position_id_, light_position[0], light_position[1], light_position[2], 1.0);
						glUniform1i(use_env_map_id_, m_bUseEnvMap);

						float cubeMapRotInvMat[9];
						lastMatrixCubeMap.inverseRotateMatrix(cubeMapRotInvMat);
						glUniformMatrix3fv(cubemap_rot_inv_mat_id_, 1, false, cubeMapRotInvMat);

						glUniform1i(tex_case_, 0);
						MtlTex mtlTex = m_pCommonData->getMtl(0);
						glUniform4f(ka_id_, mtlTex.m_fKa[0], mtlTex.m_fKa[1], mtlTex.m_fKa[2], 1.0);
						glUniform4f(kd_id_, mtlTex.m_fKd[0], mtlTex.m_fKd[1], mtlTex.m_fKd[2], 1.0);
						glUniform4f(ks_id_, mtlTex.m_fKs[0], mtlTex.m_fKs[1], mtlTex.m_fKs[2], 1.0);

						glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, teapotIdxVBO);
						int triangleCount = m_pCommonData->teapotFaceList.size();
						glDrawElements(GL_TRIANGLES, triangleCount * 3, GL_UNSIGNED_INT, 0);
					}
					else
					{
						glBindBuffer(GL_ARRAY_BUFFER, vertexVBO);
						const GLuint index_mPosition = glGetAttribLocation(prog, "MCvertex");
						glVertexAttribPointer(index_mPosition, 3, GL_FLOAT, GL_FALSE, 0, 0);
						glEnableVertexAttribArray(index_mPosition);

						glBindBuffer(GL_ARRAY_BUFFER, texCoordVBO);
						const GLuint index_mTexCoord = glGetAttribLocation(prog, "TexCoord2D0");
						glVertexAttribPointer(index_mTexCoord, 2, GL_FLOAT, GL_FALSE, 0, 0);
						glEnableVertexAttribArray(index_mTexCoord);

						glBindBuffer(GL_ARRAY_BUFFER, texCoord3DVBO);
						const GLuint index_mTexCoord3D = glGetAttribLocation(prog, "TexCoord3D0");
						glVertexAttribPointer(index_mTexCoord3D, 2, GL_FLOAT, GL_FALSE, 0, 0);
						glEnableVertexAttribArray(index_mTexCoord3D);

						glBindBuffer(GL_ARRAY_BUFFER, normalVBO);
						const GLuint index_mNormal = glGetAttribLocation(prog, "MCnormal");
						glVertexAttribPointer(index_mNormal, 3, GL_FLOAT, GL_FALSE, 0, 0);
						glEnableVertexAttribArray(index_mNormal);

#ifdef LINE
						glBindBuffer(GL_ARRAY_BUFFER, baryVBO);
						const GLuint index_mBary = glGetAttribLocation(prog, "Bary");
						glVertexAttribPointer(index_mBary, 3, GL_FLOAT, GL_FALSE, 0, 0);
						glEnableVertexAttribArray(index_mBary);

						glBindBuffer(GL_ARRAY_BUFFER, oriBaryVBO);
						const GLuint index_mOriBary = glGetAttribLocation(prog, "OriBary");
						glVertexAttribPointer(index_mOriBary, 3, GL_FLOAT, GL_FALSE, 0, 0);
						glEnableVertexAttribArray(index_mOriBary);

						glUniform1i(divide_id_, m_pCommonData->samplePointCount() - 1);
#endif

	//#ifdef TRUTH
						glBindBuffer(GL_ARRAY_BUFFER, vertexVBO_truth);
						const GLuint index_mPosition_truth = glGetAttribLocation(prog, "MCvertex_truth");
						glVertexAttribPointer(index_mPosition_truth, 3, GL_FLOAT, GL_FALSE, 0, 0);
						glEnableVertexAttribArray(index_mPosition_truth);

						glBindBuffer(GL_ARRAY_BUFFER, normalVBO_truth);
						const GLuint index_mNormal_truth = glGetAttribLocation(prog, "MCnormal_truth");
						glVertexAttribPointer(index_mNormal_truth, 3, GL_FLOAT, GL_FALSE, 0, 0);
						glEnableVertexAttribArray(index_mNormal_truth);
	//#endif

						glUniform1i(local_viewer_id_, m_bLocalViewer);
						//glUniform4f(light_position_id_, 0, 0, m_fDist, 1.0);
						glUniform4f(light_position_id_, light_position[0], light_position[1], light_position[2], 1.0);

						glActiveTexture(GL_TEXTURE0);
						glEnable(GL_TEXTURE_CUBE_MAP);
						glBindTexture(GL_TEXTURE_CUBE_MAP, texCubeName);

	//#ifdef TRUTH
						glActiveTexture(GL_TEXTURE1);
						glEnable(GL_TEXTURE_2D);
						glBindTexture(GL_TEXTURE_2D, colormap_name);
	//#endif

						glActiveTexture(GL_TEXTURE3);
						glEnable(GL_TEXTURE_3D);
						glBindTexture(GL_TEXTURE_3D, tex3DName);

						glUniform1i(use_env_map_id_, m_bUseEnvMap);
	//#ifdef TRUTH
						glUniform1i(display_truth_id_, m_bDisplayTruth);
						glUniform1i(error_texture_id_, error_texture_);
						glUniform1i(error_display_mode_id_, (int)m_eErrorDisplayMode);
	//#endif

#ifdef LINE
						glUniform1i(use_line_id_, m_bLine);
#endif

						float cubeMapRotInvMat[9];
						lastMatrixCubeMap.inverseRotateMatrix(cubeMapRotInvMat);
						glUniformMatrix3fv(cubemap_rot_inv_mat_id_, 1, false, cubeMapRotInvMat);

#ifndef MORPH
						int ttt = 0;
#endif
						for (int i = 0; i < m_pCommonData->objMtlTexListSize() + 1; ++i)
						{
							//cout << "绘制：i = " << i << endl;
							glActiveTexture(GL_TEXTURE2);
							glEnable(GL_TEXTURE_2D);
#ifdef QT_BIND
							glBindTexture(GL_TEXTURE_2D, 4);
#else
							glBindTexture(GL_TEXTURE_2D, tex2DNameList[i]);
#endif

							MtlTex mtlTex = m_pCommonData->getMtl(i);
							glUniform4f(ka_id_, mtlTex.m_fKa[0], mtlTex.m_fKa[1], mtlTex.m_fKa[2], 1.0);
							glUniform4f(kd_id_, mtlTex.m_fKd[0], mtlTex.m_fKd[1], mtlTex.m_fKd[2], 1.0);
							glUniform4f(ks_id_, mtlTex.m_fKs[0], mtlTex.m_fKs[1], mtlTex.m_fKs[2], 1.0);

							int texCase = m_nTexCase;
							//if (mtlTex.m_sTexFileName == "" && m_nTexCase == 2)
								//texCase = 0;
							glUniform1i(tex_case_, texCase);

							glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, globalIdxVBOList[i]);
							int tessTriangleCount = m_pCommonData->faceMtlCount(i);


							//cout << "tessTriangleCount = " << tessTriangleCount << endl;
							glDrawElements(GL_TRIANGLES, tessTriangleCount * 3, GL_UNSIGNED_INT, 0);


							/*glUniform4f(ka_id_, mtlTex.m_fKa[0], mtlTex.m_fKa[1], mtlTex.m_fKa[2], 1.0);
							glUniform4f(kd_id_, 1.0, 0.0, 0.0, 1.0);
							glUniform4f(ks_id_, mtlTex.m_fKs[0], mtlTex.m_fKs[1], mtlTex.m_fKs[2], 1.0);
							glDrawElements(GL_POINTS, tessTriangleCount * 3, GL_UNSIGNED_INT, 0);*/

#ifndef MORPH
							ttt += tessTriangleCount;
#endif
							//glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_INT, 0);
						}
#ifndef MORPH
						cout << "面片数 = " << ttt << endl;
						//glFinish();
						//glFlush();
						//int elapsedTime = m_calcTime.elapsed();
						//cout << "绘制\t" << elapsedTime << "\tGPU绘制" << endl;
#endif
					}
					glUseProgram(0);
				}
				else										// 使用 CPU 计算
				{
				}
			}
		}
		glActiveTexture(GL_TEXTURE0);
		glDisable(GL_TEXTURE_CUBE_MAP);

		glActiveTexture(GL_TEXTURE1);
		glDisable(GL_TEXTURE_2D);

		glActiveTexture(GL_TEXTURE2);
		glDisable(GL_TEXTURE_2D);

		glActiveTexture(GL_TEXTURE3);
		glDisable(GL_TEXTURE_3D);

		////////////////////////////////////////////////////////////////
		//GLfloat mat_diffuse_sphere[] = {0.0, 0.0, 1.0, 1.0};
		//glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse_sphere);
		//VertexCoord v0(0, 0.333, 0.999);
		//VertexCoord v1(0, -0.333, 0.999);
		//VertexCoord v01(v1 - v0);
		//VertexCoord vm(0.5 * (v0 + v1));
		//NormalCoord n0(0.301758, 0.90437, 0.301758);
		//float t0 = ((vm - v0) * v01) / (n0 * v01);
		//VertexCoord O0 = v0 + n0 * t0;
		//double rr =  (v0 - O0).norm();
		//drawSphere(rr, O0.x(), O0.y(), O0.z());

		//mat_diffuse_sphere[0] = 1.0;
		//mat_diffuse_sphere[2] = 0.0;
		//glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse_sphere);
		//NormalCoord n1(-0.301758, -0.90437, 0.301758);
		//float t1 = ((vm - v0) * v01) / (n1 * v01);
		//VertexCoord O1 = v0 + n1 * t1;
		//rr =  (v0 - O1).norm();
		//drawSphere(rr, O1.x(), O1.y(), O1.z());

		//VertexCoord center_delta = O0 - O1;
		//float t = (v0 - O0) * center_delta / (center_delta * center_delta);
		//cout << "my t = " << t << endl;

		//drawSphere(0.01, 0.666, 0.1111111, 0.999);
		////////////////////////////////////////////////////////////////

		if ((m_nDisplayMode & DISPLAY_CTRL_POINTS) != 0)		// 需要绘制控制网
		{
			glLineWidth(1.0f);
			GLfloat mat_diffuse1[] = {0.0, 0.0, 0.0, 1.0};
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse1);
			/* 控制顶点连线的绘制 */
			for (int k = 0; k < m_pCommonData->ctrlPointCount(W); ++k)
			{
				for (int j = 0; j < m_pCommonData->ctrlPointCount(V); ++j)
				{
					for (int i = 0; i < m_pCommonData->ctrlPointCount(U); ++i)
					{
						CtrlPoint ctrlPoint = m_pCommonData->getCtrlPoint(i, j, k);
						if (k < m_pCommonData->ctrlPointCount(W) - 1)
						{
							CtrlPoint ctrlPointW = m_pCommonData->getCtrlPoint(i, j, k + 1);
							glBegin(GL_LINES);
								glNormal3fv(m_fVector001);
								glVertex3f(ctrlPoint.x(), ctrlPoint.y(), ctrlPoint.z());
								glVertex3f(ctrlPointW.x(), ctrlPointW.y(), ctrlPointW.z());
							glEnd();
						}
						if (j < m_pCommonData->ctrlPointCount(V) - 1)
						{
							CtrlPoint ctrlPointV = m_pCommonData->getCtrlPoint(i, j + 1, k);
							glBegin(GL_LINES);
								glNormal3fv(m_fVector001);
								glVertex3f(ctrlPoint.x(), ctrlPoint.y(), ctrlPoint.z());
								glVertex3f(ctrlPointV.x(), ctrlPointV.y(), ctrlPointV.z());
							glEnd();
						}
						if (i < m_pCommonData->ctrlPointCount(U) - 1)
						{
							CtrlPoint ctrlPointU = m_pCommonData->getCtrlPoint(i + 1, j, k);
							glBegin(GL_LINES);
								glNormal3fv(m_fVector001);
								glVertex3f(ctrlPoint.x(), ctrlPoint.y(), ctrlPoint.z());
								glVertex3f(ctrlPointU.x(), ctrlPointU.y(), ctrlPointU.z());
							glEnd();
						}
					}
				}
			}
			/* 控制顶点的绘制 */
			drawCtrlPoint(GL_RENDER);
		}

#ifdef DRAW_TRIANGULAR_CTRL_POINTS
		if ((m_nDisplayMode & DISPLAY_TRIANGULAR_CTRL_POINTS) != 0)	// 需要绘制三角Bezier控制定点
			drawTriangularCtrlPoints();
#endif

		if ((m_nDisplayMode & DISPLAY_DIRECT_POINT) != 0)	// 需要绘制直接编辑点
		{
			drawDirectPoint(GL_RENDER);
		}
		if ((m_nDisplayMode & DISPLAY_KNOT) != 0)			// 需要绘制节点
		{
			glLineWidth(1.0f);
			GLfloat mat_diffuse1[] = {0.0, 0.0, 0.0, 1.0};
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse1);
			for (int k = 0; k < m_pCommonData->order(W) + m_pCommonData->ctrlPointCount(W) - 1; ++k)
			{
				for (int j = 0; j < m_pCommonData->order(V) + m_pCommonData->ctrlPointCount(V) - 1; ++j)
				{
					for (int i = 0; i < m_pCommonData->order(U) + m_pCommonData->ctrlPointCount(U) - 1; ++i)
					{
						double x = m_pCommonData->getKnot(U, i);
						double y = m_pCommonData->getKnot(V, j);
						double z = m_pCommonData->getKnot(W, k);
						double x1 = m_pCommonData->getKnot(U, i + 1);
						double y1 = m_pCommonData->getKnot(V, j + 1);
						double z1 = m_pCommonData->getKnot(W, k + 1);
						glBegin(GL_LINES);
							glNormal3fv(m_fVector001);
							glVertex3f(x, y, z);
							glVertex3f(x1, y, z);
							glVertex3f(x, y, z);
							glVertex3f(x, y1, z);
							glVertex3f(x, y, z);
							glVertex3f(x, y, z1);
						glEnd();
					}
				}
			}
		}

		if (m_eMouseButton == RIGHT_BUTTON)	// 需要绘制选取框
		{
			glLineWidth(1.0f);
			GLfloat mat_diffuse1[] = {1.0, 0.0, 1.0, 1.0};
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse1);

			/* 虚线 */
			glLineStipple(1, 0xFFF0);
			glEnable(GL_LINE_STIPPLE);

			double centerX, centerY, selectionWidth, selectionHeight;
			calcSelectionRange(centerX, centerY, selectionWidth, selectionHeight);
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse1);

			m_projectionMatrixStack.pushMatrix();
			m_projectionMatrixStack.loadIdentity();
			m_projectionMatrixStack.ortho2D(0.0f, width(), 0.0f, height());
			m_modelViewMatrixStack.pushMatrix();
			m_modelViewMatrixStack.loadIdentity();

			glMatrixMode(GL_PROJECTION);
			glPushMatrix();
			glLoadMatrixf(m_projectionMatrixStack.top());
			glMatrixMode(GL_MODELVIEW);
			glPushMatrix();
			glLoadMatrixf(m_modelViewMatrixStack.top());

			glBegin(GL_LINE_LOOP);
				glNormal3f(0.0, 0.0, 1.0);
				glVertex2f(centerX - selectionWidth / 2.0,
						  (double)height() - (centerY - selectionHeight / 2.0));
				glVertex2f(centerX - selectionWidth / 2.0,
						  (double)height() - (centerY + selectionHeight / 2.0));
				glVertex2f(centerX + selectionWidth / 2.0,
						  (double)height() - (centerY + selectionHeight / 2.0));
				glVertex2f(centerX + selectionWidth / 2.0,
						  (double)height() - (centerY - selectionHeight / 2.0));
			glEnd();

			m_modelViewMatrixStack.popMatrix();
			m_projectionMatrixStack.popMatrix();

			glMatrixMode(GL_MODELVIEW);
			glPopMatrix();
			glMatrixMode(GL_PROJECTION);
			glPopMatrix();

			glDisable(GL_LINE_STIPPLE);
		}
		if ((m_nDisplayMode & DISPLAY_EDIT_ICON) != 0)			// 需要绘制编辑图标
		{
			if (m_eEditMode == ROTATE)
			{
				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();
				drawCircle(RIGHTFAC, BOTTOMFAC, GL_RENDER);
				glMatrixMode(GL_PROJECTION);
				glPopMatrix();

				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();
				drawAxis(RIGHTFAC, BOTTOMFAC, GL_RENDER, false);
				glMatrixMode(GL_PROJECTION);
				glPopMatrix();
			}
			else
			{
				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();
				drawAxis(RIGHTFAC, BOTTOMFAC, GL_RENDER, true);
				glMatrixMode(GL_PROJECTION);
				glPopMatrix();
			}
		}
		if (m_bUseCubeMap)
		{
			glActiveTexture(GL_TEXTURE0);
			glEnable(GL_TEXTURE_CUBE_MAP);
			glBindTexture(GL_TEXTURE_CUBE_MAP, texCubeName);

			glActiveTexture(GL_TEXTURE2);
			glDisable(GL_TEXTURE_2D);

			glActiveTexture(GL_TEXTURE3);
			glDisable(GL_TEXTURE_3D);

			drawCubeMap();

			glActiveTexture(GL_TEXTURE0);
			glDisable(GL_TEXTURE_CUBE_MAP);
		}
	}
	/*double det = lastMatrix[0][0] * lastMatrix[1][1] * lastMatrix[2][2] +
				 lastMatrix[1][0] * lastMatrix[2][1] * lastMatrix[0][2] +
				 lastMatrix[2][0] * lastMatrix[0][1] * lastMatrix[1][2] -
				 lastMatrix[2][0] * lastMatrix[1][1] * lastMatrix[0][2] -
				 lastMatrix[1][0] * lastMatrix[0][1] * lastMatrix[2][2] -
				 lastMatrix[0][0] * lastMatrix[2][1] * lastMatrix[1][2];
	cout << "det = " << det << endl;*/
	/*double detA = lastMatrix[0] * lastMatrix[5] * lastMatrix[10] +
				  lastMatrix[4] * lastMatrix[9] * lastMatrix[2] +
				  lastMatrix[8] * lastMatrix[1] * lastMatrix[6] -
  				  lastMatrix[0] * lastMatrix[9] * lastMatrix[6] -
				  lastMatrix[4] * lastMatrix[1] * lastMatrix[10] -
				  lastMatrix[8] * lastMatrix[5] * lastMatrix[2];*/

	glEndQuery(GL_TIME_ELAPSED);
	GLuint64 elapsed_time;

	glGetQueryObjectui64v(queryID[queryFrontBuffer], GL_QUERY_RESULT, &elapsed_time);
#ifndef MORPH
	cout << "绘制\t" << elapsed_time / 1000000.0 << "\tGPU绘制" << endl;
#endif
	swapQueryBuffers();
}

void View::drawCubeMap()
{
	/* 将负责旋转cubemap的矩阵保存，不影响模型显示效果 */
	m_cubeMapMatrixStack.loadIdentity();

	m_cubeMapMatrixStack.rotate(m_fThetaCubeMap, rotateAxisCubeMap[0], rotateAxisCubeMap[1], rotateAxisCubeMap[2]);
	if (m_pCommonData->cubemap_matrix_list_.size() != 0)
	{
		//cout << "sorceIdx = " << m_pCommonData->sourceIdx << endl;
		m_cubeMapMatrixStack.multMatrix(m_pCommonData->cubemap_matrix_list_[m_pCommonData->sourceIdx]);
	}
	m_cubeMapMatrixStack.multMatrix(lastMatrixCubeMap);
	m_cubeMapMatrixStack.top(lastMatrixCubeMap);

	/* 进行实际的旋转 */
	m_cubeMapMatrixStack.pushMatrix();
	m_cubeMapMatrixStack.loadIdentity();
	m_cubeMapMatrixStack.lookAt(0.0, 0.0, 0.0,
								0.0, 0.0, -1.0,
								0.0, 1.0, 0.0);
	m_cubeMapMatrixStack.multMatrix(lastMatrixCubeMap);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadMatrixf(m_cubeMapMatrixStack.top());

	m_projectionMatrixStack.pushMatrix();
	m_projectionMatrixStack.loadIdentity();
	m_projectionMatrixStack.perspective(40.0f, (double)width() / height(), 1, 10000);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadMatrixf(m_projectionMatrixStack.top());

	m_fThetaCubeMap = 0.0;

	float unit = 5000.0f;
	glBegin(GL_QUADS);
		// posX 1265
		glTexCoord3f(1,1,-1);
		glVertex3f(unit,  unit, unit);
		glTexCoord3f(1,1,1);
		glVertex3f(unit,  unit, -unit);
		glTexCoord3f(1,-1,1);
		glVertex3f(unit, -unit, -unit);
		glTexCoord3f(1,-1,-1);
		glVertex3f(unit, -unit, unit);

		// negX 3047
		glTexCoord3f(-1,1,1);
		glVertex3f(-unit,  unit, -unit);
		glTexCoord3f(-1,1,-1);
		glVertex3f(-unit,  unit,  unit);
		glTexCoord3f(-1,-1,-1);
		glVertex3f(-unit, -unit,  unit);
		glTexCoord3f(-1,-1,1);
		glVertex3f(-unit, -unit, -unit);

		// posY 3210
		glTexCoord3f(-1,1,1);
		glVertex3f(-unit,  unit, -unit);
		glTexCoord3f(1,1,1);
		glVertex3f(unit,  unit, -unit);
		glTexCoord3f(1,1,-1);
		glVertex3f(unit, unit, unit);
		glTexCoord3f(-1,1,-1);
		glVertex3f(-unit, unit,  unit);

		// negY 7456
		glTexCoord3f(-1,-1,1);
		glVertex3f(-unit, -unit, -unit);
		glTexCoord3f(-1,-1,-1);
		glVertex3f(-unit, -unit,  unit);
		glTexCoord3f(1,-1,-1);
		glVertex3f(unit, -unit, unit);
		glTexCoord3f(1,-1,1);
		glVertex3f(unit, -unit, -unit);

		// posZ 3762
		glTexCoord3f(-1,1,1);
		glVertex3f(-unit,  unit, -unit);
		glTexCoord3f(-1,-1,1);
		glVertex3f(-unit, -unit, -unit);
		glTexCoord3f(1,-1,1);
		glVertex3f(unit, -unit, -unit);
		glTexCoord3f(1,1,1);
		glVertex3f(unit, unit, -unit);

		// negZ 0154
		glTexCoord3f(-1,1,-1);
		glVertex3f(-unit,  unit, unit);
		glTexCoord3f(1,1,-1);
		glVertex3f(unit,  unit,  unit);
		glTexCoord3f(1,-1,-1);
		glVertex3f(unit, -unit,  unit);
		glTexCoord3f(-1,-1,-1);
		glVertex3f(-unit, -unit, unit);
	glEnd();

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	//m_modelViewMatrixStack.popMatrix();
	m_cubeMapMatrixStack.popMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	m_projectionMatrixStack.popMatrix();
}

/* 绘制控制顶点 */
void View::drawCtrlPoint(GLenum mode)
{
	GLfloat mat_diffuseGreen[] = {0.0, 1.0, 0.0, 1.0};
	GLfloat mat_diffuseRed[] = {1.0, 0.0, 0.0, 1.0};
	glPointSize(7.0);
	for (int k = 0; k < m_pCommonData->ctrlPointCount(W); ++k)
	{
		if (mode == GL_SELECT)
			glLoadName(k);
		for (int j = 0; j < m_pCommonData->ctrlPointCount(V); ++j)
		{
			if (mode == GL_SELECT)
				glPushName(j);
			for (int i = 0; i < m_pCommonData->ctrlPointCount(U); ++i)
			{
				if (mode == GL_SELECT)
					glPushName(i);
				CtrlPoint ctrlPoint = m_pCommonData->getCtrlPoint(i, j, k);
				if (!ctrlPoint.selected())			// 未选中的控制顶点用绿色绘制
					glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuseGreen);
				else								// 选中的控制顶点用红色绘制
					glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuseRed);
				glBegin(GL_POINTS);
					glNormal3fv(m_fVector001);
					glVertex3f(ctrlPoint.x(), ctrlPoint.y(), ctrlPoint.z());
				glEnd();
				if (mode == GL_SELECT)
					glPopName();
			}
			if (mode == GL_SELECT)
				glPopName();
		}
	}
}

extern int index2c(int i, int j, int stride);

#ifdef DRAW_TRIANGULAR_CTRL_POINTS
void View::drawTriangularCtrlPoints()
{
	extern float *triangular_ctrl_points;
	extern int triangleNum, triangleCtrlPointNum_lower;
	extern int copy_begin_id, copy_end_id;

	float *c = triangular_ctrl_points;
	int f = triangleNum, m_ = triangleCtrlPointNum_lower;

	// 绘制三角Bezier曲面片控制顶点
	GLfloat mat_diffuse_cyan[] = {0.0, 1.0, 1.0, 1.0};
	GLfloat mat_diffuse_yellow[] = {1.0, 1.0, 0.0, 1.0};
	//glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse_cyan);
	//const int face_idx_lower = 44, face_idx_upper = 44;		// 圆柱体
	//const int face_idx_lower = 54, face_idx_upper = 55;		// 立方体
	//const int face_idx_lower = 54, face_idx_upper = 54;			// 立方体
	//const int face_idx_lower = 22, face_idx_upper = 22;			// 立方体
	//const int face_idx_lower = 22, face_idx_upper = 22;			// 立方体
	//const int face_idx_lower = 0, face_idx_upper = f - 1;		// 全部显示
	//const int face_idx_lower = copy_begin_id, face_idx_upper = copy_end_id;

	//if (draw_which_ == f)
		//draw_which_ = 0;

	if (draw_which_ % 2 == 0)		// 用于膨胀的立方体
		draw_which_ = 22;
	else
		draw_which_ = 25;

	//if (draw_which_ % 2 == 0)		// 用于立方体
		//draw_which_ = 24;
	//else
		//draw_which_ = 23;

	//if (draw_which_ > 65)
		//draw_which_ = 63;
	//else if (draw_which_ == 64)
		//draw_which_ = 65;

	//cout << "draw " << draw_which_ << ", f = " << f << endl;
	int face_idx_lower = draw_which_, face_idx_upper = draw_which_;
	//int face_idx_lower = 23, face_idx_upper = 26;
	//int face_idx_lower = 23, face_idx_upper = 23;
	for (int i = face_idx_lower; i <= face_idx_upper; ++i)
	{
		//if (i > 22 && i < 25) continue;
		//if (i == 64) continue;
		for (int j = 0; j < m_; ++j)
		{
			if (j == 4)
			{
				glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse_yellow);
				//glPointSize(18.0);
				glPointSize(5.0);
			}
			else if (j == 1)
			{
				GLfloat mat_diffuse_red[] = {1.0, 0.0, 0.0, 1.0};
				glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse_red);
				//glPointSize(18.0);
				glPointSize(5.0);
			}
			else
			{
				glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse_cyan);
				glPointSize(5.0);
			}
			glBegin(GL_POINTS);
				glNormal3fv(m_fVector001);
				glVertex3d(c[m_ * i + j], c[m_ * (i + f) + j], c[m_ * (i + f * 2) + j]);
				//cout << "cp = " << c[m_ * i + j] << ", " << c[m_ * (i + f) + j]
					 //<< ", " << c[m_ * (i + f * 2) + j] << endl;
			glEnd();
		}
	}

	// 绘制三角形左侧边
	extern int degree_lower;
	int n_ = degree_lower;
	GLfloat mat_diffuse_black[] = {0.0, 0.0, 0.0, 1.0};
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse_black);
	glLineWidth(1.0f);
	for (int i = face_idx_lower; i <= face_idx_upper; ++i)
	{
		//if (i > 22 && i < 25) continue;
		//if (i == 64) continue;
		glBegin(GL_LINE_STRIP);
			int idx = 0, gap = 0;
			for (int j = 0; j <= n_; ++j)
			{
				glNormal3fv(m_fVector001);
				glVertex3d(c[m_ * i + idx], c[m_ * (i + f) + idx], c[m_ * (i + f * 2) + idx]);
				++gap;
				idx += gap;
			}
		glEnd();
	}
	// 绘制三角形右侧边
	for (int i = face_idx_lower; i <= face_idx_upper; ++i)
	{
		//if (i > 22 && i < 25) continue;
		//if (i == 64) continue;
		glBegin(GL_LINE_STRIP);
			int idx = 0, gap = 2;
			for (int j = 0; j <= n_; ++j)
			{
				glNormal3fv(m_fVector001);
				glVertex3d(c[m_ * i + idx], c[m_ * (i + f) + idx], c[m_ * (i + f * 2) + idx]);
				idx += gap;
				++gap;
			}
		glEnd();
	}
	// 绘制三角形下侧边
	for (int i = face_idx_lower; i <= face_idx_upper; ++i)
	{
		//if (i > 22 && i < 25) continue;
		//if (i == 64) continue;
		glBegin(GL_LINE_STRIP);
			for (int j = m_ - 1; j >= m_ - degree_lower - 1; --j)
			{
				glNormal3fv(m_fVector001);
				glVertex3d(c[m_ * i + j], c[m_ * (i + f) + j], c[m_ * (i + f * 2) + j]);
			}
		glEnd();
	}
}
#endif

void View::drawDirectPoint(GLenum mode)
{
	GLfloat mat_diffuse0[] = {0.0, 0.0, 1.0, 1.0};
	GLfloat mat_diffuse1[] = {1.0, 0.0, 1.0, 1.0};
	glPointSize(10.0);
	for (int i = 0; i < m_pCommonData->directPointVectorSize(); ++i)
	{
		if (mode == GL_SELECT)
		{
			glLoadName(DIRECTPOINTSELECTBASE + i);
			//cout << "name = " << DIRECTPOINTSELECTBASE + i << endl;
		}
		DirectPoint temp = m_pCommonData->directPoint(i);
		if (temp.selected())
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse1);
		else
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse0);
		glBegin(GL_POINTS);
			glNormal3fv(m_fVector001);
			glVertex3d(temp.x(), temp.y(), temp.z());
		glEnd();
	}
}

/* 绘制用于平移和缩放控制顶点的坐标轴 */
void View::drawAxis(double rightFac, double bottomFac, GLenum mode, bool drawChip)
{
	double r_origin = m_pCommonData->length();
	double r_scale = r_origin / m_fScale;

	if (width() >= height())
	{
		GLfloat x = GLfloat(width()) / height();
		glOrtho(-x * r_scale, x * r_scale,
				-r_scale, r_scale,
				r_origin * (EYE_Z_FACTOR - 1), r_origin * (EYE_Z_FACTOR + 1));
	}
	else
	{
		GLfloat y = GLfloat(height()) / width();
		glOrtho(-r_scale, r_scale,
				-y * r_scale, y * r_scale,
				r_origin * (EYE_Z_FACTOR - 1), r_origin * (EYE_Z_FACTOR + 1));
	}
	m_modelViewMatrixStack.pushMatrix();
	m_modelViewMatrixStack.loadIdentity();
	m_modelViewMatrixStack.lookAt(0.0, 0.0, r_origin * EYE_Z_FACTOR,
								  0.0, 0.0, 0.0,
								  0.0, 1.0, 0.0);
	m_modelViewMatrixStack.translate(m_fRight * rightFac, m_fBottom * bottomFac, 0);
	m_modelViewMatrixStack.multMatrix(lastMatrix);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadMatrixf(m_modelViewMatrixStack.top());

	double shorter = m_fRight <= -m_fBottom ? m_fRight : -m_fBottom;
	int n = 20;
	double unit = shorter * 0.1;
	double rAxis = unit * 0.1;

	double R = rAxis * 2.0;
	double L = unit * 2;

	for (int axis = 3; axis > 0; --axis)
	{
		GLfloat mat_diffuse[] = {0.0f, 0.0f, 0.0f, 1.0f};
		mat_diffuse[3 - axis] = 1.0f;
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		float norm[4][3], ver[4][3];

		int axisCode = axisConvert123To124[4 - axis];
		if (mode == GL_SELECT)
			glPushName(axisCode);
		/* 该方向被选中则使用黄色绘制 */
		if ((m_nEditDirection & axisCode) != 0)
		{
			mat_diffuse[0] = 1.0f;
			mat_diffuse[1] = 1.0f;
			mat_diffuse[2] = 0.0f;
		}
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

		/* 绘制圆柱 */
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < 4; ++j)
				norm[j][0] = 0.0f;
			norm[0][1] = norm[3][1] = cos(2 * PI * i / n);
			norm[0][2] = norm[3][2] = sin(2 * PI * i / n);
			norm[1][1] = norm[2][1] = cos(2 * PI * (i + 1) / n);
			norm[1][2] = norm[2][2] = sin(2 * PI * (i + 1) / n);

			ver[0][0] = ver[1][0] = 0.0f;
			ver[0][1] = ver[3][1] = rAxis * cos(2 * PI * i / n);
			ver[0][2] = ver[3][2] = rAxis * sin(2 * PI * i / n);
			ver[1][1] = ver[2][1] = rAxis * cos(2 * PI * (i + 1) / n);
			ver[1][2] = ver[2][2] = rAxis * sin(2 * PI * (i + 1) / n);
			ver[2][0] = ver[3][0] = L;

			glBegin(GL_QUADS);
			for (int j = 0; j < 4; ++j)
			{
				glNormal3f(norm[j][axis % 3], norm[j][(axis + 1) % 3], norm[j][(axis + 2) % 3]);
				glVertex3f(ver[j][axis % 3], ver[j][(axis + 1) % 3], ver[j][(axis + 2) % 3]);
			}
			glEnd();
		}

		/* 绘制圆锥 */
		for (int i = 0; i < n; ++i)
		{
			/* 绘制侧面的一部分 */
			norm[0][0] = norm[1][0] = norm[2][0] = R * R / L;
			norm[0][1] = R * cos(2 * PI * i / n);
			norm[0][2] = R * sin(2 * PI * i / n);
			norm[1][1] = R * cos(2 * PI * (i + 1) / n);
			norm[1][2] = R * sin(2 * PI * (i + 1) / n);
			norm[2][1] = R * (norm[0][1] + norm[1][1]) / 2;
			norm[2][2] = R * (norm[0][2] + norm[1][2]) / 2;
			double length = sqrt(norm[0][0] * norm[0][0] + norm[0][1] * norm[0][1] + norm[0][2] * norm[0][2]);
			for (int j = 0; j < 2; ++j)
				for (int k = 0; k < 3; ++k)
					norm[j][k] /= length;
			length = sqrt(norm[2][0] * norm[2][0] + norm[2][1] * norm[2][1] + norm[2][2] * norm[2][2]);
			for (int k = 0; k < 3; ++k)
				norm[2][k] /= length;

			ver[0][0] = ver[1][0] = L;
			ver[0][1] = R * cos(2 * PI * i / n);
			ver[0][2] = R * sin(2 * PI * i / n);
			ver[1][1] = R * cos(2 * PI * (i + 1) / n);
			ver[1][2] = R * sin(2 * PI * (i + 1) / n);
			ver[2][0] = unit + L;
			ver[2][1] = ver[2][2] = 0.0f;

			glBegin(GL_TRIANGLES);
			for (int j = 0; j < 3; ++j)
			{
				glNormal3f(norm[j][axis % 3], norm[j][(axis + 1) % 3], norm[j][(axis + 2) % 3]);
				glVertex3f(ver[j][axis % 3], ver[j][(axis + 1) % 3], ver[j][(axis + 2) % 3]);
			}
			glEnd();

			/* 绘制圆的一部分 */
			ver[2][0] = L;
			norm[0][0] = norm[0][1] = norm[0][2] = 0.0f;
			norm[0][3 - axis] = -1.0f;
			glBegin(GL_TRIANGLES);
			for (int j = 3; j > 0; --j)
			{
				glNormal3f(norm[0][0], norm[0][1], norm[0][2]);
				glVertex3f(ver[j % 3][axis % 3], ver[j % 3][(axis + 1) % 3], ver[j % 3][(axis + 2) % 3]);
			}
			glEnd();
		}
		if (mode == GL_SELECT)
			glPopName();

		/* 绘制三个薄片 */
		if (drawChip)
		{
			if (mode == GL_SELECT)
			{
				if (axis == 3)						// x轴
					glPushName(Y_AXIS | Z_AXIS);
				else if (axis == 2)					// y轴
					glPushName(Z_AXIS | X_AXIS);
				else								// z轴
					glPushName(X_AXIS | Y_AXIS);
			}
			glDisable(GL_CULL_FACE);
			glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
			mat_diffuse[0] = 0.0f;
			mat_diffuse[1] = 0.0f;
			mat_diffuse[2] = 0.0f;
			mat_diffuse[3 - axis] = 1.0f;
			glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);

			ver[0][0] = 0.0f;  ver[0][1] = 0.0f;  ver[0][2] = 0.0f;
			ver[1][0] = 0.0f;  ver[1][1] = L;	  ver[1][2] = 0.0f;
			ver[2][0] = 0.0f;  ver[2][1] = L;	  ver[2][2] = L;
			ver[3][0] = 0.0f;  ver[3][1] = 0.0f;  ver[3][2] = L;
			norm[0][0] = 1.0f; norm[0][1] = 0.0f; norm[0][2] = 0.0f;
			norm[1][0] = 0.0f; norm[1][1] = 1.0f; norm[1][2] = 0.0f;
			norm[2][0] = 0.0f; norm[2][1] = 0.0f; norm[2][2] = 1.0f;
			glBegin(GL_QUADS);
				glNormal3f(norm[3 - axis][0], norm[3 - axis][1], norm[3 - axis][2]);
				for (int j = 0; j < 4; ++j)
					glVertex3f(ver[j][axis % 3], ver[j][(axis + 1) % 3],
							   ver[j][(axis + 2) % 3]);
			glEnd();
			glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
			glEnable(GL_CULL_FACE);
			if (mode == GL_SELECT)
				glPopName();
		}
	}

	/* x轴标记'X'，红色 */
	GLfloat mat_diffuse[] = {1.0, 0.0, 0.0, 1.0};
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glLineWidth(2.0f);
	glBegin(GL_LINES);
		glNormal3fv(m_fVector001);
		glVertex3f(L + unit, unit * -0.3f, 0.0f);
		glVertex3f(L + unit * 1.3f, unit * 0.3f, 0.0f);
		glVertex3f(L + unit, unit * 0.3f, 0.0f);
		glVertex3f(L + unit * 1.3f, unit * -0.3f, 0.0f);
	glEnd();

	/* y轴标记'Y'，绿色 */
	mat_diffuse[0] = 0.0f;
	mat_diffuse[1] = 1.0f;
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glBegin(GL_LINES);
		glNormal3fv(m_fVector001);
		glVertex3f(unit * -0.15f, L + unit * 1.6f, 0.0f);
		glVertex3f(0.0f, L + unit * 1.3f, 0.0f);
		glVertex3f(unit * 0.15f, L + unit * 1.6f, 0.0f);
		glVertex3f(0.0f, L + unit * 1.3f, 0.0f);
		glVertex3f(0.0f, L + unit * 1.3f, 0.0f);
		glVertex3f(0.0f, L + unit * 1.0f, 0.0f);
	glEnd();

	/* z轴标记'Z'，蓝色 */
	mat_diffuse[1] = 0.0f;
	mat_diffuse[2] = 1.0f;
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glBegin(GL_LINE_STRIP);
		glNormal3fv(m_fVector001);
		glVertex3f(unit * -0.15f, unit * 0.3f, L + unit * 1.1f);
		glVertex3f(unit * 0.15f, unit * 0.3f, L + unit * 1.1f);
		glVertex3f(unit * -0.15f, unit * -0.3f, L + unit * 1.1f);
		glVertex3f(unit * 0.15f, unit * -0.3f, L + unit * 1.1f);
	glEnd();

	/* 原点，白色球体 */
	GLfloat mat_diffuseWhite[] = {1.0, 1.0, 1.0, 1.0};
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuseWhite);

	int horN = 20, vertN = 20;
	R = unit * 0.3;

	if (mode == GL_SELECT)
		glPushName(X_AXIS | Y_AXIS | Z_AXIS);
	for (int i = 0; i < vertN; ++i)
	{
		for (int j = 0; j < horN; ++j)
		{
			glBegin(GL_POLYGON);
				glNormal3f(sin(PI * i / vertN) * cos(2 * PI * j / horN),
						   sin(PI * i / vertN) * sin(2 * PI * j / horN),
						   cos(PI * i / vertN));
				glVertex3f(sin(PI * i / vertN) * cos(2 * PI * j / horN) * R,
						   sin(PI * i / vertN) * sin(2 * PI * j / horN) * R,
						   cos(PI * i / vertN) * R);

				glNormal3f(sin(PI * (i + 1) / vertN) * cos(2 * PI * j / horN),
						   sin(PI * (i + 1) / vertN) * sin(2 * PI * j / horN),
						   cos(PI * (i + 1) / vertN));
				glVertex3f(sin(PI * (i + 1) / vertN) * cos(2 * PI * j / horN) * R,
						   sin(PI * (i + 1) / vertN) * sin(2 * PI * j / horN) * R,
						   cos(PI * (i + 1) / vertN) * R);

				glNormal3f(sin(PI * (i + 1) / vertN) * cos(2 * PI * (j + 1) / horN),
						   sin(PI * (i + 1) / vertN) * sin(2 * PI * (j + 1) / horN),
						   cos(PI * (i + 1) / vertN));
				glVertex3f(sin(PI * (i + 1) / vertN) * cos(2 * PI * (j + 1) / horN) * R,
						   sin(PI * (i + 1) / vertN) * sin(2 * PI * (j + 1) / horN) * R,
						   cos(PI * (i + 1) / vertN) * R);

				glNormal3f(sin(PI * i / vertN) * cos(2 * PI * (j + 1) / horN),
						   sin(PI * i / vertN) * sin(2 * PI * (j + 1) / horN),
						   cos(PI * i / vertN));
				glVertex3f(sin(PI * i / vertN) * cos(2 * PI * (j + 1) / horN) * R,
						   sin(PI * i / vertN) * sin(2 * PI * (j + 1) / horN) * R,
						   cos(PI * i / vertN) * R);
			glEnd();
		}
	}
	if (mode == GL_SELECT)
		glPopName();

	m_modelViewMatrixStack.popMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}

/* 绘制用于旋转控制顶点的三个圆 */
void View::drawCircle(double rightFac, double bottomFac, GLenum mode)
{
	double r_origin = m_pCommonData->length();
	double r_scale = r_origin / m_fScale;

	if (width() >= height())
	{
		GLfloat x = GLfloat(width()) / height();
		glOrtho(-x * r_scale, x * r_scale,
				-r_scale, r_scale,
				r_origin * (EYE_Z_FACTOR - 1), r_origin * (EYE_Z_FACTOR + 1));
	}
	else
	{
		GLfloat y = GLfloat(height()) / width();
		glOrtho(-r_scale, r_scale,
				-y * r_scale, y * r_scale,
				r_origin * (EYE_Z_FACTOR - 1), r_origin * (EYE_Z_FACTOR + 1));
	}
	m_modelViewMatrixStack.pushMatrix();
	m_modelViewMatrixStack.loadIdentity();
	m_modelViewMatrixStack.lookAt(0.0, 0.0, r_origin * EYE_Z_FACTOR,
								  0.0, 0.0, 0.0,
								  0.0, 1.0, 0.0);
	//m_modelViewMatrixStack.translate(m_fRight * rightFac, m_fBottom * bottomFac, length * EYE_Z_FACTOR / 2);
	m_modelViewMatrixStack.translate(m_fRight * rightFac, m_fBottom * bottomFac, 0);
	m_modelViewMatrixStack.multMatrix(lastMatrix);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadMatrixf(m_modelViewMatrixStack.top());

	double shorter = m_fRight <= -m_fBottom ? m_fRight : -m_fBottom;
	double unit = shorter * 0.1;

	glLineWidth(2.0f);
	for (int axis = 3; axis > 0; --axis)
	{
		int axisCode = axisConvert123To124[4 - axis];
		if (mode == GL_SELECT)
			glPushName(axisCode);
		GLfloat mat_diffuse[] = {0.0f, 0.0f, 0.0f, 1.0f};
		mat_diffuse[3 - axis] = 1.0f;
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

		if (m_nEditDirection == axisCode)
		{
			mat_diffuse[0] = 1.0f;
			mat_diffuse[1] = 1.0f;
			mat_diffuse[2] = 0.0f;
		}
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

		for (int i = 0; i < CIRCLEPOINTS; ++i)
		{
			float ver[2][3];
			ver[0][0] = ver[1][0] = 0.0f;
			ver[0][1] = unit * 2.0f * cos(2 * PI * i / CIRCLEPOINTS);
			ver[0][2] = unit * 2.0f * sin(2 * PI * i / CIRCLEPOINTS);
			ver[1][1] = unit * 2.0f * cos(2 * PI * (i + 1) / CIRCLEPOINTS);
			ver[1][2] = unit * 2.0f * sin(2 * PI * (i + 1) / CIRCLEPOINTS);
			if (mode == GL_SELECT)
				glPushName(i);
			glBegin(GL_LINES);
				glNormal3fv(m_fVector001);
				glVertex3f(ver[0][axis % 3], ver[0][(axis + 1) % 3], ver[0][(axis + 2) % 3]);
				glVertex3f(ver[1][axis % 3], ver[1][(axis + 1) % 3], ver[1][(axis + 2) % 3]);
			glEnd();
			if (mode == GL_SELECT)
				glPopName();
		}
		if (mode == GL_SELECT)
			glPopName();
	}
	m_modelViewMatrixStack.popMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}

/* 改变所有选中顶点的状态（选取缓存的具体格式参考 OpenGL 编程指南） */
void View::processHits(GLint hits)
{
	GLuint *ptr = 0, names = 0, x = 0, y = 0, z = 0;
	ptr = (GLuint *)selectBuf;
	for (int i = 0; i < hits; ++i)
	{
		names = *ptr;
		ptr += 3;
		if (names == 1)			// 选中的是直接编辑点
		{
			m_pCommonData->setDirectPointSelected(*ptr++ - DIRECTPOINTSELECTBASE);
		}
		else					// 选中的是控制顶点
		{
			for (unsigned int j = 0; j < names; ++j)
			{
				if (j == 0)
					z = *ptr;
				else if (j == 1)
					y = *ptr;
				else if (j == 2)
					x = *ptr;
				ptr++;
			}
			m_pCommonData->setCtrlPointSelected(x, y, z);	// 记录该顶点被选中
		}
		/*for (unsigned int j = 0; j < names; ++j)
		{
			if (j == 0)
				z = *ptr;
			else if (j == 1)
				y = *ptr;
			else if (j == 2)
				x = *ptr;
			ptr++;
		}
		m_pCommonData->setCtrlPointSelected(x, y, z);	// 记录该顶点被选中*/
	}
}

/* 记录选中的坐标轴 */
void View::processHitsAxis(GLint hits)
{
	//cout << "hits = " << hits << endl;
	//cout << "processAxis" << endl;
	GLuint *ptr = 0;
	ptr = (GLuint *)selectBufAxis;
	unsigned int near = 0;
	m_nEditDirection = NONE_AXIS;
	for (int i = 0; i < hits; ++i)
	{
		if (*ptr)
		{
			//cout << "nameCount = " << *ptr << endl;
			//cout << "ptr + 1 = " << *(ptr + 1) << endl;
			//cout << "ptr + 2 = " << *(ptr + 2) << endl;
			if (i == 0 || *(ptr + 1) < near)
			{
				near = *(ptr + 1);
				m_nEditDirection = *(ptr + 3);
				/*if (m_nSelectedAxis == X_SELECTED)
					m_pCommonData->setCtrlPointTranslateDirection(X_AXIS);
				else if (m_nSelectedAxis == Y_SELECTED)
					m_pCommonData->setCtrlPointTranslateDirection(Y_AXIS);
				else
					m_pCommonData->setCtrlPointTranslateDirection(Z_AXIS);*/
			}
			ptr += 4;
			//m_nSelectedAxis = *ptr;
			//cout << "selected = " << m_nSelectedAxis << endl;
		}
		else
			ptr += 3;
	}
}

/* 记录选中的旋转圆 */
void View::processHitsCircle(GLint hits)
{
	//cout << "processCircle" << endl;
	GLuint *ptr = 0;
	ptr = (GLuint *)selectBufAxis;
	unsigned int near = 0;
	m_nEditDirection = NONE_AXIS;
	for (int i = 0; i < hits; ++i)
	{
		//cout << "nameCount = " << *ptr << endl;
		//cout << "ptr + 1 = " << *(ptr + 1) << endl;
		//cout << "ptr + 2 = " << *(ptr + 2) << endl;
		if (i == 0 || *(ptr + 1) < near)
		{
			near = *(ptr + 1);
			m_nEditDirection = *(ptr + 3);
				//cout << "edit改变" << endl;
			//cout << "direction = " << m_nEditDirection << endl;
			m_nLineSeqNum = *(ptr + 4);
			//cout << "第" << m_nLineSeqNum << "段" << endl;
		}
		ptr += 5;
	}
}

void View::processHitsOriginalModel(GLint hits)
{
	GLuint *ptr = 0;
	ptr = (GLuint *)selectBufOriginalModel;
	unsigned int near = 0;
	for (int i = 0; i < hits; ++i)
	{
		if (i == 0 || *(ptr + 1) < near)
		{
			near = *(ptr + 1);
			m_nSelectedOriginalFace = *(ptr + 3);
		}
		ptr += 4;
	}
}

/* 计算选取框的范围 */
void View::calcSelectionRange(double &centerX, double &centerY, double &selectionWidth, double &selectionHeight)
{
	int smallX, smallY, bigX, bigY;	
	if (m_nCtrlSelectStartX > m_nCtrlSelectEndX)
	{
		smallX = m_nCtrlSelectEndX;
		bigX = m_nCtrlSelectStartX;
	}
	else
	{
		smallX = m_nCtrlSelectStartX;
		bigX = m_nCtrlSelectEndX;
	}
	if (m_nCtrlSelectStartY > m_nCtrlSelectEndY)
	{
		smallY = m_nCtrlSelectEndY;
		bigY = m_nCtrlSelectStartY;
	}
	else
	{
		smallY = m_nCtrlSelectStartY;
		bigY = m_nCtrlSelectEndY;
	}
	selectionWidth = (double)(bigX - smallX);
	selectionHeight = (double)(bigY - smallY);

	centerX = (double)smallX + selectionWidth / 2;
	centerY = (double)smallY + selectionHeight / 2;
}

/* 鼠标点击的消息响应函数 */
void View::mousePressEvent(QMouseEvent *event)
{
	if (!m_bLoaded)
		return;
	if (event->button() == Qt::LeftButton)				// 按下鼠标左键
	{
		m_eMouseButton = LEFT_BUTTON;					// 记录
		if (m_nEditDirection == NONE_AXIS)
		{
			if (m_bCtrlPressed)							// 如果按下了Ctrl键，表示要平移模型
			{
				m_nCursorX = event->x();
				m_nCursorY = event->y();
			}
			else if (m_bShiftPressed)					// 如果按下了Shift键，表示旋转模型和cubemap
			{
				calcSphereCoord(event->x(), event->y(), lastPos);
				calcSphereCoord(event->x(), event->y(), lastPosCubeMap);
			}
			else	// 没有按下任何键，表示只旋转模型，记录鼠标按下时的坐标，并将其转换成半球坐标
			{
				calcSphereCoord(event->x(), event->y(), lastPos);
			}
		}
		else											// 处于控制顶点编辑模式
		{
			m_nCursorX = event->x();
			m_nCursorY = event->y();
		}
	}
	else if (event->button() == Qt::RightButton)		// 按下鼠标右键，表示要生成一个选取框
	{
		m_eMouseButton = RIGHT_BUTTON;					// 记录

		/* 开始和结束坐标都记为当前鼠标坐标 */
		m_nCtrlSelectStartX = event->x();
		m_nCtrlSelectStartY = event->y();
		m_nCtrlSelectEndX = event->x();
		m_nCtrlSelectEndY = event->y();
	}
	else		// 按下滚轮，表示新建一个直接编辑点
	{
		double x = (double)(2 * event->x()) / (width() - 1) - 1;
		double y = 1 - (double)(2 * event->y()) / (height() - 1);

		const float *topMatrixPtr = m_projectionMatrixStack.top();
		double arrayP[16];
		for (int i = 0; i < 16; ++i)
			arrayP[i] = topMatrixPtr[i];
		LaGenMatDouble P(arrayP, 4, 4);
		//cout << "P:\n" << P << endl;

		topMatrixPtr = m_modelViewMatrixStack.top();
		double arrayM[16];
		for (int i = 0; i < 16; ++i)
			arrayM[i] = topMatrixPtr[i];
		LaGenMatDouble M(arrayM, 4, 4);
		//cout << "M:\n" << M << endl;

		double arrayPM[16] = {0.0};
		LaGenMatDouble PM(arrayPM, 4, 4);
		Blas_Mat_Mat_Mult(P, M, PM);
		//cout << "PM:\n" << PM << endl;

		double arrayPointXY[4] = {x, y, -1, 1};		// 投影变换后能投到(x, y)的最靠近视点的点
		LaGenMatDouble pointXY(arrayPointXY, 4, 1);
		//cout << "pointXY:\n" << pointXY << endl;

		LaGenMatDouble X(4, 1);
		LaLinearSolve(PM, X, pointXY);
		//cout << "X:\n" << X << endl;

		/* 求鼠标选中的面片 */
		GLint viewport[4];
		glGetIntegerv(GL_VIEWPORT, viewport);			// 得到视口
		glSelectBuffer(1000, selectBufOriginalModel);
		glRenderMode(GL_SELECT);

		glInitNames();
		//glPushName(0);

		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();

		gluPickMatrix((GLdouble)event->x(), (GLdouble)(viewport[3] - event->y()), 1.0, 1.0, viewport);
		glMultMatrixf(m_projectionMatrixStack.top());

		unsigned int faceCount = m_pCommonData->faceCount();
		for (unsigned int i = 0; i < faceCount; ++i)
		{
			glPushName(i);
			unsigned int totalVertex = m_pCommonData->vertexCoordIdxCount(i);
			glBegin(GL_POLYGON);
			for (unsigned int j = 0; j < totalVertex; ++j)
			{
				int index = m_pCommonData->vertexCoordIdx(i, j);
				VertexCoord v = m_pCommonData->vertexCoord(index);
				glVertex3d((GLdouble)v.x(), (GLdouble)v.y(), (GLdouble)v.z());
			}
			glEnd();
			glPopName();
		}

		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		glFlush();

		GLint hits = glRenderMode(GL_RENDER);
		//cout << "hits = " << hits << endl;
		if (hits)
		{
			processHitsOriginalModel(hits);
			//cout << "faceIdx = " << m_nSelectedOriginalFace << endl;

			int index0 = m_pCommonData->vertexCoordIdx(m_nSelectedOriginalFace, 0);
			int index1 = m_pCommonData->vertexCoordIdx(m_nSelectedOriginalFace, 1);
			int index2 = m_pCommonData->vertexCoordIdx(m_nSelectedOriginalFace, 2);
			VertexCoord v0 = m_pCommonData->vertexCoord(index0);
			VertexCoord v1 = m_pCommonData->vertexCoord(index1);
			VertexCoord v2 = m_pCommonData->vertexCoord(index2);
			//cout << "v0:(" << v0.x() << ", " << v0.y() << ", " << v0.z() << ")" << endl;
			//cout << "v1:(" << v1.x() << ", " << v1.y() << ", " << v1.z() << ")" << endl;
			//cout << "v2:(" << v2.x() << ", " << v2.y() << ", " << v2.z() << ")" << endl;

			VertexCoord v01 = v1 - v0;
			VertexCoord v02 = v2 - v0;
			//VertexCoord N = v01.cross(v02);
			VertexCoord N = cross(v01, v02);
			N.normalize();

			double d = -(N * v0);

			VertexCoord O(X(0, 0), X(1, 0), X(2, 0));
			VertexCoord V(-m_fVector001[0], -m_fVector001[1], -m_fVector001[2]);

			double t = -(N * O + d) / (N * V);

			DirectPoint result(O + V * t);
			//cout << "result:(" << result.x() << ", " << result.y() << ", " << result.z() << ")" << endl;
			m_pCommonData->directPointVectorPush(result);
			updateGL();
		}
	}
}

/* 鼠标拖动的消息响应函数 */
void View::mouseMoveEvent(QMouseEvent *event)
{
	if (!m_bLoaded)
		return;
	if (m_eMouseButton == LEFT_BUTTON)								// 鼠标左键被按下
	{
		//cout << "left" << endl;
		//cout << "editDirection = " << m_nEditDirection << endl;
		if (m_nEditDirection == NONE_AXIS)							// 不对控制定点进行编辑
		{
			if (m_bCtrlPressed)										// Ctrl键被按下，此时要平移模型
			{
				//cout << "ctrl" << endl;
				m_fDeltaX += (double)event->x() - m_nCursorX;
				m_fDeltaY -= (double)event->y() - m_nCursorY;
				m_nCursorX = event->x();
				m_nCursorY = event->y();
			}
			else if (m_bShiftPressed)								// Shift键被按下，此时要旋转模型和cubemap
			{
				//cout << "shift" << endl;
				calcSphereCoord(event->x(), event->y(), curPos);	// 将当前鼠标平面坐标转换为半球面坐标
				calcThetaAndRotateAxis(curPos, lastPos, m_fTheta, rotateAxis);
				calcSphereCoord(event->x(), event->y(), curPosCubeMap);	// 将当前鼠标平面坐标转换为半球面坐标
				calcThetaAndRotateAxis(curPosCubeMap, lastPosCubeMap, m_fThetaCubeMap, rotateAxisCubeMap);
			}
			else													// 没有按键，只旋转模型
			{
				//cout << m_modelViewMatrixStack << endl;
				calcSphereCoord(event->x(), event->y(), curPos);	// 将当前鼠标平面坐标转换为半球面坐标
				calcThetaAndRotateAxis(curPos, lastPos, m_fTheta, rotateAxis);
			}
		}
		else													// 对控制顶点进行编辑
		{
			//m_calcTime.start();
			double deltaX = (double)event->x() - m_nCursorX;
			double deltaY = (double)m_nCursorY - event->y();
			m_nCursorX = event->x();
			m_nCursorY = event->y();
			double delta = 0.0;
			if (m_eEditMode == TRANSLATE)						// 处于平移控制顶点模式
			{
				double theta;
				if (width() > height())
					theta = 2.0 * m_pCommonData->length() / (m_fScale * height());
				else
					theta = 2.0 * m_pCommonData->length() / (m_fScale * width());
				double count = 0.0;
				if (m_nEditDirection & X_AXIS)
				{
					count++;
					delta += deltaX * lastMatrix[0][0] + deltaY * lastMatrix[0][1];
				}
				if (m_nEditDirection & Y_AXIS)
				{
					count++;
					delta += deltaX * lastMatrix[1][0] + deltaY * lastMatrix[1][1];
				}
				if (m_nEditDirection & Z_AXIS)
				{
					count++;
					delta += deltaX * lastMatrix[2][0] + deltaY * lastMatrix[2][1];
				}
				delta /= sqrt(count);
				if (m_pCommonData->directPointSelectedCount())
					m_pCommonData->directPointTranslate(delta * theta, m_nEditDirection);
				else
					m_pCommonData->ctrlPointTranslate(delta * theta, m_nEditDirection);
			}
			else if (m_eEditMode == ROTATE)					// 处于旋转控制顶点模式
			{
				double cosDelta = cos(2 * PI * (m_nLineSeqNum + 1) / CIRCLEPOINTS) -
								  cos(2 * PI * m_nLineSeqNum / CIRCLEPOINTS);
				double sinDelta = sin(2 * PI * (m_nLineSeqNum + 1) / CIRCLEPOINTS) -
								  sin(2 * PI * m_nLineSeqNum / CIRCLEPOINTS);

				if (m_nEditDirection == X_AXIS)
				{
					delta = deltaX * (cosDelta * lastMatrix[1][0] + sinDelta * lastMatrix[2][0]) +
							deltaY * (cosDelta * lastMatrix[1][1] + sinDelta * lastMatrix[2][1]);
				}
				else if (m_nEditDirection == Y_AXIS)
				{
					delta = deltaX * (sinDelta * lastMatrix[0][0] + cosDelta * lastMatrix[2][0]) +
							deltaY * (sinDelta * lastMatrix[0][1] + cosDelta * lastMatrix[2][1]);
				}
				else if (m_nEditDirection == Z_AXIS)
				{
					delta = deltaX * (cosDelta * lastMatrix[0][0] + sinDelta * lastMatrix[1][0]) +
							deltaY * (cosDelta * lastMatrix[0][1] + sinDelta * lastMatrix[1][1]);
				}
				//double r = delta * 2 * PI / sqrt(width() * width() + height() * height());
				double r = delta * 2 * PI / width() * 10;
				//cout << "\tr = " << r << endl;
				//cout << "\tdirection = " << m_nEditDirection << endl;
				if (m_pCommonData->directPointSelectedCount())
					m_pCommonData->directPointRotate(r, m_nEditDirection);
				else
					m_pCommonData->ctrlPointRotate(r, m_nEditDirection);
			}
			else if (m_eEditMode == SCALE)					// 处于放缩控制顶点模式
			{
				double count = 0.0;
				if (m_nEditDirection & X_AXIS)
				{
					count++;
					delta += deltaX * lastMatrix[0][0] + deltaY * lastMatrix[0][1];
				}
				if (m_nEditDirection & Y_AXIS)
				{
					count++;
					delta += deltaX * lastMatrix[1][0] + deltaY * lastMatrix[1][1];
				}
				if (m_nEditDirection & Z_AXIS)
				{
					count++;
					delta += deltaX * lastMatrix[2][0] + deltaY * lastMatrix[2][1];
				}
				delta /= sqrt(count);
				double factor = delta / sqrt(width() * width() + height() * height());
				if (m_pCommonData->directPointSelectedCount())
					m_pCommonData->directPointScale(factor, m_nEditDirection);
				else
					m_pCommonData->ctrlPointScale(factor, m_nEditDirection);
			}
		}
		updateGL();
	}
	else if (m_eMouseButton == RIGHT_BUTTON)	// 鼠标右键被按下，此时处于生成选取框模式
	{
		m_nCtrlSelectEndX = event->x();
		m_nCtrlSelectEndY = event->y();
		updateGL();
	}
	/* 没有按下任何鼠标键，此时如果显示了编辑图标，则需要检测鼠标选中了哪个坐标轴或者旋转圆 */
	else if ((m_nDisplayMode & DISPLAY_EDIT_ICON) != 0)
	{
		GLint viewport[4];
		glGetIntegerv(GL_VIEWPORT, viewport);		// 得到视口
		glSelectBuffer(1200, selectBufAxis);
		glRenderMode(GL_SELECT);

		glInitNames();
		//glPushName(0);

		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();

		if (m_eEditMode == TRANSLATE || m_eEditMode == SCALE)
			gluPickMatrix((GLdouble)event->x(), (GLdouble)(viewport[3] - event->y()), 1.0, 1.0, viewport);
		else if (m_eEditMode == ROTATE)
			gluPickMatrix((GLdouble)event->x(), (GLdouble)(viewport[3] - event->y()), 3.0, 3.0, viewport);
		//glMultMatrixf(m_projectionMatrixStack.top());
		if (m_eEditMode == TRANSLATE || m_eEditMode == SCALE)
			drawAxis(RIGHTFAC, BOTTOMFAC, GL_SELECT, true);
		else if (m_eEditMode == ROTATE)
			drawCircle(RIGHTFAC, BOTTOMFAC, GL_SELECT);

		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		glFlush();

		GLint hits = glRenderMode(GL_RENDER);
		//cout << "hits = " << hits << endl;
		if (m_eEditMode == TRANSLATE || m_eEditMode == SCALE)
			processHitsAxis(hits);
		else
			processHitsCircle(hits);
		updateGL();
	}
}

/* 鼠标按键被放开的消息相应函数 */
void View::mouseReleaseEvent(QMouseEvent *event)
{
	if (!m_bLoaded)
		return;
	m_eMouseButton = NO_BUTTON;					// 记录此状态
	if (event->button() == Qt::RightButton)		// 放开的是右键，此时需要计算选中的控制顶点
	{
		m_nCtrlSelectEndX = event->x();
		m_nCtrlSelectEndY = event->y();

		double centerX, centerY, selectionWidth, selectionHeight;
		calcSelectionRange(centerX, centerY, selectionWidth, selectionHeight);

		if (selectionWidth > ZERO && selectionHeight > ZERO)
		{
			GLint viewport[4];
			glGetIntegerv(GL_VIEWPORT, viewport);	// 得到视口
			glSelectBuffer(MAXCTRLPOINT * MAXCTRLPOINT * MAXCTRLPOINT * 6, selectBuf);
			glRenderMode(GL_SELECT);

			glInitNames();
			glPushName(0);

			glMatrixMode(GL_PROJECTION);
			glPushMatrix();
			glLoadIdentity();

			gluPickMatrix((GLdouble)centerX, (GLdouble)(viewport[3] - centerY), selectionWidth, selectionHeight, viewport);
			glMultMatrixf(m_projectionMatrixStack.top());
			if ((m_nDisplayMode & DISPLAY_CTRL_POINTS) != 0)
				drawCtrlPoint(GL_SELECT);
			if ((m_nDisplayMode & DISPLAY_DIRECT_POINT) != 0)
				drawDirectPoint(GL_SELECT);

			glMatrixMode(GL_PROJECTION);
			glPopMatrix();
			glFlush();

			GLint hits = glRenderMode(GL_RENDER);
			processHits(hits);
		}
	}
	updateGL();
}

/* 使用鼠标指针的平面坐标计算它的半球面坐标 */
void View::calcSphereCoord(int x, int y, double *pos)
{
	double halfW, halfH, r;
	halfW = (double)width() / 2;					// 窗口宽度的一半
	halfH = (double)height() / 2;					// 窗口高度的一半
	r = sqrt(halfW * halfW + halfH * halfH);		// 窗口对角线长度的一半

	/* 防止鼠标移到窗口外时出现异常情况 */
	if (x > width())
		x = width();
	if (x < 0)
		x = 0;
	if (y > height())
		y = height();
	if (y < 0)
		y = 0;

	/* 算出x, y, z三个方向的坐标 */
	pos[0] = (double)x - halfW;
	pos[1] = halfH - y;
	pos[2] = sqrt(halfW * halfW + halfH * halfH + 1 - pos[0] * pos[0] - pos[1] * pos[1]);

	/* 单位化 */
	pos[0] /= r;
	pos[1] /= r;
	pos[2] /= r;
}

/* 计算旋转角度和旋转轴 */
void View::calcThetaAndRotateAxis(double *cPos, double *lPos, double &theta, double *rotAxis)
{
	/* 计算出和前一个坐标之间x, y, z方向的差值 */
	double dx = cPos[0] - lPos[0];
	double dy = cPos[1] - lPos[1];
	double dz = cPos[2] - lPos[2];

	double d = sqrt(dx * dx + dy * dy + dz * dz);				// 计算出鼠标在半球面上移动的距离
	/* 计算出相应的圆心角，这里使用弦长近似代替弧长 */
	theta = d * 180.0;

	/* 计算旋转轴，是两个向量的叉积 */
	rotAxis[0] = lPos[1] * cPos[2] - lPos[2] * cPos[1];
	rotAxis[1] = lPos[2] * cPos[0] - lPos[0] * cPos[2];
	rotAxis[2] = lPos[0] * cPos[1] - lPos[1] * cPos[0];

	/* 将上一个半球面坐标替换为当前半球面坐标 */
	lPos[0] = cPos[0];
	lPos[1] = cPos[1];
	lPos[2] = cPos[2];
}

/* 鼠标滚轮的消息响应函数（模型的放大缩小） */
void View::wheelEvent(QWheelEvent *event)
{
	if (event->orientation() == Qt::Vertical)
	{
		//m_fScale *= (1.0 + 0.0003 * event->delta());
		m_fScale += 0.0003 * event->delta();
		cout << "m_fScale = " << m_fScale << endl;
		setProjection();
		updateGL();
	}
}

/* 按下键的消息响应函数 */
void View::keyPressEvent(QKeyEvent *event)
{
	if (event->key() == Qt::Key_Control)
		m_bCtrlPressed = true;
	else if (event->key() == Qt::Key_Shift)
		m_bShiftPressed = true;
	else if (event->key() == Qt::Key_A)
	{
		m_bLocalViewer = !m_bLocalViewer;
		updateGL();
	}
	else if (event->key() == Qt::Key_D)
	{
		//m_fDist--;
		//cout << "dist = " << m_fDist << endl;
		light_position[2]--;
		cout << "dist = " << light_position[2] << endl;
		updateGL();
	}
	else if (event->key() == Qt::Key_E)
	{
		//m_fDist++;
		//cout << "dist = " << m_fDist << endl;
		light_position[2]++;
		cout << "dist = " << light_position[2] << endl;
		updateGL();
	}
	else if (event->key() == Qt::Key_M)
	{
		m_pCommonData->testMove(true);
		updateGL();
	}
	else if (event->key() == Qt::Key_N)
	{
		m_pCommonData->testMove(false);
		updateGL();
	}
	else if (event->key() == Qt::Key_1)
	{
		m_pCommonData->testMoveTorus(true);
		updateGL();
	}
	else if (event->key() == Qt::Key_2)
	{
		m_pCommonData->testMoveTorus(false);
		updateGL();
	}
	else if (event->key() == Qt::Key_P)
	{
		if (m_bIsOrtho)
			setOrtho(false);
		else
			setOrtho(true);
		updateGL();
	}
	else if (event->key() == Qt::Key_R)
	{
		m_pCommonData->testRotate(false);
		updateGL();
	}
	else if (event->key() == Qt::Key_Q)
	{
		m_pParent->close();
	}
	else if (event->key() == Qt::Key_Plus)
	{
		m_pCommonData->increaseCenterCtrlPointFactor();
		updateGL();
	}
	else if (event->key() == Qt::Key_Minus)
	{
		m_pCommonData->decreaseCenterCtrlPointFactor();
		updateGL();
	}
	else if (event->key() == Qt::Key_0)
	{
		m_pCommonData->cancelAllSelection();
		updateGL();
	}
	else if (event->key() == Qt::Key_V)
	{
		m_pParent->setUsePN();
	}
	else if (event->key() == Qt::Key_C)
	{
		m_pParent->changeAlgorithm();
	}
#ifdef DRAW_TRIANGULAR_CTRL_POINTS
	else if (event->key() == Qt::Key_Z)
	{
		++draw_which_;
		updateGL();
	}
#endif
//#ifdef TRUTH
	else if (event->key() == Qt::Key_T)
	{
		if (m_eErrorDisplayMode == MESH)
			m_eErrorDisplayMode = VERTEX_ERROR;
		else if (m_eErrorDisplayMode == VERTEX_ERROR)
			m_eErrorDisplayMode = NORMAL_ERROR;
		else
			m_eErrorDisplayMode = MESH;
		m_pParent->changeDisplayMode(m_eErrorDisplayMode);
		updateGL();
	}
	else if (event->key() == Qt::Key_W)
	{
		if (error_texture_)
			error_texture_ = false;
		else
			error_texture_ = true;
		m_pParent->changeErrorTexture(error_texture_);
		updateGL();
	}
	else if (event->key() == Qt::Key_O)
	{
		if (m_bDisplayTruth)
			m_bDisplayTruth = false;
		else
			m_bDisplayTruth = true;
		m_pParent->changeDisplayTruth(m_bDisplayTruth);
		updateGL();
	}
//#endif

#ifdef LINE
	else if (event->key() == Qt::Key_L)
	{
		if (m_bLine)
			m_bLine = false;
		else
			m_bLine = true;
		updateGL();
	}
#endif
}

/* 放开键的消息响应函数 */
void View::keyReleaseEvent(QKeyEvent *event)
{
	if (event->key() == Qt::Key_Control)
		m_bCtrlPressed = false;
	else if (event->key() == Qt::Key_Shift)
		m_bShiftPressed = false;
}

/* 将顶点编辑模式设为平移模式 */
void View::setModeCtrlPointTrans(bool checked)
{
	if (checked)
	{
		m_eEditMode = TRANSLATE;
		if (!m_pCommonData->loadingEdit())
			updateGL();
	}
}

/* 将顶点编辑模式设为旋转模式 */
void View::setModeCtrlPointRotate(bool checked)
{
	if (checked)
	{
		m_eEditMode = ROTATE;
		if (!m_pCommonData->loadingEdit())
			updateGL();
	}
}

/* 将顶点编辑模式设为缩放模式 */
void View::setModeCtrlPointScale(bool checked)
{
	if (checked)
	{
		m_eEditMode = SCALE;
		if (!m_pCommonData->loadingEdit())
			updateGL();
	}
}

/* 设置当前显示模型 */
void View::setDisplayModel()
{
	m_nDisplayMode ^= DISPLAY_MODEL;
	if (!m_pCommonData->loadingEdit())
		updateGL();
}

/* 设置当前显示控制网 */
void View::setDisplayCtrlPoint()
{
	m_nDisplayMode ^= DISPLAY_CTRL_POINTS;
	if (!m_pCommonData->loadingEdit())
		updateGL();
}

/* 设置当前显示直接编辑点 */
void View::setDisplayDirectPoint()
{
	m_nDisplayMode ^= DISPLAY_DIRECT_POINT;
	if (!m_pCommonData->loadingEdit())
		updateGL();
}

void View::setDisplayEditIcon()
{
	m_nDisplayMode ^= DISPLAY_EDIT_ICON;
	if (!m_pCommonData->loadingEdit())
		updateGL();
}

void View::setDisplayTriangularCP()
{
	m_nDisplayMode ^= DISPLAY_TRIANGULAR_CTRL_POINTS;
	if (!m_pCommonData->loadingEdit())
		updateGL();
}

/* 设置当前显示节点 */
void View::setDisplayKnot()
{
	m_nDisplayMode ^= DISPLAY_KNOT;
	if (!m_pCommonData->loadingEdit())
		updateGL();
}

/* 设置绘制模式为GL_FLAT */
void View::setRenderFlat(bool checked)
{
	if (checked)
		m_bIsFlat = true;
	if (!m_pCommonData->loadingEdit())
		updateGL();
}

/* 设置绘制模式为GL_SMOOTH */
void View::setRenderSmooth(bool checked)
{
	if (checked)
		m_bIsFlat = false;
	if (!m_pCommonData->loadingEdit())
		updateGL();
}

void View::setOrtho(bool is_ortho)
{
	if (is_ortho)
		m_bIsOrtho = true;
	else
		m_bIsOrtho = false;
	setProjection();
	updateGL();
}


void View::changeTex(const QString &texName)
{
	string choice(texName.toUtf8().constData());
	//cout << "texChanged: " << choice << endl;
	if (choice == string("无纹理"))
	{
		m_nTexCase = 0;
	}
	else if (choice == string("自带2D纹理"))
	{
		m_nTexCase = 2;
		load2DTexture();
	}
	else
	{
		m_nTexCase = 3;
		load3DTexture(choice);
	}
	if (m_bLoaded && !m_pCommonData->loadingEdit())
		updateGL();
}

void View::changeCubeMap(const QString &cubeMapName)
{
	string choice(cubeMapName.toUtf8().constData());
	//cout << "cubeMapChanged: " << choice << endl;
	if (choice == string("无纹理"))
	{
		m_bUseCubeMap = false;
	}
	else
	{
		m_bUseCubeMap = true;
		loadCubeMap(QString(choice.c_str()));
	}
	if (m_bLoaded && !m_pCommonData->loadingEdit())
		updateGL();
}

void View::setEnvMap(bool checked)
{
	if (checked)
		m_bUseEnvMap = true;
	if (!m_pCommonData->loadingEdit())
		updateGL();
}

void View::setNoEnvMap(bool checked)
{
	if (checked)
		m_bUseEnvMap = false;
	if (!m_pCommonData->loadingEdit())
		updateGL();
}

GLint View::getUniLoc(GLuint program, const GLchar *name)
{
	GLint loc;
	loc = glGetUniformLocation(program, name);
	if (loc == -1)
		cout << "uniform 变量 \"" << name << "\" 未定义" << endl;
	//printOpenGLError(); // Check for OpenGL errors

	return loc;
}

void View::printShaderInfoLog(GLuint shader, string shader_name)
{
	int infologLen = 0;
	int charsWritten = 0;
	GLchar *infoLog;
	glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infologLen);
	if (infologLen > 0)
	{
		infoLog = new GLchar[infologLen];
		if (infoLog == NULL)
		{
			cout << "ERROR: 无法创建 shader InfoLog 缓存" << endl;
			exit(1);
		}
		glGetShaderInfoLog(shader, infologLen, &charsWritten, infoLog);
		cout << shader_name << " 编译信息:\n" << infoLog << "\n" << endl;
		delete []infoLog;
	}
}

void View::printProgramInfoLog(GLuint program)
{
	int infologLen = 0;
	int charsWritten = 0;
	GLchar *infoLog;
	glGetProgramiv(program, GL_INFO_LOG_LENGTH, &infologLen);
	if (infologLen > 0)
	{
		infoLog = new GLchar[infologLen];
		if (infoLog == NULL)
		{
			cout << "ERROR: 无法创建 program InfoLog 缓存" << endl;
			exit(1);
		}
		glGetProgramInfoLog(program, infologLen, &charsWritten, infoLog);
		cout << "程序链接信息:\n" << infoLog << "\n" << endl;
		delete []infoLog;
	}
}

int View::readShaderSource(const char *fileName, string &vertexShader, string &fragmentShader)
{
	string vsFileName = string("../") + fileName + ".vs";
	string fsFileName = string("../") + fileName + ".fs";

	std::ifstream inFile(vsFileName.c_str());
	string line;
	while(getline(inFile, line))
	{
		line += '\n';
		vertexShader += line;
	}
	inFile.close();

	inFile.open(fsFileName.c_str());
	while(getline(inFile, line))
	{
		line += '\n';
		fragmentShader += line;
	}
	inFile.close();
 
	return 1;
}

/* 安装shader */
bool View::installShaders(const GLchar *vertexShader, const GLchar *fragmentShader)
{
	vs = glCreateShader(GL_VERTEX_SHADER);
	fs = glCreateShader(GL_FRAGMENT_SHADER);

	glShaderSource(vs, 1, &vertexShader, NULL);
	glShaderSource(fs, 1, &fragmentShader, NULL);

	GLint vert_status, frag_status;

	glCompileShader(vs);
	glGetShaderiv(vs, GL_COMPILE_STATUS, &vert_status);
	if (vert_status == GL_FALSE)
		printShaderInfoLog(vs, "vertex shader");

	glCompileShader(fs);
	glGetShaderiv(fs, GL_COMPILE_STATUS, &frag_status);
	if (frag_status == GL_FALSE)
		printShaderInfoLog(fs, "fragment shader");

	if ((vert_status == GL_FALSE) || (frag_status == GL_FALSE))
		return 0;

	prog = glCreateProgram();
	glAttachShader(prog, vs);
	glAttachShader(prog, fs);

	glLinkProgram(prog);
	GLint link_status;
	glGetProgramiv(prog, GL_LINK_STATUS, &link_status);
	if (link_status == GL_FALSE)
	{
		cout << "程序链接失败！" << endl;
		printProgramInfoLog(prog);
		return false;
	}

	glUseProgram(prog);

	glUniform1f(getUniLoc(prog, "PointLightSource0.constantAttenuation"), 1.0);
	glUniform1f(getUniLoc(prog, "PointLightSource0.linearAttenuation"), 0.0);
	glUniform1f(getUniLoc(prog, "PointLightSource0.quadraticAttenuation"), 0.0);
	glUniform4f(getUniLoc(prog, "PointLightSource0.ambient"), 0.5, 0.5, 0.5, 1.0);
	glUniform4f(getUniLoc(prog, "PointLightSource0.diffuse"), 0.8, 0.8, 0.8, 1.0);
	glUniform4f(getUniLoc(prog, "PointLightSource0.specular"), 1.0, 1.0, 1.0, 1.0);
	//glUniform4f(getUniLoc(prog, "PointLightSource0.specular"), 0.4, 0.4, 0.4, 1.0);

	glUniform1i(getUniLoc(prog, "TexSamplerCube"), 0);
//#ifdef TRUTH
	glUniform1i(getUniLoc(prog, "TexSamplerColormap"), 1);
//#endif
	glUniform1i(getUniLoc(prog, "TexSampler2D"), 2);
	glUniform1i(getUniLoc(prog, "TexSampler3D"), 3);

	tex_case_ = getUniLoc(prog, "TexCase");

	MVMatrix_id_ = getUniLoc(prog, "MVMatrix");
	PMatrix_id_ = getUniLoc(prog, "PMatrix");
	NormalMatrix_id_ = getUniLoc(prog, "NormalMatrix");
	cubemap_rot_inv_mat_id_ = getUniLoc(prog, "CubeMapRotInvMat");

	local_viewer_id_ = getUniLoc(prog, "LocalViewer");
	light_position_id_ = getUniLoc(prog, "PointLightSource0.position");
	use_env_map_id_ = getUniLoc(prog, "UseEnvMap");

//#ifdef TRUTH
	display_truth_id_ = getUniLoc(prog, "DisplayTruth");
	error_texture_id_ = getUniLoc(prog, "ErrorDisplayMode");
	error_display_mode_id_ = getUniLoc(prog, "ErrorDisplayMode");
//#endif

#ifdef LINE
	divide_id_ = getUniLoc(prog, "Divide");
	use_line_id_ = getUniLoc(prog, "UseLine");
#endif

	ka_id_ = getUniLoc(prog, "Mtl.ka");
	kd_id_ = getUniLoc(prog, "Mtl.kd");
	ks_id_ = getUniLoc(prog, "Mtl.ks");

	//min_vertex_id_ = getUniLoc(prog, "min_vertex");
	//delta_vertex_inverse_id_ = getUniLoc(prog, "delta_vertex_inverse");

	glUseProgram(0);

	return true;
}
