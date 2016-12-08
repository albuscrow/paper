#ifndef __VIEW_H__
#define __VIEW_H__

#define GL_GLEXT_PROTOTYPES

#include <QtOpenGL>
#include <QGLWidget>
#include <GL/glu.h>
#include <GL/glext.h>
#include <string>
#include <fstream>
#include "common_data.h"
#include "matrix_stack.h"

/*
 * 鼠标键
 */
enum MouseButton
{
	NO_BUTTON,
	LEFT_BUTTON,
	RIGHT_BUTTON
};

/*
 * 显示模式
 */
enum DisplayMode
{
	DISPLAY_MODEL = 1,				// 显示模型
	DISPLAY_CTRL_POINTS = 2,		// 显示控制网
	DISPLAY_KNOT = 4,				// 显示节点
	DISPLAY_DIRECT_POINT = 8,		// 显示直接编辑点
	DISPLAY_EDIT_ICON = 16,			// 显示编辑图标
	DISPLAY_TRIANGULAR_CTRL_POINTS = 32	// 显示三角Bezier曲面片控制顶点
};

/*
 * 编辑模式
 */
enum EditMode
{
	TRANSLATE,
	ROTATE,
	SCALE
};

#ifdef TRUTH
/*
 * 误差显示模式
 */
enum ErrorDisplayMode
{
	MESH,
	VERTEX_ERROR,
	NORMAL_ERROR
};
#endif

class Widget;

/*
 * 显示类
 * 调用 OpenGL 绘制模型，显示结果
 */
class View : public QGLWidget
{
	Q_OBJECT
private:
	QOpenGLTimerQuery		qt_query;
	//std::ofstream fout;	// 测试用输出文件，测试完成后删除
	static const int EYE_Z_FACTOR;			// 视点离原点的距离因子
	static const int MAXCTRLPOINT = 15;		// 每个方向上最大控制顶点数
	static const double ZERO;				// 很小的正数
	static const double PI;					// 圆周率
	static const int CIRCLEPOINTS = 40;		// 圆圈的采样数
	static const int DIRECTPOINTSELECTBASE = 10000;
	static const double RIGHTFAC;
	static const double BOTTOMFAC;

	Widget *m_pParent;
	CommonData *m_pCommonData;
	matrix_stack::ModelViewMatrixStack m_modelViewMatrixStack;
	matrix_stack::ProjectionMatrixStack m_projectionMatrixStack;
	matrix_stack::ModelViewMatrixStack m_cubeMapMatrixStack;
	bool m_bLoaded;					// obj文件是否已经载入
	bool m_bIsFlat;					// 绘制模式（FLAT 或 SMOOTH）
	bool m_bIsOrtho;				// 是否使用平行投影
	bool m_bCtrlPressed;			// Ctrl是否按下
	bool m_bShiftPressed;			// Shift是否按下
	MouseButton m_eMouseButton;		// 按下的鼠标按键
	int m_nDisplayMode;				// 物体、控制顶点等是否显示
	int m_nTexCase;					// 使用纹理情况
	bool m_bUseCubeMap;				// 是否使用cubemap
	std::string cubeMapName;
	bool m_bUseEnvMap;				// 是否使用环境映射
#ifdef TRUTH
	ErrorDisplayMode m_eErrorDisplayMode;
	bool m_bDisplayTruth;
	bool error_texture_;
	GLuint colormap_name;
#endif
	GLuint *tex2DNameList;
	GLuint tex3DName;
	GLuint texCubeName;
	matrix_stack::Matrix4x4 lastMatrix;			// 用于鼠标旋转模型
	matrix_stack::Matrix4x4 lastMatrixCubeMap;	// 用于鼠标旋转模型
	double lastPos[3], curPos[3];				// 用于鼠标旋转模型
	double lastPosCubeMap[3], curPosCubeMap[3]; // 用于鼠标旋转cubemap
	float light_position[4];		// 光源的位置
	double rotateAxis[3];			// 旋转轴
	double m_fTheta;				// 旋转角度
	double rotateAxisCubeMap[3];
	double m_fThetaCubeMap;
	double m_fScale;				// 缩放比例
	double m_fDeltaX, m_fDeltaY;	// 移动距离
	int m_nCursorX, m_nCursorY;		// 鼠标指针位置
	GLfloat m_fVector001[3];		// 垂直屏幕向外的法向量
	double m_fRight, m_fBottom;		// 屏幕最右端和最下端的坐标
	EditMode m_eEditMode;			// 控制顶点的编辑模式
	int m_nLineSeqNum;		// 旋转控制顶点时鼠标选择的线段（拼接成圆）的编号
	int m_nSelectedOriginalFace;

	/* 选取控制顶点的缓存 */
	GLuint selectBuf[MAXCTRLPOINT * MAXCTRLPOINT * MAXCTRLPOINT * 6];
	GLuint selectBufAxis[1200];
	GLuint selectBufOriginalModel[1000];
	int m_nEditDirection;					// 控制顶点的编辑方向
	int axisConvert123To124[4];				// 1,2,3转化为1,2,4

	GLint m_nCtrlSelectStartX, m_nCtrlSelectStartY;		// 选取框开始点
	GLint m_nCtrlSelectEndX, m_nCtrlSelectEndY;			// 选取框结束点

	QTime m_calcTime;

	//GLuint vs, tcs, tes, fs, prog;
	GLuint prog;
	GLint MVMatrix_id_, PMatrix_id_, NormalMatrix_id_;
	GLint local_viewer_id_, light_position_id_, use_env_map_id_;
#ifdef DRAW_TRIANGULAR_CTRL_POINTS
	GLint draw_which_;
#endif
#ifdef TRUTH
	GLint display_truth_id_, error_texture_id_, error_display_mode_id_;
#endif
//#ifdef LINE
	//GLint divide_id_, use_line_id_;
//#endif
#ifdef LINE
	GLint use_line_id_;
#endif
	GLint cubemap_rot_inv_mat_id_, ka_id_, kd_id_, ks_id_, tex_case_;
	//GLint min_vertex_id_, delta_vertex_inverse_id_;
	GLuint texCoordVBO;

	bool m_bLocalViewer;
	bool m_bLine;
	//float m_fDist;

	// two query buffers: front and back
	// the array to store the two sets of queries.
	unsigned int queryID[2];
	unsigned int queryBackBuffer, queryFrontBuffer;

	GLuint64 startTime, stopTime;

	//std::vector<GLuint> globalIdxVBOList;
public:
	View(QWidget *parent, CommonData *commonData);
	~View();
	void genTextures();
	void genBuffers();
	void delTextures();
	void delBuffers();
	void loadObj(const QString &texName);
	int displayMode() const	{return m_nDisplayMode;}
	bool isFlat() const		{return m_bIsFlat;}
	void setBufferData();
	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void mouseReleaseEvent(QMouseEvent *event);
	EditMode getEditMode() { return m_eEditMode; }

private:
	void swapQueryBuffers();
	QSize sizeHint() const {return QSize(900, 900);}
	void setProjection();
	void initializeGL();
	void resizeGL(int width, int height);
	void paintGL();
	void drawCtrlPoint(GLenum mode);
#ifdef DRAW_TRIANGULAR_CTRL_POINTS
	void drawTriangularCtrlPoints();
#endif
	void drawDirectPoint(GLenum mode);
	void drawAxis(double rightFac, double bottomFac, GLenum mode, bool drawChip);
	void drawCircle(double rightFac, double bottomFac, GLenum mode);
	void drawCubeMap();
	void processHits(GLint hits);
	void processHitsAxis(GLint hits);
	void processHitsCircle(GLint hits);
	void processHitsOriginalModel(GLint hits);
	void calcSelectionRange(double &centerX, double &centerY, double &width, double &height);
	void calcSphereCoord(int x, int y, double *lastPos);
	void calcThetaAndRotateAxis(double *cPos, double *lPos, double &theta, double *rotAxis);
	void wheelEvent(QWheelEvent *event);
	void keyPressEvent(QKeyEvent *event);
	void keyReleaseEvent(QKeyEvent *event);

	GLint getUniLoc(GLuint program, const GLchar *name);
	void printShaderInfoLog(GLuint shader, std::string shade_name);
	void printProgramInfoLog(GLuint program);
#ifdef LINE
	int readShaderSource(const char *fileName, std::string &, std::string &, std::string &, std::string &, std::string &);
	bool installShaders(const GLchar *, const GLchar *, const GLchar *, const GLchar *, const GLchar *);
#else
	int readShaderSource(const char *fileName, std::string &, std::string &, std::string &, std::string &);
	bool installShaders(const GLchar *, const GLchar *, const GLchar *, const GLchar *);
#endif

	void load2DTexture();
	void loadColormap();
	void load3DTexture(const std::string &textureName);
	void loadCubeMap(const QString &textureName);
public slots:
	void setModeCtrlPointTrans(bool checked);
	void setModeCtrlPointRotate(bool checked);
	void setModeCtrlPointScale(bool checked);
	void setDisplayModel();
	void setDisplayCtrlPoint();
	void setDisplayKnot();
	void setDisplayDirectPoint();
	void setDisplayEditIcon();
	void setDisplayTriangularCP();
	void setRenderFlat(bool checked);
	void setRenderSmooth(bool checked);
	void changeTex(const QString &texName);
	void changeCubeMap(const QString &cubeMapName);
	void setEnvMap(bool checked);
	void setNoEnvMap(bool checked);
	void setOrtho(bool is_ortho);
};

#endif
