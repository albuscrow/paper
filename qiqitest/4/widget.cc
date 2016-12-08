#include "widget.h"
#include "matrix_stack.h"
#include <string>
#include <fstream>
#include <sstream>
#include <QString>

using namespace std;

void regGLBuffer();

Widget::Widget()
{
	fout.open("target");

	m_pCommonData = new CommonData();
	m_pView = new View(this, m_pCommonData);
	setWindowTitle(tr("AFFD"));

	/* 不响应键盘事件，以免影响View类相应Ctrl键按下的事件 */
	setFocusPolicy(Qt::StrongFocus);				// 对键盘事件的响应方式
	//setFocusPolicy(Qt::NoFocus);

	/* 控制面板组合框 */
	QGroupBox *ctrlGroupBox = new QGroupBox(this);
	ctrlGroupBox->setTitle(tr("控制面板"));

	/* 载入OBJ按钮 */
	loadObjButton = new QPushButton(ctrlGroupBox);
	loadObjButton->setText(tr("载入obj文件"));
	loadObjButton->setFocusPolicy(Qt::NoFocus);

	/* 设置B样条体属性按钮 */
	bSplineVolumeSetButton = new QPushButton(ctrlGroupBox);
	bSplineVolumeSetButton->setText(tr("设置B样条体"));
	bSplineVolumeSetButton->setEnabled(false);

	loadEditButton = new QPushButton(ctrlGroupBox);
	loadEditButton->setText(tr("载入edit文件"));

	saveEditButton = new QPushButton(ctrlGroupBox);
	saveEditButton->setText(tr("保存edit文件"));
	saveEditButton->setEnabled(false);

	loadTargetButton = new QPushButton(ctrlGroupBox);
	loadTargetButton->setText(tr("载入动画序列"));
	loadTargetButton->setEnabled(false);

	saveTargetButton = new QPushButton(ctrlGroupBox);
	saveTargetButton->setText(tr("保存控制顶点"));
	saveTargetButton->setEnabled(false);

	/* 变形按钮 */
	morphButton = new QPushButton(ctrlGroupBox);
	morphButton->setText(tr("开始变形"));
	morphButton->setEnabled(false);

	stopMorphButton = new QPushButton(ctrlGroupBox);
	stopMorphButton->setText(tr("停止变形"));
	stopMorphButton->setEnabled(false);

	morphTimer = new QTimer(this);
	connect(morphTimer, SIGNAL(timeout()), this, SLOT(morph()));

	QGridLayout *buttonLayout = new QGridLayout;
	buttonLayout->addWidget(loadObjButton, 0, 0, 1, 1);
	buttonLayout->addWidget(bSplineVolumeSetButton, 0, 1, 1, 1);
	buttonLayout->addWidget(loadEditButton, 1, 0, 1, 1);
	buttonLayout->addWidget(saveEditButton, 1, 1, 1, 1);
	buttonLayout->addWidget(loadTargetButton, 2, 0, 1, 1);
	buttonLayout->addWidget(saveTargetButton, 2, 1, 1, 1);
	buttonLayout->addWidget(morphButton, 3, 0, 1, 1);
	buttonLayout->addWidget(stopMorphButton, 3, 1, 1, 1);

	/* 算法组合框 */
	QGroupBox *algorithmGroupBox = new QGroupBox(ctrlGroupBox);
	algorithmGroupBox->setTitle(tr("算法"));

	FFDRButton = new QRadioButton(algorithmGroupBox);
	FFDRButton->setText(tr("FFD"));

	AFFDRButton = new QRadioButton(algorithmGroupBox);
	AFFDRButton->setText(tr("AFFD"));

	if (m_pCommonData->algorithm() == FFD)
		FFDRButton->setChecked(true);
	else
		AFFDRButton->setChecked(true);

	QHBoxLayout *algorithmLayout = new QHBoxLayout(algorithmGroupBox);
	algorithmLayout->addWidget(FFDRButton);
	algorithmLayout->addWidget(AFFDRButton);

	/* 绘制模式组合框 */
	QGroupBox *renderGroupBox = new QGroupBox(ctrlGroupBox);
	renderGroupBox->setTitle(tr("绘制模式"));

	renderFlatRButton = new QRadioButton(renderGroupBox);
	renderFlatRButton->setText(tr("FLAT"));

	renderSmoothRButton = new QRadioButton(renderGroupBox);
	renderSmoothRButton->setText(tr("SMOOTH"));

	if (m_pView->isFlat())
		renderFlatRButton->setChecked(true);
	else
		renderSmoothRButton->setChecked(true);

	QHBoxLayout *renderLayout = new QHBoxLayout(renderGroupBox);
	renderLayout->addWidget(renderFlatRButton);
	renderLayout->addWidget(renderSmoothRButton);

	/* 计算设备组合框 */
	QGroupBox *calcDeviceGroupBox = new QGroupBox(ctrlGroupBox);
	calcDeviceGroupBox->setTitle(tr("计算设备"));

	calcCpuRButton = new QRadioButton(calcDeviceGroupBox);
	calcCpuRButton->setText(tr("CPU"));

	calcGpuRButton = new QRadioButton(calcDeviceGroupBox);
	calcGpuRButton->setText(tr("GPU"));

	if (m_pCommonData->isGPU())
		calcGpuRButton->setChecked(true);
	else
		calcCpuRButton->setChecked(true);

	QHBoxLayout *calcDeviceLayout = new QHBoxLayout(calcDeviceGroupBox);
	calcDeviceLayout->addWidget(calcCpuRButton);
	calcDeviceLayout->addWidget(calcGpuRButton);

	/* 显示模式组合框 */
	QGroupBox *displayModeGroupBox = new QGroupBox(ctrlGroupBox);
	displayModeGroupBox->setTitle(tr("显示"));

	displayModelCheckBox = new QCheckBox(displayModeGroupBox);
	displayModelCheckBox->setText(tr("模型"));
	if ((m_pView->displayMode() & DISPLAY_MODEL) != 0)
		displayModelCheckBox->setChecked(true);

	displayCtrlPointCheckBox = new QCheckBox(displayModeGroupBox);
	displayCtrlPointCheckBox->setText(tr("控制顶点"));
	if ((m_pView->displayMode() & DISPLAY_CTRL_POINTS) != 0)
		displayCtrlPointCheckBox->setChecked(true);

	displayKnotCheckBox = new QCheckBox(displayModeGroupBox);
	displayKnotCheckBox->setText(tr("节点"));
	if ((m_pView->displayMode() & DISPLAY_KNOT) != 0)
		displayKnotCheckBox->setChecked(true);

	displayDirectPointCheckBox = new QCheckBox(displayModeGroupBox);
	displayDirectPointCheckBox->setText(tr("直接编辑点"));
	if ((m_pView->displayMode() & DISPLAY_DIRECT_POINT) != 0)
		displayDirectPointCheckBox->setChecked(true);

	displayEditIconCheckBox = new QCheckBox(displayModeGroupBox);
	displayEditIconCheckBox->setText(tr("编辑图标"));
	if ((m_pView->displayMode() & DISPLAY_EDIT_ICON) != 0)
		displayEditIconCheckBox->setChecked(true);

#ifdef DRAW_TRIANGULAR_CTRL_POINTS
	display_triangular_cp_checkBox = new QCheckBox(displayModeGroupBox);
	display_triangular_cp_checkBox->setText(tr("结果控制顶点"));
	if ((m_pView->displayMode() & DISPLAY_TRIANGULAR_CTRL_POINTS) != 0)
		display_triangular_cp_checkBox->setChecked(true);
#endif

	QGridLayout *displayModeLayout = new QGridLayout(displayModeGroupBox);
	displayModeLayout->addWidget(displayModelCheckBox, 0, 0);
	displayModeLayout->addWidget(displayCtrlPointCheckBox, 0, 1);
	displayModeLayout->addWidget(displayKnotCheckBox, 1, 0);
	displayModeLayout->addWidget(displayDirectPointCheckBox, 1, 1);
	displayModeLayout->addWidget(displayEditIconCheckBox, 0, 2);
#ifdef DRAW_TRIANGULAR_CTRL_POINTS
	displayModeLayout->addWidget(display_triangular_cp_checkBox, 1, 2);
#endif

	/* 采样点密度组合框 */
	QGroupBox *sampleGroupBox = new QGroupBox(ctrlGroupBox);
	sampleGroupBox->setTitle(tr("采样点数目"));

	sampleSlider = new QSlider(Qt::Horizontal);
	sampleSlider->setRange(2, 30);
	sampleSlider->setValue(m_pCommonData->samplePointCount());

	sampleSpinBox = new QSpinBox;
	sampleSpinBox->setRange(2, 30);
	sampleSpinBox->setValue(m_pCommonData->samplePointCount());

	QHBoxLayout *sampleLayout = new QHBoxLayout(sampleGroupBox);
	sampleLayout->addWidget(sampleSlider);
	sampleLayout->addWidget(sampleSpinBox);

	/* 编辑模式组合框 */
	QGroupBox *modeGroupBox = new QGroupBox(ctrlGroupBox);
	modeGroupBox->setTitle(tr("控制顶点编辑方式"));

	ctrlPointTransRButton = new QRadioButton(modeGroupBox);
	ctrlPointTransRButton->setText(tr("平移"));

	ctrlPointRotateRButton = new QRadioButton(modeGroupBox);
	ctrlPointRotateRButton->setText(tr("旋转"));

	ctrlPointScaleRButton = new QRadioButton(modeGroupBox);
	ctrlPointScaleRButton->setText(tr("放缩"));
 
	if (m_pView->getEditMode() == TRANSLATE)
		ctrlPointTransRButton->setChecked(true);
	else if (m_pView->getEditMode() == ROTATE)
		ctrlPointRotateRButton->setChecked(true);
	else
		ctrlPointScaleRButton->setChecked(true);

	QHBoxLayout *modeLayout = new QHBoxLayout(modeGroupBox);
	modeLayout->addWidget(ctrlPointTransRButton);
	modeLayout->addWidget(ctrlPointRotateRButton);
	modeLayout->addWidget(ctrlPointScaleRButton);

	/* 纹理选择组合框 */
	QGroupBox *texGroupBox = new QGroupBox(ctrlGroupBox);
	texGroupBox->setTitle(tr("纹理"));

	texComboBox = new QComboBox();
	ifstream fin("/home/cym/program/textures/3D_textures/filelist");
	string texName;
	QStringList texNameList;
	while(fin >> texName)
		texNameList << QString(tr(texName.c_str()));
	fin.close();
	fin.clear();
	texComboBox->addItem(tr("无纹理"));
	texComboBox->insertSeparator(1);
	texComboBox->addItem(tr("自带2D纹理"));
	texComboBox->insertSeparator(3);
	texComboBox->addItems(texNameList);

	texComboBox->setCurrentIndex(2);

	QVBoxLayout *texLayout = new QVBoxLayout(texGroupBox);
	texLayout->addWidget(texComboBox);

	/* cubemap选择下拉菜单 */
	QGroupBox *cubeMapGroupBox = new QGroupBox(ctrlGroupBox);
	cubeMapGroupBox->setTitle(tr("CubeMap"));

	cubeMapComboBox = new QComboBox();
	fin.open("/home/cym/program/textures/CubeMaps/filelist");
	texNameList.clear();
	while(fin >> texName)
		texNameList << QString(tr(texName.c_str()));
	fin.close();
	fin.clear();
	cubeMapComboBox->addItem(tr("无纹理"));
	cubeMapComboBox->insertSeparator(1);
	cubeMapComboBox->addItems(texNameList);

	QVBoxLayout *cubeMapLayout = new QVBoxLayout(cubeMapGroupBox);
	cubeMapLayout->addWidget(cubeMapComboBox);

	/* 是否使用环境映照组合框 */
	QGroupBox *envMapGroupBox = new QGroupBox(ctrlGroupBox);
	envMapGroupBox->setTitle(tr("是否使用环境映照"));

	envMapRButton = new QRadioButton(envMapGroupBox);
	envMapRButton->setText(tr("是"));

	noEnvMapRButton = new QRadioButton(envMapGroupBox);
	noEnvMapRButton->setText(tr("否"));
	noEnvMapRButton->setChecked(true);

	QHBoxLayout *envMapLayout = new QHBoxLayout(envMapGroupBox);
	envMapLayout->addWidget(envMapRButton);
	envMapLayout->addWidget(noEnvMapRButton);

	/* 是否调整切割点组合框 */
	QGroupBox *adjust_split_points_groupbox = new QGroupBox(ctrlGroupBox);
	adjust_split_points_groupbox->setTitle(tr("是否调整切割点"));

	adjust_split_points_rbutton = new QRadioButton(adjust_split_points_groupbox);
	adjust_split_points_rbutton->setText(tr("是"));

	no_adjust_split_points_rbutton = new QRadioButton(adjust_split_points_groupbox);
	no_adjust_split_points_rbutton->setText(tr("否"));
	no_adjust_split_points_rbutton->setChecked(true);
	if (m_pCommonData->adjustSplitPoints())
		adjust_split_points_rbutton->setChecked(true);

	QHBoxLayout *adjust_split_points_layout = new QHBoxLayout(adjust_split_points_groupbox);
	adjust_split_points_layout->addWidget(adjust_split_points_rbutton);
	adjust_split_points_layout->addWidget(no_adjust_split_points_rbutton);

	/* 是否调整侧影轮廓线组合框 */
	QGroupBox *adjust_silhouette_groupbox = new QGroupBox(ctrlGroupBox);
	adjust_silhouette_groupbox->setTitle(tr("是否调整侧影轮廓线"));

	adjust_silhouette_rbutton = new QRadioButton(adjust_silhouette_groupbox);
	adjust_silhouette_rbutton->setText(tr("是"));

	no_adjust_silhouette_rbutton = new QRadioButton(adjust_silhouette_groupbox);
	no_adjust_silhouette_rbutton->setText(tr("否"));
	no_adjust_silhouette_rbutton->setChecked(true);
	if (m_pCommonData->adjustSilhouette())
		adjust_silhouette_rbutton->setChecked(true);

	QHBoxLayout *adjust_silhouette_layout = new QHBoxLayout(adjust_silhouette_groupbox);
	adjust_silhouette_layout->addWidget(adjust_silhouette_rbutton);
	adjust_silhouette_layout->addWidget(no_adjust_silhouette_rbutton);

	if (m_pCommonData->getAlgorithmType() == CYM)
	{
		if (m_pCommonData->usePN())
			algorithm_label = new QLabel(tr("CYM+PN"));
		else
			algorithm_label = new QLabel(tr("CYM"));
	}
	else if (m_pCommonData->getAlgorithmType() == PN_CUTTING)
		algorithm_label = new QLabel(tr("PN_CUTTING"));
	else
		algorithm_label = new QLabel(tr("PN_NO_CUTTING"));
#ifdef TRUTH
	/* 状态显示组合框 */
	QGroupBox *statusGroupBox = new QGroupBox(ctrlGroupBox);
	statusGroupBox->setTitle("状态");

	truth_label = new QLabel(tr("近似模型"));
	error_texture_label = new QLabel(tr("误差映射：纹理"));
	error_display_mode_label = new QLabel(tr("显示：模型"));

	QGridLayout *status_layout = new QGridLayout(statusGroupBox);
	status_layout->addWidget(truth_label, 0, 0, 1, 1);
	status_layout->addWidget(error_texture_label, 0, 1, 1, 1);
	status_layout->addWidget(algorithm_label, 1, 0, 1, 1);
	status_layout->addWidget(error_display_mode_label, 1, 1, 1, 1);
#endif

	/* 布局 */
	QHBoxLayout *tex2Layout = new QHBoxLayout;
	tex2Layout->addWidget(texGroupBox);
	tex2Layout->addWidget(cubeMapGroupBox);
	QVBoxLayout *ctrlLayout = new QVBoxLayout(ctrlGroupBox);
	ctrlLayout->addLayout(buttonLayout);
	ctrlLayout->addWidget(algorithmGroupBox);
	ctrlLayout->addWidget(renderGroupBox);
	ctrlLayout->addWidget(calcDeviceGroupBox);
	ctrlLayout->addWidget(displayModeGroupBox);
	ctrlLayout->addWidget(sampleGroupBox);
	ctrlLayout->addWidget(modeGroupBox);
	//ctrlLayout->addWidget(texGroupBox);
	//ctrlLayout->addWidget(cubeMapGroupBox);
	ctrlLayout->addLayout(tex2Layout);
	ctrlLayout->addWidget(envMapGroupBox);
	ctrlLayout->addWidget(adjust_split_points_groupbox);
	ctrlLayout->addWidget(adjust_silhouette_groupbox);
#ifdef TRUTH
	ctrlLayout->addWidget(statusGroupBox);
#else
	ctrlLayout->addWidget(algorithm_label);
#endif

	QVBoxLayout *rightLayout = new QVBoxLayout;
	rightLayout->addWidget(ctrlGroupBox);
	rightLayout->addStretch();

	QHBoxLayout *mainLayout = new QHBoxLayout(this);
	mainLayout->addWidget(m_pView);
	mainLayout->addLayout(rightLayout);

	//ctrlGroupBox->setFixedWidth(ctrlGroupBox->sizeHint().width());
	ctrlGroupBox->setFixedWidth(300);
	ctrlGroupBox->setFixedHeight(ctrlGroupBox->sizeHint().height());

	/* 连接 */
	connect(loadObjButton, SIGNAL(clicked()), this, SLOT(loadObj()));
	connect(bSplineVolumeSetButton, SIGNAL(clicked()), this, SLOT(setBSplineVolume()));
	connect(loadEditButton, SIGNAL(clicked()), this, SLOT(loadEdit()));
	connect(saveEditButton, SIGNAL(clicked()), this, SLOT(saveEdit()));
	connect(loadTargetButton, SIGNAL(clicked()), this, SLOT(loadTarget()));
	connect(saveTargetButton, SIGNAL(clicked()), this, SLOT(saveTarget()));
	connect(morphButton, SIGNAL(clicked()), this, SLOT(startTimer()));
	connect(stopMorphButton, SIGNAL(clicked()), this, SLOT(stopTimer()));

	connect(FFDRButton, SIGNAL(clicked(bool)), this, SLOT(setAlgorithmFFD(bool)));
	connect(AFFDRButton, SIGNAL(clicked(bool)), this, SLOT(setAlgorithmAFFD(bool)));

	connect(renderFlatRButton, SIGNAL(clicked(bool)), m_pView, SLOT(setRenderFlat(bool)));
	connect(renderSmoothRButton, SIGNAL(clicked(bool)), m_pView, SLOT(setRenderSmooth(bool)));

	connect(calcCpuRButton, SIGNAL(clicked(bool)), this, SLOT(setCalcDeviceCPU(bool)));
	connect(calcGpuRButton, SIGNAL(clicked(bool)), this, SLOT(setCalcDeviceGPU(bool)));

	connect(displayModelCheckBox, SIGNAL(clicked()), m_pView, SLOT(setDisplayModel()));
	connect(displayCtrlPointCheckBox, SIGNAL(clicked()), m_pView, SLOT(setDisplayCtrlPoint()));
	connect(displayKnotCheckBox, SIGNAL(clicked()), m_pView, SLOT(setDisplayKnot()));
	connect(displayDirectPointCheckBox, SIGNAL(clicked()), m_pView, SLOT(setDisplayDirectPoint()));
	connect(displayEditIconCheckBox, SIGNAL(clicked()), m_pView, SLOT(setDisplayEditIcon()));
#ifdef DRAW_TRIANGULAR_CTRL_POINTS
	connect(display_triangular_cp_checkBox, SIGNAL(clicked()), m_pView, SLOT(setDisplayTriangularCP()));
#endif

	connect(sampleSlider, SIGNAL(valueChanged(int)), sampleSpinBox, SLOT(setValue(int)));
	connect(sampleSlider, SIGNAL(valueChanged(int)), this, SLOT(setSamplePointCount(int)));
	connect(sampleSpinBox, SIGNAL(valueChanged(int)), sampleSlider, SLOT(setValue(int)));
	connect(sampleSpinBox, SIGNAL(valueChanged(int)), this, SLOT(setSamplePointCount(int)));

	connect(ctrlPointTransRButton, SIGNAL(clicked(bool)), m_pView, SLOT(setModeCtrlPointTrans(bool)));
	connect(ctrlPointRotateRButton, SIGNAL(clicked(bool)), m_pView, SLOT(setModeCtrlPointRotate(bool)));
	connect(ctrlPointScaleRButton, SIGNAL(clicked(bool)), m_pView, SLOT(setModeCtrlPointScale(bool)));

	connect(texComboBox, SIGNAL(currentIndexChanged(const QString &)), m_pView, SLOT(changeTex(const QString &)));
	connect(cubeMapComboBox, SIGNAL(currentIndexChanged(const QString &)), m_pView, SLOT(changeCubeMap(const QString &)));

	connect(envMapRButton, SIGNAL(clicked(bool)), m_pView, SLOT(setEnvMap(bool)));
	connect(noEnvMapRButton, SIGNAL(clicked(bool)), m_pView, SLOT(setNoEnvMap(bool)));

	connect(adjust_split_points_rbutton, SIGNAL(clicked(bool)), this, SLOT(setAdjustSplitPoints(bool)));
	connect(no_adjust_split_points_rbutton, SIGNAL(clicked(bool)), this, SLOT(setNoAdjustSplitPoints(bool)));

	connect(adjust_silhouette_rbutton, SIGNAL(clicked(bool)), this, SLOT(setAdjustSilhouette(bool)));
	connect(no_adjust_silhouette_rbutton, SIGNAL(clicked(bool)), this, SLOT(setNoAdjustSilhouette(bool)));
}

Widget::~Widget()
{
	fout.close();
	fout.clear();

	delete m_pView;

	delete loadObjButton;
	delete bSplineVolumeSetButton;

	delete loadEditButton;
	delete saveEditButton;

	delete loadTargetButton;
	delete saveTargetButton;

	delete morphButton;
	delete stopMorphButton;

	delete morphTimer;

	delete FFDRButton;
	delete AFFDRButton;

	delete renderFlatRButton;
	delete renderSmoothRButton;

	delete calcCpuRButton;
	delete calcGpuRButton;

	delete sampleSlider;
	delete sampleSpinBox;

	delete ctrlPointTransRButton;
	delete ctrlPointRotateRButton;
	delete ctrlPointScaleRButton;

	delete texComboBox;
	delete cubeMapComboBox;

	delete envMapRButton;
	delete noEnvMapRButton;

	delete adjust_split_points_rbutton;
	delete no_adjust_split_points_rbutton;

	delete adjust_silhouette_rbutton;
	delete no_adjust_silhouette_rbutton;

	delete algorithm_label;
#ifdef TRUTH
	delete truth_label;
	delete error_texture_label;
	delete error_display_mode_label;
#endif
}

void Widget::closeEvent(QCloseEvent *event)
{
	delete m_pCommonData;			// closeEvent函数在窗口关闭时会调用
									// 而且早于上面的~Widget函数。
									// 这里是调用cudaGraphicsUnregisterResource的合适地方，
									// 不会返回错误，而且系统也不会卡
	event->accept();
}

void Widget::morph()
{
	m_pCommonData->morph();
	m_pView->updateGL();
}

void Widget::startTimer()
{
	morphTimer->start(30);
	//morphTimer->start(16.6666);
	//morphTimer->start(33.333333);
}

void Widget::stopTimer()
{
	morphTimer->stop();
}

void printCudaError(const char *file, const char *function, int line);

void Widget::loadObj()
{
//#define DEFAULT_OPEN
#ifdef DEFAULT_OPEN
	//QString file_name(tr("/home/cym/program/make_obj/cylinder/cylinder.obj"));
	QString file_name(tr("/home/cym/program/OBJ/simple/01.Cube.obj"));
	//QString file_name(tr("/home/cym/program/make_obj/cube/cube.obj"));
	//QString file_name(tr("/home/cym/program/OBJ/simple/test.obj"));
	//QString file_name(tr("/home/cym/program/OBJ/simple/knotbox.obj"));
	//QString file_name(tr("/home/cym/program/OBJ/simple/Tetrahedra.obj"));
	//QString file_name(tr("/home/cym/program/OBJ/simple/test.obj"));
	//QString file_name(tr("/home/cym/program/OBJ/simple/torus.obj"));
	//QString file_name(tr("/home/cym/program/OBJ/coke2/coke.obj"));
	//QString file_name(tr("/home/cym/program/OBJ/algorithm2/snail/snail.obj"));
	//QString file_name(tr("/home/cym/program/make_obj/drum/drum_cym.obj"));
	//QString file_name(tr("/home/cym/program/make_obj/cone/cone_cym.obj"));
	//QString file_name(tr("/home/cym/program/make_obj/double_sphere/double_sphere_cym.obj"));
	//QString file_name(tr("/home/cym/program/make_obj/bishop/biship_cym_area_average_normal.obj"));

	currentFileName = string(file_name.toLocal8Bit().data());
	bSplineVolumeSetButton->setEnabled(true);
	saveEditButton->setEnabled(true);
	loadTargetButton->setEnabled(true);
	saveTargetButton->setEnabled(true);
	m_pCommonData->loadObj(currentFileName.c_str());

	printCudaError(__FILE__, __FUNCTION__, __LINE__);
	m_pCommonData->preCalc();
	printCudaError(__FILE__, __FUNCTION__, __LINE__);
	//m_pView->delGlobalIdxVBOList();
	m_pView->delTextures();
	m_pView->delBuffers();
	printCudaError(__FILE__, __FUNCTION__, __LINE__);
	//m_pView->genGlobalIdxVBOList();
	m_pView->genBuffers();
	m_pView->genTextures();
	printCudaError(__FILE__, __FUNCTION__, __LINE__);
	m_pView->setBufferData();
	printCudaError(__FILE__, __FUNCTION__, __LINE__);
	regGLBuffer();
	m_pCommonData->execute();

	m_pView->loadObj(texComboBox->currentText());
	m_pView->updateGL();
#else
	QString qCurrentFileName = QFileDialog::getOpenFileName(this, tr("打开obj文件"),
				//QString(tr("/home/cym/program/OBJ/simple/01.Cube.obj")),
				//QString(tr("/home/cym/program/make_obj/cube/cube.obj")),
				//QString(tr("/home/cym/program/OBJ/simple/03.Cube2.obj")),
				//QString(tr("/home/cym/program/OBJ/simple/torus.obj")),
				//QString(tr("/home/cym/program/OBJ/simple")),
				//QString(tr("/home/cym/program/OBJ/simple/knotbox.obj")),
				//QString(tr("/home/cym/program/OBJ/simple/test.obj")),
				//QString(tr("/home/cym/program/OBJ/simple/Tetrahedra.obj")),
				//QString(tr("/home/cym/program/OBJ/simple/drum_cym.obj")),
				//QString(tr("/home/cym/program/make_obj/cylinder/cylinder.obj")),
				//QString(tr("/home/cym/program/make_obj/drum/drum_cym.obj")),
				//QString(tr("/home/cym/program/make_obj/double_sphere/double_sphere_cym.obj")),



				//QString(tr("/home/cym/program/make_obj/bishop/biship_cym_area_average_normal.obj")),
				//QString(tr("/home/cym/program/OBJ/兔/rabbit_cym.obj")),
				//QString(tr("/home/cym/program/make_obj/greek_vase/greek_vase_cym_area_average_normal.obj")),
				QString(tr("/home/cym/program/make_obj/greek_vase_no_tex/greek_vase_cym_area_average_normal.obj")),
				//QString(tr("/home/cym/program/OBJ/rabbit_real/rabbit.obj")),
				//QString(tr("/home/cym/program/OBJ/ship/ship.obj")),



				//QString(tr("/home/cym/E/model/OBJ/兔子.obj")),
				//QString(tr("/home/cym/E/3/")),
				//QString(tr("/home/cym/program/OBJ/simple/03.Cube4.obj")),
				//QString(tr("/home/cym/program/OBJ/simple/03.Cube8.obj")),
				//QString(tr("/home/cym/program/OBJ/test/right.obj")),
				//QString(tr("/home/cym/program/OBJ/")),
				//QString(tr("/home/cym/program/OBJ/box2/49.obj")),
				//QString(tr("/home/cym/program/OBJ/coke/coke_can.obj")),
				//QString(tr("/home/cym/program/OBJ/algorithm2/")),
				//QString(tr("/home/cym/program/OBJ/algorithm2/amphora/amphora.obj")),
				//QString(tr("/home/cym/program/OBJ/algorithm2/spaceship/spaceship.obj")),
				//QString(tr("/home/cym/program/OBJ/algorithm2/snail/snail.obj")),

				tr("obj文件 (*.obj)"));
				//tr("obj文件 (*.obj)(*.obj)"));
				//tr("obj文件( *.obj *.mtl)"));
	if (!qCurrentFileName.isEmpty())
	{
		//说明QByteArray ba = fileName.toLocal8Bit();
		//说明const char *c_str = ba.data();
		currentFileName = string(qCurrentFileName.toLocal8Bit().data());
		bSplineVolumeSetButton->setEnabled(true);
		saveEditButton->setEnabled(true);
		loadTargetButton->setEnabled(true);
		saveTargetButton->setEnabled(true);
		m_pCommonData->loadObj(currentFileName.c_str());

		printCudaError(__FILE__, __FUNCTION__, __LINE__);
		m_pCommonData->preCalc();
		printCudaError(__FILE__, __FUNCTION__, __LINE__);
		//m_pView->delGlobalIdxVBOList();
		m_pView->delTextures();
		m_pView->delBuffers();
		printCudaError(__FILE__, __FUNCTION__, __LINE__);
		//m_pView->genGlobalIdxVBOList();
		m_pView->genBuffers();
		m_pView->genTextures();
		printCudaError(__FILE__, __FUNCTION__, __LINE__);
		m_pView->setBufferData();
		printCudaError(__FILE__, __FUNCTION__, __LINE__);
		regGLBuffer();
		m_pCommonData->execute();

		m_pView->loadObj(texComboBox->currentText());
		m_pView->updateGL();
	}
#endif
}

/* 载入 edit 文件 */
void Widget::loadEdit()
{
	QString fileName = QFileDialog::getOpenFileName(this, tr("打开edit文件"),
		//作废QString(QDir::currentPath() + tr("/OBJ/")),
		//QString(tr("/home/cym/program/OBJ/chairIV.edit")),
		//QString(tr("/home/cym/program/OBJ/algorithm2/snail/0.edit")),



		//QString(tr("/home/cym/program/make_obj/bishop/step0.edit")),
		//QString(tr("/home/cym/program/OBJ/兔/1.edit")),
		//QString(tr("/home/cym/program/make_obj/greek_vase_no_tex/1.edit")),
		//QString(tr("/home/cym/program/OBJ/rabbit_real/2.edit")),
		QString(tr("/home/cym/program/OBJ/ship/b.edit")),
		//QString(tr("/home/cym/program/make_obj/greek_vase/begin.edit")),



													tr("edit文件 (*.edit)"));
	if (!fileName.isEmpty())
	{
		m_pCommonData->setLoadingEdit(true);
		QByteArray ba = fileName.toLocal8Bit();
		const char *c_str = ba.data();
		ifstream fin(c_str);

		vector<CtrlPoint> ctrlPointVector;
		vector<DirectPoint> directPointVector;
		string line;
		//cout << "begin!" << endl;
		while(getline(fin, line))
		{
			istringstream iss(line);
			string type;
			iss >> type;
			if (type == "文件名")
				currentFileName = parseFileName(iss);
			else if (type == "算法")
				parseAlgorithm(iss);
			//else if (type == "绘制模式")
				//parseRender(iss);
			//else if (type == "计算设备")
				//parseDevice(iss);
			else if (type == "显示")
				parseDisplay(iss);
			else if (type == "采样点数目")
				parseSamplePointCount(iss);
			else if (type == "控制顶点编辑方式")
				parseCtrlPointEditMode(iss);
			else if (type == "纹理")
				parseTex(iss);
			else if (type == "CubeMap")
				parseCubeMap(iss);
			else if (type == "是否使用环境映照")
				parseEnvMap(iss);
			else if (type == "B样条体阶数")
				parseBSplineOrder(iss);
			else if (type == "B样条体顶点数")
				parseBSplineCtrlPointCount(iss);
			else if (type == "控制顶点")
				parseCtrlPoint(iss, fin, ctrlPointVector);
			else if (type == "直接编辑点")
				parseDirectPoint(iss, fin, directPointVector);
		}
		fin.close();
		fin.clear();
		//cout << "finish!" << endl;
		bSplineVolumeSetButton->setEnabled(true);
		saveEditButton->setEnabled(true);
		loadTargetButton->setEnabled(true);
		saveTargetButton->setEnabled(true);
		m_pCommonData->loadObj(currentFileName.c_str());

		m_pCommonData->preCalc();
		//m_pView->delGlobalIdxVBOList();
		m_pView->delTextures();
		m_pView->delBuffers();
		//m_pView->genGlobalIdxVBOList();
		m_pView->genBuffers();
		m_pView->genTextures();
		m_pView->setBufferData();
		regGLBuffer();
		m_pCommonData->execute();

		m_pView->loadObj(texComboBox->currentText());
		m_pCommonData->setAllCtrlPoint(ctrlPointVector);
		//m_pCommonData->setAllDirectPoint(directPointVector);
		m_pCommonData->setLoadingEdit(false);
		m_pView->updateGL();
	}
}

void Widget::loadTarget()
{
	QString fileName = QFileDialog::getOpenFileName(this, tr("打开target文件"),
													//作废QString(QDir::currentPath() + tr("/OBJ/box1")),
													//QString(tr("/home/cym/program/OBJ/algorithm2/bluebird/output.target")),
													//QString(tr("/home/cym/program/OBJ/algorithm2/amphora/output.target")),
													//QString(tr("/home/cym/program/OBJ/algorithm2/snail/output.target")),
													//QString(tr("/home/cym/program/OBJ/ship/make_animation/output.target")),
													QString(tr("/home/cym/program/make_obj/greek_vase/make_animation/output.target")),
													//QString(tr("/home/cym/program/make_obj/bishop/make_animation/output.target")),
													tr("target文件 (*.target)"));
	if (!fileName.isEmpty())
	{
		QByteArray ba = fileName.toLocal8Bit();
		const char *c_str = ba.data();
		ifstream fin(c_str);

		vector<vector<CtrlPoint> > animationList;
		vector<int> animationTimeList;
		vector<matrix_stack::Matrix4x4> cubemap_matrix_list;
		string line;
		while(getline(fin, line))
		{
			istringstream iss(line);
			string type;
			int time;
			iss >> type >> time;
			if (type == "控制顶点")
			{
				//cout << "ctrl, time = " << time << endl;
				vector<CtrlPoint> emptyVector;
				parseCtrlPoint(iss, fin, emptyVector);
				animationList.push_back(emptyVector);
				animationTimeList.push_back(time);
			}
			else if (type == "cubemap")
			{
				parseCubemap(fin, cubemap_matrix_list);
			}
			else if (type == "路径")
			{
				//cout << "path, time = " << time << endl;
				parsePath(iss, animationList, animationTimeList, time);
			}
		}
		fin.close();
		fin.clear();

		m_pCommonData->setTargetCtrlPoint(animationList, animationTimeList, cubemap_matrix_list);
		morphButton->setEnabled(true);
		stopMorphButton->setEnabled(true);
	}
}

void Widget::saveTarget()
{
	cout << "saveTarget" << endl;
	int countU = m_pCommonData->ctrlPointCount(U);
	int countV = m_pCommonData->ctrlPointCount(V);
	int countW = m_pCommonData->ctrlPointCount(W);
	fout << countU * countV * countW << endl;
	for (int i = 0; i < countU; ++i)
		for (int j = 0; j < countV; ++j)
			for (int k = 0; k < countW; ++k)
			{
				CtrlPoint cp = m_pCommonData->getCtrlPoint(i, j, k);
				fout << cp.x() << ' ' << cp.y() << ' ' << cp.z() << endl;
			}
}

string Widget::parseFileName(istringstream &iss)
{
	string fileName;
	iss >> fileName;
	return fileName;
}

void Widget::parseAlgorithm(istringstream &iss)
{
	string algorithm;
	iss >> algorithm;
	if (algorithm == "FFD")
	{
		FFDRButton->click();
	}
	else
	{
		AFFDRButton->click();
	}
}

void Widget::parseRender(istringstream &iss)
{
	string render;
	iss >> render;
	if (render == "FLAT")
	{
		renderFlatRButton->click();
	}
	else
	{
		renderSmoothRButton->click();
	}
}

void Widget::parseDevice(istringstream &iss)
{
	string device;
	iss >> device;
	if (device == "CPU")
	{
		calcCpuRButton->click();
	}
	else
	{
		calcGpuRButton->click();
	}
}

void Widget::parseDisplay(istringstream &iss)
{
	string display;
	/* 将显示模式和复选框清零 */
	int displayMode = m_pView->displayMode();
	if (displayMode & DISPLAY_MODEL)
		m_pView->setDisplayModel();
	if (displayMode & DISPLAY_CTRL_POINTS)
		m_pView->setDisplayCtrlPoint();
	if (displayMode & DISPLAY_KNOT)
		m_pView->setDisplayKnot();
	if (displayMode & DISPLAY_DIRECT_POINT)
		m_pView->setDisplayDirectPoint();
	if (displayMode & DISPLAY_EDIT_ICON)
		m_pView->setDisplayEditIcon();
	if (displayMode & DISPLAY_TRIANGULAR_CTRL_POINTS)
		m_pView->setDisplayTriangularCP();
	displayModelCheckBox->setChecked(false);
	displayCtrlPointCheckBox->setChecked(false);
	displayKnotCheckBox->setChecked(false);
	displayDirectPointCheckBox->setChecked(false);
	displayEditIconCheckBox->setChecked(false);
#ifdef DRAW_TRIANGULAR_CTRL_POINTS
	display_triangular_cp_checkBox->setChecked(false);
#endif

	while(iss >> display)
	{
		if (display == "Model")
		{
			displayModelCheckBox->click();
		}
		else if (display == "CtrlPoint")
		{
			displayCtrlPointCheckBox->click();
		}
		else if (display == "Knot")
		{
			displayKnotCheckBox->click();
		}
		else if (display == "DirectPoint")
		{
			displayDirectPointCheckBox->click();
		}
		else if (display == "EditIcon")
		{
			displayEditIconCheckBox->click();
		}
#ifdef DRAW_TRIANGULAR_CTRL_POINTS
		else if (display == "TriangularCP")
		{
			display_triangular_cp_checkBox->click();
		}
#endif
	}
}

void Widget::parseSamplePointCount(istringstream &iss)
{
	int samplePointCount;
	iss >> samplePointCount;
	sampleSlider->setValue(samplePointCount);
	sampleSpinBox->setValue(samplePointCount);
}

void Widget::parseCtrlPointEditMode(istringstream &iss)
{
	string ctrlPointEditMode;
	iss >> ctrlPointEditMode;
	if (ctrlPointEditMode == "Trans")
	{
		ctrlPointTransRButton->click();
	}
	else if (ctrlPointEditMode == "Rotate")
	{
		ctrlPointRotateRButton->click();
	}
	else if (ctrlPointEditMode == "Scale")
	{
		ctrlPointScaleRButton->click();
	}
}

void Widget::parseTex(istringstream &iss)
{
	int texIdx;
	iss >> texIdx;
	texComboBox->setCurrentIndex(texIdx);
}

void Widget::parseCubeMap(istringstream &iss)
{
	int cubeMapIdx;
	iss >> cubeMapIdx;
	cubeMapComboBox->setCurrentIndex(cubeMapIdx);
}

void Widget::parseEnvMap(istringstream &iss)
{
	string useEnvMap;
	iss >> useEnvMap;
	if (useEnvMap == "Yes")
	{
		envMapRButton->click();
	}
	else
	{
		noEnvMapRButton->click();
	}
}

void Widget::parseBSplineOrder(istringstream &iss)
{
	int orderU, orderV, orderW;
	iss >> orderU >> orderV >> orderW;
	m_pCommonData->setOrder(orderU, U);
	m_pCommonData->setOrder(orderV, V);
	m_pCommonData->setOrder(orderW, W);
}

void Widget::parseBSplineCtrlPointCount(istringstream &iss)
{
	int ctrlPointCountU, ctrlPointCountV, ctrlPointCountW;
	iss >> ctrlPointCountU >> ctrlPointCountV >> ctrlPointCountW;
	m_pCommonData->setCtrlPointCount(ctrlPointCountU, U);
	m_pCommonData->setCtrlPointCount(ctrlPointCountV, V);
	m_pCommonData->setCtrlPointCount(ctrlPointCountW, W);
}

void Widget::parseCtrlPoint(istringstream &iss, ifstream &fin, vector<CtrlPoint> &ctrlPointVector)
{
	int ctrlPointCount;
	iss >> ctrlPointCount;
	for (int i = 0; i < ctrlPointCount; ++i)
	{
		string line;
		getline(fin, line);
		double x, y, z;
		istringstream ctrlPointStream(line);
		ctrlPointStream >> x >> y >> z;
		ctrlPointVector.push_back(CtrlPoint(x, y, z));
	}
}

void Widget::parseCubemap(ifstream &fin, vector<matrix_stack::Matrix4x4> &cubemap_matrix_list)
{
	string line;
	getline(fin, line);
	istringstream line_stream(line);
	float matrix[16];
	for (int i = 0; i < 16; ++i)
	{
		line_stream >> matrix[i];
	}
	cubemap_matrix_list.push_back(matrix_stack::Matrix4x4(matrix));
}

double Widget::polynomial(double *coord, int degree, double x)
{
	double y = coord[degree];
	for (int i = degree - 1; i >= 0; --i)
		y = y * x + coord[i];
	return y;
}

void Widget::parsePath(istringstream &iss, vector<vector<CtrlPoint> > &animationList,
					   vector<int> &animationTimeList, int time)
{
	int degree;
	iss >> degree;
	cout << "degree = " << degree << endl;
	double coord[100];
	for (int i = 0; i <= degree; ++i)
	{
		iss >> coord[i];
	}
	double xBegin, xEnd;
	int sampleCount;
	iss >> xBegin >> xEnd >> sampleCount;
	cout << "(" << xBegin << ", " << xEnd << "), sampleCount = " << sampleCount << endl;
	int ctrlPointCountU = m_pCommonData->ctrlPointCount(U);
	int ctrlPointCountV = m_pCommonData->ctrlPointCount(V);
	int ctrlPointCountW = m_pCommonData->ctrlPointCount(W);
	int ctrlPointCount = ctrlPointCountU * ctrlPointCountV * ctrlPointCountW;
	cout << "ctrlPointCount = " << ctrlPointCount << endl;
	int lastIdx = animationList.size() - 1;
	cout << "lastIdx = " << lastIdx << endl;
	for (int i = 0; i < sampleCount; ++i)
	{
		//cout << "i = " << i << endl;
		double x = xBegin * (1 - (double)i / (sampleCount - 1)) + xEnd * (double)i / (sampleCount - 1);
		double y = polynomial(coord, degree, x);
		//cout << "x = " << x << ", y = " << y << endl;
		if (lastIdx == -1)
		{
			cout << "target文件中，定义路径之前必须定义至少一个控制顶点序列！" << endl;
			return;
		}
		vector<CtrlPoint> tempList;
		for (int j = 0; j < ctrlPointCount; ++j)
		{
			CtrlPoint cp = animationList[lastIdx][j];
			cp += CtrlPoint(x, y, 0.0);
			tempList.push_back(cp);
		}
		animationList.push_back(tempList);
		animationTimeList.push_back(time);
	}
}

#ifdef TRUTH
void Widget::changeDisplayTruth(bool truth)
{
	if (truth)
		truth_label->setText(tr("精确模型\t"));
	else
		truth_label->setText(tr("近似模型\t"));
}

void Widget::changeErrorTexture(bool error_texture)
{
	if (error_texture)
		error_texture_label->setText(tr("误差映射：纹理"));
	else
		error_texture_label->setText(tr("误差映射：计算"));
}

void Widget::changeDisplayMode(ErrorDisplayMode mode)
{
	if (mode == MESH)
		error_display_mode_label->setText(tr("显示：模型"));
	else if (mode == VERTEX_ERROR)
		error_display_mode_label->setText(tr("显示：顶点误差"));
	else
		error_display_mode_label->setText(tr("显示：法向误差"));
}
#endif

/* 按下键的消息响应函数 */
void Widget::keyPressEvent(QKeyEvent *event)
{
	if (event->key() == Qt::Key_Q)
		close();
}

void Widget::parseDirectPoint(istringstream &iss, ifstream &fin, vector<DirectPoint> &directPointVector)
{
	int directPointCount;
	iss >> directPointCount;
	for (int i = 0; i < directPointCount; ++i)
	{
		string line;
		getline(fin, line);
		double u, v, w, x, y, z;
		istringstream directPointStream(line);
		directPointStream >> u >> v >> w >> x >> y >> z;
		directPointVector.push_back(DirectPoint(u, v, w, x, y, z));
	}
}

void Widget::saveEdit()
{
	QString fileName = QFileDialog::getSaveFileName(this, tr("保存edit文件"),
													//作废QString(QDir::currentPath() + tr("/OBJ/")),
													QString(tr("~/program/OBJ/")),
													tr("edit文件 (*.edit)"));
	if (!fileName.isEmpty())
	{
		if (!fileName.endsWith(".edit"))
			fileName += ".edit";
		bSplineVolumeSetButton->setEnabled(true);
		QByteArray ba = fileName.toLocal8Bit();
		const char *c_str = ba.data();
		ofstream fout(c_str);

		fout << "文件名\t\t\t\t" << currentFileName << endl;

		fout << "算法\t\t\t\t";
		if (FFDRButton->isChecked())
			fout << "FFD" << endl;
		else
			fout << "AFFD" << endl;

		fout << "绘制模式\t\t\t";
		if (renderFlatRButton->isChecked())
			fout << "FLAT" << endl;
		else
			fout << "SMOOTH" << endl;

		fout << "计算设备\t\t\t";
		if (calcCpuRButton->isChecked())
			fout << "CPU" << endl;
		else
			fout << "GPU" << endl;

		fout << "显示\t\t\t\t";
		if (displayModelCheckBox->isChecked())
			fout << "Model ";
		if (displayCtrlPointCheckBox->isChecked())
			fout << "CtrlPoint ";
		if (displayKnotCheckBox->isChecked())
			fout << "Knot ";
		if (displayDirectPointCheckBox->isChecked())
			fout << "DirectPoint ";
		if (displayEditIconCheckBox->isChecked())
			fout << "EditIcon";
#ifdef DRAW_TRIANGULAR_CTRL_POINTS
		if (display_triangular_cp_checkBox->isChecked())
			fout << "TriangularCP";
#endif
		fout << endl;

		fout << "采样点数目\t\t\t" << sampleSlider->value() << endl;

		fout << "控制顶点编辑方式\t";
		if (ctrlPointTransRButton->isChecked())
			fout << "Trans" << endl;
		else if (ctrlPointRotateRButton->isChecked())
			fout << "Rotate" << endl;
		else
			fout << "Scale" << endl;

		fout << "纹理\t\t\t\t";
		fout << texComboBox->currentIndex() << endl;

		fout << "CubeMap\t\t\t\t";
		fout << cubeMapComboBox->currentIndex() << endl;

		fout << "是否使用环境映照\t";
		if (envMapRButton->isChecked())
			fout << "Yes" << endl;
		else
			fout << "No" << endl;

		m_pCommonData->saveEdit(fout);

		fout.close();
		fout.clear();
	}
}

void Widget::setBSplineVolume()					// 按下“设置B样条体”按钮后的响应函数
{
	bSplineVolumeSetDialog = new BSplineVolumeSetDialog(m_pCommonData, this);
	if (bSplineVolumeSetDialog->exec())
	{
		m_pCommonData->resetFirstLoad();
		m_pCommonData->preCalc();
		m_pView->delBuffers();
		m_pView->genBuffers();
		m_pView->setBufferData();
		regGLBuffer();
		m_pCommonData->execute();
		m_pView->updateGL();
	}
	delete bSplineVolumeSetDialog;
}

void Widget::setUsePN()
{
	m_pCommonData->setUsePN();
	if (m_pCommonData->getAlgorithmType() == CYM)
	{
		if (m_pCommonData->usePN())
			algorithm_label->setText(tr("CYM+PN"));
		else
			algorithm_label->setText(tr("CYM"));
	}
	m_pView->updateGL();
}

void Widget::changeAlgorithm()
{
	m_pCommonData->setAlgorithmType();
	if (m_pCommonData->getAlgorithmType() == CYM)
	{
		if (m_pCommonData->usePN())
			algorithm_label->setText(tr("CYM+PN"));
		else
			algorithm_label->setText(tr("CYM"));
	}
	else if (m_pCommonData->getAlgorithmType() == PN_CUTTING)
		algorithm_label->setText(tr("PN_CUTTING"));
	else
		algorithm_label->setText(tr("PN_NO_CUTTING"));
	cout << windowTitle().toUtf8().constData() << endl;
	m_pCommonData->resetFirstLoad();
	m_pCommonData->preCalc(false);
	m_pView->delBuffers();
	m_pView->genBuffers();
	m_pView->setBufferData();
	regGLBuffer();
	m_pCommonData->execute();
	m_pView->updateGL();
}

void Widget::setAlgorithmFFD(bool checked)		// 将当前编辑模式设置为Obj编辑模式（模型的平移、旋转、缩放和控制顶点选取操作）
{
	if (checked)
		m_pCommonData->setAlgorithm(FFD);
	if (!m_pCommonData->loadingEdit())
		m_pView->updateGL();
}

void Widget::setAlgorithmAFFD(bool checked)		// 将当前编辑模式设置为控制顶点编辑模式（控制顶点的平移、旋转、缩放操作）
{
	if (checked)
		m_pCommonData->setAlgorithm(AFFD);
	if (!m_pCommonData->loadingEdit())
		m_pView->updateGL();
}

void Widget::setCalcDeviceCPU(bool checked)
{
	if (checked)
	{
		m_pCommonData->setCalcDevice(false);
		m_pCommonData->calcFinalResult();
		if (!m_pCommonData->loadingEdit())
			m_pView->updateGL();
	}
}

void Widget::setCalcDeviceGPU(bool checked)
{
	if (checked)
	{
		m_pCommonData->setCalcDevice(true);
		m_pCommonData->calcFinalResult();
		if (!m_pCommonData->loadingEdit())
			m_pView->updateGL();
	}
}

void Widget::setSamplePointCount(int count)
{
	cout << "sampleCount = " << count << endl;
	m_pCommonData->setSamplePointCount(count);
	m_pCommonData->newSamplePointTesslate();

	//if (m_pCommonData->loaded())
	if (bSplineVolumeSetButton->isEnabled())
	{
		m_pView->delBuffers();
		m_pView->genBuffers();
		m_pView->setBufferData();
		regGLBuffer();
		m_pCommonData->callTesslateD();
	}
	if (!m_pCommonData->loadingEdit())
		m_pView->updateGL();
}

/********************/

void Widget::setAdjustSplitPoints(bool checked)
{
	if (checked)
		m_pCommonData->setAdjustSplitPoints(true);
	m_pCommonData->preCalc(false);
	m_pCommonData->execute();
	m_pView->updateGL();
}

void Widget::setNoAdjustSplitPoints(bool checked)
{
	if (checked)
		m_pCommonData->setAdjustSplitPoints(false);
	m_pCommonData->preCalc(false);
	m_pCommonData->execute();
	m_pView->updateGL();
}

/********************/

void Widget::setAdjustSilhouette(bool checked)
{
	if (checked)
		m_pCommonData->setAdjustSilhouette(true);
	m_pCommonData->execute();
	m_pView->updateGL();
}

void Widget::setNoAdjustSilhouette(bool checked)
{
	if (checked)
		m_pCommonData->setAdjustSilhouette(false);
	//if (!m_pCommonData->loadingEdit())

	//m_pCommonData->preCalc();
	//m_pView->delTextures();
	//m_pView->delBuffers();
	//m_pView->genBuffers();
	//m_pView->genTextures();
	//m_pView->setBufferData();
	//regGLBuffer();
	m_pCommonData->execute();

	//m_pView->loadObj(texComboBox->currentText());
	m_pView->updateGL();
}
