/*
 * 程序的主框架、主界面类
 */

#ifndef __WIDGET_H__
#define __WIDGET_H__

#include <QWidget>
#include <fstream>
#include "view.h"
#include "b_spline_volume_set_dialog.h"
#include "common_data.h"

class Widget : public QWidget
{
	Q_OBJECT
public:
	Widget();
	~Widget();
	void closeEvent(QCloseEvent *event);
	void setUsePN();
	void changeAlgorithm();
#ifdef TRUTH
	void changeDisplayTruth(bool truth);
	void changeErrorTexture(bool use_texture);
	void changeDisplayMode(ErrorDisplayMode mode);
#endif
private:
	View *m_pView;
	CommonData *m_pCommonData;
	std::ofstream fout;

	BSplineVolumeSetDialog *bSplineVolumeSetDialog;
	QPushButton *loadObjButton;
	QPushButton *bSplineVolumeSetButton;

	QPushButton *loadEditButton;
	QPushButton *saveEditButton;

	QPushButton *loadTargetButton;
	QPushButton *saveTargetButton;

	QPushButton *morphButton;
	QPushButton *stopMorphButton;

	QTimer *morphTimer;

	QRadioButton *FFDRButton;
	QRadioButton *AFFDRButton;

	QRadioButton *renderFlatRButton;
	QRadioButton *renderSmoothRButton;

	QRadioButton *calcCpuRButton;
	QRadioButton *calcGpuRButton;

	QCheckBox *displayModelCheckBox;
	QCheckBox *displayCtrlPointCheckBox;
	QCheckBox *displayKnotCheckBox;
	QCheckBox *displayDirectPointCheckBox;
	QCheckBox *displayEditIconCheckBox;
#ifdef DRAW_TRIANGULAR_CTRL_POINTS
	QCheckBox *display_triangular_cp_checkBox;
#endif

	QSlider *sampleSlider;
	QSpinBox *sampleSpinBox;

	QRadioButton *ctrlPointTransRButton;
	QRadioButton *ctrlPointRotateRButton;
	QRadioButton *ctrlPointScaleRButton;

	QComboBox *texComboBox;
	QComboBox *cubeMapComboBox;

	QRadioButton *envMapRButton;
	QRadioButton *noEnvMapRButton;

	QRadioButton *adjust_split_points_rbutton;
	QRadioButton *no_adjust_split_points_rbutton;

	QRadioButton *adjust_silhouette_rbutton;
	QRadioButton *no_adjust_silhouette_rbutton;

	QLabel *algorithm_label;
#ifdef TRUTH
	QLabel *truth_label;
	QLabel *error_texture_label;
	QLabel *error_display_mode_label;
#endif

	std::string currentFileName;
	std::string parseFileName(std::istringstream &);
	void parseAlgorithm(std::istringstream &);
	void parseRender(std::istringstream &);
	void parseDevice(std::istringstream &);
	void parseDisplay(std::istringstream &);
	void parseSamplePointCount(std::istringstream &);
	void parseCtrlPointEditMode(std::istringstream &);
	void parseTex(std::istringstream &);
	void parseCubeMap(std::istringstream &);
	void parseEnvMap(std::istringstream &);
	void parseBSplineOrder(std::istringstream &);
	void parseBSplineCtrlPointCount(std::istringstream &);
	void parseCtrlPoint(std::istringstream &, std::ifstream &, std::vector<CtrlPoint> &);
	void parseCubemap(std::ifstream &, std::vector<matrix_stack::Matrix4x4> &);
	void parsePath(std::istringstream &, std::vector<std::vector<CtrlPoint> > &,
				   std::vector<int> &, int);
	void parseDirectPoint(std::istringstream &, std::ifstream &, std::vector<DirectPoint> &);
	double polynomial(double *coord, int degree, double x);
	void keyPressEvent(QKeyEvent *event);
private slots:
	void morph();
	void startTimer();
	void stopTimer();
	void loadObj();
	void loadEdit();
	void loadTarget();
	void saveTarget();
	void saveEdit();
	void setBSplineVolume();
	void setAlgorithmFFD(bool checked);
	void setAlgorithmAFFD(bool checked);
	void setCalcDeviceCPU(bool checked);
	void setCalcDeviceGPU(bool checked);
	void setSamplePointCount(int count);
	void setAdjustSplitPoints(bool cheched);
	void setNoAdjustSplitPoints(bool cheched);
	void setAdjustSilhouette(bool cheched);
	void setNoAdjustSilhouette(bool cheched);
};

#endif
