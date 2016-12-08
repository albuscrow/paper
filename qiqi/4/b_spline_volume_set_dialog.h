/*
 * 定义了一个设置B样条体参数的对话框
 */

#ifndef __BSPLINEVOLUMESETDIALOG_H__
#define __BSPLINEVOLUMESETDIALOG_H__

#include <QSpinBox>
#include <QPushButton>
#include <QDialog>
#include "common_data.h"

class BSplineVolumeSetDialog : public QDialog
{
	Q_OBJECT
private:
	QSpinBox *orderSpinBox[3];				// 三个方向的阶数
	QSpinBox *ctrlPointSpinBox[3];			// 三个方向的控制顶点数
	QPushButton *cancelButton, *okButton;
	CommonData *m_pCommonData;
public:
	BSplineVolumeSetDialog(CommonData *commonData, QWidget *parent = 0);
private slots:
	void orderUChanged(int orderU);
	void orderVChanged(int orderV);
	void orderWChanged(int orderW);
	void ctrlPointUChanged(int ctrlPointU);
	void ctrlPointVChanged(int ctrlPointV);
	void ctrlPointWChanged(int ctrlPointW);
	void okToChange();
};

#endif
