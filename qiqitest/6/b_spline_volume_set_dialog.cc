#include "b_spline_volume_set_dialog.h"
#include <QGroupBox>
#include <QLabel>
#include <QHBoxLayout>
#include <QVBoxLayout>

BSplineVolumeSetDialog::BSplineVolumeSetDialog(CommonData *commonData, QWidget *parent)
: QDialog(parent)
{
	m_pCommonData = commonData;
	setWindowTitle(tr("设置B样条体"));

	/* 阶数框 */
	QGroupBox *orderGroupBox = new QGroupBox(this);
	orderGroupBox->setTitle(tr("阶数"));

	QLabel *uOrderLabel = new QLabel(tr("u:"));
	QLabel *vOrderLabel = new QLabel(tr("v:"));
	QLabel *wOrderLabel = new QLabel(tr("w:"));

	for (int i = 0; i < 3; ++i)
	{
		orderSpinBox[i] = new QSpinBox(orderGroupBox);
		orderSpinBox[i]->setRange(2, 4);							// 2-4阶
		orderSpinBox[i]->setValue(m_pCommonData->order(i));
	}

	QHBoxLayout *orderLayout = new QHBoxLayout(orderGroupBox);
	orderLayout->addWidget(uOrderLabel);
	orderLayout->addWidget(orderSpinBox[0]);
	orderLayout->addWidget(vOrderLabel);
	orderLayout->addWidget(orderSpinBox[1]);
	orderLayout->addWidget(wOrderLabel);
	orderLayout->addWidget(orderSpinBox[2]);

	/* 控制顶点数框 */
	QGroupBox *ctrlPointGroupBox = new QGroupBox(this);
	ctrlPointGroupBox->setTitle(tr("控制顶点数"));

	QLabel *ctrlPointULabel = new QLabel(tr("u:"));
	QLabel *ctrlPointVLabel = new QLabel(tr("v:"));
	QLabel *ctrlPointWLabel = new QLabel(tr("w:"));

	for (int i = 0; i < 3; ++i)
	{
		ctrlPointSpinBox[i] = new QSpinBox(ctrlPointGroupBox);
		ctrlPointSpinBox[i]->setRange(2, 15);						// 2-15个控制顶点
		ctrlPointSpinBox[i]->setValue(m_pCommonData->ctrlPointCount(i));
	}

	QHBoxLayout *ctrlPointLayout = new QHBoxLayout(ctrlPointGroupBox);
	ctrlPointLayout->addWidget(ctrlPointULabel);
	ctrlPointLayout->addWidget(ctrlPointSpinBox[0]);
	ctrlPointLayout->addWidget(ctrlPointVLabel);
	ctrlPointLayout->addWidget(ctrlPointSpinBox[1]);
	ctrlPointLayout->addWidget(ctrlPointWLabel);
	ctrlPointLayout->addWidget(ctrlPointSpinBox[2]);

	/* 取消、确定按钮 */
	cancelButton = new QPushButton;
	cancelButton->setText(tr("取消"));

	okButton = new QPushButton;
	okButton->setText(tr("确定"));

	QHBoxLayout *cancelOKLayout = new QHBoxLayout;
	cancelOKLayout->addWidget(cancelButton);
	cancelOKLayout->addWidget(okButton);

	QVBoxLayout *mainLayout = new QVBoxLayout(this);
	mainLayout->addWidget(orderGroupBox);
	mainLayout->addWidget(ctrlPointGroupBox);
	mainLayout->addLayout(cancelOKLayout);

	/* 连接 */
	connect(orderSpinBox[0], SIGNAL(valueChanged(int)), this, SLOT(orderUChanged(int)));
	connect(orderSpinBox[1], SIGNAL(valueChanged(int)), this, SLOT(orderVChanged(int)));
	connect(orderSpinBox[2], SIGNAL(valueChanged(int)), this, SLOT(orderWChanged(int)));
	connect(ctrlPointSpinBox[0], SIGNAL(valueChanged(int)), this, SLOT(ctrlPointUChanged(int)));
	connect(ctrlPointSpinBox[1], SIGNAL(valueChanged(int)), this, SLOT(ctrlPointVChanged(int)));
	connect(ctrlPointSpinBox[2], SIGNAL(valueChanged(int)), this, SLOT(ctrlPointWChanged(int)));
	connect(cancelButton, SIGNAL(clicked()), this, SLOT(reject()));
	connect(cancelButton, SIGNAL(clicked()), this, SLOT(close()));
	connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
	connect(okButton, SIGNAL(clicked()), this, SLOT(okToChange()));
}

void BSplineVolumeSetDialog::orderUChanged(int orderU)
{
	if (orderU > ctrlPointSpinBox[U]->value())
		ctrlPointSpinBox[U]->setValue(orderU);
}

void BSplineVolumeSetDialog::orderVChanged(int orderV)
{
	if (orderV > ctrlPointSpinBox[V]->value())
		ctrlPointSpinBox[V]->setValue(orderV);
}

void BSplineVolumeSetDialog::orderWChanged(int orderW)
{
	if (orderW > ctrlPointSpinBox[W]->value())
		ctrlPointSpinBox[W]->setValue(orderW);
}

void BSplineVolumeSetDialog::ctrlPointUChanged(int ctrlPointU)
{
	if (ctrlPointU < orderSpinBox[U]->value())
		orderSpinBox[U]->setValue(ctrlPointU);
}

void BSplineVolumeSetDialog::ctrlPointVChanged(int ctrlPointV)
{
	if (ctrlPointV < orderSpinBox[V]->value())
		orderSpinBox[V]->setValue(ctrlPointV);
}

void BSplineVolumeSetDialog::ctrlPointWChanged(int ctrlPointW)
{
	if (ctrlPointW < orderSpinBox[W]->value())
		orderSpinBox[W]->setValue(ctrlPointW);
}

void BSplineVolumeSetDialog::okToChange()
{
	for (int i = 0; i < 3; ++i)
	{
		m_pCommonData->setOrder(orderSpinBox[i]->value(), i);
		m_pCommonData->setCtrlPointCount(ctrlPointSpinBox[i]->value(), i);
	}
	close();
}
