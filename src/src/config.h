/*
    Surface viewer - A 3D Mathematical surface inspection & renderer tool
    Copyright (C) 2011  B.J. Conijn <b.j.conijn@student.tue.nl>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CONFIG_H
#define CONFIG_H

#include <QtCore/QObject>
#include <QtCore/QVariant>
#include <QtGui/QWidget>
#include <QtGui/QVBoxLayout>
#include <QtGui/QVector3D>
#include <QPushButton>

//Basically fixes the broken Qt macro, which interferes with C++11 user-defined literals.
#undef QLOCATION
#define QLOCATION "\0" __FILE__ ":" QTOSTRING(__LINE__)

class QSpinBox;
namespace KTextEditor {
class Document;
class View;
}

class QComboBox;
class QTextEdit;
class QColorDialog;
class QPushButton;
class QCheckBox;
class QLabel;
class QDoubleSpinBox;
class Config;

class Config : public QWidget {
    Q_OBJECT;

public:
    explicit Config(QWidget *parent = 0, Qt::WindowFlags f = 0);
    virtual void done();

    QVBoxLayout *vbox;
};

class Vector3DItem : public QObject {
	Q_OBJECT;
private:
	QVector3D value;
	QDoubleSpinBox *x, *y, *z;
	bool changing; // Used to prevent a value changed event when setting the value.
public:
	explicit Vector3DItem(Config * parent, QString name);
public slots:
	virtual void setVector(QVector3D v);
	virtual void setRange(double min, double max);
	virtual void setSingleStep(double v);
private slots:
	virtual void setX(double x);
	virtual void setY(double y);
	virtual void setZ(double z);
signals:
	void vectorChanged(QVector3D);
};

class RangeItem : public QObject {
	Q_OBJECT;
private:
	double cmin, cmax;
	QDoubleSpinBox *smin;
	QDoubleSpinBox *smax;
	bool changing; // Used to prevent a value changed event when setting the value.
public:
	explicit RangeItem(Config * parent, QString name);
public slots:
	virtual void setRange(double min, double max);
	virtual void setSingleStep(double v);
	virtual void setDecimals(int v);
	virtual void setValue(double vmin, double vmax);
private slots:
	virtual void setMin(double v);
	virtual void setMax(double v);
signals:
	void valueChanged(double,double);
};

class DoubleItem : public QObject {
	Q_OBJECT;
private:
	QDoubleSpinBox *spinner;
	bool changing; // Used to prevent a value changed event when setting the value.
	double * var_ptr;
public:
	explicit DoubleItem(Config * parent, QString name, double * var=NULL);
public slots:
	virtual void setRange(double min, double max);
	virtual void setSingleStep(double v);
	virtual void setDecimals(int v);
	virtual void setValue(double v);
private slots:
	virtual void setValueProxy(double v);
signals:
	void valueChanged(double);
};

class IntegerItem : public QObject {
    Q_OBJECT;
private:
    QSpinBox *spinner;
    bool changing; // Used to prevent a value changed event when setting the value.
    int * var_ptr;
public:
    explicit IntegerItem(Config * parent, QString name, int * var=NULL);
public slots:
    virtual void setRange(int min, int max);
    virtual void setSingleStep(int v);
    virtual void setValue(int v);
private slots:
    virtual void setValueProxy(int v);
signals:
    void valueChanged(int);
};

class LabelItem : public QObject {
	Q_OBJECT;
private:
	QLabel * label;
public:
	explicit LabelItem(Config *parent, QString name);
public slots:
	virtual void setValue(QString v);
	virtual void setValue(double v);
	virtual void setValue(int v);
};

class BoolItem : public QObject {
	Q_OBJECT;
private:
	QCheckBox * cbox;
	bool changing; // Used to prevent a value changed event when setting the value.
    bool * var_ptr;
public:
	explicit BoolItem(Config * parent, QString name, bool * var=NULL);
public slots:
	virtual void setValue(bool v);
private slots:
	virtual void setValueProxy(bool v);
signals:
	void valueChanged(bool);
};

class ColorItem : public QObject {
	Q_OBJECT;
private:
	QPushButton * button;
	QColorDialog * col;
	bool alpha;
public:
	explicit ColorItem(Config * parent, QString name, bool alpha=false);
public slots:
	virtual void setValue(QColor v);
private slots:
	virtual void setValueProxy(QColor v);
signals:
	void valueChanged(QColor);
};

class TextBoxItem : public QObject {
	Q_OBJECT;
private:
#ifdef KDE4_FOUND
	KTextEditor::Document * edit;
	KTextEditor::View * view;
	QVariant def_color;
#else
	QTextEdit * edit;
#endif
	QTimer * timer;
	bool changing; // Used to prevent a value changed event when setting the value.
public:
	explicit TextBoxItem(Config * parent, QString name);
	void setLanguage(QString lang);
public slots:
	virtual void setValue(QString v);
private slots:
	void changed();
	void send();
signals:
	void valueChanged(QString);
};

class OptionItem : public QObject {
	Q_OBJECT;
private:
	QComboBox * list;
	bool changing; // Used to prevent a value changed event when setting the value.
    int * var_ptr;
public:
	explicit OptionItem(Config * parent, QString name, int * var = NULL);
	int addItem(QString text);
	int addSeparator();
	int other;
public slots:
	virtual void setValue(int i);
	virtual void setToOther();
private slots:
	virtual void setValueProxy(int i);
signals:
	void valueChanged(int);
};

class ButtonItem : public QPushButton {
    Q_OBJECT;
public:
    explicit ButtonItem(Config * parent, QString caption);
public slots:
    void enable();
    void disable();
};

#endif // CONFIG_H
