/****************************************************************************
** Meta object code from reading C++ file 'config.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.5)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "config.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'config.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.5. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_Config[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       0,    0, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

       0        // eod
};

static const char qt_meta_stringdata_Config[] = {
    "Config\0"
};

void Config::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    Q_UNUSED(_o);
    Q_UNUSED(_id);
    Q_UNUSED(_c);
    Q_UNUSED(_a);
}

const QMetaObjectExtraData Config::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject Config::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_Config,
      qt_meta_data_Config, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &Config::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *Config::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *Config::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_Config))
        return static_cast<void*>(const_cast< Config*>(this));
    return QWidget::qt_metacast(_clname);
}

int Config::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    return _id;
}
static const uint qt_meta_data_Vector3DItem[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       7,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      14,   13,   13,   13, 0x05,

 // slots: signature, parameters, type, tag, flags
      41,   39,   13,   13, 0x0a,
      70,   62,   13,   13, 0x0a,
      94,   39,   13,   13, 0x0a,
     118,  116,   13,   13, 0x08,
     133,  131,   13,   13, 0x08,
     148,  146,   13,   13, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_Vector3DItem[] = {
    "Vector3DItem\0\0vectorChanged(QVector3D)\0"
    "v\0setVector(QVector3D)\0min,max\0"
    "setRange(double,double)\0setSingleStep(double)\0"
    "x\0setX(double)\0y\0setY(double)\0z\0"
    "setZ(double)\0"
};

void Vector3DItem::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        Vector3DItem *_t = static_cast<Vector3DItem *>(_o);
        switch (_id) {
        case 0: _t->vectorChanged((*reinterpret_cast< QVector3D(*)>(_a[1]))); break;
        case 1: _t->setVector((*reinterpret_cast< QVector3D(*)>(_a[1]))); break;
        case 2: _t->setRange((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2]))); break;
        case 3: _t->setSingleStep((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 4: _t->setX((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 5: _t->setY((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 6: _t->setZ((*reinterpret_cast< double(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData Vector3DItem::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject Vector3DItem::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_Vector3DItem,
      qt_meta_data_Vector3DItem, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &Vector3DItem::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *Vector3DItem::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *Vector3DItem::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_Vector3DItem))
        return static_cast<void*>(const_cast< Vector3DItem*>(this));
    return QObject::qt_metacast(_clname);
}

int Vector3DItem::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 7)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 7;
    }
    return _id;
}

// SIGNAL 0
void Vector3DItem::vectorChanged(QVector3D _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
static const uint qt_meta_data_RangeItem[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       7,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      13,   11,   10,   10, 0x05,

 // slots: signature, parameters, type, tag, flags
      49,   41,   10,   10, 0x0a,
      75,   73,   10,   10, 0x0a,
      97,   73,   10,   10, 0x0a,
     124,  114,   10,   10, 0x0a,
     148,   73,   10,   10, 0x08,
     163,   73,   10,   10, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_RangeItem[] = {
    "RangeItem\0\0,\0valueChanged(double,double)\0"
    "min,max\0setRange(double,double)\0v\0"
    "setSingleStep(double)\0setDecimals(int)\0"
    "vmin,vmax\0setValue(double,double)\0"
    "setMin(double)\0setMax(double)\0"
};

void RangeItem::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        RangeItem *_t = static_cast<RangeItem *>(_o);
        switch (_id) {
        case 0: _t->valueChanged((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2]))); break;
        case 1: _t->setRange((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2]))); break;
        case 2: _t->setSingleStep((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 3: _t->setDecimals((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 4: _t->setValue((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2]))); break;
        case 5: _t->setMin((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 6: _t->setMax((*reinterpret_cast< double(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData RangeItem::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject RangeItem::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_RangeItem,
      qt_meta_data_RangeItem, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &RangeItem::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *RangeItem::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *RangeItem::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_RangeItem))
        return static_cast<void*>(const_cast< RangeItem*>(this));
    return QObject::qt_metacast(_clname);
}

int RangeItem::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 7)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 7;
    }
    return _id;
}

// SIGNAL 0
void RangeItem::valueChanged(double _t1, double _t2)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
static const uint qt_meta_data_DoubleItem[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       6,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      12,   11,   11,   11, 0x05,

 // slots: signature, parameters, type, tag, flags
      41,   33,   11,   11, 0x0a,
      67,   65,   11,   11, 0x0a,
      89,   65,   11,   11, 0x0a,
     106,   65,   11,   11, 0x0a,
     123,   65,   11,   11, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_DoubleItem[] = {
    "DoubleItem\0\0valueChanged(double)\0"
    "min,max\0setRange(double,double)\0v\0"
    "setSingleStep(double)\0setDecimals(int)\0"
    "setValue(double)\0setValueProxy(double)\0"
};

void DoubleItem::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        DoubleItem *_t = static_cast<DoubleItem *>(_o);
        switch (_id) {
        case 0: _t->valueChanged((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 1: _t->setRange((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2]))); break;
        case 2: _t->setSingleStep((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 3: _t->setDecimals((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 4: _t->setValue((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 5: _t->setValueProxy((*reinterpret_cast< double(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData DoubleItem::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject DoubleItem::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_DoubleItem,
      qt_meta_data_DoubleItem, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &DoubleItem::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *DoubleItem::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *DoubleItem::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_DoubleItem))
        return static_cast<void*>(const_cast< DoubleItem*>(this));
    return QObject::qt_metacast(_clname);
}

int DoubleItem::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 6)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 6;
    }
    return _id;
}

// SIGNAL 0
void DoubleItem::valueChanged(double _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
static const uint qt_meta_data_IntegerItem[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       5,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      13,   12,   12,   12, 0x05,

 // slots: signature, parameters, type, tag, flags
      39,   31,   12,   12, 0x0a,
      59,   57,   12,   12, 0x0a,
      78,   57,   12,   12, 0x0a,
      92,   57,   12,   12, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_IntegerItem[] = {
    "IntegerItem\0\0valueChanged(int)\0min,max\0"
    "setRange(int,int)\0v\0setSingleStep(int)\0"
    "setValue(int)\0setValueProxy(int)\0"
};

void IntegerItem::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        IntegerItem *_t = static_cast<IntegerItem *>(_o);
        switch (_id) {
        case 0: _t->valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: _t->setRange((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 2: _t->setSingleStep((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 3: _t->setValue((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 4: _t->setValueProxy((*reinterpret_cast< int(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData IntegerItem::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject IntegerItem::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_IntegerItem,
      qt_meta_data_IntegerItem, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &IntegerItem::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *IntegerItem::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *IntegerItem::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_IntegerItem))
        return static_cast<void*>(const_cast< IntegerItem*>(this));
    return QObject::qt_metacast(_clname);
}

int IntegerItem::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 5)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 5;
    }
    return _id;
}

// SIGNAL 0
void IntegerItem::valueChanged(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
static const uint qt_meta_data_LabelItem[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      13,   11,   10,   10, 0x0a,
      31,   11,   10,   10, 0x0a,
      48,   11,   10,   10, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_LabelItem[] = {
    "LabelItem\0\0v\0setValue(QString)\0"
    "setValue(double)\0setValue(int)\0"
};

void LabelItem::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        LabelItem *_t = static_cast<LabelItem *>(_o);
        switch (_id) {
        case 0: _t->setValue((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 1: _t->setValue((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 2: _t->setValue((*reinterpret_cast< int(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData LabelItem::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject LabelItem::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_LabelItem,
      qt_meta_data_LabelItem, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &LabelItem::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *LabelItem::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *LabelItem::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_LabelItem))
        return static_cast<void*>(const_cast< LabelItem*>(this));
    return QObject::qt_metacast(_clname);
}

int LabelItem::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 3)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 3;
    }
    return _id;
}
static const uint qt_meta_data_BoolItem[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      10,    9,    9,    9, 0x05,

 // slots: signature, parameters, type, tag, flags
      31,   29,    9,    9, 0x0a,
      46,   29,    9,    9, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_BoolItem[] = {
    "BoolItem\0\0valueChanged(bool)\0v\0"
    "setValue(bool)\0setValueProxy(bool)\0"
};

void BoolItem::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        BoolItem *_t = static_cast<BoolItem *>(_o);
        switch (_id) {
        case 0: _t->valueChanged((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 1: _t->setValue((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 2: _t->setValueProxy((*reinterpret_cast< bool(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData BoolItem::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject BoolItem::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_BoolItem,
      qt_meta_data_BoolItem, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &BoolItem::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *BoolItem::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *BoolItem::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_BoolItem))
        return static_cast<void*>(const_cast< BoolItem*>(this));
    return QObject::qt_metacast(_clname);
}

int BoolItem::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 3)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 3;
    }
    return _id;
}

// SIGNAL 0
void BoolItem::valueChanged(bool _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
static const uint qt_meta_data_ColorItem[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      11,   10,   10,   10, 0x05,

 // slots: signature, parameters, type, tag, flags
      34,   32,   10,   10, 0x0a,
      51,   32,   10,   10, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_ColorItem[] = {
    "ColorItem\0\0valueChanged(QColor)\0v\0"
    "setValue(QColor)\0setValueProxy(QColor)\0"
};

void ColorItem::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        ColorItem *_t = static_cast<ColorItem *>(_o);
        switch (_id) {
        case 0: _t->valueChanged((*reinterpret_cast< QColor(*)>(_a[1]))); break;
        case 1: _t->setValue((*reinterpret_cast< QColor(*)>(_a[1]))); break;
        case 2: _t->setValueProxy((*reinterpret_cast< QColor(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData ColorItem::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject ColorItem::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_ColorItem,
      qt_meta_data_ColorItem, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &ColorItem::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *ColorItem::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *ColorItem::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_ColorItem))
        return static_cast<void*>(const_cast< ColorItem*>(this));
    return QObject::qt_metacast(_clname);
}

int ColorItem::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 3)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 3;
    }
    return _id;
}

// SIGNAL 0
void ColorItem::valueChanged(QColor _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
static const uint qt_meta_data_TextBoxItem[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      13,   12,   12,   12, 0x05,

 // slots: signature, parameters, type, tag, flags
      37,   35,   12,   12, 0x0a,
      55,   12,   12,   12, 0x08,
      65,   12,   12,   12, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_TextBoxItem[] = {
    "TextBoxItem\0\0valueChanged(QString)\0v\0"
    "setValue(QString)\0changed()\0send()\0"
};

void TextBoxItem::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        TextBoxItem *_t = static_cast<TextBoxItem *>(_o);
        switch (_id) {
        case 0: _t->valueChanged((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 1: _t->setValue((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 2: _t->changed(); break;
        case 3: _t->send(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData TextBoxItem::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject TextBoxItem::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_TextBoxItem,
      qt_meta_data_TextBoxItem, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &TextBoxItem::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *TextBoxItem::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *TextBoxItem::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_TextBoxItem))
        return static_cast<void*>(const_cast< TextBoxItem*>(this));
    return QObject::qt_metacast(_clname);
}

int TextBoxItem::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 4)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 4;
    }
    return _id;
}

// SIGNAL 0
void TextBoxItem::valueChanged(QString _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
static const uint qt_meta_data_OptionItem[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      12,   11,   11,   11, 0x05,

 // slots: signature, parameters, type, tag, flags
      32,   30,   11,   11, 0x0a,
      46,   11,   11,   11, 0x0a,
      59,   30,   11,   11, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_OptionItem[] = {
    "OptionItem\0\0valueChanged(int)\0i\0"
    "setValue(int)\0setToOther()\0"
    "setValueProxy(int)\0"
};

void OptionItem::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        OptionItem *_t = static_cast<OptionItem *>(_o);
        switch (_id) {
        case 0: _t->valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: _t->setValue((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: _t->setToOther(); break;
        case 3: _t->setValueProxy((*reinterpret_cast< int(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData OptionItem::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject OptionItem::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_OptionItem,
      qt_meta_data_OptionItem, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &OptionItem::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *OptionItem::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *OptionItem::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_OptionItem))
        return static_cast<void*>(const_cast< OptionItem*>(this));
    return QObject::qt_metacast(_clname);
}

int OptionItem::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 4)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 4;
    }
    return _id;
}

// SIGNAL 0
void OptionItem::valueChanged(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
static const uint qt_meta_data_ButtonItem[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      12,   11,   11,   11, 0x0a,
      21,   11,   11,   11, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_ButtonItem[] = {
    "ButtonItem\0\0enable()\0disable()\0"
};

void ButtonItem::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        ButtonItem *_t = static_cast<ButtonItem *>(_o);
        switch (_id) {
        case 0: _t->enable(); break;
        case 1: _t->disable(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObjectExtraData ButtonItem::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject ButtonItem::staticMetaObject = {
    { &QPushButton::staticMetaObject, qt_meta_stringdata_ButtonItem,
      qt_meta_data_ButtonItem, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &ButtonItem::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *ButtonItem::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *ButtonItem::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_ButtonItem))
        return static_cast<void*>(const_cast< ButtonItem*>(this));
    return QPushButton::qt_metacast(_clname);
}

int ButtonItem::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QPushButton::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 2)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 2;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
