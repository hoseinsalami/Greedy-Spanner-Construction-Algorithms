/****************************************************************************
** Meta object code from reading C++ file 'window.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.5)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "window.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'window.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.5. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_Window[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       8,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       4,       // signalCount

 // signals: signature, parameters, type, tag, flags
       8,    7,    7,    7, 0x05,
      20,    7,    7,    7, 0x05,
      35,    7,    7,    7, 0x05,
      54,    7,    7,    7, 0x05,

 // slots: signature, parameters, type, tag, flags
      71,    7,    7,    7, 0x08,
      81,    7,    7,    7, 0x08,
      92,    7,    7,    7, 0x08,
     100,    7,    7,    7, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_Window[] = {
    "Window\0\0solved(int)\0newTask(Task*)\0"
    "quitWorkerThread()\0stopWorkerTask()\0"
    "aboutQt()\0generate()\0solve()\0"
    "createWorker()\0"
};

void Window::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        Window *_t = static_cast<Window *>(_o);
        switch (_id) {
        case 0: _t->solved((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: _t->newTask((*reinterpret_cast< Task*(*)>(_a[1]))); break;
        case 2: _t->quitWorkerThread(); break;
        case 3: _t->stopWorkerTask(); break;
        case 4: _t->aboutQt(); break;
        case 5: _t->generate(); break;
        case 6: _t->solve(); break;
        case 7: _t->createWorker(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData Window::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject Window::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_Window,
      qt_meta_data_Window, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &Window::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *Window::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *Window::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_Window))
        return static_cast<void*>(const_cast< Window*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int Window::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 8)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 8;
    }
    return _id;
}

// SIGNAL 0
void Window::solved(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void Window::newTask(Task * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void Window::quitWorkerThread()
{
    QMetaObject::activate(this, &staticMetaObject, 2, 0);
}

// SIGNAL 3
void Window::stopWorkerTask()
{
    QMetaObject::activate(this, &staticMetaObject, 3, 0);
}
QT_END_MOC_NAMESPACE
