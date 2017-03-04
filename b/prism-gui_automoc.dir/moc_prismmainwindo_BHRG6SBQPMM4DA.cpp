/****************************************************************************
** Meta object code from reading C++ file 'prismmainwindow.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.8.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../Qt/prismmainwindow.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'prismmainwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.8.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_PRISMMainWindow_t {
    QByteArrayData data[12];
    char stringdata0[221];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_PRISMMainWindow_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_PRISMMainWindow_t qt_meta_stringdata_PRISMMainWindow = {
    {
QT_MOC_LITERAL(0, 0, 15), // "PRISMMainWindow"
QT_MOC_LITERAL(1, 16, 22), // "setInterpolationFactor"
QT_MOC_LITERAL(2, 39, 0), // ""
QT_MOC_LITERAL(3, 40, 29), // "setFilenameAtoms_fromLineEdit"
QT_MOC_LITERAL(4, 70, 27), // "setFilenameAtoms_fromDialog"
QT_MOC_LITERAL(5, 98, 30), // "setFilenameOutput_fromLineEdit"
QT_MOC_LITERAL(6, 129, 28), // "setFilenameOutput_fromDialog"
QT_MOC_LITERAL(7, 158, 16), // "setFilenameAtoms"
QT_MOC_LITERAL(8, 175, 11), // "std::string"
QT_MOC_LITERAL(9, 187, 8), // "filename"
QT_MOC_LITERAL(10, 196, 17), // "setFilenameOutput"
QT_MOC_LITERAL(11, 214, 6) // "launch"

    },
    "PRISMMainWindow\0setInterpolationFactor\0"
    "\0setFilenameAtoms_fromLineEdit\0"
    "setFilenameAtoms_fromDialog\0"
    "setFilenameOutput_fromLineEdit\0"
    "setFilenameOutput_fromDialog\0"
    "setFilenameAtoms\0std::string\0filename\0"
    "setFilenameOutput\0launch"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_PRISMMainWindow[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       8,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   54,    2, 0x0a /* Public */,
       3,    0,   55,    2, 0x0a /* Public */,
       4,    0,   56,    2, 0x0a /* Public */,
       5,    0,   57,    2, 0x0a /* Public */,
       6,    0,   58,    2, 0x0a /* Public */,
       7,    1,   59,    2, 0x0a /* Public */,
      10,    1,   62,    2, 0x0a /* Public */,
      11,    0,   65,    2, 0x0a /* Public */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 8,    9,
    QMetaType::Void, 0x80000000 | 8,    9,
    QMetaType::Void,

       0        // eod
};

void PRISMMainWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        PRISMMainWindow *_t = static_cast<PRISMMainWindow *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->setInterpolationFactor(); break;
        case 1: _t->setFilenameAtoms_fromLineEdit(); break;
        case 2: _t->setFilenameAtoms_fromDialog(); break;
        case 3: _t->setFilenameOutput_fromLineEdit(); break;
        case 4: _t->setFilenameOutput_fromDialog(); break;
        case 5: _t->setFilenameAtoms((*reinterpret_cast< const std::string(*)>(_a[1]))); break;
        case 6: _t->setFilenameOutput((*reinterpret_cast< const std::string(*)>(_a[1]))); break;
        case 7: _t->launch(); break;
        default: ;
        }
    }
}

const QMetaObject PRISMMainWindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_PRISMMainWindow.data,
      qt_meta_data_PRISMMainWindow,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *PRISMMainWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *PRISMMainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_PRISMMainWindow.stringdata0))
        return static_cast<void*>(const_cast< PRISMMainWindow*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int PRISMMainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 8)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 8;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 8)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 8;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
