/****************************************************************************
** Meta object code from reading C++ file 'prism_progressbar.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.8.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "prism_progressbar.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'prism_progressbar.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.8.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_prism_progressbar_t {
    QByteArrayData data[12];
    char stringdata0[185];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_prism_progressbar_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_prism_progressbar_t qt_meta_stringdata_prism_progressbar = {
    {
QT_MOC_LITERAL(0, 0, 17), // "prism_progressbar"
QT_MOC_LITERAL(1, 18, 24), // "updateDescriptionMessage"
QT_MOC_LITERAL(2, 43, 0), // ""
QT_MOC_LITERAL(3, 44, 3), // "str"
QT_MOC_LITERAL(4, 48, 16), // "updateCalcStatus"
QT_MOC_LITERAL(5, 65, 17), // "updateProgressBar"
QT_MOC_LITERAL(6, 83, 5), // "value"
QT_MOC_LITERAL(7, 89, 16), // "setStepPotential"
QT_MOC_LITERAL(8, 106, 27), // "update_calculatingPotential"
QT_MOC_LITERAL(9, 134, 23), // "updateCalcStatusMessage"
QT_MOC_LITERAL(10, 158, 17), // "updateDescription"
QT_MOC_LITERAL(11, 176, 8) // "setTitle"

    },
    "prism_progressbar\0updateDescriptionMessage\0"
    "\0str\0updateCalcStatus\0updateProgressBar\0"
    "value\0setStepPotential\0"
    "update_calculatingPotential\0"
    "updateCalcStatusMessage\0updateDescription\0"
    "setTitle"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_prism_progressbar[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       8,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       3,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   54,    2, 0x06 /* Public */,
       4,    1,   57,    2, 0x06 /* Public */,
       5,    1,   60,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       7,    0,   63,    2, 0x0a /* Public */,
       8,    2,   64,    2, 0x0a /* Public */,
       9,    1,   69,    2, 0x0a /* Public */,
      10,    1,   72,    2, 0x0a /* Public */,
      11,    1,   75,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void, QMetaType::QString,    3,
    QMetaType::Void, QMetaType::QString,    3,
    QMetaType::Void, QMetaType::Int,    6,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void, QMetaType::Long, QMetaType::Long,    2,    2,
    QMetaType::Void, QMetaType::QString,    3,
    QMetaType::Void, QMetaType::QString,    3,
    QMetaType::Void, QMetaType::QString,    3,

       0        // eod
};

void prism_progressbar::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        prism_progressbar *_t = static_cast<prism_progressbar *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->updateDescriptionMessage((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 1: _t->updateCalcStatus((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 2: _t->updateProgressBar((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 3: _t->setStepPotential(); break;
        case 4: _t->update_calculatingPotential((*reinterpret_cast< long(*)>(_a[1])),(*reinterpret_cast< long(*)>(_a[2]))); break;
        case 5: _t->updateCalcStatusMessage((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 6: _t->updateDescription((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 7: _t->setTitle((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (prism_progressbar::*_t)(const QString );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&prism_progressbar::updateDescriptionMessage)) {
                *result = 0;
                return;
            }
        }
        {
            typedef void (prism_progressbar::*_t)(const QString );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&prism_progressbar::updateCalcStatus)) {
                *result = 1;
                return;
            }
        }
        {
            typedef void (prism_progressbar::*_t)(int );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&prism_progressbar::updateProgressBar)) {
                *result = 2;
                return;
            }
        }
    }
}

const QMetaObject prism_progressbar::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_prism_progressbar.data,
      qt_meta_data_prism_progressbar,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *prism_progressbar::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *prism_progressbar::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_prism_progressbar.stringdata0))
        return static_cast<void*>(const_cast< prism_progressbar*>(this));
    return QDialog::qt_metacast(_clname);
}

int prism_progressbar::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
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

// SIGNAL 0
void prism_progressbar::updateDescriptionMessage(const QString _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void prism_progressbar::updateCalcStatus(const QString _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void prism_progressbar::updateProgressBar(int _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
