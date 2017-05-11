/****************************************************************************
** Meta object code from reading C++ file 'prism_qthreads.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.9.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "prism_qthreads.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'prism_qthreads.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.9.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_PotentialThread_t {
    QByteArrayData data[1];
    char stringdata0[16];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_PotentialThread_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_PotentialThread_t qt_meta_stringdata_PotentialThread = {
    {
QT_MOC_LITERAL(0, 0, 15) // "PotentialThread"

    },
    "PotentialThread"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_PotentialThread[] = {

 // content:
       7,       // revision
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

void PotentialThread::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    Q_UNUSED(_o);
    Q_UNUSED(_id);
    Q_UNUSED(_c);
    Q_UNUSED(_a);
}

const QMetaObject PotentialThread::staticMetaObject = {
    { &QThread::staticMetaObject, qt_meta_stringdata_PotentialThread.data,
      qt_meta_data_PotentialThread,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *PotentialThread::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *PotentialThread::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_PotentialThread.stringdata0))
        return static_cast<void*>(const_cast< PotentialThread*>(this));
    return QThread::qt_metacast(_clname);
}

int PotentialThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QThread::qt_metacall(_c, _id, _a);
    return _id;
}
struct qt_meta_stringdata_SMatrixThread_t {
    QByteArrayData data[4];
    char stringdata0[54];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_SMatrixThread_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_SMatrixThread_t qt_meta_stringdata_SMatrixThread = {
    {
QT_MOC_LITERAL(0, 0, 13), // "SMatrixThread"
QT_MOC_LITERAL(1, 14, 19), // "potentialCalculated"
QT_MOC_LITERAL(2, 34, 0), // ""
QT_MOC_LITERAL(3, 35, 18) // "ScompactCalculated"

    },
    "SMatrixThread\0potentialCalculated\0\0"
    "ScompactCalculated"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_SMatrixThread[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       2,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    0,   24,    2, 0x06 /* Public */,
       3,    0,   25,    2, 0x06 /* Public */,

 // signals: parameters
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void SMatrixThread::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        SMatrixThread *_t = static_cast<SMatrixThread *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->potentialCalculated(); break;
        case 1: _t->ScompactCalculated(); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (SMatrixThread::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&SMatrixThread::potentialCalculated)) {
                *result = 0;
                return;
            }
        }
        {
            typedef void (SMatrixThread::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&SMatrixThread::ScompactCalculated)) {
                *result = 1;
                return;
            }
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject SMatrixThread::staticMetaObject = {
    { &QThread::staticMetaObject, qt_meta_stringdata_SMatrixThread.data,
      qt_meta_data_SMatrixThread,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *SMatrixThread::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *SMatrixThread::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_SMatrixThread.stringdata0))
        return static_cast<void*>(const_cast< SMatrixThread*>(this));
    return QThread::qt_metacast(_clname);
}

int SMatrixThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QThread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 2)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 2;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 2)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 2;
    }
    return _id;
}

// SIGNAL 0
void SMatrixThread::potentialCalculated()
{
    QMetaObject::activate(this, &staticMetaObject, 0, nullptr);
}

// SIGNAL 1
void SMatrixThread::ScompactCalculated()
{
    QMetaObject::activate(this, &staticMetaObject, 1, nullptr);
}
struct qt_meta_stringdata_ProbeThread_t {
    QByteArrayData data[16];
    char stringdata0[274];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_ProbeThread_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_ProbeThread_t qt_meta_stringdata_ProbeThread = {
    {
QT_MOC_LITERAL(0, 0, 11), // "ProbeThread"
QT_MOC_LITERAL(1, 12, 19), // "potentialCalculated"
QT_MOC_LITERAL(2, 32, 0), // ""
QT_MOC_LITERAL(3, 33, 18), // "ScompactCalculated"
QT_MOC_LITERAL(4, 52, 18), // "signalProbeK_PRISM"
QT_MOC_LITERAL(5, 71, 37), // "PRISM::Array2D<PRISM_FLOAT_PR..."
QT_MOC_LITERAL(6, 109, 18), // "signalProbeR_PRISM"
QT_MOC_LITERAL(7, 128, 23), // "signalProbeK_Multislice"
QT_MOC_LITERAL(8, 152, 23), // "signalProbeR_Multislice"
QT_MOC_LITERAL(9, 176, 17), // "signalProbe_diffK"
QT_MOC_LITERAL(10, 194, 17), // "signalProbe_diffR"
QT_MOC_LITERAL(11, 212, 18), // "signal_pearsonReal"
QT_MOC_LITERAL(12, 231, 3), // "str"
QT_MOC_LITERAL(13, 235, 15), // "signal_pearsonK"
QT_MOC_LITERAL(14, 251, 12), // "signal_RReal"
QT_MOC_LITERAL(15, 264, 9) // "signal_RK"

    },
    "ProbeThread\0potentialCalculated\0\0"
    "ScompactCalculated\0signalProbeK_PRISM\0"
    "PRISM::Array2D<PRISM_FLOAT_PRECISION>\0"
    "signalProbeR_PRISM\0signalProbeK_Multislice\0"
    "signalProbeR_Multislice\0signalProbe_diffK\0"
    "signalProbe_diffR\0signal_pearsonReal\0"
    "str\0signal_pearsonK\0signal_RReal\0"
    "signal_RK"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_ProbeThread[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      12,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
      12,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    0,   74,    2, 0x06 /* Public */,
       3,    0,   75,    2, 0x06 /* Public */,
       4,    1,   76,    2, 0x06 /* Public */,
       6,    1,   79,    2, 0x06 /* Public */,
       7,    1,   82,    2, 0x06 /* Public */,
       8,    1,   85,    2, 0x06 /* Public */,
       9,    2,   88,    2, 0x06 /* Public */,
      10,    2,   93,    2, 0x06 /* Public */,
      11,    1,   98,    2, 0x06 /* Public */,
      13,    1,  101,    2, 0x06 /* Public */,
      14,    1,  104,    2, 0x06 /* Public */,
      15,    1,  107,    2, 0x06 /* Public */,

 // signals: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 5,    2,
    QMetaType::Void, 0x80000000 | 5,    2,
    QMetaType::Void, 0x80000000 | 5,    2,
    QMetaType::Void, 0x80000000 | 5,    2,
    QMetaType::Void, 0x80000000 | 5, 0x80000000 | 5,    2,    2,
    QMetaType::Void, 0x80000000 | 5, 0x80000000 | 5,    2,    2,
    QMetaType::Void, QMetaType::QString,   12,
    QMetaType::Void, QMetaType::QString,   12,
    QMetaType::Void, QMetaType::QString,   12,
    QMetaType::Void, QMetaType::QString,   12,

       0        // eod
};

void ProbeThread::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        ProbeThread *_t = static_cast<ProbeThread *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->potentialCalculated(); break;
        case 1: _t->ScompactCalculated(); break;
        case 2: _t->signalProbeK_PRISM((*reinterpret_cast< PRISM::Array2D<PRISM_FLOAT_PRECISION>(*)>(_a[1]))); break;
        case 3: _t->signalProbeR_PRISM((*reinterpret_cast< PRISM::Array2D<PRISM_FLOAT_PRECISION>(*)>(_a[1]))); break;
        case 4: _t->signalProbeK_Multislice((*reinterpret_cast< PRISM::Array2D<PRISM_FLOAT_PRECISION>(*)>(_a[1]))); break;
        case 5: _t->signalProbeR_Multislice((*reinterpret_cast< PRISM::Array2D<PRISM_FLOAT_PRECISION>(*)>(_a[1]))); break;
        case 6: _t->signalProbe_diffK((*reinterpret_cast< PRISM::Array2D<PRISM_FLOAT_PRECISION>(*)>(_a[1])),(*reinterpret_cast< PRISM::Array2D<PRISM_FLOAT_PRECISION>(*)>(_a[2]))); break;
        case 7: _t->signalProbe_diffR((*reinterpret_cast< PRISM::Array2D<PRISM_FLOAT_PRECISION>(*)>(_a[1])),(*reinterpret_cast< PRISM::Array2D<PRISM_FLOAT_PRECISION>(*)>(_a[2]))); break;
        case 8: _t->signal_pearsonReal((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 9: _t->signal_pearsonK((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 10: _t->signal_RReal((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 11: _t->signal_RK((*reinterpret_cast< QString(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (ProbeThread::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&ProbeThread::potentialCalculated)) {
                *result = 0;
                return;
            }
        }
        {
            typedef void (ProbeThread::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&ProbeThread::ScompactCalculated)) {
                *result = 1;
                return;
            }
        }
        {
            typedef void (ProbeThread::*_t)(PRISM::Array2D<PRISM_FLOAT_PRECISION> );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&ProbeThread::signalProbeK_PRISM)) {
                *result = 2;
                return;
            }
        }
        {
            typedef void (ProbeThread::*_t)(PRISM::Array2D<PRISM_FLOAT_PRECISION> );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&ProbeThread::signalProbeR_PRISM)) {
                *result = 3;
                return;
            }
        }
        {
            typedef void (ProbeThread::*_t)(PRISM::Array2D<PRISM_FLOAT_PRECISION> );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&ProbeThread::signalProbeK_Multislice)) {
                *result = 4;
                return;
            }
        }
        {
            typedef void (ProbeThread::*_t)(PRISM::Array2D<PRISM_FLOAT_PRECISION> );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&ProbeThread::signalProbeR_Multislice)) {
                *result = 5;
                return;
            }
        }
        {
            typedef void (ProbeThread::*_t)(PRISM::Array2D<PRISM_FLOAT_PRECISION> , PRISM::Array2D<PRISM_FLOAT_PRECISION> );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&ProbeThread::signalProbe_diffK)) {
                *result = 6;
                return;
            }
        }
        {
            typedef void (ProbeThread::*_t)(PRISM::Array2D<PRISM_FLOAT_PRECISION> , PRISM::Array2D<PRISM_FLOAT_PRECISION> );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&ProbeThread::signalProbe_diffR)) {
                *result = 7;
                return;
            }
        }
        {
            typedef void (ProbeThread::*_t)(QString );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&ProbeThread::signal_pearsonReal)) {
                *result = 8;
                return;
            }
        }
        {
            typedef void (ProbeThread::*_t)(QString );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&ProbeThread::signal_pearsonK)) {
                *result = 9;
                return;
            }
        }
        {
            typedef void (ProbeThread::*_t)(QString );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&ProbeThread::signal_RReal)) {
                *result = 10;
                return;
            }
        }
        {
            typedef void (ProbeThread::*_t)(QString );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&ProbeThread::signal_RK)) {
                *result = 11;
                return;
            }
        }
    }
}

const QMetaObject ProbeThread::staticMetaObject = {
    { &QThread::staticMetaObject, qt_meta_stringdata_ProbeThread.data,
      qt_meta_data_ProbeThread,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *ProbeThread::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *ProbeThread::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_ProbeThread.stringdata0))
        return static_cast<void*>(const_cast< ProbeThread*>(this));
    return QThread::qt_metacast(_clname);
}

int ProbeThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QThread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 12)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 12;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 12)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 12;
    }
    return _id;
}

// SIGNAL 0
void ProbeThread::potentialCalculated()
{
    QMetaObject::activate(this, &staticMetaObject, 0, nullptr);
}

// SIGNAL 1
void ProbeThread::ScompactCalculated()
{
    QMetaObject::activate(this, &staticMetaObject, 1, nullptr);
}

// SIGNAL 2
void ProbeThread::signalProbeK_PRISM(PRISM::Array2D<PRISM_FLOAT_PRECISION> _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void ProbeThread::signalProbeR_PRISM(PRISM::Array2D<PRISM_FLOAT_PRECISION> _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}

// SIGNAL 4
void ProbeThread::signalProbeK_Multislice(PRISM::Array2D<PRISM_FLOAT_PRECISION> _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 4, _a);
}

// SIGNAL 5
void ProbeThread::signalProbeR_Multislice(PRISM::Array2D<PRISM_FLOAT_PRECISION> _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 5, _a);
}

// SIGNAL 6
void ProbeThread::signalProbe_diffK(PRISM::Array2D<PRISM_FLOAT_PRECISION> _t1, PRISM::Array2D<PRISM_FLOAT_PRECISION> _t2)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 6, _a);
}

// SIGNAL 7
void ProbeThread::signalProbe_diffR(PRISM::Array2D<PRISM_FLOAT_PRECISION> _t1, PRISM::Array2D<PRISM_FLOAT_PRECISION> _t2)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 7, _a);
}

// SIGNAL 8
void ProbeThread::signal_pearsonReal(QString _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 8, _a);
}

// SIGNAL 9
void ProbeThread::signal_pearsonK(QString _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 9, _a);
}

// SIGNAL 10
void ProbeThread::signal_RReal(QString _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 10, _a);
}

// SIGNAL 11
void ProbeThread::signal_RK(QString _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 11, _a);
}
struct qt_meta_stringdata_FullPRISMCalcThread_t {
    QByteArrayData data[7];
    char stringdata0[93];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_FullPRISMCalcThread_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_FullPRISMCalcThread_t qt_meta_stringdata_FullPRISMCalcThread = {
    {
QT_MOC_LITERAL(0, 0, 19), // "FullPRISMCalcThread"
QT_MOC_LITERAL(1, 20, 19), // "potentialCalculated"
QT_MOC_LITERAL(2, 40, 0), // ""
QT_MOC_LITERAL(3, 41, 18), // "ScompactCalculated"
QT_MOC_LITERAL(4, 60, 16), // "outputCalculated"
QT_MOC_LITERAL(5, 77, 11), // "signalTitle"
QT_MOC_LITERAL(6, 89, 3) // "str"

    },
    "FullPRISMCalcThread\0potentialCalculated\0"
    "\0ScompactCalculated\0outputCalculated\0"
    "signalTitle\0str"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_FullPRISMCalcThread[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       4,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    0,   34,    2, 0x06 /* Public */,
       3,    0,   35,    2, 0x06 /* Public */,
       4,    0,   36,    2, 0x06 /* Public */,
       5,    1,   37,    2, 0x06 /* Public */,

 // signals: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QString,    6,

       0        // eod
};

void FullPRISMCalcThread::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        FullPRISMCalcThread *_t = static_cast<FullPRISMCalcThread *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->potentialCalculated(); break;
        case 1: _t->ScompactCalculated(); break;
        case 2: _t->outputCalculated(); break;
        case 3: _t->signalTitle((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (FullPRISMCalcThread::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&FullPRISMCalcThread::potentialCalculated)) {
                *result = 0;
                return;
            }
        }
        {
            typedef void (FullPRISMCalcThread::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&FullPRISMCalcThread::ScompactCalculated)) {
                *result = 1;
                return;
            }
        }
        {
            typedef void (FullPRISMCalcThread::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&FullPRISMCalcThread::outputCalculated)) {
                *result = 2;
                return;
            }
        }
        {
            typedef void (FullPRISMCalcThread::*_t)(const QString );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&FullPRISMCalcThread::signalTitle)) {
                *result = 3;
                return;
            }
        }
    }
}

const QMetaObject FullPRISMCalcThread::staticMetaObject = {
    { &QThread::staticMetaObject, qt_meta_stringdata_FullPRISMCalcThread.data,
      qt_meta_data_FullPRISMCalcThread,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *FullPRISMCalcThread::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *FullPRISMCalcThread::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_FullPRISMCalcThread.stringdata0))
        return static_cast<void*>(const_cast< FullPRISMCalcThread*>(this));
    return QThread::qt_metacast(_clname);
}

int FullPRISMCalcThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QThread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 4)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 4;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 4)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 4;
    }
    return _id;
}

// SIGNAL 0
void FullPRISMCalcThread::potentialCalculated()
{
    QMetaObject::activate(this, &staticMetaObject, 0, nullptr);
}

// SIGNAL 1
void FullPRISMCalcThread::ScompactCalculated()
{
    QMetaObject::activate(this, &staticMetaObject, 1, nullptr);
}

// SIGNAL 2
void FullPRISMCalcThread::outputCalculated()
{
    QMetaObject::activate(this, &staticMetaObject, 2, nullptr);
}

// SIGNAL 3
void FullPRISMCalcThread::signalTitle(const QString _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}
struct qt_meta_stringdata_FullMultisliceCalcThread_t {
    QByteArrayData data[7];
    char stringdata0[98];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_FullMultisliceCalcThread_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_FullMultisliceCalcThread_t qt_meta_stringdata_FullMultisliceCalcThread = {
    {
QT_MOC_LITERAL(0, 0, 24), // "FullMultisliceCalcThread"
QT_MOC_LITERAL(1, 25, 19), // "potentialCalculated"
QT_MOC_LITERAL(2, 45, 0), // ""
QT_MOC_LITERAL(3, 46, 18), // "ScompactCalculated"
QT_MOC_LITERAL(4, 65, 16), // "outputCalculated"
QT_MOC_LITERAL(5, 82, 11), // "signalTitle"
QT_MOC_LITERAL(6, 94, 3) // "str"

    },
    "FullMultisliceCalcThread\0potentialCalculated\0"
    "\0ScompactCalculated\0outputCalculated\0"
    "signalTitle\0str"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_FullMultisliceCalcThread[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       4,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    0,   34,    2, 0x06 /* Public */,
       3,    0,   35,    2, 0x06 /* Public */,
       4,    0,   36,    2, 0x06 /* Public */,
       5,    1,   37,    2, 0x06 /* Public */,

 // signals: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QString,    6,

       0        // eod
};

void FullMultisliceCalcThread::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        FullMultisliceCalcThread *_t = static_cast<FullMultisliceCalcThread *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->potentialCalculated(); break;
        case 1: _t->ScompactCalculated(); break;
        case 2: _t->outputCalculated(); break;
        case 3: _t->signalTitle((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (FullMultisliceCalcThread::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&FullMultisliceCalcThread::potentialCalculated)) {
                *result = 0;
                return;
            }
        }
        {
            typedef void (FullMultisliceCalcThread::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&FullMultisliceCalcThread::ScompactCalculated)) {
                *result = 1;
                return;
            }
        }
        {
            typedef void (FullMultisliceCalcThread::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&FullMultisliceCalcThread::outputCalculated)) {
                *result = 2;
                return;
            }
        }
        {
            typedef void (FullMultisliceCalcThread::*_t)(const QString );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&FullMultisliceCalcThread::signalTitle)) {
                *result = 3;
                return;
            }
        }
    }
}

const QMetaObject FullMultisliceCalcThread::staticMetaObject = {
    { &QThread::staticMetaObject, qt_meta_stringdata_FullMultisliceCalcThread.data,
      qt_meta_data_FullMultisliceCalcThread,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *FullMultisliceCalcThread::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *FullMultisliceCalcThread::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_FullMultisliceCalcThread.stringdata0))
        return static_cast<void*>(const_cast< FullMultisliceCalcThread*>(this));
    return QThread::qt_metacast(_clname);
}

int FullMultisliceCalcThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QThread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 4)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 4;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 4)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 4;
    }
    return _id;
}

// SIGNAL 0
void FullMultisliceCalcThread::potentialCalculated()
{
    QMetaObject::activate(this, &staticMetaObject, 0, nullptr);
}

// SIGNAL 1
void FullMultisliceCalcThread::ScompactCalculated()
{
    QMetaObject::activate(this, &staticMetaObject, 1, nullptr);
}

// SIGNAL 2
void FullMultisliceCalcThread::outputCalculated()
{
    QMetaObject::activate(this, &staticMetaObject, 2, nullptr);
}

// SIGNAL 3
void FullMultisliceCalcThread::signalTitle(const QString _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
