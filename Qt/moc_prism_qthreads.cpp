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
struct qt_meta_stringdata_PRISMThread_t {
    QByteArrayData data[5];
    char stringdata0[80];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_PRISMThread_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_PRISMThread_t qt_meta_stringdata_PRISMThread = {
    {
QT_MOC_LITERAL(0, 0, 11), // "PRISMThread"
QT_MOC_LITERAL(1, 12, 19), // "potentialCalculated"
QT_MOC_LITERAL(2, 32, 0), // ""
QT_MOC_LITERAL(3, 33, 29), // "signalErrorReadingAtomsDialog"
QT_MOC_LITERAL(4, 63, 16) // "outputCalculated"

    },
    "PRISMThread\0potentialCalculated\0\0"
    "signalErrorReadingAtomsDialog\0"
    "outputCalculated"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_PRISMThread[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       3,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    0,   29,    2, 0x06 /* Public */,
       3,    0,   30,    2, 0x06 /* Public */,
       4,    0,   31,    2, 0x06 /* Public */,

 // signals: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void PRISMThread::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        PRISMThread *_t = static_cast<PRISMThread *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->potentialCalculated(); break;
        case 1: _t->signalErrorReadingAtomsDialog(); break;
        case 2: _t->outputCalculated(); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (PRISMThread::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&PRISMThread::potentialCalculated)) {
                *result = 0;
                return;
            }
        }
        {
            typedef void (PRISMThread::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&PRISMThread::signalErrorReadingAtomsDialog)) {
                *result = 1;
                return;
            }
        }
        {
            typedef void (PRISMThread::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&PRISMThread::outputCalculated)) {
                *result = 2;
                return;
            }
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject PRISMThread::staticMetaObject = {
    { &QThread::staticMetaObject, qt_meta_stringdata_PRISMThread.data,
      qt_meta_data_PRISMThread,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *PRISMThread::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *PRISMThread::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_PRISMThread.stringdata0))
        return static_cast<void*>(const_cast< PRISMThread*>(this));
    return QThread::qt_metacast(_clname);
}

int PRISMThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QThread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 3)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 3;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 3)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 3;
    }
    return _id;
}

// SIGNAL 0
void PRISMThread::potentialCalculated()
{
    QMetaObject::activate(this, &staticMetaObject, 0, nullptr);
}

// SIGNAL 1
void PRISMThread::signalErrorReadingAtomsDialog()
{
    QMetaObject::activate(this, &staticMetaObject, 1, nullptr);
}

// SIGNAL 2
void PRISMThread::outputCalculated()
{
    QMetaObject::activate(this, &staticMetaObject, 2, nullptr);
}
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
    { &PRISMThread::staticMetaObject, qt_meta_stringdata_PotentialThread.data,
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
    return PRISMThread::qt_metacast(_clname);
}

int PotentialThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = PRISMThread::qt_metacall(_c, _id, _a);
    return _id;
}
struct qt_meta_stringdata_ProbeThread_t {
    QByteArrayData data[14];
    char stringdata0[243];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_ProbeThread_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_ProbeThread_t qt_meta_stringdata_ProbeThread = {
    {
QT_MOC_LITERAL(0, 0, 11), // "ProbeThread"
QT_MOC_LITERAL(1, 12, 18), // "signalProbeK_PRISM"
QT_MOC_LITERAL(2, 31, 0), // ""
QT_MOC_LITERAL(3, 32, 45), // "Prismatic::Array2D<PRISMATIC_..."
QT_MOC_LITERAL(4, 78, 18), // "signalProbeR_PRISM"
QT_MOC_LITERAL(5, 97, 23), // "signalProbeK_Multislice"
QT_MOC_LITERAL(6, 121, 23), // "signalProbeR_Multislice"
QT_MOC_LITERAL(7, 145, 17), // "signalProbe_diffK"
QT_MOC_LITERAL(8, 163, 17), // "signalProbe_diffR"
QT_MOC_LITERAL(9, 181, 18), // "signal_pearsonReal"
QT_MOC_LITERAL(10, 200, 3), // "str"
QT_MOC_LITERAL(11, 204, 15), // "signal_pearsonK"
QT_MOC_LITERAL(12, 220, 12), // "signal_RReal"
QT_MOC_LITERAL(13, 233, 9) // "signal_RK"

    },
    "ProbeThread\0signalProbeK_PRISM\0\0"
    "Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>\0"
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
      10,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
      10,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   64,    2, 0x06 /* Public */,
       4,    1,   67,    2, 0x06 /* Public */,
       5,    1,   70,    2, 0x06 /* Public */,
       6,    1,   73,    2, 0x06 /* Public */,
       7,    2,   76,    2, 0x06 /* Public */,
       8,    2,   81,    2, 0x06 /* Public */,
       9,    1,   86,    2, 0x06 /* Public */,
      11,    1,   89,    2, 0x06 /* Public */,
      12,    1,   92,    2, 0x06 /* Public */,
      13,    1,   95,    2, 0x06 /* Public */,

 // signals: parameters
    QMetaType::Void, 0x80000000 | 3,    2,
    QMetaType::Void, 0x80000000 | 3,    2,
    QMetaType::Void, 0x80000000 | 3,    2,
    QMetaType::Void, 0x80000000 | 3,    2,
    QMetaType::Void, 0x80000000 | 3, 0x80000000 | 3,    2,    2,
    QMetaType::Void, 0x80000000 | 3, 0x80000000 | 3,    2,    2,
    QMetaType::Void, QMetaType::QString,   10,
    QMetaType::Void, QMetaType::QString,   10,
    QMetaType::Void, QMetaType::QString,   10,
    QMetaType::Void, QMetaType::QString,   10,

       0        // eod
};

void ProbeThread::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        ProbeThread *_t = static_cast<ProbeThread *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->signalProbeK_PRISM((*reinterpret_cast< Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>(*)>(_a[1]))); break;
        case 1: _t->signalProbeR_PRISM((*reinterpret_cast< Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>(*)>(_a[1]))); break;
        case 2: _t->signalProbeK_Multislice((*reinterpret_cast< Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>(*)>(_a[1]))); break;
        case 3: _t->signalProbeR_Multislice((*reinterpret_cast< Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>(*)>(_a[1]))); break;
        case 4: _t->signalProbe_diffK((*reinterpret_cast< Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>(*)>(_a[1])),(*reinterpret_cast< Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>(*)>(_a[2]))); break;
        case 5: _t->signalProbe_diffR((*reinterpret_cast< Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>(*)>(_a[1])),(*reinterpret_cast< Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>(*)>(_a[2]))); break;
        case 6: _t->signal_pearsonReal((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 7: _t->signal_pearsonK((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 8: _t->signal_RReal((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 9: _t->signal_RK((*reinterpret_cast< QString(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (ProbeThread::*_t)(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&ProbeThread::signalProbeK_PRISM)) {
                *result = 0;
                return;
            }
        }
        {
            typedef void (ProbeThread::*_t)(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&ProbeThread::signalProbeR_PRISM)) {
                *result = 1;
                return;
            }
        }
        {
            typedef void (ProbeThread::*_t)(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&ProbeThread::signalProbeK_Multislice)) {
                *result = 2;
                return;
            }
        }
        {
            typedef void (ProbeThread::*_t)(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&ProbeThread::signalProbeR_Multislice)) {
                *result = 3;
                return;
            }
        }
        {
            typedef void (ProbeThread::*_t)(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> , Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&ProbeThread::signalProbe_diffK)) {
                *result = 4;
                return;
            }
        }
        {
            typedef void (ProbeThread::*_t)(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> , Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&ProbeThread::signalProbe_diffR)) {
                *result = 5;
                return;
            }
        }
        {
            typedef void (ProbeThread::*_t)(QString );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&ProbeThread::signal_pearsonReal)) {
                *result = 6;
                return;
            }
        }
        {
            typedef void (ProbeThread::*_t)(QString );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&ProbeThread::signal_pearsonK)) {
                *result = 7;
                return;
            }
        }
        {
            typedef void (ProbeThread::*_t)(QString );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&ProbeThread::signal_RReal)) {
                *result = 8;
                return;
            }
        }
        {
            typedef void (ProbeThread::*_t)(QString );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&ProbeThread::signal_RK)) {
                *result = 9;
                return;
            }
        }
    }
}

const QMetaObject ProbeThread::staticMetaObject = {
    { &PRISMThread::staticMetaObject, qt_meta_stringdata_ProbeThread.data,
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
    return PRISMThread::qt_metacast(_clname);
}

int ProbeThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = PRISMThread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 10)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 10;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 10)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 10;
    }
    return _id;
}

// SIGNAL 0
void ProbeThread::signalProbeK_PRISM(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void ProbeThread::signalProbeR_PRISM(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void ProbeThread::signalProbeK_Multislice(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void ProbeThread::signalProbeR_Multislice(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}

// SIGNAL 4
void ProbeThread::signalProbe_diffK(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> _t1, Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> _t2)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 4, _a);
}

// SIGNAL 5
void ProbeThread::signalProbe_diffR(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> _t1, Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> _t2)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 5, _a);
}

// SIGNAL 6
void ProbeThread::signal_pearsonReal(QString _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 6, _a);
}

// SIGNAL 7
void ProbeThread::signal_pearsonK(QString _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 7, _a);
}

// SIGNAL 8
void ProbeThread::signal_RReal(QString _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 8, _a);
}

// SIGNAL 9
void ProbeThread::signal_RK(QString _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 9, _a);
}
struct qt_meta_stringdata_FullPRISMCalcThread_t {
    QByteArrayData data[4];
    char stringdata0[37];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_FullPRISMCalcThread_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_FullPRISMCalcThread_t qt_meta_stringdata_FullPRISMCalcThread = {
    {
QT_MOC_LITERAL(0, 0, 19), // "FullPRISMCalcThread"
QT_MOC_LITERAL(1, 20, 11), // "signalTitle"
QT_MOC_LITERAL(2, 32, 0), // ""
QT_MOC_LITERAL(3, 33, 3) // "str"

    },
    "FullPRISMCalcThread\0signalTitle\0\0str"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_FullPRISMCalcThread[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   19,    2, 0x06 /* Public */,

 // signals: parameters
    QMetaType::Void, QMetaType::QString,    3,

       0        // eod
};

void FullPRISMCalcThread::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        FullPRISMCalcThread *_t = static_cast<FullPRISMCalcThread *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->signalTitle((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (FullPRISMCalcThread::*_t)(const QString );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&FullPRISMCalcThread::signalTitle)) {
                *result = 0;
                return;
            }
        }
    }
}

const QMetaObject FullPRISMCalcThread::staticMetaObject = {
    { &PRISMThread::staticMetaObject, qt_meta_stringdata_FullPRISMCalcThread.data,
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
    return PRISMThread::qt_metacast(_clname);
}

int FullPRISMCalcThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = PRISMThread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 1)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 1;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 1)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 1;
    }
    return _id;
}

// SIGNAL 0
void FullPRISMCalcThread::signalTitle(const QString _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
struct qt_meta_stringdata_FullMultisliceCalcThread_t {
    QByteArrayData data[4];
    char stringdata0[42];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_FullMultisliceCalcThread_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_FullMultisliceCalcThread_t qt_meta_stringdata_FullMultisliceCalcThread = {
    {
QT_MOC_LITERAL(0, 0, 24), // "FullMultisliceCalcThread"
QT_MOC_LITERAL(1, 25, 11), // "signalTitle"
QT_MOC_LITERAL(2, 37, 0), // ""
QT_MOC_LITERAL(3, 38, 3) // "str"

    },
    "FullMultisliceCalcThread\0signalTitle\0"
    "\0str"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_FullMultisliceCalcThread[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   19,    2, 0x06 /* Public */,

 // signals: parameters
    QMetaType::Void, QMetaType::QString,    3,

       0        // eod
};

void FullMultisliceCalcThread::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        FullMultisliceCalcThread *_t = static_cast<FullMultisliceCalcThread *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->signalTitle((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (FullMultisliceCalcThread::*_t)(const QString );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&FullMultisliceCalcThread::signalTitle)) {
                *result = 0;
                return;
            }
        }
    }
}

const QMetaObject FullMultisliceCalcThread::staticMetaObject = {
    { &PRISMThread::staticMetaObject, qt_meta_stringdata_FullMultisliceCalcThread.data,
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
    return PRISMThread::qt_metacast(_clname);
}

int FullMultisliceCalcThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = PRISMThread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 1)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 1;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 1)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 1;
    }
    return _id;
}

// SIGNAL 0
void FullMultisliceCalcThread::signalTitle(const QString _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
