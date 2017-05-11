/****************************************************************************
** Meta object code from reading C++ file 'prismmainwindow.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.9.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "prismmainwindow.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'prismmainwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.9.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_PRISMMainWindow_t {
    QByteArrayData data[87];
    char stringdata0[1829];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_PRISMMainWindow_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_PRISMMainWindow_t qt_meta_stringdata_PRISMMainWindow = {
    {
QT_MOC_LITERAL(0, 0, 15), // "PRISMMainWindow"
QT_MOC_LITERAL(1, 16, 23), // "setInterpolationFactorX"
QT_MOC_LITERAL(2, 40, 0), // ""
QT_MOC_LITERAL(3, 41, 23), // "setInterpolationFactorY"
QT_MOC_LITERAL(4, 65, 27), // "setFilenameAtoms_fromDialog"
QT_MOC_LITERAL(5, 93, 30), // "setFilenameOutput_fromLineEdit"
QT_MOC_LITERAL(6, 124, 28), // "setFilenameOutput_fromDialog"
QT_MOC_LITERAL(7, 153, 10), // "setNumGPUs"
QT_MOC_LITERAL(8, 164, 7), // "numGPUs"
QT_MOC_LITERAL(9, 172, 13), // "setNumThreads"
QT_MOC_LITERAL(10, 186, 10), // "numThreads"
QT_MOC_LITERAL(11, 197, 8), // "setNumFP"
QT_MOC_LITERAL(12, 206, 5), // "numFP"
QT_MOC_LITERAL(13, 212, 25), // "setPixelSize_fromLineEdit"
QT_MOC_LITERAL(14, 238, 24), // "setPotBound_fromLineEdit"
QT_MOC_LITERAL(15, 263, 30), // "setprobeSemiangle_fromLineEdit"
QT_MOC_LITERAL(16, 294, 30), // "setSliceThickness_fromLineEdit"
QT_MOC_LITERAL(17, 325, 24), // "setCellDimX_fromLineEdit"
QT_MOC_LITERAL(18, 350, 24), // "setCellDimY_fromLineEdit"
QT_MOC_LITERAL(19, 375, 24), // "setCellDimZ_fromLineEdit"
QT_MOC_LITERAL(20, 400, 21), // "setTileX_fromLineEdit"
QT_MOC_LITERAL(21, 422, 21), // "setTileY_fromLineEdit"
QT_MOC_LITERAL(22, 444, 21), // "setTileZ_fromLineEdit"
QT_MOC_LITERAL(23, 466, 18), // "setE0_fromLineEdit"
QT_MOC_LITERAL(24, 485, 27), // "setprobe_stepX_fromLineEdit"
QT_MOC_LITERAL(25, 513, 27), // "setprobe_stepY_fromLineEdit"
QT_MOC_LITERAL(26, 541, 13), // "setAlgo_PRISM"
QT_MOC_LITERAL(27, 555, 18), // "setAlgo_Multislice"
QT_MOC_LITERAL(28, 574, 18), // "calculatePotential"
QT_MOC_LITERAL(29, 593, 16), // "calculateSMatrix"
QT_MOC_LITERAL(30, 610, 12), // "calculateAll"
QT_MOC_LITERAL(31, 623, 14), // "calculateProbe"
QT_MOC_LITERAL(32, 638, 20), // "updatePotentialImage"
QT_MOC_LITERAL(33, 659, 22), // "updatePotentialDisplay"
QT_MOC_LITERAL(34, 682, 25), // "updatePotentialFloatImage"
QT_MOC_LITERAL(35, 708, 17), // "updateOutputImage"
QT_MOC_LITERAL(36, 726, 19), // "updateOutputDisplay"
QT_MOC_LITERAL(37, 746, 22), // "updateOutputFloatImage"
QT_MOC_LITERAL(38, 769, 27), // "updateSliders_fromLineEdits"
QT_MOC_LITERAL(39, 797, 31), // "updateSliders_fromLineEdits_ang"
QT_MOC_LITERAL(40, 829, 20), // "updateContrastPotMin"
QT_MOC_LITERAL(41, 850, 20), // "updateContrastPotMax"
QT_MOC_LITERAL(42, 871, 20), // "updateContrastAngMin"
QT_MOC_LITERAL(43, 892, 20), // "updateContrastAngMax"
QT_MOC_LITERAL(44, 913, 26), // "updateSlider_lineEdits_min"
QT_MOC_LITERAL(45, 940, 26), // "updateSlider_lineEdits_max"
QT_MOC_LITERAL(46, 967, 30), // "updateSlider_lineEdits_max_ang"
QT_MOC_LITERAL(47, 998, 3), // "val"
QT_MOC_LITERAL(48, 1002, 30), // "updateSlider_lineEdits_min_ang"
QT_MOC_LITERAL(49, 1033, 14), // "updateAlphaMax"
QT_MOC_LITERAL(50, 1048, 11), // "resizeEvent"
QT_MOC_LITERAL(51, 1060, 13), // "QResizeEvent*"
QT_MOC_LITERAL(52, 1074, 5), // "event"
QT_MOC_LITERAL(53, 1080, 12), // "redrawImages"
QT_MOC_LITERAL(54, 1093, 22), // "saveCurrentOutputImage"
QT_MOC_LITERAL(55, 1116, 19), // "toggleStreamingMode"
QT_MOC_LITERAL(56, 1136, 28), // "toggleSaveProjectedPotential"
QT_MOC_LITERAL(57, 1165, 19), // "enableOutputWidgets"
QT_MOC_LITERAL(58, 1185, 29), // "setprobe_defocus_fromLineEdit"
QT_MOC_LITERAL(59, 1215, 26), // "setRandomSeed_fromLineEdit"
QT_MOC_LITERAL(60, 1242, 24), // "setprobe_C3_fromLineEdit"
QT_MOC_LITERAL(61, 1267, 24), // "setprobe_C5_fromLineEdit"
QT_MOC_LITERAL(62, 1292, 35), // "setDetector_angle_step_fromLi..."
QT_MOC_LITERAL(63, 1328, 27), // "setprobe_Xtilt_fromLineEdit"
QT_MOC_LITERAL(64, 1356, 27), // "setprobe_Ytilt_fromLineEdit"
QT_MOC_LITERAL(65, 1384, 14), // "toggle3DOutput"
QT_MOC_LITERAL(66, 1399, 14), // "toggle4DOutput"
QT_MOC_LITERAL(67, 1414, 20), // "toggleThermalEffects"
QT_MOC_LITERAL(68, 1435, 31), // "setscan_WindowXMin_fromLineEdit"
QT_MOC_LITERAL(69, 1467, 31), // "setscan_WindowXMax_fromLineEdit"
QT_MOC_LITERAL(70, 1499, 31), // "setscan_WindowYMin_fromLineEdit"
QT_MOC_LITERAL(71, 1531, 31), // "setscan_WindowYMax_fromLineEdit"
QT_MOC_LITERAL(72, 1563, 16), // "resetCalculation"
QT_MOC_LITERAL(73, 1580, 13), // "newRandomSeed"
QT_MOC_LITERAL(74, 1594, 18), // "updateProbeK_PRISM"
QT_MOC_LITERAL(75, 1613, 37), // "PRISM::Array2D<PRISM_FLOAT_PR..."
QT_MOC_LITERAL(76, 1651, 18), // "updateProbeR_PRISM"
QT_MOC_LITERAL(77, 1670, 23), // "updateProbeK_Multislice"
QT_MOC_LITERAL(78, 1694, 23), // "updateProbeR_Multislice"
QT_MOC_LITERAL(79, 1718, 17), // "updateProbe_diffR"
QT_MOC_LITERAL(80, 1736, 12), // "arr_contrast"
QT_MOC_LITERAL(81, 1749, 17), // "updateProbe_diffK"
QT_MOC_LITERAL(82, 1767, 18), // "update_pearsonReal"
QT_MOC_LITERAL(83, 1786, 3), // "str"
QT_MOC_LITERAL(84, 1790, 15), // "update_pearsonK"
QT_MOC_LITERAL(85, 1806, 12), // "update_RReal"
QT_MOC_LITERAL(86, 1819, 9) // "update_RK"

    },
    "PRISMMainWindow\0setInterpolationFactorX\0"
    "\0setInterpolationFactorY\0"
    "setFilenameAtoms_fromDialog\0"
    "setFilenameOutput_fromLineEdit\0"
    "setFilenameOutput_fromDialog\0setNumGPUs\0"
    "numGPUs\0setNumThreads\0numThreads\0"
    "setNumFP\0numFP\0setPixelSize_fromLineEdit\0"
    "setPotBound_fromLineEdit\0"
    "setprobeSemiangle_fromLineEdit\0"
    "setSliceThickness_fromLineEdit\0"
    "setCellDimX_fromLineEdit\0"
    "setCellDimY_fromLineEdit\0"
    "setCellDimZ_fromLineEdit\0setTileX_fromLineEdit\0"
    "setTileY_fromLineEdit\0setTileZ_fromLineEdit\0"
    "setE0_fromLineEdit\0setprobe_stepX_fromLineEdit\0"
    "setprobe_stepY_fromLineEdit\0setAlgo_PRISM\0"
    "setAlgo_Multislice\0calculatePotential\0"
    "calculateSMatrix\0calculateAll\0"
    "calculateProbe\0updatePotentialImage\0"
    "updatePotentialDisplay\0updatePotentialFloatImage\0"
    "updateOutputImage\0updateOutputDisplay\0"
    "updateOutputFloatImage\0"
    "updateSliders_fromLineEdits\0"
    "updateSliders_fromLineEdits_ang\0"
    "updateContrastPotMin\0updateContrastPotMax\0"
    "updateContrastAngMin\0updateContrastAngMax\0"
    "updateSlider_lineEdits_min\0"
    "updateSlider_lineEdits_max\0"
    "updateSlider_lineEdits_max_ang\0val\0"
    "updateSlider_lineEdits_min_ang\0"
    "updateAlphaMax\0resizeEvent\0QResizeEvent*\0"
    "event\0redrawImages\0saveCurrentOutputImage\0"
    "toggleStreamingMode\0toggleSaveProjectedPotential\0"
    "enableOutputWidgets\0setprobe_defocus_fromLineEdit\0"
    "setRandomSeed_fromLineEdit\0"
    "setprobe_C3_fromLineEdit\0"
    "setprobe_C5_fromLineEdit\0"
    "setDetector_angle_step_fromLineEdit\0"
    "setprobe_Xtilt_fromLineEdit\0"
    "setprobe_Ytilt_fromLineEdit\0toggle3DOutput\0"
    "toggle4DOutput\0toggleThermalEffects\0"
    "setscan_WindowXMin_fromLineEdit\0"
    "setscan_WindowXMax_fromLineEdit\0"
    "setscan_WindowYMin_fromLineEdit\0"
    "setscan_WindowYMax_fromLineEdit\0"
    "resetCalculation\0newRandomSeed\0"
    "updateProbeK_PRISM\0"
    "PRISM::Array2D<PRISM_FLOAT_PRECISION>\0"
    "updateProbeR_PRISM\0updateProbeK_Multislice\0"
    "updateProbeR_Multislice\0updateProbe_diffR\0"
    "arr_contrast\0updateProbe_diffK\0"
    "update_pearsonReal\0str\0update_pearsonK\0"
    "update_RReal\0update_RK"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_PRISMMainWindow[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      76,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,  394,    2, 0x0a /* Public */,
       3,    0,  395,    2, 0x0a /* Public */,
       4,    0,  396,    2, 0x0a /* Public */,
       5,    0,  397,    2, 0x0a /* Public */,
       6,    0,  398,    2, 0x0a /* Public */,
       7,    1,  399,    2, 0x0a /* Public */,
       9,    1,  402,    2, 0x0a /* Public */,
      11,    1,  405,    2, 0x0a /* Public */,
      13,    0,  408,    2, 0x0a /* Public */,
      14,    0,  409,    2, 0x0a /* Public */,
      15,    0,  410,    2, 0x0a /* Public */,
      16,    0,  411,    2, 0x0a /* Public */,
      17,    0,  412,    2, 0x0a /* Public */,
      18,    0,  413,    2, 0x0a /* Public */,
      19,    0,  414,    2, 0x0a /* Public */,
      20,    0,  415,    2, 0x0a /* Public */,
      21,    0,  416,    2, 0x0a /* Public */,
      22,    0,  417,    2, 0x0a /* Public */,
      23,    0,  418,    2, 0x0a /* Public */,
      24,    0,  419,    2, 0x0a /* Public */,
      25,    0,  420,    2, 0x0a /* Public */,
      26,    0,  421,    2, 0x0a /* Public */,
      27,    0,  422,    2, 0x0a /* Public */,
      28,    0,  423,    2, 0x0a /* Public */,
      29,    0,  424,    2, 0x0a /* Public */,
      30,    0,  425,    2, 0x0a /* Public */,
      31,    0,  426,    2, 0x0a /* Public */,
      32,    0,  427,    2, 0x0a /* Public */,
      33,    0,  428,    2, 0x0a /* Public */,
      34,    0,  429,    2, 0x0a /* Public */,
      35,    0,  430,    2, 0x0a /* Public */,
      36,    0,  431,    2, 0x0a /* Public */,
      37,    0,  432,    2, 0x0a /* Public */,
      38,    0,  433,    2, 0x0a /* Public */,
      39,    0,  434,    2, 0x0a /* Public */,
      40,    0,  435,    2, 0x0a /* Public */,
      41,    0,  436,    2, 0x0a /* Public */,
      42,    0,  437,    2, 0x0a /* Public */,
      43,    0,  438,    2, 0x0a /* Public */,
      44,    1,  439,    2, 0x0a /* Public */,
      45,    1,  442,    2, 0x0a /* Public */,
      46,    1,  445,    2, 0x0a /* Public */,
      48,    1,  448,    2, 0x0a /* Public */,
      49,    0,  451,    2, 0x0a /* Public */,
      50,    1,  452,    2, 0x0a /* Public */,
      53,    0,  455,    2, 0x0a /* Public */,
      54,    0,  456,    2, 0x0a /* Public */,
      55,    0,  457,    2, 0x0a /* Public */,
      56,    0,  458,    2, 0x0a /* Public */,
      57,    0,  459,    2, 0x0a /* Public */,
      58,    0,  460,    2, 0x0a /* Public */,
      59,    0,  461,    2, 0x0a /* Public */,
      60,    0,  462,    2, 0x0a /* Public */,
      61,    0,  463,    2, 0x0a /* Public */,
      62,    0,  464,    2, 0x0a /* Public */,
      63,    0,  465,    2, 0x0a /* Public */,
      64,    0,  466,    2, 0x0a /* Public */,
      65,    0,  467,    2, 0x0a /* Public */,
      66,    0,  468,    2, 0x0a /* Public */,
      67,    0,  469,    2, 0x0a /* Public */,
      68,    0,  470,    2, 0x0a /* Public */,
      69,    0,  471,    2, 0x0a /* Public */,
      70,    0,  472,    2, 0x0a /* Public */,
      71,    0,  473,    2, 0x0a /* Public */,
      72,    0,  474,    2, 0x0a /* Public */,
      73,    0,  475,    2, 0x0a /* Public */,
      74,    1,  476,    2, 0x0a /* Public */,
      76,    1,  479,    2, 0x0a /* Public */,
      77,    1,  482,    2, 0x0a /* Public */,
      78,    1,  485,    2, 0x0a /* Public */,
      79,    2,  488,    2, 0x0a /* Public */,
      81,    2,  493,    2, 0x0a /* Public */,
      82,    1,  498,    2, 0x0a /* Public */,
      84,    1,  501,    2, 0x0a /* Public */,
      85,    1,  504,    2, 0x0a /* Public */,
      86,    1,  507,    2, 0x0a /* Public */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,    8,
    QMetaType::Void, QMetaType::Int,   10,
    QMetaType::Void, QMetaType::Int,   12,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,   47,
    QMetaType::Void, QMetaType::Int,   47,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 51,   52,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 75,    2,
    QMetaType::Void, 0x80000000 | 75,    2,
    QMetaType::Void, 0x80000000 | 75,    2,
    QMetaType::Void, 0x80000000 | 75,    2,
    QMetaType::Void, 0x80000000 | 75, 0x80000000 | 75,    2,   80,
    QMetaType::Void, 0x80000000 | 75, 0x80000000 | 75,    2,   80,
    QMetaType::Void, QMetaType::QString,   83,
    QMetaType::Void, QMetaType::QString,   83,
    QMetaType::Void, QMetaType::QString,   83,
    QMetaType::Void, QMetaType::QString,   83,

       0        // eod
};

void PRISMMainWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        PRISMMainWindow *_t = static_cast<PRISMMainWindow *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->setInterpolationFactorX(); break;
        case 1: _t->setInterpolationFactorY(); break;
        case 2: _t->setFilenameAtoms_fromDialog(); break;
        case 3: _t->setFilenameOutput_fromLineEdit(); break;
        case 4: _t->setFilenameOutput_fromDialog(); break;
        case 5: _t->setNumGPUs((*reinterpret_cast< const int(*)>(_a[1]))); break;
        case 6: _t->setNumThreads((*reinterpret_cast< const int(*)>(_a[1]))); break;
        case 7: _t->setNumFP((*reinterpret_cast< const int(*)>(_a[1]))); break;
        case 8: _t->setPixelSize_fromLineEdit(); break;
        case 9: _t->setPotBound_fromLineEdit(); break;
        case 10: _t->setprobeSemiangle_fromLineEdit(); break;
        case 11: _t->setSliceThickness_fromLineEdit(); break;
        case 12: _t->setCellDimX_fromLineEdit(); break;
        case 13: _t->setCellDimY_fromLineEdit(); break;
        case 14: _t->setCellDimZ_fromLineEdit(); break;
        case 15: _t->setTileX_fromLineEdit(); break;
        case 16: _t->setTileY_fromLineEdit(); break;
        case 17: _t->setTileZ_fromLineEdit(); break;
        case 18: _t->setE0_fromLineEdit(); break;
        case 19: _t->setprobe_stepX_fromLineEdit(); break;
        case 20: _t->setprobe_stepY_fromLineEdit(); break;
        case 21: _t->setAlgo_PRISM(); break;
        case 22: _t->setAlgo_Multislice(); break;
        case 23: _t->calculatePotential(); break;
        case 24: _t->calculateSMatrix(); break;
        case 25: _t->calculateAll(); break;
        case 26: _t->calculateProbe(); break;
        case 27: _t->updatePotentialImage(); break;
        case 28: _t->updatePotentialDisplay(); break;
        case 29: _t->updatePotentialFloatImage(); break;
        case 30: _t->updateOutputImage(); break;
        case 31: _t->updateOutputDisplay(); break;
        case 32: _t->updateOutputFloatImage(); break;
        case 33: _t->updateSliders_fromLineEdits(); break;
        case 34: _t->updateSliders_fromLineEdits_ang(); break;
        case 35: _t->updateContrastPotMin(); break;
        case 36: _t->updateContrastPotMax(); break;
        case 37: _t->updateContrastAngMin(); break;
        case 38: _t->updateContrastAngMax(); break;
        case 39: _t->updateSlider_lineEdits_min((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 40: _t->updateSlider_lineEdits_max((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 41: _t->updateSlider_lineEdits_max_ang((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 42: _t->updateSlider_lineEdits_min_ang((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 43: _t->updateAlphaMax(); break;
        case 44: _t->resizeEvent((*reinterpret_cast< QResizeEvent*(*)>(_a[1]))); break;
        case 45: _t->redrawImages(); break;
        case 46: _t->saveCurrentOutputImage(); break;
        case 47: _t->toggleStreamingMode(); break;
        case 48: _t->toggleSaveProjectedPotential(); break;
        case 49: _t->enableOutputWidgets(); break;
        case 50: _t->setprobe_defocus_fromLineEdit(); break;
        case 51: _t->setRandomSeed_fromLineEdit(); break;
        case 52: _t->setprobe_C3_fromLineEdit(); break;
        case 53: _t->setprobe_C5_fromLineEdit(); break;
        case 54: _t->setDetector_angle_step_fromLineEdit(); break;
        case 55: _t->setprobe_Xtilt_fromLineEdit(); break;
        case 56: _t->setprobe_Ytilt_fromLineEdit(); break;
        case 57: _t->toggle3DOutput(); break;
        case 58: _t->toggle4DOutput(); break;
        case 59: _t->toggleThermalEffects(); break;
        case 60: _t->setscan_WindowXMin_fromLineEdit(); break;
        case 61: _t->setscan_WindowXMax_fromLineEdit(); break;
        case 62: _t->setscan_WindowYMin_fromLineEdit(); break;
        case 63: _t->setscan_WindowYMax_fromLineEdit(); break;
        case 64: _t->resetCalculation(); break;
        case 65: _t->newRandomSeed(); break;
        case 66: _t->updateProbeK_PRISM((*reinterpret_cast< PRISM::Array2D<PRISM_FLOAT_PRECISION>(*)>(_a[1]))); break;
        case 67: _t->updateProbeR_PRISM((*reinterpret_cast< PRISM::Array2D<PRISM_FLOAT_PRECISION>(*)>(_a[1]))); break;
        case 68: _t->updateProbeK_Multislice((*reinterpret_cast< PRISM::Array2D<PRISM_FLOAT_PRECISION>(*)>(_a[1]))); break;
        case 69: _t->updateProbeR_Multislice((*reinterpret_cast< PRISM::Array2D<PRISM_FLOAT_PRECISION>(*)>(_a[1]))); break;
        case 70: _t->updateProbe_diffR((*reinterpret_cast< PRISM::Array2D<PRISM_FLOAT_PRECISION>(*)>(_a[1])),(*reinterpret_cast< PRISM::Array2D<PRISM_FLOAT_PRECISION>(*)>(_a[2]))); break;
        case 71: _t->updateProbe_diffK((*reinterpret_cast< PRISM::Array2D<PRISM_FLOAT_PRECISION>(*)>(_a[1])),(*reinterpret_cast< PRISM::Array2D<PRISM_FLOAT_PRECISION>(*)>(_a[2]))); break;
        case 72: _t->update_pearsonReal((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 73: _t->update_pearsonK((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 74: _t->update_RReal((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 75: _t->update_RK((*reinterpret_cast< QString(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObject PRISMMainWindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_PRISMMainWindow.data,
      qt_meta_data_PRISMMainWindow,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *PRISMMainWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *PRISMMainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
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
        if (_id < 76)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 76;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 76)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 76;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
