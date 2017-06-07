/****************************************************************************
** Meta object code from reading C++ file 'prismmainwindow.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.8.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "prismmainwindow.h"
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
    QByteArrayData data[116];
    char stringdata0[2615];
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
QT_MOC_LITERAL(13, 212, 13), // "setNumStreams"
QT_MOC_LITERAL(14, 226, 26), // "setPixelSizeX_fromLineEdit"
QT_MOC_LITERAL(15, 253, 26), // "setPixelSizeY_fromLineEdit"
QT_MOC_LITERAL(16, 280, 24), // "setPotBound_fromLineEdit"
QT_MOC_LITERAL(17, 305, 30), // "setprobeSemiangle_fromLineEdit"
QT_MOC_LITERAL(18, 336, 30), // "setSliceThickness_fromLineEdit"
QT_MOC_LITERAL(19, 367, 24), // "setCellDimX_fromLineEdit"
QT_MOC_LITERAL(20, 392, 24), // "setCellDimY_fromLineEdit"
QT_MOC_LITERAL(21, 417, 24), // "setCellDimZ_fromLineEdit"
QT_MOC_LITERAL(22, 442, 21), // "setTileX_fromLineEdit"
QT_MOC_LITERAL(23, 464, 21), // "setTileY_fromLineEdit"
QT_MOC_LITERAL(24, 486, 21), // "setTileZ_fromLineEdit"
QT_MOC_LITERAL(25, 508, 24), // "setBatchGPU_fromLineEdit"
QT_MOC_LITERAL(26, 533, 24), // "setBatchCPU_fromLineEdit"
QT_MOC_LITERAL(27, 558, 18), // "setE0_fromLineEdit"
QT_MOC_LITERAL(28, 577, 27), // "setprobe_stepX_fromLineEdit"
QT_MOC_LITERAL(29, 605, 27), // "setprobe_stepY_fromLineEdit"
QT_MOC_LITERAL(30, 633, 13), // "setAlgo_PRISM"
QT_MOC_LITERAL(31, 647, 18), // "setAlgo_Multislice"
QT_MOC_LITERAL(32, 666, 18), // "calculatePotential"
QT_MOC_LITERAL(33, 685, 12), // "calculateAll"
QT_MOC_LITERAL(34, 698, 14), // "calculateProbe"
QT_MOC_LITERAL(35, 713, 20), // "updatePotentialImage"
QT_MOC_LITERAL(36, 734, 22), // "updatePotentialDisplay"
QT_MOC_LITERAL(37, 757, 25), // "updatePotentialFloatImage"
QT_MOC_LITERAL(38, 783, 17), // "updateOutputImage"
QT_MOC_LITERAL(39, 801, 19), // "updateOutputDisplay"
QT_MOC_LITERAL(40, 821, 22), // "updateOutputFloatImage"
QT_MOC_LITERAL(41, 844, 27), // "updateSliders_fromLineEdits"
QT_MOC_LITERAL(42, 872, 31), // "updateSliders_fromLineEdits_ang"
QT_MOC_LITERAL(43, 904, 20), // "updateContrastPotMin"
QT_MOC_LITERAL(44, 925, 20), // "updateContrastPotMax"
QT_MOC_LITERAL(45, 946, 20), // "updateContrastAngMin"
QT_MOC_LITERAL(46, 967, 20), // "updateContrastAngMax"
QT_MOC_LITERAL(47, 988, 26), // "updateSlider_lineEdits_min"
QT_MOC_LITERAL(48, 1015, 26), // "updateSlider_lineEdits_max"
QT_MOC_LITERAL(49, 1042, 30), // "updateSlider_lineEdits_max_ang"
QT_MOC_LITERAL(50, 1073, 3), // "val"
QT_MOC_LITERAL(51, 1077, 30), // "updateSlider_lineEdits_min_ang"
QT_MOC_LITERAL(52, 1108, 14), // "updateAlphaMax"
QT_MOC_LITERAL(53, 1123, 11), // "resizeEvent"
QT_MOC_LITERAL(54, 1135, 13), // "QResizeEvent*"
QT_MOC_LITERAL(55, 1149, 5), // "event"
QT_MOC_LITERAL(56, 1155, 12), // "redrawImages"
QT_MOC_LITERAL(57, 1168, 22), // "saveCurrentOutputImage"
QT_MOC_LITERAL(58, 1191, 16), // "setStreamingMode"
QT_MOC_LITERAL(59, 1208, 28), // "toggleSaveProjectedPotential"
QT_MOC_LITERAL(60, 1237, 19), // "enableOutputWidgets"
QT_MOC_LITERAL(61, 1257, 29), // "setprobe_defocus_fromLineEdit"
QT_MOC_LITERAL(62, 1287, 26), // "setRandomSeed_fromLineEdit"
QT_MOC_LITERAL(63, 1314, 24), // "setprobe_C3_fromLineEdit"
QT_MOC_LITERAL(64, 1339, 24), // "setprobe_C5_fromLineEdit"
QT_MOC_LITERAL(65, 1364, 35), // "setDetector_angle_step_fromLi..."
QT_MOC_LITERAL(66, 1400, 27), // "setprobe_Xtilt_fromLineEdit"
QT_MOC_LITERAL(67, 1428, 27), // "setprobe_Ytilt_fromLineEdit"
QT_MOC_LITERAL(68, 1456, 14), // "toggle3DOutput"
QT_MOC_LITERAL(69, 1471, 14), // "toggle4DOutput"
QT_MOC_LITERAL(70, 1486, 20), // "toggleThermalEffects"
QT_MOC_LITERAL(71, 1507, 31), // "setscan_WindowXMin_fromLineEdit"
QT_MOC_LITERAL(72, 1539, 31), // "setscan_WindowXMax_fromLineEdit"
QT_MOC_LITERAL(73, 1571, 31), // "setscan_WindowYMin_fromLineEdit"
QT_MOC_LITERAL(74, 1603, 31), // "setscan_WindowYMax_fromLineEdit"
QT_MOC_LITERAL(75, 1635, 16), // "resetCalculation"
QT_MOC_LITERAL(76, 1652, 13), // "newRandomSeed"
QT_MOC_LITERAL(77, 1666, 18), // "updateProbeK_PRISM"
QT_MOC_LITERAL(78, 1685, 37), // "PRISM::Array2D<PRISM_FLOAT_PR..."
QT_MOC_LITERAL(79, 1723, 18), // "updateProbeR_PRISM"
QT_MOC_LITERAL(80, 1742, 23), // "updateProbeK_Multislice"
QT_MOC_LITERAL(81, 1766, 23), // "updateProbeR_Multislice"
QT_MOC_LITERAL(82, 1790, 17), // "updateProbe_diffR"
QT_MOC_LITERAL(83, 1808, 12), // "arr_contrast"
QT_MOC_LITERAL(84, 1821, 17), // "updateProbe_diffK"
QT_MOC_LITERAL(85, 1839, 18), // "update_pearsonReal"
QT_MOC_LITERAL(86, 1858, 3), // "str"
QT_MOC_LITERAL(87, 1862, 15), // "update_pearsonK"
QT_MOC_LITERAL(88, 1878, 12), // "update_RReal"
QT_MOC_LITERAL(89, 1891, 9), // "update_RK"
QT_MOC_LITERAL(90, 1901, 17), // "potentialReceived"
QT_MOC_LITERAL(91, 1919, 37), // "PRISM::Array3D<PRISM_FLOAT_PR..."
QT_MOC_LITERAL(92, 1957, 14), // "outputReceived"
QT_MOC_LITERAL(93, 1972, 30), // "displayErrorReadingAtomsDialog"
QT_MOC_LITERAL(94, 2003, 25), // "setscan_WindowYMin_edited"
QT_MOC_LITERAL(95, 2029, 25), // "setscan_WindowYMax_edited"
QT_MOC_LITERAL(96, 2055, 20), // "setinterpYSet_edited"
QT_MOC_LITERAL(97, 2076, 23), // "setpixelSizeYSet_edited"
QT_MOC_LITERAL(98, 2100, 23), // "setprobeStepYSet_edited"
QT_MOC_LITERAL(99, 2124, 23), // "setprobeTiltYSet_edited"
QT_MOC_LITERAL(100, 2148, 34), // "checkInput_lineEdit_scanWindo..."
QT_MOC_LITERAL(101, 2183, 34), // "checkInput_lineEdit_scanWindo..."
QT_MOC_LITERAL(102, 2218, 34), // "checkInput_lineEdit_scanWindo..."
QT_MOC_LITERAL(103, 2253, 34), // "checkInput_lineEdit_scanWindo..."
QT_MOC_LITERAL(104, 2288, 28), // "checkInput_lineEdit_cellDimX"
QT_MOC_LITERAL(105, 2317, 28), // "checkInput_lineEdit_cellDimY"
QT_MOC_LITERAL(106, 2346, 28), // "checkInput_lineEdit_cellDimZ"
QT_MOC_LITERAL(107, 2375, 25), // "checkInput_lineEdit_tileX"
QT_MOC_LITERAL(108, 2401, 25), // "checkInput_lineEdit_tileY"
QT_MOC_LITERAL(109, 2427, 25), // "checkInput_lineEdit_tileZ"
QT_MOC_LITERAL(110, 2453, 30), // "checkInput_lineEdit_pixelSizeX"
QT_MOC_LITERAL(111, 2484, 30), // "checkInput_lineEdit_pixelSizeY"
QT_MOC_LITERAL(112, 2515, 34), // "checkInput_lineEdit_interpFac..."
QT_MOC_LITERAL(113, 2550, 34), // "checkInput_lineEdit_interpFac..."
QT_MOC_LITERAL(114, 2585, 18), // "userHasSetCellDims"
QT_MOC_LITERAL(115, 2604, 10) // "resetLinks"

    },
    "PRISMMainWindow\0setInterpolationFactorX\0"
    "\0setInterpolationFactorY\0"
    "setFilenameAtoms_fromDialog\0"
    "setFilenameOutput_fromLineEdit\0"
    "setFilenameOutput_fromDialog\0setNumGPUs\0"
    "numGPUs\0setNumThreads\0numThreads\0"
    "setNumFP\0numFP\0setNumStreams\0"
    "setPixelSizeX_fromLineEdit\0"
    "setPixelSizeY_fromLineEdit\0"
    "setPotBound_fromLineEdit\0"
    "setprobeSemiangle_fromLineEdit\0"
    "setSliceThickness_fromLineEdit\0"
    "setCellDimX_fromLineEdit\0"
    "setCellDimY_fromLineEdit\0"
    "setCellDimZ_fromLineEdit\0setTileX_fromLineEdit\0"
    "setTileY_fromLineEdit\0setTileZ_fromLineEdit\0"
    "setBatchGPU_fromLineEdit\0"
    "setBatchCPU_fromLineEdit\0setE0_fromLineEdit\0"
    "setprobe_stepX_fromLineEdit\0"
    "setprobe_stepY_fromLineEdit\0setAlgo_PRISM\0"
    "setAlgo_Multislice\0calculatePotential\0"
    "calculateAll\0calculateProbe\0"
    "updatePotentialImage\0updatePotentialDisplay\0"
    "updatePotentialFloatImage\0updateOutputImage\0"
    "updateOutputDisplay\0updateOutputFloatImage\0"
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
    "setStreamingMode\0toggleSaveProjectedPotential\0"
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
    "update_RReal\0update_RK\0potentialReceived\0"
    "PRISM::Array3D<PRISM_FLOAT_PRECISION>\0"
    "outputReceived\0displayErrorReadingAtomsDialog\0"
    "setscan_WindowYMin_edited\0"
    "setscan_WindowYMax_edited\0"
    "setinterpYSet_edited\0setpixelSizeYSet_edited\0"
    "setprobeStepYSet_edited\0setprobeTiltYSet_edited\0"
    "checkInput_lineEdit_scanWindowXMin\0"
    "checkInput_lineEdit_scanWindowXMax\0"
    "checkInput_lineEdit_scanWindowYMin\0"
    "checkInput_lineEdit_scanWindowYMax\0"
    "checkInput_lineEdit_cellDimX\0"
    "checkInput_lineEdit_cellDimY\0"
    "checkInput_lineEdit_cellDimZ\0"
    "checkInput_lineEdit_tileX\0"
    "checkInput_lineEdit_tileY\0"
    "checkInput_lineEdit_tileZ\0"
    "checkInput_lineEdit_pixelSizeX\0"
    "checkInput_lineEdit_pixelSizeY\0"
    "checkInput_lineEdit_interpFactor_x\0"
    "checkInput_lineEdit_interpFactor_y\0"
    "userHasSetCellDims\0resetLinks"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_PRISMMainWindow[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
     104,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,  534,    2, 0x0a /* Public */,
       3,    0,  535,    2, 0x0a /* Public */,
       4,    0,  536,    2, 0x0a /* Public */,
       5,    0,  537,    2, 0x0a /* Public */,
       6,    0,  538,    2, 0x0a /* Public */,
       7,    1,  539,    2, 0x0a /* Public */,
       9,    1,  542,    2, 0x0a /* Public */,
      11,    1,  545,    2, 0x0a /* Public */,
      13,    1,  548,    2, 0x0a /* Public */,
      14,    0,  551,    2, 0x0a /* Public */,
      15,    0,  552,    2, 0x0a /* Public */,
      16,    0,  553,    2, 0x0a /* Public */,
      17,    0,  554,    2, 0x0a /* Public */,
      18,    0,  555,    2, 0x0a /* Public */,
      19,    0,  556,    2, 0x0a /* Public */,
      20,    0,  557,    2, 0x0a /* Public */,
      21,    0,  558,    2, 0x0a /* Public */,
      22,    0,  559,    2, 0x0a /* Public */,
      23,    0,  560,    2, 0x0a /* Public */,
      24,    0,  561,    2, 0x0a /* Public */,
      25,    0,  562,    2, 0x0a /* Public */,
      26,    0,  563,    2, 0x0a /* Public */,
      27,    0,  564,    2, 0x0a /* Public */,
      28,    0,  565,    2, 0x0a /* Public */,
      29,    0,  566,    2, 0x0a /* Public */,
      30,    0,  567,    2, 0x0a /* Public */,
      31,    0,  568,    2, 0x0a /* Public */,
      32,    0,  569,    2, 0x0a /* Public */,
      33,    0,  570,    2, 0x0a /* Public */,
      34,    0,  571,    2, 0x0a /* Public */,
      35,    0,  572,    2, 0x0a /* Public */,
      36,    0,  573,    2, 0x0a /* Public */,
      37,    0,  574,    2, 0x0a /* Public */,
      38,    0,  575,    2, 0x0a /* Public */,
      39,    0,  576,    2, 0x0a /* Public */,
      40,    0,  577,    2, 0x0a /* Public */,
      41,    0,  578,    2, 0x0a /* Public */,
      42,    0,  579,    2, 0x0a /* Public */,
      43,    0,  580,    2, 0x0a /* Public */,
      44,    0,  581,    2, 0x0a /* Public */,
      45,    0,  582,    2, 0x0a /* Public */,
      46,    0,  583,    2, 0x0a /* Public */,
      47,    1,  584,    2, 0x0a /* Public */,
      48,    1,  587,    2, 0x0a /* Public */,
      49,    1,  590,    2, 0x0a /* Public */,
      51,    1,  593,    2, 0x0a /* Public */,
      52,    0,  596,    2, 0x0a /* Public */,
      53,    1,  597,    2, 0x0a /* Public */,
      56,    0,  600,    2, 0x0a /* Public */,
      57,    0,  601,    2, 0x0a /* Public */,
      58,    1,  602,    2, 0x0a /* Public */,
      59,    0,  605,    2, 0x0a /* Public */,
      60,    0,  606,    2, 0x0a /* Public */,
      61,    0,  607,    2, 0x0a /* Public */,
      62,    0,  608,    2, 0x0a /* Public */,
      63,    0,  609,    2, 0x0a /* Public */,
      64,    0,  610,    2, 0x0a /* Public */,
      65,    0,  611,    2, 0x0a /* Public */,
      66,    0,  612,    2, 0x0a /* Public */,
      67,    0,  613,    2, 0x0a /* Public */,
      68,    0,  614,    2, 0x0a /* Public */,
      69,    0,  615,    2, 0x0a /* Public */,
      70,    0,  616,    2, 0x0a /* Public */,
      71,    0,  617,    2, 0x0a /* Public */,
      72,    0,  618,    2, 0x0a /* Public */,
      73,    0,  619,    2, 0x0a /* Public */,
      74,    0,  620,    2, 0x0a /* Public */,
      75,    0,  621,    2, 0x0a /* Public */,
      76,    0,  622,    2, 0x0a /* Public */,
      77,    1,  623,    2, 0x0a /* Public */,
      79,    1,  626,    2, 0x0a /* Public */,
      80,    1,  629,    2, 0x0a /* Public */,
      81,    1,  632,    2, 0x0a /* Public */,
      82,    2,  635,    2, 0x0a /* Public */,
      84,    2,  640,    2, 0x0a /* Public */,
      85,    1,  645,    2, 0x0a /* Public */,
      87,    1,  648,    2, 0x0a /* Public */,
      88,    1,  651,    2, 0x0a /* Public */,
      89,    1,  654,    2, 0x0a /* Public */,
      90,    1,  657,    2, 0x0a /* Public */,
      92,    1,  660,    2, 0x0a /* Public */,
      93,    0,  663,    2, 0x0a /* Public */,
      94,    0,  664,    2, 0x0a /* Public */,
      95,    0,  665,    2, 0x0a /* Public */,
      96,    0,  666,    2, 0x0a /* Public */,
      97,    0,  667,    2, 0x0a /* Public */,
      98,    0,  668,    2, 0x0a /* Public */,
      99,    0,  669,    2, 0x0a /* Public */,
     100,    0,  670,    2, 0x0a /* Public */,
     101,    0,  671,    2, 0x0a /* Public */,
     102,    0,  672,    2, 0x0a /* Public */,
     103,    0,  673,    2, 0x0a /* Public */,
     104,    0,  674,    2, 0x0a /* Public */,
     105,    0,  675,    2, 0x0a /* Public */,
     106,    0,  676,    2, 0x0a /* Public */,
     107,    0,  677,    2, 0x0a /* Public */,
     108,    0,  678,    2, 0x0a /* Public */,
     109,    0,  679,    2, 0x0a /* Public */,
     110,    0,  680,    2, 0x0a /* Public */,
     111,    0,  681,    2, 0x0a /* Public */,
     112,    0,  682,    2, 0x0a /* Public */,
     113,    0,  683,    2, 0x0a /* Public */,
     114,    0,  684,    2, 0x0a /* Public */,
     115,    0,  685,    2, 0x0a /* Public */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,    8,
    QMetaType::Void, QMetaType::Int,   10,
    QMetaType::Void, QMetaType::Int,   12,
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
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,   50,
    QMetaType::Void, QMetaType::Int,   50,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 54,   55,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,    2,
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
    QMetaType::Void, 0x80000000 | 78,    2,
    QMetaType::Void, 0x80000000 | 78,    2,
    QMetaType::Void, 0x80000000 | 78,    2,
    QMetaType::Void, 0x80000000 | 78,    2,
    QMetaType::Void, 0x80000000 | 78, 0x80000000 | 78,    2,   83,
    QMetaType::Void, 0x80000000 | 78, 0x80000000 | 78,    2,   83,
    QMetaType::Void, QMetaType::QString,   86,
    QMetaType::Void, QMetaType::QString,   86,
    QMetaType::Void, QMetaType::QString,   86,
    QMetaType::Void, QMetaType::QString,   86,
    QMetaType::Void, 0x80000000 | 91,    2,
    QMetaType::Void, 0x80000000 | 91,    2,
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
        case 8: _t->setNumStreams((*reinterpret_cast< const int(*)>(_a[1]))); break;
        case 9: _t->setPixelSizeX_fromLineEdit(); break;
        case 10: _t->setPixelSizeY_fromLineEdit(); break;
        case 11: _t->setPotBound_fromLineEdit(); break;
        case 12: _t->setprobeSemiangle_fromLineEdit(); break;
        case 13: _t->setSliceThickness_fromLineEdit(); break;
        case 14: _t->setCellDimX_fromLineEdit(); break;
        case 15: _t->setCellDimY_fromLineEdit(); break;
        case 16: _t->setCellDimZ_fromLineEdit(); break;
        case 17: _t->setTileX_fromLineEdit(); break;
        case 18: _t->setTileY_fromLineEdit(); break;
        case 19: _t->setTileZ_fromLineEdit(); break;
        case 20: _t->setBatchGPU_fromLineEdit(); break;
        case 21: _t->setBatchCPU_fromLineEdit(); break;
        case 22: _t->setE0_fromLineEdit(); break;
        case 23: _t->setprobe_stepX_fromLineEdit(); break;
        case 24: _t->setprobe_stepY_fromLineEdit(); break;
        case 25: _t->setAlgo_PRISM(); break;
        case 26: _t->setAlgo_Multislice(); break;
        case 27: _t->calculatePotential(); break;
        case 28: _t->calculateAll(); break;
        case 29: _t->calculateProbe(); break;
        case 30: _t->updatePotentialImage(); break;
        case 31: _t->updatePotentialDisplay(); break;
        case 32: _t->updatePotentialFloatImage(); break;
        case 33: _t->updateOutputImage(); break;
        case 34: _t->updateOutputDisplay(); break;
        case 35: _t->updateOutputFloatImage(); break;
        case 36: _t->updateSliders_fromLineEdits(); break;
        case 37: _t->updateSliders_fromLineEdits_ang(); break;
        case 38: _t->updateContrastPotMin(); break;
        case 39: _t->updateContrastPotMax(); break;
        case 40: _t->updateContrastAngMin(); break;
        case 41: _t->updateContrastAngMax(); break;
        case 42: _t->updateSlider_lineEdits_min((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 43: _t->updateSlider_lineEdits_max((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 44: _t->updateSlider_lineEdits_max_ang((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 45: _t->updateSlider_lineEdits_min_ang((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 46: _t->updateAlphaMax(); break;
        case 47: _t->resizeEvent((*reinterpret_cast< QResizeEvent*(*)>(_a[1]))); break;
        case 48: _t->redrawImages(); break;
        case 49: _t->saveCurrentOutputImage(); break;
        case 50: _t->setStreamingMode((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 51: _t->toggleSaveProjectedPotential(); break;
        case 52: _t->enableOutputWidgets(); break;
        case 53: _t->setprobe_defocus_fromLineEdit(); break;
        case 54: _t->setRandomSeed_fromLineEdit(); break;
        case 55: _t->setprobe_C3_fromLineEdit(); break;
        case 56: _t->setprobe_C5_fromLineEdit(); break;
        case 57: _t->setDetector_angle_step_fromLineEdit(); break;
        case 58: _t->setprobe_Xtilt_fromLineEdit(); break;
        case 59: _t->setprobe_Ytilt_fromLineEdit(); break;
        case 60: _t->toggle3DOutput(); break;
        case 61: _t->toggle4DOutput(); break;
        case 62: _t->toggleThermalEffects(); break;
        case 63: _t->setscan_WindowXMin_fromLineEdit(); break;
        case 64: _t->setscan_WindowXMax_fromLineEdit(); break;
        case 65: _t->setscan_WindowYMin_fromLineEdit(); break;
        case 66: _t->setscan_WindowYMax_fromLineEdit(); break;
        case 67: _t->resetCalculation(); break;
        case 68: _t->newRandomSeed(); break;
        case 69: _t->updateProbeK_PRISM((*reinterpret_cast< PRISM::Array2D<PRISM_FLOAT_PRECISION>(*)>(_a[1]))); break;
        case 70: _t->updateProbeR_PRISM((*reinterpret_cast< PRISM::Array2D<PRISM_FLOAT_PRECISION>(*)>(_a[1]))); break;
        case 71: _t->updateProbeK_Multislice((*reinterpret_cast< PRISM::Array2D<PRISM_FLOAT_PRECISION>(*)>(_a[1]))); break;
        case 72: _t->updateProbeR_Multislice((*reinterpret_cast< PRISM::Array2D<PRISM_FLOAT_PRECISION>(*)>(_a[1]))); break;
        case 73: _t->updateProbe_diffR((*reinterpret_cast< PRISM::Array2D<PRISM_FLOAT_PRECISION>(*)>(_a[1])),(*reinterpret_cast< PRISM::Array2D<PRISM_FLOAT_PRECISION>(*)>(_a[2]))); break;
        case 74: _t->updateProbe_diffK((*reinterpret_cast< PRISM::Array2D<PRISM_FLOAT_PRECISION>(*)>(_a[1])),(*reinterpret_cast< PRISM::Array2D<PRISM_FLOAT_PRECISION>(*)>(_a[2]))); break;
        case 75: _t->update_pearsonReal((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 76: _t->update_pearsonK((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 77: _t->update_RReal((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 78: _t->update_RK((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 79: _t->potentialReceived((*reinterpret_cast< PRISM::Array3D<PRISM_FLOAT_PRECISION>(*)>(_a[1]))); break;
        case 80: _t->outputReceived((*reinterpret_cast< PRISM::Array3D<PRISM_FLOAT_PRECISION>(*)>(_a[1]))); break;
        case 81: _t->displayErrorReadingAtomsDialog(); break;
        case 82: _t->setscan_WindowYMin_edited(); break;
        case 83: _t->setscan_WindowYMax_edited(); break;
        case 84: _t->setinterpYSet_edited(); break;
        case 85: _t->setpixelSizeYSet_edited(); break;
        case 86: _t->setprobeStepYSet_edited(); break;
        case 87: _t->setprobeTiltYSet_edited(); break;
        case 88: _t->checkInput_lineEdit_scanWindowXMin(); break;
        case 89: _t->checkInput_lineEdit_scanWindowXMax(); break;
        case 90: _t->checkInput_lineEdit_scanWindowYMin(); break;
        case 91: _t->checkInput_lineEdit_scanWindowYMax(); break;
        case 92: _t->checkInput_lineEdit_cellDimX(); break;
        case 93: _t->checkInput_lineEdit_cellDimY(); break;
        case 94: _t->checkInput_lineEdit_cellDimZ(); break;
        case 95: _t->checkInput_lineEdit_tileX(); break;
        case 96: _t->checkInput_lineEdit_tileY(); break;
        case 97: _t->checkInput_lineEdit_tileZ(); break;
        case 98: _t->checkInput_lineEdit_pixelSizeX(); break;
        case 99: _t->checkInput_lineEdit_pixelSizeY(); break;
        case 100: _t->checkInput_lineEdit_interpFactor_x(); break;
        case 101: _t->checkInput_lineEdit_interpFactor_y(); break;
        case 102: _t->userHasSetCellDims(); break;
        case 103: _t->resetLinks(); break;
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
        if (_id < 104)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 104;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 104)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 104;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
