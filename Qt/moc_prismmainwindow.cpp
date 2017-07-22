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
    QByteArrayData data[121];
    char stringdata0[2731];
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
QT_MOC_LITERAL(28, 577, 26), // "setprobeStepX_fromLineEdit"
QT_MOC_LITERAL(29, 604, 26), // "setprobeStepY_fromLineEdit"
QT_MOC_LITERAL(30, 631, 13), // "setAlgo_PRISM"
QT_MOC_LITERAL(31, 645, 18), // "setAlgo_Multislice"
QT_MOC_LITERAL(32, 664, 18), // "calculatePotential"
QT_MOC_LITERAL(33, 683, 12), // "calculateAll"
QT_MOC_LITERAL(34, 696, 14), // "calculateProbe"
QT_MOC_LITERAL(35, 711, 20), // "updatePotentialImage"
QT_MOC_LITERAL(36, 732, 22), // "updatePotentialDisplay"
QT_MOC_LITERAL(37, 755, 25), // "updatePotentialFloatImage"
QT_MOC_LITERAL(38, 781, 17), // "updateOutputImage"
QT_MOC_LITERAL(39, 799, 19), // "updateOutputDisplay"
QT_MOC_LITERAL(40, 819, 22), // "updateOutputFloatImage"
QT_MOC_LITERAL(41, 842, 27), // "updateSliders_fromLineEdits"
QT_MOC_LITERAL(42, 870, 31), // "updateSliders_fromLineEdits_ang"
QT_MOC_LITERAL(43, 902, 20), // "updateContrastPotMin"
QT_MOC_LITERAL(44, 923, 20), // "updateContrastPotMax"
QT_MOC_LITERAL(45, 944, 20), // "updateContrastAngMin"
QT_MOC_LITERAL(46, 965, 20), // "updateContrastAngMax"
QT_MOC_LITERAL(47, 986, 26), // "updateSlider_lineEdits_min"
QT_MOC_LITERAL(48, 1013, 26), // "updateSlider_lineEdits_max"
QT_MOC_LITERAL(49, 1040, 30), // "updateSlider_lineEdits_max_ang"
QT_MOC_LITERAL(50, 1071, 3), // "val"
QT_MOC_LITERAL(51, 1075, 30), // "updateSlider_lineEdits_min_ang"
QT_MOC_LITERAL(52, 1106, 14), // "updateAlphaMax"
QT_MOC_LITERAL(53, 1121, 11), // "resizeEvent"
QT_MOC_LITERAL(54, 1133, 13), // "QResizeEvent*"
QT_MOC_LITERAL(55, 1147, 5), // "event"
QT_MOC_LITERAL(56, 1153, 12), // "redrawImages"
QT_MOC_LITERAL(57, 1166, 22), // "saveCurrentOutputImage"
QT_MOC_LITERAL(58, 1189, 16), // "setStreamingMode"
QT_MOC_LITERAL(59, 1206, 28), // "toggleSaveProjectedPotential"
QT_MOC_LITERAL(60, 1235, 19), // "enableOutputWidgets"
QT_MOC_LITERAL(61, 1255, 29), // "setprobe_defocus_fromLineEdit"
QT_MOC_LITERAL(62, 1285, 26), // "setRandomSeed_fromLineEdit"
QT_MOC_LITERAL(63, 1312, 24), // "setprobe_C3_fromLineEdit"
QT_MOC_LITERAL(64, 1337, 24), // "setprobe_C5_fromLineEdit"
QT_MOC_LITERAL(65, 1362, 33), // "setdetectorAngleStep_fromLine..."
QT_MOC_LITERAL(66, 1396, 27), // "setprobe_Xtilt_fromLineEdit"
QT_MOC_LITERAL(67, 1424, 27), // "setprobe_Ytilt_fromLineEdit"
QT_MOC_LITERAL(68, 1452, 14), // "toggle3DOutput"
QT_MOC_LITERAL(69, 1467, 14), // "toggle4DOutput"
QT_MOC_LITERAL(70, 1482, 20), // "toggleThermalEffects"
QT_MOC_LITERAL(71, 1503, 15), // "toggleOccupancy"
QT_MOC_LITERAL(72, 1519, 31), // "setscan_WindowXMin_fromLineEdit"
QT_MOC_LITERAL(73, 1551, 31), // "setscan_WindowXMax_fromLineEdit"
QT_MOC_LITERAL(74, 1583, 31), // "setscan_WindowYMin_fromLineEdit"
QT_MOC_LITERAL(75, 1615, 31), // "setscan_WindowYMax_fromLineEdit"
QT_MOC_LITERAL(76, 1647, 16), // "resetCalculation"
QT_MOC_LITERAL(77, 1664, 13), // "newRandomSeed"
QT_MOC_LITERAL(78, 1678, 18), // "updateProbeK_PRISM"
QT_MOC_LITERAL(79, 1697, 45), // "Prismatic::Array2D<PRISMATIC_..."
QT_MOC_LITERAL(80, 1743, 18), // "updateProbeR_PRISM"
QT_MOC_LITERAL(81, 1762, 23), // "updateProbeK_Multislice"
QT_MOC_LITERAL(82, 1786, 23), // "updateProbeR_Multislice"
QT_MOC_LITERAL(83, 1810, 17), // "updateProbe_diffR"
QT_MOC_LITERAL(84, 1828, 12), // "arr_contrast"
QT_MOC_LITERAL(85, 1841, 17), // "updateProbe_diffK"
QT_MOC_LITERAL(86, 1859, 18), // "update_pearsonReal"
QT_MOC_LITERAL(87, 1878, 3), // "str"
QT_MOC_LITERAL(88, 1882, 15), // "update_pearsonK"
QT_MOC_LITERAL(89, 1898, 12), // "update_RReal"
QT_MOC_LITERAL(90, 1911, 9), // "update_RK"
QT_MOC_LITERAL(91, 1921, 17), // "potentialReceived"
QT_MOC_LITERAL(92, 1939, 45), // "Prismatic::Array3D<PRISMATIC_..."
QT_MOC_LITERAL(93, 1985, 14), // "outputReceived"
QT_MOC_LITERAL(94, 2000, 30), // "displayErrorReadingAtomsDialog"
QT_MOC_LITERAL(95, 2031, 25), // "setscan_WindowYMin_edited"
QT_MOC_LITERAL(96, 2057, 25), // "setscan_WindowYMax_edited"
QT_MOC_LITERAL(97, 2083, 20), // "setinterpYSet_edited"
QT_MOC_LITERAL(98, 2104, 23), // "setpixelSizeYSet_edited"
QT_MOC_LITERAL(99, 2128, 23), // "setprobeStepYSet_edited"
QT_MOC_LITERAL(100, 2152, 23), // "setprobeTiltYSet_edited"
QT_MOC_LITERAL(101, 2176, 34), // "checkInput_lineEdit_scanWindo..."
QT_MOC_LITERAL(102, 2211, 34), // "checkInput_lineEdit_scanWindo..."
QT_MOC_LITERAL(103, 2246, 34), // "checkInput_lineEdit_scanWindo..."
QT_MOC_LITERAL(104, 2281, 34), // "checkInput_lineEdit_scanWindo..."
QT_MOC_LITERAL(105, 2316, 28), // "checkInput_lineEdit_cellDimX"
QT_MOC_LITERAL(106, 2345, 28), // "checkInput_lineEdit_cellDimY"
QT_MOC_LITERAL(107, 2374, 28), // "checkInput_lineEdit_cellDimZ"
QT_MOC_LITERAL(108, 2403, 25), // "checkInput_lineEdit_tileX"
QT_MOC_LITERAL(109, 2429, 25), // "checkInput_lineEdit_tileY"
QT_MOC_LITERAL(110, 2455, 25), // "checkInput_lineEdit_tileZ"
QT_MOC_LITERAL(111, 2481, 30), // "checkInput_lineEdit_pixelSizeX"
QT_MOC_LITERAL(112, 2512, 30), // "checkInput_lineEdit_pixelSizeY"
QT_MOC_LITERAL(113, 2543, 34), // "checkInput_lineEdit_interpFac..."
QT_MOC_LITERAL(114, 2578, 34), // "checkInput_lineEdit_interpFac..."
QT_MOC_LITERAL(115, 2613, 18), // "userHasSetCellDims"
QT_MOC_LITERAL(116, 2632, 10), // "resetLinks"
QT_MOC_LITERAL(117, 2643, 24), // "moveBothPotentialSliders"
QT_MOC_LITERAL(118, 2668, 27), // "updateSlider_PotentialCombo"
QT_MOC_LITERAL(119, 2696, 19), // "openSaveAtomsDialog"
QT_MOC_LITERAL(120, 2716, 14) // "saveAtomCoords"

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
    "setprobeStepX_fromLineEdit\0"
    "setprobeStepY_fromLineEdit\0setAlgo_PRISM\0"
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
    "setdetectorAngleStep_fromLineEdit\0"
    "setprobe_Xtilt_fromLineEdit\0"
    "setprobe_Ytilt_fromLineEdit\0toggle3DOutput\0"
    "toggle4DOutput\0toggleThermalEffects\0"
    "toggleOccupancy\0setscan_WindowXMin_fromLineEdit\0"
    "setscan_WindowXMax_fromLineEdit\0"
    "setscan_WindowYMin_fromLineEdit\0"
    "setscan_WindowYMax_fromLineEdit\0"
    "resetCalculation\0newRandomSeed\0"
    "updateProbeK_PRISM\0"
    "Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>\0"
    "updateProbeR_PRISM\0updateProbeK_Multislice\0"
    "updateProbeR_Multislice\0updateProbe_diffR\0"
    "arr_contrast\0updateProbe_diffK\0"
    "update_pearsonReal\0str\0update_pearsonK\0"
    "update_RReal\0update_RK\0potentialReceived\0"
    "Prismatic::Array3D<PRISMATIC_FLOAT_PRECISION>\0"
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
    "userHasSetCellDims\0resetLinks\0"
    "moveBothPotentialSliders\0"
    "updateSlider_PotentialCombo\0"
    "openSaveAtomsDialog\0saveAtomCoords"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_PRISMMainWindow[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
     109,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,  559,    2, 0x0a /* Public */,
       3,    0,  560,    2, 0x0a /* Public */,
       4,    0,  561,    2, 0x0a /* Public */,
       5,    0,  562,    2, 0x0a /* Public */,
       6,    0,  563,    2, 0x0a /* Public */,
       7,    1,  564,    2, 0x0a /* Public */,
       9,    1,  567,    2, 0x0a /* Public */,
      11,    1,  570,    2, 0x0a /* Public */,
      13,    1,  573,    2, 0x0a /* Public */,
      14,    0,  576,    2, 0x0a /* Public */,
      15,    0,  577,    2, 0x0a /* Public */,
      16,    0,  578,    2, 0x0a /* Public */,
      17,    0,  579,    2, 0x0a /* Public */,
      18,    0,  580,    2, 0x0a /* Public */,
      19,    0,  581,    2, 0x0a /* Public */,
      20,    0,  582,    2, 0x0a /* Public */,
      21,    0,  583,    2, 0x0a /* Public */,
      22,    0,  584,    2, 0x0a /* Public */,
      23,    0,  585,    2, 0x0a /* Public */,
      24,    0,  586,    2, 0x0a /* Public */,
      25,    0,  587,    2, 0x0a /* Public */,
      26,    0,  588,    2, 0x0a /* Public */,
      27,    0,  589,    2, 0x0a /* Public */,
      28,    0,  590,    2, 0x0a /* Public */,
      29,    0,  591,    2, 0x0a /* Public */,
      30,    0,  592,    2, 0x0a /* Public */,
      31,    0,  593,    2, 0x0a /* Public */,
      32,    0,  594,    2, 0x0a /* Public */,
      33,    0,  595,    2, 0x0a /* Public */,
      34,    0,  596,    2, 0x0a /* Public */,
      35,    0,  597,    2, 0x0a /* Public */,
      36,    0,  598,    2, 0x0a /* Public */,
      37,    0,  599,    2, 0x0a /* Public */,
      38,    0,  600,    2, 0x0a /* Public */,
      39,    0,  601,    2, 0x0a /* Public */,
      40,    0,  602,    2, 0x0a /* Public */,
      41,    0,  603,    2, 0x0a /* Public */,
      42,    0,  604,    2, 0x0a /* Public */,
      43,    0,  605,    2, 0x0a /* Public */,
      44,    0,  606,    2, 0x0a /* Public */,
      45,    0,  607,    2, 0x0a /* Public */,
      46,    0,  608,    2, 0x0a /* Public */,
      47,    1,  609,    2, 0x0a /* Public */,
      48,    1,  612,    2, 0x0a /* Public */,
      49,    1,  615,    2, 0x0a /* Public */,
      51,    1,  618,    2, 0x0a /* Public */,
      52,    0,  621,    2, 0x0a /* Public */,
      53,    1,  622,    2, 0x0a /* Public */,
      56,    0,  625,    2, 0x0a /* Public */,
      57,    0,  626,    2, 0x0a /* Public */,
      58,    1,  627,    2, 0x0a /* Public */,
      59,    0,  630,    2, 0x0a /* Public */,
      60,    0,  631,    2, 0x0a /* Public */,
      61,    0,  632,    2, 0x0a /* Public */,
      62,    0,  633,    2, 0x0a /* Public */,
      63,    0,  634,    2, 0x0a /* Public */,
      64,    0,  635,    2, 0x0a /* Public */,
      65,    0,  636,    2, 0x0a /* Public */,
      66,    0,  637,    2, 0x0a /* Public */,
      67,    0,  638,    2, 0x0a /* Public */,
      68,    0,  639,    2, 0x0a /* Public */,
      69,    0,  640,    2, 0x0a /* Public */,
      70,    0,  641,    2, 0x0a /* Public */,
      71,    0,  642,    2, 0x0a /* Public */,
      72,    0,  643,    2, 0x0a /* Public */,
      73,    0,  644,    2, 0x0a /* Public */,
      74,    0,  645,    2, 0x0a /* Public */,
      75,    0,  646,    2, 0x0a /* Public */,
      76,    0,  647,    2, 0x0a /* Public */,
      77,    0,  648,    2, 0x0a /* Public */,
      78,    1,  649,    2, 0x0a /* Public */,
      80,    1,  652,    2, 0x0a /* Public */,
      81,    1,  655,    2, 0x0a /* Public */,
      82,    1,  658,    2, 0x0a /* Public */,
      83,    2,  661,    2, 0x0a /* Public */,
      85,    2,  666,    2, 0x0a /* Public */,
      86,    1,  671,    2, 0x0a /* Public */,
      88,    1,  674,    2, 0x0a /* Public */,
      89,    1,  677,    2, 0x0a /* Public */,
      90,    1,  680,    2, 0x0a /* Public */,
      91,    1,  683,    2, 0x0a /* Public */,
      93,    1,  686,    2, 0x0a /* Public */,
      94,    0,  689,    2, 0x0a /* Public */,
      95,    0,  690,    2, 0x0a /* Public */,
      96,    0,  691,    2, 0x0a /* Public */,
      97,    0,  692,    2, 0x0a /* Public */,
      98,    0,  693,    2, 0x0a /* Public */,
      99,    0,  694,    2, 0x0a /* Public */,
     100,    0,  695,    2, 0x0a /* Public */,
     101,    0,  696,    2, 0x0a /* Public */,
     102,    0,  697,    2, 0x0a /* Public */,
     103,    0,  698,    2, 0x0a /* Public */,
     104,    0,  699,    2, 0x0a /* Public */,
     105,    0,  700,    2, 0x0a /* Public */,
     106,    0,  701,    2, 0x0a /* Public */,
     107,    0,  702,    2, 0x0a /* Public */,
     108,    0,  703,    2, 0x0a /* Public */,
     109,    0,  704,    2, 0x0a /* Public */,
     110,    0,  705,    2, 0x0a /* Public */,
     111,    0,  706,    2, 0x0a /* Public */,
     112,    0,  707,    2, 0x0a /* Public */,
     113,    0,  708,    2, 0x0a /* Public */,
     114,    0,  709,    2, 0x0a /* Public */,
     115,    0,  710,    2, 0x0a /* Public */,
     116,    0,  711,    2, 0x0a /* Public */,
     117,    1,  712,    2, 0x0a /* Public */,
     118,    1,  715,    2, 0x0a /* Public */,
     119,    0,  718,    2, 0x0a /* Public */,
     120,    2,  719,    2, 0x0a /* Public */,

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
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 79,    2,
    QMetaType::Void, 0x80000000 | 79,    2,
    QMetaType::Void, 0x80000000 | 79,    2,
    QMetaType::Void, 0x80000000 | 79,    2,
    QMetaType::Void, 0x80000000 | 79, 0x80000000 | 79,    2,   84,
    QMetaType::Void, 0x80000000 | 79, 0x80000000 | 79,    2,   84,
    QMetaType::Void, QMetaType::QString,   87,
    QMetaType::Void, QMetaType::QString,   87,
    QMetaType::Void, QMetaType::QString,   87,
    QMetaType::Void, QMetaType::QString,   87,
    QMetaType::Void, 0x80000000 | 92,    2,
    QMetaType::Void, 0x80000000 | 92,    2,
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
    QMetaType::Void,
    QMetaType::Void, QMetaType::QString, QMetaType::QString,    2,    2,

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
        case 23: _t->setprobeStepX_fromLineEdit(); break;
        case 24: _t->setprobeStepY_fromLineEdit(); break;
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
        case 57: _t->setdetectorAngleStep_fromLineEdit(); break;
        case 58: _t->setprobe_Xtilt_fromLineEdit(); break;
        case 59: _t->setprobe_Ytilt_fromLineEdit(); break;
        case 60: _t->toggle3DOutput(); break;
        case 61: _t->toggle4DOutput(); break;
        case 62: _t->toggleThermalEffects(); break;
        case 63: _t->toggleOccupancy(); break;
        case 64: _t->setscan_WindowXMin_fromLineEdit(); break;
        case 65: _t->setscan_WindowXMax_fromLineEdit(); break;
        case 66: _t->setscan_WindowYMin_fromLineEdit(); break;
        case 67: _t->setscan_WindowYMax_fromLineEdit(); break;
        case 68: _t->resetCalculation(); break;
        case 69: _t->newRandomSeed(); break;
        case 70: _t->updateProbeK_PRISM((*reinterpret_cast< Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>(*)>(_a[1]))); break;
        case 71: _t->updateProbeR_PRISM((*reinterpret_cast< Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>(*)>(_a[1]))); break;
        case 72: _t->updateProbeK_Multislice((*reinterpret_cast< Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>(*)>(_a[1]))); break;
        case 73: _t->updateProbeR_Multislice((*reinterpret_cast< Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>(*)>(_a[1]))); break;
        case 74: _t->updateProbe_diffR((*reinterpret_cast< Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>(*)>(_a[1])),(*reinterpret_cast< Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>(*)>(_a[2]))); break;
        case 75: _t->updateProbe_diffK((*reinterpret_cast< Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>(*)>(_a[1])),(*reinterpret_cast< Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>(*)>(_a[2]))); break;
        case 76: _t->update_pearsonReal((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 77: _t->update_pearsonK((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 78: _t->update_RReal((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 79: _t->update_RK((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 80: _t->potentialReceived((*reinterpret_cast< Prismatic::Array3D<PRISMATIC_FLOAT_PRECISION>(*)>(_a[1]))); break;
        case 81: _t->outputReceived((*reinterpret_cast< Prismatic::Array3D<PRISMATIC_FLOAT_PRECISION>(*)>(_a[1]))); break;
        case 82: _t->displayErrorReadingAtomsDialog(); break;
        case 83: _t->setscan_WindowYMin_edited(); break;
        case 84: _t->setscan_WindowYMax_edited(); break;
        case 85: _t->setinterpYSet_edited(); break;
        case 86: _t->setpixelSizeYSet_edited(); break;
        case 87: _t->setprobeStepYSet_edited(); break;
        case 88: _t->setprobeTiltYSet_edited(); break;
        case 89: _t->checkInput_lineEdit_scanWindowXMin(); break;
        case 90: _t->checkInput_lineEdit_scanWindowXMax(); break;
        case 91: _t->checkInput_lineEdit_scanWindowYMin(); break;
        case 92: _t->checkInput_lineEdit_scanWindowYMax(); break;
        case 93: _t->checkInput_lineEdit_cellDimX(); break;
        case 94: _t->checkInput_lineEdit_cellDimY(); break;
        case 95: _t->checkInput_lineEdit_cellDimZ(); break;
        case 96: _t->checkInput_lineEdit_tileX(); break;
        case 97: _t->checkInput_lineEdit_tileY(); break;
        case 98: _t->checkInput_lineEdit_tileZ(); break;
        case 99: _t->checkInput_lineEdit_pixelSizeX(); break;
        case 100: _t->checkInput_lineEdit_pixelSizeY(); break;
        case 101: _t->checkInput_lineEdit_interpFactor_x(); break;
        case 102: _t->checkInput_lineEdit_interpFactor_y(); break;
        case 103: _t->userHasSetCellDims(); break;
        case 104: _t->resetLinks(); break;
        case 105: _t->moveBothPotentialSliders((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 106: _t->updateSlider_PotentialCombo((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 107: _t->openSaveAtomsDialog(); break;
        case 108: _t->saveAtomCoords((*reinterpret_cast< QString(*)>(_a[1])),(*reinterpret_cast< QString(*)>(_a[2]))); break;
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
        if (_id < 109)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 109;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 109)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 109;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
