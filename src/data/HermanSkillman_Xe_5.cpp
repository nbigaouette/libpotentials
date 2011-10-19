
#include <cstring> // memcpy()

#include <Memory.hpp>

#include "data/HermanSkillman.hpp"

// **************************************************************
fdouble * Load_HermanSkillman_Xe_5()
{
    const int N = 440;
    fdouble * HS_Xe_5 = (fdouble *) calloc_and_check(2*N, sizeof(fdouble));

    // Data is on the stack
    const fdouble stack_HS_Xe_5[2*N] = {
        5.8557995E-04, -184432.5,
        1.1711599E-03, -91948.73,
        1.7567398E-03, -61119.87,
        2.3423198E-03, -45704.85,
        2.9278998E-03, -36455.60,
        3.5134798E-03, -30289.22,
        4.0990598E-03, -25884.69,
        4.6846396E-03, -22581.37,
        5.2702194E-03, -20012.23,
        5.8557992E-03, -17957.04,
        6.4413790E-03, -16275.70,
        7.0269587E-03, -14874.75,
        7.6125385E-03, -13689.51,
        8.1981188E-03, -12673.77,
        8.7836990E-03, -11793.65,
        9.3692793E-03, -11023.72,
        9.9548595E-03, -10344.56,
        1.0540440E-02, -9741.035,
        1.1126020E-02, -9201.210,
        1.1711600E-02, -8715.533,
        1.2297180E-02, -8276.272,
        1.2882761E-02, -7877.100,
        1.3468341E-02, -7512.787,
        1.4053921E-02, -7178.973,
        1.4639501E-02, -6872.000,
        1.5225082E-02, -6588.768,
        1.5810661E-02, -6326.641,
        1.6396241E-02, -6083.356,
        1.6981822E-02, -5856.965,
        1.7567402E-02, -5645.778,
        1.8152982E-02, -5448.320,
        1.8738562E-02, -5263.300,
        1.9324142E-02, -5089.588,
        1.9909723E-02, -4926.183,
        2.0495303E-02, -4772.203,
        2.1080883E-02, -4626.859,
        2.1666463E-02, -4489.453,
        2.2252044E-02, -4359.355,
        2.2837624E-02, -4236.003,
        2.3423204E-02, -4118.887,
        2.4594365E-02, -3912.135,
        2.5765525E-02, -3714.234,
        2.6936686E-02, -3533.770,
        2.8107846E-02, -3368.560,
        2.9279007E-02, -3216.767,
        3.0450167E-02, -3076.844,
        3.1621326E-02, -2947.463,
        3.2792486E-02, -2827.492,
        3.3963647E-02, -2715.953,
        3.5134807E-02, -2612.004,
        3.6305968E-02, -2514.904,
        3.7477128E-02, -2424.010,
        3.8648289E-02, -2338.756,
        3.9819449E-02, -2258.645,
        4.0990610E-02, -2183.228,
        4.2161770E-02, -2112.113,
        4.3332931E-02, -2044.948,
        4.4504091E-02, -1981.422,
        4.5675252E-02, -1921.250,
        4.6846412E-02, -1864.179,
        4.8017573E-02, -1809.979,
        4.9188733E-02, -1758.445,
        5.0359894E-02, -1709.386,
        5.1531054E-02, -1662.634,
        5.2702215E-02, -1618.032,
        5.3873375E-02, -1575.440,
        5.5044536E-02, -1534.727,
        5.6215696E-02, -1495.773,
        5.7386857E-02, -1458.470,
        5.8558017E-02, -1422.718,
        5.9729178E-02, -1388.425,
        6.0900338E-02, -1355.504,
        6.2071498E-02, -1323.878,
        6.3242659E-02, -1293.474,
        6.4413816E-02, -1264.224,
        6.5584973E-02, -1236.066,
        6.6756129E-02, -1208.942,
        6.7927286E-02, -1182.798,
        6.9098443E-02, -1157.583,
        7.0269600E-02, -1133.250,
        7.2611921E-02, -1091.856,
        7.4954242E-02, -1048.478,
        7.7296562E-02, -1007.873,
        7.9638883E-02, -969.7921,
        8.1981204E-02, -934.0171,
        8.4323525E-02, -900.3552,
        8.6665846E-02, -868.6290,
        8.9008167E-02, -838.6837,
        9.1350488E-02, -810.3798,
        9.3692809E-02, -783.5936,
        9.6035130E-02, -758.2095,
        9.8377451E-02, -734.1251,
        0.1007198, -711.2478,
        0.1030621, -689.4952,
        0.1054044, -668.7879,
        0.1077467, -649.0558,
        0.1100891, -630.2352,
        0.1124314, -612.2687,
        0.1147737, -595.1010,
        0.1171160, -578.6827,
        0.1194583, -562.9687,
        0.1218007, -547.9179,
        0.1241430, -533.4905,
        0.1264853, -519.6508,
        0.1288276, -506.3660,
        0.1311699, -493.6061,
        0.1335122, -481.3416,
        0.1358546, -469.5459,
        0.1381969, -458.1947,
        0.1405392, -447.2654,
        0.1428815, -436.7359,
        0.1452238, -426.5864,
        0.1475661, -416.7981,
        0.1499084, -407.3542,
        0.1522508, -398.2372,
        0.1545931, -389.4317,
        0.1569354, -380.9233,
        0.1592777, -372.6988,
        0.1616200, -364.7445,
        0.1639623, -357.0484,
        0.1686470, -344.7434,
        0.1733316, -330.8681,
        0.1780162, -317.8281,
        0.1827009, -305.5561,
        0.1873855, -293.9914,
        0.1920702, -283.0816,
        0.1967548, -272.7741,
        0.2014395, -263.0247,
        0.2061241, -253.7932,
        0.2108087, -245.0448,
        0.2154934, -236.7438,
        0.2201780, -228.8604,
        0.2248627, -221.3669,
        0.2295473, -214.2393,
        0.2342319, -207.4526,
        0.2389166, -200.9857,
        0.2436012, -194.8187,
        0.2482859, -188.9347,
        0.2529705, -183.3153,
        0.2576551, -177.9452,
        0.2623398, -172.8096,
        0.2670244, -167.8965,
        0.2717090, -163.1916,
        0.2763937, -158.6837,
        0.2810783, -154.3616,
        0.2857629, -150.2164,
        0.2904475, -146.2369,
        0.2951322, -142.4147,
        0.2998168, -138.7412,
        0.3045014, -135.2095,
        0.3091860, -131.8112,
        0.3138707, -128.5399,
        0.3185553, -125.3887,
        0.3232399, -122.3528,
        0.3279245, -119.4251,
        0.3326092, -116.6008,
        0.3372938, -113.8746,
        0.3419784, -111.2427,
        0.3466631, -108.6996,
        0.3513477, -106.2416,
        0.3607170, -102.5154,
        0.3700863, -98.08318,
        0.3794555, -93.92526,
        0.3888248, -90.01948,
        0.3981941, -86.34568,
        0.4075634, -82.88750,
        0.4169327, -79.62665,
        0.4263020, -76.54897,
        0.4356712, -73.64133,
        0.4450405, -70.89301,
        0.4544098, -68.29151,
        0.4637791, -65.82713,
        0.4731484, -63.49052,
        0.4825177, -61.27431,
        0.4918869, -59.16891,
        0.5012562, -57.16723,
        0.5106255, -55.26238,
        0.5199947, -53.44901,
        0.5293640, -51.72002,
        0.5387332, -50.07031,
        0.5481025, -48.49488,
        0.5574718, -46.98991,
        0.5668410, -45.55040,
        0.5762103, -44.17274,
        0.5855795, -42.85338,
        0.5949488, -41.58955,
        0.6043180, -40.37777,
        0.6136873, -39.21547,
        0.6230565, -38.10005,
        0.6324258, -37.02950,
        0.6417950, -36.00122,
        0.6511643, -35.01334,
        0.6605335, -34.06386,
        0.6699028, -33.15123,
        0.6792721, -32.27349,
        0.6886413, -31.42919,
        0.6980106, -30.61674,
        0.7073798, -29.83499,
        0.7167491, -29.08228,
        0.7261183, -28.35745,
        0.7448569, -27.31132,
        0.7635955, -26.01567,
        0.7823340, -24.81100,
        0.8010726, -23.68981,
        0.8198112, -22.64550,
        0.8385497, -21.67221,
        0.8572883, -20.76363,
        0.8760269, -19.91441,
        0.8947654, -19.12005,
        0.9135040, -18.37655,
        0.9322426, -17.67940,
        0.9509811, -17.02489,
        0.9697197, -16.40994,
        0.9884583, -15.83186,
        1.007197, -15.28740,
        1.025935, -14.77403,
        1.044674, -14.28959,
        1.063412, -13.83229,
        1.082151, -13.39980,
        1.100889, -12.99036,
        1.119628, -12.60250,
        1.138366, -12.23495,
        1.157105, -11.88609,
        1.175843, -11.55466,
        1.194582, -11.23964,
        1.213320, -10.94013,
        1.232059, -10.65498,
        1.250797, -10.38323,
        1.269536, -10.12418,
        1.288274, -9.877151,
        1.307013, -9.641299,
        1.325751, -9.415898,
        1.344490, -9.200445,
        1.363228, -8.994408,
        1.381967, -8.797173,
        1.400705, -8.608162,
        1.419444, -8.427022,
        1.438182, -8.253345,
        1.456921, -8.086639,
        1.475659, -7.926423,
        1.513137, -7.676244,
        1.550614, -7.391720,
        1.588091, -7.127681,
        1.625568, -6.882062,
        1.663045, -6.653044,
        1.700522, -6.439140,
        1.737999, -6.238696,
        1.775477, -6.050460,
        1.812954, -5.873331,
        1.850431, -5.706548,
        1.887908, -5.548968,
        1.925385, -5.399880,
        1.962862, -5.258623,
        2.000339, -5.124498,
        2.037816, -4.997105,
        2.075293, -4.876118,
        2.112770, -4.760716,
        2.150247, -4.650729,
        2.187724, -4.570960,
        2.225201, -4.493975,
        2.262678, -4.419541,
        2.300155, -4.347532,
        2.337632, -4.277833,
        2.375109, -4.210332,
        2.412586, -4.144929,
        2.450063, -4.081527,
        2.487540, -4.020035,
        2.525017, -3.960369,
        2.562495, -3.902447,
        2.599972, -3.846196,
        2.637449, -3.791543,
        2.674926, -3.738422,
        2.712403, -3.686768,
        2.749880, -3.636523,
        2.787357, -3.587629,
        2.824834, -3.540031,
        2.862311, -3.493681,
        2.899788, -3.448528,
        2.937265, -3.404528,
        2.974742, -3.361636,
        3.049696, -3.279015,
        3.124650, -3.200358,
        3.199605, -3.125386,
        3.274559, -3.053847,
        3.349513, -2.985509,
        3.424467, -2.920162,
        3.499422, -2.857615,
        3.574376, -2.797691,
        3.649330, -2.740229,
        3.724284, -2.685080,
        3.799239, -2.632106,
        3.874193, -2.581183,
        3.949147, -2.532192,
        4.024101, -2.485027,
        4.099055, -2.439586,
        4.174009, -2.395778,
        4.248963, -2.353515,
        4.323917, -2.312718,
        4.398871, -2.273310,
        4.473825, -2.235224,
        4.548779, -2.198392,
        4.623734, -2.162754,
        4.698688, -2.128254,
        4.773642, -2.094837,
        4.848596, -2.062453,
        4.923550, -2.031055,
        4.998504, -2.000599,
        5.073458, -1.971042,
        5.148412, -1.942347,
        5.223366, -1.914474,
        5.298320, -1.887391,
        5.373274, -1.861063,
        5.448228, -1.835459,
        5.523182, -1.810551,
        5.598136, -1.786309,
        5.673090, -1.762708,
        5.748044, -1.739722,
        5.822998, -1.717328,
        5.897952, -1.695504,
        5.972906, -1.674227,
        6.122815, -1.633236,
        6.272723, -1.594204,
        6.422632, -1.556994,
        6.572540, -1.521482,
        6.722449, -1.487553,
        6.872357, -1.455105,
        7.022266, -1.424042,
        7.172174, -1.394277,
        7.322083, -1.365732,
        7.471992, -1.338331,
        7.621900, -1.312009,
        7.771809, -1.286702,
        7.921717, -1.262353,
        8.071626, -1.238908,
        8.221534, -1.216318,
        8.371442, -1.194537,
        8.521350, -1.173523,
        8.671258, -1.153235,
        8.821166, -1.133637,
        8.971074, -1.114694,
        9.120982, -1.096373,
        9.270890, -1.078645,
        9.420798, -1.061481,
        9.570706, -1.044855,
        9.720614, -1.028742,
        9.870522, -1.013118,
        10.02043, -0.9979611,
        10.17034, -0.9832515,
        10.32025, -0.9689691,
        10.47015, -0.9550957,
        10.62006, -0.9416140,
        10.76997, -0.9285076,
        10.91988, -0.9157611,
        11.06979, -0.9033598,
        11.21970, -0.8912898,
        11.36960, -0.8795382,
        11.51951, -0.8680924,
        11.66942, -0.8569407,
        11.81933, -0.8460718,
        11.96924, -0.8354753,
        12.26905, -0.8150589,
        12.56887, -0.7956165,
        12.86869, -0.7770801,
        13.16850, -0.7593877,
        13.46832, -0.7424831,
        13.76814, -0.7263146,
        14.06796, -0.7108354,
        14.36777, -0.6960021,
        14.66759, -0.6817753,
        14.96741, -0.6681184,
        15.26722, -0.6549979,
        15.56704, -0.6423829,
        15.86686, -0.6302445,
        16.16667, -0.6185564,
        16.46649, -0.6072940,
        16.76631, -0.5964343,
        17.06612, -0.5859562,
        17.36594, -0.5758399,
        17.66575, -0.5660670,
        17.96557, -0.5566202,
        18.26539, -0.5474836,
        18.56520, -0.5386421,
        18.86502, -0.5300816,
        19.16483, -0.5217890,
        19.46465, -0.5137518,
        19.76447, -0.5059585,
        20.06428, -0.4983981,
        20.36410, -0.4910603,
        20.66392, -0.4839354,
        20.96373, -0.4770143,
        21.26355, -0.4702884,
        21.56336, -0.4637495,
        21.86318, -0.4573900,
        22.16300, -0.4512025,
        22.46281, -0.4451802,
        22.76263, -0.4393166,
        23.06244, -0.4336054,
        23.36226, -0.4280408,
        23.66208, -0.4226172,
        23.96189, -0.4173293,
        24.56153, -0.4071408,
        25.16116, -0.3974379,
        25.76080, -0.3881868,
        26.36043, -0.3793565,
        26.96006, -0.3709190,
        27.55970, -0.3628487,
        28.15933, -0.3551221,
        28.75897, -0.3477176,
        29.35860, -0.3406157,
        29.95823, -0.3337981,
        30.55787, -0.3272479,
        31.15750, -0.3209500,
        31.75714, -0.3148898,
        32.35677, -0.3090543,
        32.95640, -0.3034312,
        33.55603, -0.2980090,
        34.15567, -0.2927772,
        34.75530, -0.2877259,
        35.35493, -0.2828460,
        35.95456, -0.2781288,
        36.55420, -0.2735664,
        37.15383, -0.2691513,
        37.75346, -0.2648764,
        38.35309, -0.2607352,
        38.95272, -0.2567215,
        39.55236, -0.2528294,
        40.15199, -0.2490537,
        40.75162, -0.2453890,
        41.35125, -0.2418306,
        41.95089, -0.2383740,
        42.55052, -0.2350148,
        43.15015, -0.2317489,
        43.74978, -0.2285726,
        44.34941, -0.2254821,
        44.94905, -0.2224741,
        45.54868, -0.2195453,
        46.14831, -0.2166927,
        46.74794, -0.2139132,
        47.34758, -0.2112041,
        47.94721, -0.2085627
    };

    // Copy to heap
    memcpy(HS_Xe_5, stack_HS_Xe_5, 2*N*sizeof(fdouble));

    return HS_Xe_5;
}

// ********** End of file ***************************************
