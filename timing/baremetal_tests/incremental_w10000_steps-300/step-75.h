#pragma once

const bool step75_is_reconstruct = true;

const int step75_factor339_height = 7;
const int step75_factor339_width = 3;
int step75_factor339_ridx[] = {21, 22, 23, 24, 25, 26, 27, };
float step75_factor339_data[] = {
1.0000000, 0.0000000, 0.0000000, -0.9994170, -0.0341407, 0.0878746, -0.0000000, 
0.0000000, 1.0000000, 0.0000000, 0.0341407, -0.9994170, 1.0057700, -0.0000000, 
0.0000000, 0.0000000, 1.0000000, -0.0000000, -0.0000000, -1.0000000, -0.0000000, 
};

const int step75_factor339_num_blks = 1;
int step75_factor339_A_blk_start[] = {0, };
int step75_factor339_B_blk_start[] = {21, };
int step75_factor339_blk_width[] = {6, };

const int step75_factor340_height = 4;
const int step75_factor340_width = 3;
int step75_factor340_ridx[] = {3, 4, 5, 9, };
float step75_factor340_data[] = {
0.0010000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0010000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0010000, 0.0000000, 
};

const int step75_factor340_num_blks = 1;
int step75_factor340_A_blk_start[] = {0, };
int step75_factor340_B_blk_start[] = {3, };
int step75_factor340_blk_width[] = {3, };

const int step75_factor341_height = 7;
const int step75_factor341_width = 3;
int step75_factor341_ridx[] = {0, 1, 2, 3, 4, 5, 9, };
float step75_factor341_data[] = {
1.0000000, 0.0000000, 0.0000000, -0.9996110, -0.0278942, 0.0613971, -0.0000000, 
0.0000000, 1.0000000, 0.0000000, 0.0278942, -0.9996110, 0.9262350, -0.0000000, 
0.0000000, 0.0000000, 1.0000000, -0.0000000, -0.0000000, -1.0000000, -0.0000000, 
};

const int step75_factor341_num_blks = 1;
int step75_factor341_A_blk_start[] = {0, };
int step75_factor341_B_blk_start[] = {0, };
int step75_factor341_blk_width[] = {6, };

const int step75_factor342_height = 4;
const int step75_factor342_width = 3;
int step75_factor342_ridx[] = {6, 7, 8, 9, };
float step75_factor342_data[] = {
0.0010000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0010000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0010000, 0.0000000, 
};

const int step75_factor342_num_blks = 1;
int step75_factor342_A_blk_start[] = {0, };
int step75_factor342_B_blk_start[] = {6, };
int step75_factor342_blk_width[] = {3, };

const int step75_factor343_height = 7;
const int step75_factor343_width = 3;
int step75_factor343_ridx[] = {3, 4, 5, 6, 7, 8, 9, };
float step75_factor343_data[] = {
1.0000000, 0.0000000, 0.0000000, -0.9999950, 0.0032549, 0.0167794, -0.0000000, 
0.0000000, 1.0000000, 0.0000000, -0.0032549, -0.9999950, 1.0084100, -0.0000000, 
0.0000000, 0.0000000, 1.0000000, -0.0000000, -0.0000000, -1.0000000, -0.0000000, 
};

const int step75_factor343_num_blks = 1;
int step75_factor343_A_blk_start[] = {0, };
int step75_factor343_B_blk_start[] = {3, };
int step75_factor343_blk_width[] = {6, };

const int step75_factor338_height = 4;
const int step75_factor338_width = 3;
int step75_factor338_ridx[] = {0, 1, 2, 9, };
float step75_factor338_data[] = {
0.0010000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0010000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0010000, 0.0000000, 
};

const int step75_factor338_num_blks = 1;
int step75_factor338_A_blk_start[] = {0, };
int step75_factor338_B_blk_start[] = {0, };
int step75_factor338_blk_width[] = {3, };

const int step75_node7_num_factors = 1;
const bool step75_node7_marked = false;
const bool step75_node7_fixed = true;
int step75_node7_factor_height[] = {step75_factor339_height, };
int step75_node7_factor_width[] = {step75_factor339_width, };
int* step75_node7_factor_ridx[] = {step75_factor339_ridx, };
float* step75_node7_factor_data[] = {step75_factor339_data, };
int step75_node7_factor_num_blks[] = {step75_factor339_num_blks, };
int* step75_node7_factor_A_blk_start[] = {step75_factor339_A_blk_start, };
int* step75_node7_factor_B_blk_start[] = {step75_factor339_B_blk_start, };
int* step75_node7_factor_blk_width[] = {step75_factor339_blk_width, };
const int step75_node7_parent = 8;
const int step75_node7_height = 28;
const int step75_node7_width = 24;
float step75_node7_data[] = {
1.08971, 0.00294391, 0.00600652, 0.0203832, -0.917449, 0.974457, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.15203e-08, 
0.0, 1.0285, -0.180397, 0.971989, 0.0242222, 0.0275821, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.63373e-07, 
0.0, 0.0, 1.32374, 0.132368, 0.0074639, -0.756095, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.06517e-07, 
0.0, 0.0, 0.0, 1.01848, -0.00572537, 0.106262, -0.981845, 0.0043028, -0.0161944, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.93645e-07, 
0.0, 0.0, 0.0, 0.0, 1.07592, -0.149937, -0.00929781, -0.929402, 0.944931, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.01226e-07, 
0.0, 0.0, 0.0, 0.0, 0.0, 1.2541, 0.082082, -0.111481, -0.683037, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.42229e-08, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.01447, 0.00466633, 0.060118, -0.985712, 0.00650749, -0.0324992, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.24838e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.06007, -0.202784, -0.00188857, -0.943344, 0.95001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.10777e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.27658, 0.0461202, -0.150156, -0.630904, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -7.71316e-13, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.01304, 0.0114094, 0.0248542, -0.986233, 0.0421132, -0.0310218, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.24881e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.04278, -0.1973, -0.0301214, -0.958565, 0.986816, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.84086e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.29376, 0.0143528, -0.146991, -0.621853, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.25275e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.01303, 0.0145798, -0.00437611, -0.987104, 0.00793999, 0.00991942, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.04018e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.02838, -0.168409, 0.00617318, -0.972486, 0.908176, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.82801e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.29202, -0.0025387, -0.126733, -0.655572, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.20867e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.01271, 0.0133495, -0.0148521, -0.0111543, 0.987386, -0.938718, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.99292e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.01881, -0.131287, -0.981331, -0.0240254, 0.126191, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.48234e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.26512, -0.101968, 0.0090984, -0.788366, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.19327e-13, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.01315, -0.0114845, -0.0713722, -0.986569, 0.02993, -0.117564, 0.0, 0.0, 0.0, -8.41517e-13, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.01207, -0.0154507, -0.0411572, -0.987284, 0.974003, 0.0, 0.0, 0.0, -5.41513e-13, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.18035, -0.0601935, -0.0111137, -0.841565, 0.0, 0.0, 0.0, -1.80886e-13, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.01063, -0.011651, -0.0370384, -0.988908, -0.0337817, 0.0869506, -2.90337e-13, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.01199, -0.0344973, 0.022351, -0.987969, 0.994855, -2.52234e-13, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.14685, -0.0312652, -0.0308091, -0.839219, -6.28493e-14, 
};
const int step75_node7_num_blks = 1;
int step75_node7_A_blk_start[] = {0, };
int step75_node7_B_blk_start[] = {0, };
int step75_node7_blk_width[] = {3, };
const float step75_node7_H_correct_cond = 93.68343029354715;
float step75_node7_H_correct_data[] = {
1.1874678841, 0.0032080081661, 0.006545364909199999, 0.022211776872, -0.99975334979, 1.06187553747, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.614218611299999e-08, 
0.0032080081661, 1.0578209166060881, -0.1855206318457068, 0.999750692806312, 0.022211645414409998, 0.031236903556869998, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.82430801626373e-07, 
0.006545364909199999, -0.1855206318457068, 1.7848667434915102, -6.512145359839107e-07, -2.499487999928939e-08, -0.9999958279340598, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.1639772211355983e-08, 
0.022211776872, 0.999750692806312, -6.512145359839107e-07, 2.0000008887872402, -2.7823399999264606e-08, 0.05481498851929999, -0.9999894956000001, 0.004382315744, -0.016493672512000002, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.6086929013196e-07, 
-0.99975334979, 0.022211645414409998, -2.499487999928939e-08, -2.7823399999264606e-08, 1.9999917186386869, -1.06091852482782, -0.004382273827550001, -0.9999868349620359, 1.016762880451928, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.018844219585e-08, 
1.06187553747, 0.031236903556869998, -0.9999958279340598, 0.05481498851929999, -1.06091852482782, 3.1285463927274098, 3.085479700089263e-07, 6.497075999958792e-07, -0.9999976703798, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.5497322576599994e-08, 
0.0, 0.0, 0.0, -0.9999894956000001, -0.004382273827550001, 3.085479700089263e-07, 1.999992888919796, -1.1103280000170917e-08, 0.012037466192889992, -0.99997525264, 0.0066016533803000006, -0.032969463424, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.9018405990179863e-07, 
0.0, 0.0, 0.0, 0.004382315744, -0.9999868349620359, 6.497075999958792e-07, -1.1103280000170917e-08, 2.000004784588509, -1.01682949918238, -0.00660167387686, -0.9999803079841884, 1.006925448708064, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -9.915950949496085e-08, 
0.0, 0.0, 0.0, -0.016493672512000002, 1.016762880451928, -0.9999976703798, 0.012037466192889992, -1.01682949918238, 3.0340884177013603, 6.267887999829507e-08, 1.4049981996430222e-07, -1.0000000430656, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.703085820626124e-08, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.99997525264, -0.00660167387686, 6.267887999829507e-08, 2.000008828088685, 1.4019999997193301e-08, 0.0263215711519, -0.99909347832, 0.042662356128, -0.031426324272, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -6.5539507759543e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0066016533803000006, -0.9999803079841884, 1.4049981996430222e-07, 1.4019999997193301e-08, 2.0000073769064604, -1.007120623125528, -0.0426623202822, -0.9990919243559201, 1.0286780483550801, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0893728365282e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.032969463424, 1.006925448708064, -1.0000000430656, 0.0263215711519, -1.007120623125528, 3.01497501417428, -1.4805999997665542e-09, 4.882354399847233e-07, -0.99999835610156, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.1581860916600024e-13, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.99909347832, -0.0426623202822, -1.4805999997665542e-09, 1.9999986127948004, -6.541540000182016e-08, -0.012488019034700003, -0.99996596512, 0.0080434480697, 0.0100486700426, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.990287712600038e-14, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.042662356128, -0.9990919243559201, 4.882354399847233e-07, -6.541540000182016e-08, 2.00000472988828, -1.029079162213338, -0.0080434040508, -0.9999693892137981, 0.934094658039716, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.915209809940005e-13, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.031426324272, 1.0286780483550801, -0.99999835610156, -0.012488019034700003, -1.029079162213338, 3.059165745559972, 5.440820000020947e-09, -9.221556389120129e-07, -1.0000005558970562, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0908202394778e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.99996596512, -0.0080434040508, 5.440820000020947e-09, 2.0000004040650023, -1.7802339998768105e-08, -0.0175617347946, -0.011296071153000001, 0.99993567606, -0.9506491057800001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.9903903697628e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0080434480697, -0.9999693892137981, -9.221556389120129e-07, -1.7802339998768105e-08, 2.00000534217645, -0.933981855743344, -0.9999387404378499, -0.011296208366999998, 0.116033236769, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.4055151669082e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0100486700426, 0.934094658039716, -1.0000005558970562, -0.0175617347946, -0.933981855743344, 2.8726421646965465, -8.838397002936333e-08, 3.48872000004739e-08, -1.0000028981292, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.1312834080164007e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.011296071153000001, -0.9999387404378499, -8.838397002936333e-08, 2.0000053454934896, 3.3213999971572873e-09, -0.10528703817559998, -0.99954238235, 0.0303235795, -0.1191099666, 0.0, 0.0, 0.0, -3.32950205107e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.99993567606, -0.011296208366999998, 3.48872000004739e-08, 3.3213999971572873e-09, 2.00000869236397, -0.9518991855319, -0.0303237157235, -0.999544248965, 0.987109379968, 0.0, 0.0, 0.0, -2.5601702937026996e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.9506491057800001, 0.116033236769, -1.0000028981292, -0.10528703817559998, -0.9518991855319, 2.9171954395243302, 1.0980684000336336e-07, 3.157799999609588e-09, -0.9999994745813, 0.0, 0.0, 0.0, 1.5507013782545e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.99954238235, -0.0303237157235, 1.0980684000336336e-07, 2.00000856121509, -4.2754249998393086e-08, 0.0891223862799, -0.9994200920399999, -0.034140799470999995, 0.08787488487799999, 5.699666231476001e-13, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0303235795, -0.999544248965, 3.157799999609588e-09, -4.2754249998393086e-08, 2.00000852178469, -0.9902627556600999, 0.034140755598, -0.9994211577233, 1.0057702500094, 2.5957526034720004e-13, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.1191099666, 0.987109379968, -0.9999994745813, 0.0891223862799, -0.9902627556600999, 2.9885616166118494, 2.6294900001292658e-08, 6.676598000447319e-08, -0.99999863264454, -3.28899565157e-13, 
};
float step75_node7_M_correct_data[] = {
1.08971, 0.00294391, 0.00600652, 0.0203832, -0.917449, 0.974457, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.15203e-08, 
0.0, 1.0285, -0.180397, 0.971989, 0.0242222, 0.0275821, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.63373e-07, 
0.0, 0.0, 1.32374, 0.132368, 0.0074639, -0.756095, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.06517e-07, 
0.0, 0.0, 0.0, 1.01848, -0.00572537, 0.106262, -0.981845, 0.0043028, -0.0161944, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.93645e-07, 
0.0, 0.0, 0.0, 0.0, 1.07592, -0.149937, -0.00929781, -0.929402, 0.944931, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.01226e-07, 
0.0, 0.0, 0.0, 0.0, 0.0, 1.2541, 0.082082, -0.111481, -0.683037, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.42229e-08, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.01447, 0.00466633, 0.060118, -0.985712, 0.00650749, -0.0324992, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.24838e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.06007, -0.202784, -0.00188857, -0.943344, 0.95001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.10777e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.27658, 0.0461202, -0.150156, -0.630904, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -7.71316e-13, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.01304, 0.0114094, 0.0248542, -0.986233, 0.0421132, -0.0310218, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.24881e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.04278, -0.1973, -0.0301214, -0.958565, 0.986816, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.84086e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.29376, 0.0143528, -0.146991, -0.621853, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.25275e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.01303, 0.0145798, -0.00437611, -0.987104, 0.00793999, 0.00991942, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.04018e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.02838, -0.168409, 0.00617318, -0.972486, 0.908176, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.82801e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.29202, -0.0025387, -0.126733, -0.655572, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.20867e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.01271, 0.0133495, -0.0148521, -0.0111543, 0.987386, -0.938718, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.99292e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.01881, -0.131287, -0.981331, -0.0240254, 0.126191, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.48234e-12, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.26512, -0.101968, 0.0090984, -0.788366, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.19327e-13, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.01315, -0.0114845, -0.0713722, -0.986569, 0.02993, -0.117564, 0.0, 0.0, 0.0, -8.41517e-13, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.01207, -0.0154507, -0.0411572, -0.987284, 0.974003, 0.0, 0.0, 0.0, -5.41513e-13, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.18035, -0.0601935, -0.0111137, -0.841565, 0.0, 0.0, 0.0, -1.80886e-13, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.01063, -0.011651, -0.0370384, -0.988908, -0.0337817, 0.0869506, -2.90337e-13, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.01199, -0.0344973, 0.022351, -0.987969, 0.994855, -2.52234e-13, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.14685, -0.0312652, -0.0308091, -0.839219, -6.28493e-14, 
};


const int step75_node8_num_factors = 5;
const bool step75_node8_marked = true;
const bool step75_node8_fixed = false;
int step75_node8_factor_height[] = {step75_factor338_height, step75_factor340_height, step75_factor341_height, step75_factor342_height, step75_factor343_height, };
int step75_node8_factor_width[] = {step75_factor338_width, step75_factor340_width, step75_factor341_width, step75_factor342_width, step75_factor343_width, };
int* step75_node8_factor_ridx[] = {step75_factor338_ridx, step75_factor340_ridx, step75_factor341_ridx, step75_factor342_ridx, step75_factor343_ridx, };
float* step75_node8_factor_data[] = {step75_factor338_data, step75_factor340_data, step75_factor341_data, step75_factor342_data, step75_factor343_data, };
int step75_node8_factor_num_blks[] = {step75_factor338_num_blks, step75_factor340_num_blks, step75_factor341_num_blks, step75_factor342_num_blks, step75_factor343_num_blks, };
int* step75_node8_factor_A_blk_start[] = {step75_factor338_A_blk_start, step75_factor340_A_blk_start, step75_factor341_A_blk_start, step75_factor342_A_blk_start, step75_factor343_A_blk_start, };
int* step75_node8_factor_B_blk_start[] = {step75_factor338_B_blk_start, step75_factor340_B_blk_start, step75_factor341_B_blk_start, step75_factor342_B_blk_start, step75_factor343_B_blk_start, };
int* step75_node8_factor_blk_width[] = {step75_factor338_blk_width, step75_factor340_blk_width, step75_factor341_blk_width, step75_factor342_blk_width, step75_factor343_blk_width, };
const int step75_node8_parent = 9;
const int step75_node8_height = 10;
const int step75_node8_width = 9;
float step75_node8_data[] = {
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
};
const int step75_node8_num_blks = 0;
int step75_node8_A_blk_start[] = {};
int step75_node8_B_blk_start[] = {};
int step75_node8_blk_width[] = {};
const float step75_node8_H_correct_cond = 8896.73039185979;
float step75_node8_H_correct_data[] = {
1.0205848576, -0.012288155264000001, -0.015974116928, -0.9996112649600001, -0.027894140736000002, 0.061397032928, 0.0, 0.0, 0.0, -1.6548034272e-13, 
-0.012288155264000001, 1.02182416156496, -0.04821200044208, 0.0278941584964, -0.99960959515496, 0.92623294793908, 0.0, 0.0, 0.0, -1.903731507492e-13, 
-0.015974116928, -0.04821200044208, 1.31769353579045, 9.524765999625454e-08, 2.60315179996461e-07, -0.99999660178794, 0.0, 0.0, 0.0, -3.013782442e-14, 
-0.9996112649600001, 0.0278941584964, 9.524765999625454e-08, 2.0000077757622, -8.225477999800251e-08, -0.035536611995230005, -0.99999778711, 0.0032549572326, 0.016779386058, 1.9185771611079976e-14, 
-0.027894140736000002, -0.99960959515496, 2.60315179996461e-07, -8.225477999800251e-08, 1.99999296448602, -0.9275870809693754, -0.0032548659125000014, -0.9999904019259901, 1.0084058156283, 6.993763598739996e-15, 
0.061397032928, 0.92623294793908, -0.99999660178794, -0.035536611995230005, -0.9275870809693754, 2.8616853962180326, -1.1749369998972317e-09, -8.521754637664819e-08, -1.0000019037608954, -1.81107780599601e-13, 
0.0, 0.0, 0.0, -0.99999778711, -0.0032548659125000014, -1.1749369998972317e-09, 1.0000007346299296, -3.910552399876119e-08, -0.0200615299896, 1.4049640288328201e-13, 
0.0, 0.0, 0.0, 0.0032549572326, -0.9999904019259901, -8.521754637664819e-08, -3.910552399876119e-08, 1.0000009740862539, -1.0083503806754521, 1.8652966247064495e-13, 
0.0, 0.0, 0.0, 0.016779386058, 1.0084058156283, -1.0000019037608954, -0.0200615299896, -1.0083503806754521, 2.01717450376588, -1.6285328859803313e-13, 
};
float step75_node8_M_correct_data[] = {
1.01024, -0.0121636, -0.0158122, -0.989479, -0.0276114, 0.0607747, 0.0, 0.0, 0.0, -1.63803e-13, 
0.0, 1.01078, -0.0478881, 0.0156894, -0.989281, 0.917086, 0.0, 0.0, 0.0, -1.90314e-13, 
0.0, 0.0, 1.1468, -0.0129878, -0.0416909, -0.832855, 0.0, 0.0, 0.0, -3.64856e-14, 
0.0, 0.0, 0.0, 1.01021, -0.0122165, -0.000600673, -0.989891, 0.00322206, 0.0166098, -1.38963e-13, 
0.0, 0.0, 0.0, 0.0, 1.00929, -0.0528917, -0.0152066, -0.990747, 0.999325, -1.87282e-13, 
0.0, 0.0, 0.0, 0.0, 0.0, 1.14913, -0.00121736, -0.0456001, -0.82422, -3.21937e-14, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.14101, -0.0846178, 0.0749828, 3.62145e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0957676, -0.517558, -9.46258e-17, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.255775, -1.47282e-18, 
};


const int step75_node9_num_factors = 0;
const bool step75_node9_marked = true;
const bool step75_node9_fixed = false;
int step75_node9_factor_height[] = {};
int step75_node9_factor_width[] = {};
int* step75_node9_factor_ridx[] = {};
float* step75_node9_factor_data[] = {};
int step75_node9_factor_num_blks[] = {};
int* step75_node9_factor_A_blk_start[] = {};
int* step75_node9_factor_B_blk_start[] = {};
int* step75_node9_factor_blk_width[] = {};
const int step75_node9_parent = -1;
const int step75_node9_height = 1;
const int step75_node9_width = 1;
float step75_node9_data[] = {
0, 
};
const int step75_node9_num_blks = 0;
int step75_node9_A_blk_start[] = {};
int step75_node9_B_blk_start[] = {};
int step75_node9_blk_width[] = {};
const float step75_node9_H_correct_cond = 1.0;
float step75_node9_H_correct_data[] = {
1.4352998416e-50, 
};
float step75_node9_M_correct_data[] = {
-1.19804e-25, 
};


const int step75_nnodes = 10;
bool step75_node_marked[] = {false, false, false, false, false, false, false, step75_node7_marked, step75_node8_marked, step75_node9_marked, };
bool step75_node_fixed[] = {false, false, false, false, false, false, false, step75_node7_fixed, step75_node8_fixed, step75_node9_fixed, };
int step75_node_num_factors[] = {0, 0, 0, 0, 0, 0, 0, step75_node7_num_factors, step75_node8_num_factors, step75_node9_num_factors, };
int* step75_node_factor_height[] = {0, 0, 0, 0, 0, 0, 0, step75_node7_factor_height, step75_node8_factor_height, step75_node9_factor_height, };
int* step75_node_factor_width[] = {0, 0, 0, 0, 0, 0, 0, step75_node7_factor_width, step75_node8_factor_width, step75_node9_factor_width, };
int** step75_node_factor_ridx[] = {0, 0, 0, 0, 0, 0, 0, step75_node7_factor_ridx, step75_node8_factor_ridx, step75_node9_factor_ridx, };
float** step75_node_factor_data[] = {0, 0, 0, 0, 0, 0, 0, step75_node7_factor_data, step75_node8_factor_data, step75_node9_factor_data, };
int* step75_node_factor_num_blks[] = {0, 0, 0, 0, 0, 0, 0, step75_node7_factor_num_blks, step75_node8_factor_num_blks, step75_node9_factor_num_blks, };
int** step75_node_factor_A_blk_start[] = {0, 0, 0, 0, 0, 0, 0, step75_node7_factor_A_blk_start, step75_node8_factor_A_blk_start, step75_node9_factor_A_blk_start, };
int** step75_node_factor_B_blk_start[] = {0, 0, 0, 0, 0, 0, 0, step75_node7_factor_B_blk_start, step75_node8_factor_B_blk_start, step75_node9_factor_B_blk_start, };
int** step75_node_factor_blk_width[] = {0, 0, 0, 0, 0, 0, 0, step75_node7_factor_blk_width, step75_node8_factor_blk_width, step75_node9_factor_blk_width, };
int step75_node_parent[] = {0, 0, 0, 0, 0, 0, 0, step75_node7_parent, step75_node8_parent, step75_node9_parent, };
int step75_node_height[] = {0, 0, 0, 0, 0, 0, 0, step75_node7_height, step75_node8_height, step75_node9_height, };
int step75_node_width[] = {0, 0, 0, 0, 0, 0, 0, step75_node7_width, step75_node8_width, step75_node9_width, };
float* step75_node_data[] = {0, 0, 0, 0, 0, 0, 0, step75_node7_data, step75_node8_data, step75_node9_data, };
int step75_node_num_blks[] = {0, 0, 0, 0, 0, 0, 0, step75_node7_num_blks, step75_node8_num_blks, step75_node9_num_blks, };
int* step75_node_A_blk_start[] = {0, 0, 0, 0, 0, 0, 0, step75_node7_A_blk_start, step75_node8_A_blk_start, step75_node9_A_blk_start, };
int* step75_node_B_blk_start[] = {0, 0, 0, 0, 0, 0, 0, step75_node7_B_blk_start, step75_node8_B_blk_start, step75_node9_B_blk_start, };
int* step75_node_blk_width[] = {0, 0, 0, 0, 0, 0, 0, step75_node7_blk_width, step75_node8_blk_width, step75_node9_blk_width, };
float step75_node_H_correct_cond[] = {0, 0, 0, 0, 0, 0, 0, step75_node7_H_correct_cond, step75_node8_H_correct_cond, step75_node9_H_correct_cond, };
float* step75_node_H_correct_data[] = {0, 0, 0, 0, 0, 0, 0, step75_node7_H_correct_data, step75_node8_H_correct_data, step75_node9_H_correct_data, };
float* step75_node_M_correct_data[] = {0, 0, 0, 0, 0, 0, 0, step75_node7_M_correct_data, step75_node8_M_correct_data, step75_node9_M_correct_data, };


