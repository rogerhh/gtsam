#pragma once

const bool step83_is_reconstruct = true;

const int step83_factor359_height = 7;
const int step83_factor359_width = 3;
int step83_factor359_ridx[] = {3, 4, 5, 6, 7, 8, 9, };
float step83_factor359_data[] = {
1.0000000, 0.0000000, 0.0000000, -0.9999850, -0.0055218, -0.0533144, 0.0000000, 
0.0000000, 1.0000000, 0.0000000, 0.0055218, -0.9999850, 0.9372260, -0.0000000, 
0.0000000, 0.0000000, 1.0000000, -0.0000000, -0.0000000, -1.0000000, -0.0000000, 
};

const int step83_factor359_num_blks = 1;
int step83_factor359_A_blk_start[] = {0, };
int step83_factor359_B_blk_start[] = {3, };
int step83_factor359_blk_width[] = {6, };

const int step83_factor354_height = 4;
const int step83_factor354_width = 3;
int step83_factor354_ridx[] = {0, 1, 2, 9, };
float step83_factor354_data[] = {
0.0010000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0010000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0010000, 0.0000000, 
};

const int step83_factor354_num_blks = 1;
int step83_factor354_A_blk_start[] = {0, };
int step83_factor354_B_blk_start[] = {0, };
int step83_factor354_blk_width[] = {3, };

const int step83_factor355_height = 7;
const int step83_factor355_width = 3;
int step83_factor355_ridx[] = {21, 22, 23, 24, 25, 26, 27, };
float step83_factor355_data[] = {
1.0000000, 0.0000000, 0.0000000, 0.0002037, 1.0000000, -0.9591730, 0.0000000, 
0.0000000, 1.0000000, 0.0000000, -1.0000000, 0.0002037, -0.0398352, 0.0000000, 
0.0000000, 0.0000000, 1.0000000, -0.0000000, -0.0000000, -1.0000000, 0.0000000, 
};

const int step83_factor355_num_blks = 1;
int step83_factor355_A_blk_start[] = {0, };
int step83_factor355_B_blk_start[] = {21, };
int step83_factor355_blk_width[] = {6, };

const int step83_factor356_height = 4;
const int step83_factor356_width = 3;
int step83_factor356_ridx[] = {3, 4, 5, 9, };
float step83_factor356_data[] = {
0.0010000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0010000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0010000, 0.0000000, 
};

const int step83_factor356_num_blks = 1;
int step83_factor356_A_blk_start[] = {0, };
int step83_factor356_B_blk_start[] = {3, };
int step83_factor356_blk_width[] = {3, };

const int step83_factor357_height = 7;
const int step83_factor357_width = 3;
int step83_factor357_ridx[] = {0, 1, 2, 3, 4, 5, 9, };
float step83_factor357_data[] = {
1.0000000, 0.0000000, 0.0000000, -0.9999910, -0.0042715, 0.0622255, -0.0000000, 
0.0000000, 1.0000000, 0.0000000, 0.0042715, -0.9999910, 1.0416500, 0.0000000, 
0.0000000, 0.0000000, 1.0000000, -0.0000000, -0.0000000, -1.0000000, -0.0000000, 
};

const int step83_factor357_num_blks = 1;
int step83_factor357_A_blk_start[] = {0, };
int step83_factor357_B_blk_start[] = {0, };
int step83_factor357_blk_width[] = {6, };

const int step83_factor358_height = 4;
const int step83_factor358_width = 3;
int step83_factor358_ridx[] = {6, 7, 8, 9, };
float step83_factor358_data[] = {
0.0010000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0010000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0010000, 0.0000000, 
};

const int step83_factor358_num_blks = 1;
int step83_factor358_A_blk_start[] = {0, };
int step83_factor358_B_blk_start[] = {6, };
int step83_factor358_blk_width[] = {3, };

const int step83_node8_num_factors = 1;
const bool step83_node8_marked = false;
const bool step83_node8_fixed = true;
int step83_node8_factor_height[] = {step83_factor355_height, };
int step83_node8_factor_width[] = {step83_factor355_width, };
int* step83_node8_factor_ridx[] = {step83_factor355_ridx, };
float* step83_node8_factor_data[] = {step83_factor355_data, };
int step83_node8_factor_num_blks[] = {step83_factor355_num_blks, };
int* step83_node8_factor_A_blk_start[] = {step83_factor355_A_blk_start, };
int* step83_node8_factor_B_blk_start[] = {step83_factor355_B_blk_start, };
int* step83_node8_factor_blk_width[] = {step83_factor355_blk_width, };
const int step83_node8_parent = 9;
const int step83_node8_height = 28;
const int step83_node8_width = 24;
float step83_node8_data[] = {
1.01024, -0.0121636, -0.0158122, -0.989479, -0.0276114, 0.0607747, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.63803e-13, 
0.0, 1.01078, -0.0478881, 0.0156894, -0.989281, 0.917086, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.90314e-13, 
0.0, 0.0, 1.1468, -0.0129878, -0.0416909, -0.832855, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.64856e-14, 
0.0, 0.0, 0.0, 1.01021, -0.0122165, -0.000600698, -0.989891, 0.00322206, 0.0166098, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.38963e-13, 
0.0, 0.0, 0.0, 0.0, 1.00929, -0.0528922, -0.0152066, -0.990747, 0.999325, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.87282e-13, 
0.0, 0.0, 0.0, 0.0, 0.0, 1.14913, -0.00121739, -0.0456005, -0.824219, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.21938e-14, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.00989, -0.011815, 0.0104697, -0.00860113, 0.990167, -0.903063, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.23236e-15, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.00806, -0.0553401, -0.992065, 0.00298854, -0.0602893, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.34652e-15, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.15574, -0.0474249, -0.00882668, -0.859952, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.05303e-17, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.00672, 0.0109889, -0.0499998, -0.99286, 0.0303817, -0.0890087, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -6.66014e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.00964, -0.0239741, -0.0194876, -0.990324, 1.07499, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.65712e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.12805, -0.0444217, -0.0197004, -0.867583, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.72805e-17, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.00592, 0.0099318, -0.0492734, -0.993865, -0.0222367, 0.0927671, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.77928e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.00889, -0.0355963, 0.0319553, -0.990726, 1.0173, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.43051e-15, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.12426, -0.0425467, -0.0323429, -0.853196, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.13916e-17, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.00469, 0.00814438, -0.0467066, -0.995234, -0.0140124, 0.024117, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.84485e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.00839, -0.045991, 0.0219991, -0.991465, 1.06156, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -7.7213e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.135, -0.0400636, -0.0407515, -0.837049, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -5.66438e-17, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0037, 0.00621002, -0.0418861, -0.995927, -0.0276201, 0.0967563, 0.0, 0.0, 0.0, 2.32682e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.00752, -0.0515262, 0.033654, -0.991982, 1.01737, 0.0, 0.0, 0.0, -7.01533e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.14646, -0.0348739, -0.0455927, -0.822994, 0.0, 0.0, 0.0, -1.00235e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.00289, 0.00427418, -0.0351205, 0.000203087, 0.997121, -0.956412, 6.74071e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.00654, -0.0530632, -0.993508, -0.00403185, -0.0355152, 1.01087e-15, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.15567, -0.0456114, 0.0301172, -0.895997, -7.25471e-18, 
};
const int step83_node8_num_blks = 1;
int step83_node8_A_blk_start[] = {0, };
int step83_node8_B_blk_start[] = {0, };
int step83_node8_blk_width[] = {3, };
const float step83_node8_H_correct_cond = 204.86463078069744;
float step83_node8_H_correct_data[] = {
1.0205848576, -0.012288155264000001, -0.015974116928, -0.9996112649600001, -0.027894140736000002, 0.061397032928, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.6548034272e-13, 
-0.012288155264000001, 1.02182416156496, -0.04821200044208, 0.0278941584964, -0.99960959515496, 0.92623294793908, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.903731507492e-13, 
-0.015974116928, -0.04821200044208, 1.31769353579045, 9.524765999625454e-08, 2.60315179996461e-07, -0.99999660178794, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.013782442e-14, 
-0.9996112649600001, 0.0278941584964, 9.524765999625454e-08, 2.0000077757622, -8.225477999800251e-08, -0.03553663725048, -0.99999778711, 0.0032549572326, 0.016779386058, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.9185771611079976e-14, 
-0.027894140736000002, -0.99960959515496, 2.60315179996461e-07, -8.225477999800251e-08, 1.99999296448602, -0.9275875853089629, -0.0032548659125000014, -0.9999904019259901, 1.0084058156283, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.993763598739996e-15, 
0.061397032928, 0.92623294793908, -0.99999660178794, -0.03553663725048, -0.9275875853089629, 2.8616854491400168, -3.2982620000688727e-09, -4.957659787996791e-08, -1.0000012547086403, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.8110779839752602e-13, 
0.0, 0.0, 0.0, -0.99999778711, -0.0032548659125000014, -3.2982620000688727e-09, 1.9999947267029718, 6.837743500223197e-08, -0.020061585775390002, -0.0086861951757, 0.99995975163, -0.9119942930699999, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.4168988594478202e-13, 
0.0, 0.0, 0.0, 0.0032549572326, -0.9999904019259901, -4.957659787996791e-08, 6.837743500223197e-08, 1.9999939631048935, -1.0083497702048119, -0.9999594215490499, -0.008686195472600002, -0.05010554241299999, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.8519745262252002e-13, 
0.0, 0.0, 0.0, 0.016779386058, 1.0084058156283, -1.0000012547086403, -0.020061585775390002, -1.0083497702048119, 3.0171703899281397, 7.112973900179487e-08, 1.8194246000553225e-08, -0.9999993072801701, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.6286529934317803e-13, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0086861951757, -0.9999594215490499, 7.112973900179487e-08, 2.000001223202287, -3.119947800004458e-08, 0.05802560561449, -0.9995320192, 0.030585865024000004, -0.089606838464, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.557197085776698e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.99995975163, -0.008686195472600002, 1.8194246000553225e-08, -3.119947800004458e-08, 2.000011215063364, -0.911527490512482, -0.030585899718000003, -0.9995368618968701, 1.08437479589657, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.982173992063004e-15, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.9119942930699999, -0.05010554241299999, -0.9999993072801701, 0.05802560561449, -0.911527490512482, 2.83424656393834, 1.0041415999348273e-07, 1.1464740000275693e-08, -0.99999850371074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0186128551854e-15, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.9995320192, -0.030585899718000003, 1.0041415999348273e-07, 1.9999990799846499, -1.5164919999881347e-08, 0.0563986159811, -0.9997486807999999, -0.022368341264, 0.093316281232, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5302297338156499e-15, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.030585865024000004, -0.9995368618968701, 1.1464740000275693e-08, -1.5164919999881347e-08, 2.0000104511822903, -1.08660302490871, 0.022368514210000003, -0.9997544045970601, 1.02726514128378, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.2127050040592003e-15, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.089606838464, 1.08437479589657, -0.99999850371074, 0.0563986159811, -1.08660302490871, 3.18388182278594, 6.430360999815197e-08, 2.8973579998504822e-08, -0.99999720137514, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.291181370471e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.9997486807999999, 0.022368514210000003, 6.430360999815197e-08, 2.0000009972039803, -7.896669996313448e-09, -0.07031482685230002, -0.99990164746, -0.014078118156, 0.02423010873, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.80779766456072e-15, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.022368341264, -0.9997544045970601, 2.8973579998504822e-08, -7.896669996313448e-09, 1.999995264108884, -1.0290908218530779, 0.014078108564079999, -0.999897513660312, 1.07066290641246, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.112176096784599e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.093316281232, 1.02726514128378, -0.99999720137514, -0.07031482685230002, -1.0290908218530779, 3.0639701178229704, 4.973630000198819e-08, -1.4123160006673533e-08, -0.9999992440322001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.3795346110538001e-15, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.99990164746, 0.014078108564079999, 4.973630000198819e-08, 1.99999345720177, -7.191050000333322e-08, -0.009154576035600003, -0.9996119299, -0.02772229437, 0.09711429831000001, 0.0, 0.0, 0.0, 1.0990956573526799e-15, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.014078118156, -0.999897513660312, -1.4123160006673533e-08, -7.191050000333322e-08, 1.9999949930794105, -1.0709003106700221, 0.02772235149146, -0.999613226013402, 1.025621480958126, 0.0, 0.0, 0.0, 7.487837959334005e-17, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.02423010873, 1.07066290641246, -0.9999992440322001, -0.009154576035600003, -1.0709003106700221, 3.1469222179496503, -9.621410000011176e-08, -2.4564299000180287e-07, -1.00000365539143, 0.0, 0.0, 0.0, -8.820940413544001e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.9996119299, 0.02772235149146, -9.621410000011176e-08, 2.00000772204621, -3.1254570003160266e-08, -0.0686446293985, 0.00020367392143, 1.0000026796900001, -0.9591760306800001, 4.2417097276050003e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.02772229437, -0.999613226013402, -2.4564299000180287e-07, -3.1254570003160266e-08, 2.0000108927559723, -1.0279229721445202, -1.0000046742896063, 0.0002036563367800004, -0.03983534645015999, 1.71441358316908e-15, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.09711429831000001, 1.025621480958126, -1.00000365539143, -0.0686446293985, -1.0279229721445202, 3.06134492414018, -1.454493835048821e-07, 9.930641999952032e-08, -1.00000263518336, -6.944101334686001e-16, 
};
float step83_node8_M_correct_data[] = {
1.01024, -0.0121636, -0.0158122, -0.989479, -0.0276114, 0.0607747, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.63803e-13, 
0.0, 1.01078, -0.0478881, 0.0156894, -0.989281, 0.917086, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.90314e-13, 
0.0, 0.0, 1.1468, -0.0129878, -0.0416909, -0.832855, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.64856e-14, 
0.0, 0.0, 0.0, 1.01021, -0.0122165, -0.000600698, -0.989891, 0.00322206, 0.0166098, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.38963e-13, 
0.0, 0.0, 0.0, 0.0, 1.00929, -0.0528922, -0.0152066, -0.990747, 0.999325, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.87282e-13, 
0.0, 0.0, 0.0, 0.0, 0.0, 1.14913, -0.00121739, -0.0456005, -0.824219, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.21938e-14, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.00989, -0.011815, 0.0104697, -0.00860113, 0.990167, -0.903063, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.23236e-15, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.00806, -0.0553401, -0.992065, 0.00298854, -0.0602893, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.34652e-15, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.15574, -0.0474249, -0.00882668, -0.859952, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.05303e-17, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.00672, 0.0109889, -0.0499998, -0.99286, 0.0303817, -0.0890087, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -6.66014e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.00964, -0.0239741, -0.0194876, -0.990324, 1.07499, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.65712e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.12805, -0.0444217, -0.0197004, -0.867583, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.72805e-17, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.00592, 0.0099318, -0.0492734, -0.993865, -0.0222367, 0.0927671, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.77928e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.00889, -0.0355963, 0.0319553, -0.990726, 1.0173, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.43051e-15, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.12426, -0.0425467, -0.0323429, -0.853196, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.13916e-17, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.00469, 0.00814438, -0.0467066, -0.995234, -0.0140124, 0.024117, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.84485e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.00839, -0.045991, 0.0219991, -0.991465, 1.06156, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -7.7213e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.135, -0.0400636, -0.0407515, -0.837049, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -5.66438e-17, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0037, 0.00621002, -0.0418861, -0.995927, -0.0276201, 0.0967563, 0.0, 0.0, 0.0, 2.32682e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.00752, -0.0515262, 0.033654, -0.991982, 1.01737, 0.0, 0.0, 0.0, -7.01533e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.14646, -0.0348739, -0.0455927, -0.822994, 0.0, 0.0, 0.0, -1.00235e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.00289, 0.00427418, -0.0351205, 0.000203087, 0.997121, -0.956412, 6.74071e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.00654, -0.0530632, -0.993508, -0.00403185, -0.0355152, 1.01087e-15, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.15567, -0.0456114, 0.0301172, -0.895997, -7.25471e-18, 
};


const int step83_node9_num_factors = 5;
const bool step83_node9_marked = true;
const bool step83_node9_fixed = false;
int step83_node9_factor_height[] = {step83_factor354_height, step83_factor356_height, step83_factor357_height, step83_factor358_height, step83_factor359_height, };
int step83_node9_factor_width[] = {step83_factor354_width, step83_factor356_width, step83_factor357_width, step83_factor358_width, step83_factor359_width, };
int* step83_node9_factor_ridx[] = {step83_factor354_ridx, step83_factor356_ridx, step83_factor357_ridx, step83_factor358_ridx, step83_factor359_ridx, };
float* step83_node9_factor_data[] = {step83_factor354_data, step83_factor356_data, step83_factor357_data, step83_factor358_data, step83_factor359_data, };
int step83_node9_factor_num_blks[] = {step83_factor354_num_blks, step83_factor356_num_blks, step83_factor357_num_blks, step83_factor358_num_blks, step83_factor359_num_blks, };
int* step83_node9_factor_A_blk_start[] = {step83_factor354_A_blk_start, step83_factor356_A_blk_start, step83_factor357_A_blk_start, step83_factor358_A_blk_start, step83_factor359_A_blk_start, };
int* step83_node9_factor_B_blk_start[] = {step83_factor354_B_blk_start, step83_factor356_B_blk_start, step83_factor357_B_blk_start, step83_factor358_B_blk_start, step83_factor359_B_blk_start, };
int* step83_node9_factor_blk_width[] = {step83_factor354_blk_width, step83_factor356_blk_width, step83_factor357_blk_width, step83_factor358_blk_width, step83_factor359_blk_width, };
const int step83_node9_parent = 10;
const int step83_node9_height = 10;
const int step83_node9_width = 9;
float step83_node9_data[] = {
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
const int step83_node9_num_blks = 0;
int step83_node9_A_blk_start[] = {};
int step83_node9_B_blk_start[] = {};
int step83_node9_blk_width[] = {};
const float step83_node9_H_correct_cond = 5579.562784623637;
float step83_node9_H_correct_data[] = {
1.0108693763999999, -0.0028345001724, -0.036318384492, -0.99999374826, -0.0042714665448, 0.062225745426, 0.0, 0.0, 0.0, -2.1726724032e-16, 
-0.0028345001724, 1.0048337561014085, 0.021319248580372, 0.00427145271686, -0.9999942435329832, 1.0416603029284341, 0.0, 0.0, 0.0, 7.482936266212e-17, 
-0.036318384492, 0.021319248580372, 1.2028154236197202, -1.4617191200376537e-07, 1.1320794400245176e-07, -1.0000040127587801, 0.0, 0.0, 0.0, -2.4694715185520002e-17, 
-0.99999374826, 0.00427145271686, -1.4617191200376537e-07, 1.9999971744162062, 4.181700000399326e-09, -0.05777552514650001, -0.9999837049099999, -0.0055218222023999996, -0.05331437062299999, 4.44412734724424e-16, 
-0.0042714665448, -0.9999942435329832, 1.1320794400245176e-07, 4.181700000399326e-09, 2.0000030221846665, -1.041912920057898, 0.00552184464971, -0.9999854746755256, 0.937226992182063, -2.9742349535219e-16, 
0.062225745426, 1.0416603029284341, -1.0000040127587801, -0.05777552514650001, -1.041912920057898, 3.0889164995774903, 4.4786839994187764e-08, 3.132782399871471e-08, -0.99999775430682, 6.42686847024e-17, 
0.0, 0.0, 0.0, -0.9999837049099999, 0.00552184464971, 4.4786839994187764e-08, 1.00000170194088, -7.566160000512264e-09, 0.058488861389100005, -2.2326717120768e-16, 
0.0, 0.0, 0.0, -0.0055218222023999996, -0.9999854746755256, 3.132782399871471e-08, -7.566160000512264e-09, 1.0000002962160783, -0.936917052971032, 2.208150132858e-16, 
0.0, 0.0, 0.0, -0.05331437062299999, 0.937226992182063, -0.99999775430682, 0.058488861389100005, -0.936917052971032, 1.88123538526625, -2.1124414144408601e-16, 
};
float step83_node9_M_correct_data[] = {
1.00542, -0.00281922, -0.0361226, -0.994603, -0.00424844, 0.0618903, 0.0, 0.0, 0.0, -2.16096e-16, 
0.0, 1.00241, 0.0211664, 0.00146392, -0.997602, 1.03933, 0.0, 0.0, 0.0, 7.40417e-17, 
0.0, 0.0, 1.09593, -0.0328112, 0.0191274, -0.930504, 0.0, 0.0, 0.0, -3.10858e-17, 
0.0, 0.0, 0.0, 1.00483, -0.00212723, -0.0281358, -0.995177, -0.00549528, -0.0530581, 2.27257e-16, 
0.0, 0.0, 0.0, 0.0, 1.0022, 0.0128976, 0.0033974, -0.997802, 0.935057, -2.22909e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 1.06681, -0.0262876, 0.0119184, -0.950076, -1.77794e-17, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0944556, -0.0186918, -0.23784, 3.37058e-17, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.062208, 0.0429342, 7.82427e-18, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.207444, 1.69595e-19, 
};


const int step83_node10_num_factors = 0;
const bool step83_node10_marked = true;
const bool step83_node10_fixed = false;
int step83_node10_factor_height[] = {};
int step83_node10_factor_width[] = {};
int* step83_node10_factor_ridx[] = {};
float* step83_node10_factor_data[] = {};
int step83_node10_factor_num_blks[] = {};
int* step83_node10_factor_A_blk_start[] = {};
int* step83_node10_factor_B_blk_start[] = {};
int* step83_node10_factor_blk_width[] = {};
const int step83_node10_parent = -1;
const int step83_node10_height = 1;
const int step83_node10_width = 1;
float step83_node10_data[] = {
0, 
};
const int step83_node10_num_blks = 0;
int step83_node10_A_blk_start[] = {};
int step83_node10_B_blk_start[] = {};
int step83_node10_blk_width[] = {};
const float step83_node10_H_correct_cond = 1.0;
float step83_node10_H_correct_data[] = {
2.4334128035999997e-62, 
};
float step83_node10_M_correct_data[] = {
-1.55994e-31, 
};


const int step83_nnodes = 11;
bool step83_node_marked[] = {false, false, false, false, false, false, false, false, step83_node8_marked, step83_node9_marked, step83_node10_marked, };
bool step83_node_fixed[] = {false, false, false, false, false, false, false, false, step83_node8_fixed, step83_node9_fixed, step83_node10_fixed, };
int step83_node_num_factors[] = {0, 0, 0, 0, 0, 0, 0, 0, step83_node8_num_factors, step83_node9_num_factors, step83_node10_num_factors, };
int* step83_node_factor_height[] = {0, 0, 0, 0, 0, 0, 0, 0, step83_node8_factor_height, step83_node9_factor_height, step83_node10_factor_height, };
int* step83_node_factor_width[] = {0, 0, 0, 0, 0, 0, 0, 0, step83_node8_factor_width, step83_node9_factor_width, step83_node10_factor_width, };
int** step83_node_factor_ridx[] = {0, 0, 0, 0, 0, 0, 0, 0, step83_node8_factor_ridx, step83_node9_factor_ridx, step83_node10_factor_ridx, };
float** step83_node_factor_data[] = {0, 0, 0, 0, 0, 0, 0, 0, step83_node8_factor_data, step83_node9_factor_data, step83_node10_factor_data, };
int* step83_node_factor_num_blks[] = {0, 0, 0, 0, 0, 0, 0, 0, step83_node8_factor_num_blks, step83_node9_factor_num_blks, step83_node10_factor_num_blks, };
int** step83_node_factor_A_blk_start[] = {0, 0, 0, 0, 0, 0, 0, 0, step83_node8_factor_A_blk_start, step83_node9_factor_A_blk_start, step83_node10_factor_A_blk_start, };
int** step83_node_factor_B_blk_start[] = {0, 0, 0, 0, 0, 0, 0, 0, step83_node8_factor_B_blk_start, step83_node9_factor_B_blk_start, step83_node10_factor_B_blk_start, };
int** step83_node_factor_blk_width[] = {0, 0, 0, 0, 0, 0, 0, 0, step83_node8_factor_blk_width, step83_node9_factor_blk_width, step83_node10_factor_blk_width, };
int step83_node_parent[] = {0, 0, 0, 0, 0, 0, 0, 0, step83_node8_parent, step83_node9_parent, step83_node10_parent, };
int step83_node_height[] = {0, 0, 0, 0, 0, 0, 0, 0, step83_node8_height, step83_node9_height, step83_node10_height, };
int step83_node_width[] = {0, 0, 0, 0, 0, 0, 0, 0, step83_node8_width, step83_node9_width, step83_node10_width, };
float* step83_node_data[] = {0, 0, 0, 0, 0, 0, 0, 0, step83_node8_data, step83_node9_data, step83_node10_data, };
int step83_node_num_blks[] = {0, 0, 0, 0, 0, 0, 0, 0, step83_node8_num_blks, step83_node9_num_blks, step83_node10_num_blks, };
int* step83_node_A_blk_start[] = {0, 0, 0, 0, 0, 0, 0, 0, step83_node8_A_blk_start, step83_node9_A_blk_start, step83_node10_A_blk_start, };
int* step83_node_B_blk_start[] = {0, 0, 0, 0, 0, 0, 0, 0, step83_node8_B_blk_start, step83_node9_B_blk_start, step83_node10_B_blk_start, };
int* step83_node_blk_width[] = {0, 0, 0, 0, 0, 0, 0, 0, step83_node8_blk_width, step83_node9_blk_width, step83_node10_blk_width, };
float step83_node_H_correct_cond[] = {0, 0, 0, 0, 0, 0, 0, 0, step83_node8_H_correct_cond, step83_node9_H_correct_cond, step83_node10_H_correct_cond, };
float* step83_node_H_correct_data[] = {0, 0, 0, 0, 0, 0, 0, 0, step83_node8_H_correct_data, step83_node9_H_correct_data, step83_node10_H_correct_data, };
float* step83_node_M_correct_data[] = {0, 0, 0, 0, 0, 0, 0, 0, step83_node8_M_correct_data, step83_node9_M_correct_data, step83_node10_M_correct_data, };


