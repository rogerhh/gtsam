#pragma once

const bool step8_is_reconstruct = true;

const int step8_factor0_height = 4;
const int step8_factor0_width = 3;
int step8_factor0_ridx[] = {0, 1, 2, 27, };
float step8_factor0_data[] = {
0.0010000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0010000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0010000, 0.0000000, 
};

const int step8_factor0_num_blks = 1;
int step8_factor0_A_blk_start[] = {0, };
int step8_factor0_B_blk_start[] = {0, };
int step8_factor0_blk_width[] = {3, };

const int step8_factor1_height = 4;
const int step8_factor1_width = 3;
int step8_factor1_ridx[] = {3, 4, 5, 27, };
float step8_factor1_data[] = {
0.0010000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0010000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0010000, 0.0000000, 
};

const int step8_factor1_num_blks = 1;
int step8_factor1_A_blk_start[] = {0, };
int step8_factor1_B_blk_start[] = {3, };
int step8_factor1_blk_width[] = {3, };

const int step8_factor2_height = 4;
const int step8_factor2_width = 3;
int step8_factor2_ridx[] = {6, 7, 8, 27, };
float step8_factor2_data[] = {
0.0010000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0010000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0010000, 0.0000000, 
};

const int step8_factor2_num_blks = 1;
int step8_factor2_A_blk_start[] = {0, };
int step8_factor2_B_blk_start[] = {6, };
int step8_factor2_blk_width[] = {3, };

const int step8_factor3_height = 4;
const int step8_factor3_width = 3;
int step8_factor3_ridx[] = {0, 1, 2, 27, };
float step8_factor3_data[] = {
1.0000000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 1.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 1.0000000, 0.0000000, 
};

const int step8_factor3_num_blks = 1;
int step8_factor3_A_blk_start[] = {0, };
int step8_factor3_B_blk_start[] = {0, };
int step8_factor3_blk_width[] = {3, };

const int step8_factor4_height = 7;
const int step8_factor4_width = 3;
int step8_factor4_ridx[] = {0, 1, 2, 3, 4, 5, 27, };
float step8_factor4_data[] = {
1.0000000, 0.0000000, 0.0000000, -0.9999670, 0.0081837, 0.0335822, -0.0000000, 
0.0000000, 1.0000000, 0.0000000, -0.0081837, -0.9999670, 0.9990980, -0.0000000, 
0.0000000, 0.0000000, 1.0000000, -0.0000000, -0.0000000, -1.0000000, 0.0000000, 
};

const int step8_factor4_num_blks = 1;
int step8_factor4_A_blk_start[] = {0, };
int step8_factor4_B_blk_start[] = {0, };
int step8_factor4_blk_width[] = {6, };

const int step8_factor5_height = 7;
const int step8_factor5_width = 3;
int step8_factor5_ridx[] = {3, 4, 5, 6, 7, 8, 27, };
float step8_factor5_data[] = {
1.0000000, 0.0000000, 0.0000000, -0.9999840, 0.0056968, 0.0178761, -0.0000000, 
0.0000000, 1.0000000, 0.0000000, -0.0056968, -0.9999840, 1.0034800, 0.0000000, 
0.0000000, 0.0000000, 1.0000000, -0.0000000, -0.0000000, -1.0000000, -0.0000000, 
};

const int step8_factor5_num_blks = 1;
int step8_factor5_A_blk_start[] = {0, };
int step8_factor5_B_blk_start[] = {3, };
int step8_factor5_blk_width[] = {6, };

const int step8_factor6_height = 4;
const int step8_factor6_width = 3;
int step8_factor6_ridx[] = {9, 10, 11, 27, };
float step8_factor6_data[] = {
0.0010000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0010000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0010000, 0.0000000, 
};

const int step8_factor6_num_blks = 1;
int step8_factor6_A_blk_start[] = {0, };
int step8_factor6_B_blk_start[] = {9, };
int step8_factor6_blk_width[] = {3, };

const int step8_factor7_height = 7;
const int step8_factor7_width = 3;
int step8_factor7_ridx[] = {6, 7, 8, 9, 10, 11, 27, };
float step8_factor7_data[] = {
1.0000000, 0.0000000, 0.0000000, -0.9994120, 0.0342925, 0.0169251, 0.0000000, 
0.0000000, 1.0000000, 0.0000000, -0.0342925, -0.9994120, 0.9733340, -0.0000000, 
0.0000000, 0.0000000, 1.0000000, -0.0000000, -0.0000000, -1.0000000, 0.0000000, 
};

const int step8_factor7_num_blks = 1;
int step8_factor7_A_blk_start[] = {0, };
int step8_factor7_B_blk_start[] = {6, };
int step8_factor7_blk_width[] = {6, };

const int step8_factor8_height = 4;
const int step8_factor8_width = 3;
int step8_factor8_ridx[] = {12, 13, 14, 27, };
float step8_factor8_data[] = {
0.0010000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0010000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0010000, 0.0000000, 
};

const int step8_factor8_num_blks = 1;
int step8_factor8_A_blk_start[] = {0, };
int step8_factor8_B_blk_start[] = {12, };
int step8_factor8_blk_width[] = {3, };

const int step8_factor9_height = 7;
const int step8_factor9_width = 3;
int step8_factor9_ridx[] = {9, 10, 11, 12, 13, 14, 27, };
float step8_factor9_data[] = {
1.0000000, 0.0000000, 0.0000000, -0.9999320, 0.0116626, -0.0030360, -0.0000000, 
0.0000000, 1.0000000, 0.0000000, -0.0116626, -0.9999320, 1.0380501, -0.0000000, 
0.0000000, 0.0000000, 1.0000000, -0.0000000, -0.0000000, -1.0000000, 0.0000000, 
};

const int step8_factor9_num_blks = 1;
int step8_factor9_A_blk_start[] = {0, };
int step8_factor9_B_blk_start[] = {9, };
int step8_factor9_blk_width[] = {6, };

const int step8_factor10_height = 4;
const int step8_factor10_width = 3;
int step8_factor10_ridx[] = {15, 16, 17, 27, };
float step8_factor10_data[] = {
0.0010000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0010000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0010000, 0.0000000, 
};

const int step8_factor10_num_blks = 1;
int step8_factor10_A_blk_start[] = {0, };
int step8_factor10_B_blk_start[] = {15, };
int step8_factor10_blk_width[] = {3, };

const int step8_factor11_height = 7;
const int step8_factor11_width = 3;
int step8_factor11_ridx[] = {12, 13, 14, 15, 16, 17, 27, };
float step8_factor11_data[] = {
1.0000000, 0.0000000, 0.0000000, -0.9999580, -0.0091544, 0.0613274, 0.0000000, 
0.0000000, 1.0000000, 0.0000000, 0.0091544, -0.9999580, 0.9927050, 0.0000000, 
0.0000000, 0.0000000, 1.0000000, -0.0000000, -0.0000000, -1.0000000, -0.0000000, 
};

const int step8_factor11_num_blks = 1;
int step8_factor11_A_blk_start[] = {0, };
int step8_factor11_B_blk_start[] = {12, };
int step8_factor11_blk_width[] = {6, };

const int step8_factor12_height = 4;
const int step8_factor12_width = 3;
int step8_factor12_ridx[] = {18, 19, 20, 27, };
float step8_factor12_data[] = {
0.0010000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0010000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0010000, 0.0000000, 
};

const int step8_factor12_num_blks = 1;
int step8_factor12_A_blk_start[] = {0, };
int step8_factor12_B_blk_start[] = {18, };
int step8_factor12_blk_width[] = {3, };

const int step8_factor13_height = 7;
const int step8_factor13_width = 3;
int step8_factor13_ridx[] = {15, 16, 17, 18, 19, 20, 27, };
float step8_factor13_data[] = {
1.0000000, 0.0000000, 0.0000000, -0.0182853, 0.9998330, -1.1123500, 0.0000000, 
0.0000000, 1.0000000, 0.0000000, -0.9998330, -0.0182853, 0.0291199, -0.0000000, 
0.0000000, 0.0000000, 1.0000000, -0.0000000, -0.0000000, -1.0000000, -0.0000000, 
};

const int step8_factor13_num_blks = 1;
int step8_factor13_A_blk_start[] = {0, };
int step8_factor13_B_blk_start[] = {15, };
int step8_factor13_blk_width[] = {6, };

const int step8_factor14_height = 4;
const int step8_factor14_width = 3;
int step8_factor14_ridx[] = {21, 22, 23, 27, };
float step8_factor14_data[] = {
0.0010000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0010000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0010000, 0.0000000, 
};

const int step8_factor14_num_blks = 1;
int step8_factor14_A_blk_start[] = {0, };
int step8_factor14_B_blk_start[] = {21, };
int step8_factor14_blk_width[] = {3, };

const int step8_factor15_height = 7;
const int step8_factor15_width = 3;
int step8_factor15_ridx[] = {18, 19, 20, 21, 22, 23, 27, };
float step8_factor15_data[] = {
1.0000000, 0.0000000, 0.0000000, -0.9999300, 0.0118554, 0.0343397, 0.0000000, 
0.0000000, 1.0000000, 0.0000000, -0.0118554, -0.9999300, 1.0272900, 0.0000000, 
0.0000000, 0.0000000, 1.0000000, -0.0000000, -0.0000000, -1.0000000, 0.0000000, 
};

const int step8_factor15_num_blks = 1;
int step8_factor15_A_blk_start[] = {0, };
int step8_factor15_B_blk_start[] = {18, };
int step8_factor15_blk_width[] = {6, };

const int step8_factor16_height = 4;
const int step8_factor16_width = 3;
int step8_factor16_ridx[] = {24, 25, 26, 27, };
float step8_factor16_data[] = {
0.0010000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0010000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0010000, 0.0000000, 
};

const int step8_factor16_num_blks = 1;
int step8_factor16_A_blk_start[] = {0, };
int step8_factor16_B_blk_start[] = {24, };
int step8_factor16_blk_width[] = {3, };

const int step8_factor17_height = 7;
const int step8_factor17_width = 3;
int step8_factor17_ridx[] = {21, 22, 23, 24, 25, 26, 27, };
float step8_factor17_data[] = {
1.0000000, 0.0000000, 0.0000000, -0.9994110, 0.0343158, -0.0023019, -0.0000000, 
0.0000000, 1.0000000, 0.0000000, -0.0343158, -0.9994110, 0.9178200, 0.0000000, 
0.0000000, 0.0000000, 1.0000000, -0.0000000, -0.0000000, -1.0000000, -0.0000000, 
};

const int step8_factor17_num_blks = 1;
int step8_factor17_A_blk_start[] = {0, };
int step8_factor17_B_blk_start[] = {21, };
int step8_factor17_blk_width[] = {6, };

int step8_factor_max_num_blks = 1;

int step8_factor_max_height = 7;

const int step8_num_factors = 18;
const int step8_factor_height[] = {step8_factor0_height, step8_factor1_height, step8_factor2_height, step8_factor3_height, step8_factor4_height, step8_factor5_height, step8_factor6_height, step8_factor7_height, step8_factor8_height, step8_factor9_height, step8_factor10_height, step8_factor11_height, step8_factor12_height, step8_factor13_height, step8_factor14_height, step8_factor15_height, step8_factor16_height, step8_factor17_height, };
const int step8_factor_width[] = {step8_factor0_width, step8_factor1_width, step8_factor2_width, step8_factor3_width, step8_factor4_width, step8_factor5_width, step8_factor6_width, step8_factor7_width, step8_factor8_width, step8_factor9_width, step8_factor10_width, step8_factor11_width, step8_factor12_width, step8_factor13_width, step8_factor14_width, step8_factor15_width, step8_factor16_width, step8_factor17_width, };
int* step8_factor_ridx[] = {step8_factor0_ridx, step8_factor1_ridx, step8_factor2_ridx, step8_factor3_ridx, step8_factor4_ridx, step8_factor5_ridx, step8_factor6_ridx, step8_factor7_ridx, step8_factor8_ridx, step8_factor9_ridx, step8_factor10_ridx, step8_factor11_ridx, step8_factor12_ridx, step8_factor13_ridx, step8_factor14_ridx, step8_factor15_ridx, step8_factor16_ridx, step8_factor17_ridx, };
float* step8_factor_data[] = {step8_factor0_data, step8_factor1_data, step8_factor2_data, step8_factor3_data, step8_factor4_data, step8_factor5_data, step8_factor6_data, step8_factor7_data, step8_factor8_data, step8_factor9_data, step8_factor10_data, step8_factor11_data, step8_factor12_data, step8_factor13_data, step8_factor14_data, step8_factor15_data, step8_factor16_data, step8_factor17_data, };
int step8_factor_num_blks[] = {step8_factor0_num_blks, step8_factor1_num_blks, step8_factor2_num_blks, step8_factor3_num_blks, step8_factor4_num_blks, step8_factor5_num_blks, step8_factor6_num_blks, step8_factor7_num_blks, step8_factor8_num_blks, step8_factor9_num_blks, step8_factor10_num_blks, step8_factor11_num_blks, step8_factor12_num_blks, step8_factor13_num_blks, step8_factor14_num_blks, step8_factor15_num_blks, step8_factor16_num_blks, step8_factor17_num_blks, };
int* step8_factor_A_blk_start[] = {step8_factor0_A_blk_start, step8_factor1_A_blk_start, step8_factor2_A_blk_start, step8_factor3_A_blk_start, step8_factor4_A_blk_start, step8_factor5_A_blk_start, step8_factor6_A_blk_start, step8_factor7_A_blk_start, step8_factor8_A_blk_start, step8_factor9_A_blk_start, step8_factor10_A_blk_start, step8_factor11_A_blk_start, step8_factor12_A_blk_start, step8_factor13_A_blk_start, step8_factor14_A_blk_start, step8_factor15_A_blk_start, step8_factor16_A_blk_start, step8_factor17_A_blk_start, };
int* step8_factor_B_blk_start[] = {step8_factor0_B_blk_start, step8_factor1_B_blk_start, step8_factor2_B_blk_start, step8_factor3_B_blk_start, step8_factor4_B_blk_start, step8_factor5_B_blk_start, step8_factor6_B_blk_start, step8_factor7_B_blk_start, step8_factor8_B_blk_start, step8_factor9_B_blk_start, step8_factor10_B_blk_start, step8_factor11_B_blk_start, step8_factor12_B_blk_start, step8_factor13_B_blk_start, step8_factor14_B_blk_start, step8_factor15_B_blk_start, step8_factor16_B_blk_start, step8_factor17_B_blk_start, };
int* step8_factor_blk_width[] = {step8_factor0_blk_width, step8_factor1_blk_width, step8_factor2_blk_width, step8_factor3_blk_width, step8_factor4_blk_width, step8_factor5_blk_width, step8_factor6_blk_width, step8_factor7_blk_width, step8_factor8_blk_width, step8_factor9_blk_width, step8_factor10_blk_width, step8_factor11_blk_width, step8_factor12_blk_width, step8_factor13_blk_width, step8_factor14_blk_width, step8_factor15_blk_width, step8_factor16_blk_width, step8_factor17_blk_width, };

const int step8_vio_node_height = 28;
const int step8_vio_node_width = 27;
const int step8_vio_marginal_width = 0;

const float step8_H_correct_cond = 3403.412109375;
float step8_H_correct_data[] = {
2.0000010, 0.0000000, 0.0000000, -0.9999670, 0.0081837, 0.0335822, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.0000000, 
0.0000000, 2.0000010, 0.0000000, -0.0081837, -0.9999670, 0.9990980, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 2.0000010, 0.0000000, 0.0000000, -1.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
-0.9999670, -0.0081837, 0.0000000, 2.0000019, 0.0000000, -0.0417574, -0.9999840, 0.0056968, 0.0178761, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.0000000, 
0.0081837, -0.9999670, 0.0000000, 0.0000000, 2.0000019, -0.9987902, -0.0056968, -0.9999840, 1.0034800, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
0.0335822, 0.9990980, -1.0000000, -0.0417574, -0.9987902, 2.9993255, 0.0000000, 0.0000000, -1.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.0000000, 
0.0000000, 0.0000000, 0.0000000, -0.9999840, -0.0056968, 0.0000000, 2.0000014, 0.0000000, -0.0235924, -0.9994120, 0.0342925, 0.0169251, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0056968, -0.9999840, 0.0000000, 0.0000000, 2.0000014, -1.0033621, -0.0342925, -0.9994120, 0.9733340, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0178761, 1.0034800, -1.0000000, -0.0235924, -1.0033621, 3.0072925, 0.0000000, 0.0000000, -1.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.9994120, -0.0342925, 0.0000000, 2.0000014, 0.0000000, -0.0502932, -0.9999320, 0.0116626, -0.0030360, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0342925, -0.9994120, 0.0000000, 0.0000000, 2.0000014, -0.9721813, -0.0116626, -0.9999320, 1.0380501, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0169251, 0.9733340, -1.0000000, -0.0502932, -0.9721813, 2.9476666, 0.0000000, 0.0000000, -1.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.9999320, -0.0116626, 0.0000000, 2.0000010, -0.0000000, -0.0090705, -0.9999580, -0.0091544, 0.0613274, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0116626, -0.9999320, 0.0000000, -0.0000000, 2.0000010, -1.0380148, 0.0091544, -0.9999580, 0.9927050, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.0030360, 1.0380501, -1.0000000, -0.0090705, -1.0380148, 3.0775580, 0.0000000, 0.0000000, -1.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.9999580, 0.0091544, 0.0000000, 2.0000010, -0.0000000, -0.0522372, -0.0182853, 0.9998330, -1.1123500, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.0091544, -0.9999580, 0.0000000, -0.0000000, 2.0000010, -0.9932247, -0.9998330, -0.0182853, 0.0291199, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0613274, 0.9927050, -1.0000000, -0.0522372, -0.9932247, 2.9892254, 0.0000000, 0.0000000, -1.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.0182853, -0.9998330, 0.0000000, 2.0000014, 0.0000000, -0.0087754, -0.9999300, 0.0118554, 0.0343397, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.9998330, -0.0182853, 0.0000000, 0.0000000, 2.0000014, -1.1126966, -0.0118554, -0.9999300, 1.0272900, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -1.1123500, 0.0291199, -1.0000000, -0.0087754, -1.1126966, 3.2381713, 0.0000000, 0.0000000, -1.0000000, 0.0000000, 0.0000000, 0.0000000, -0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.9999300, -0.0118554, 0.0000000, 2.0000014, -0.0000000, -0.0465162, -0.9994110, 0.0343158, -0.0023019, -0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0118554, -0.9999300, 0.0000000, -0.0000000, 2.0000014, -1.0268110, -0.0343158, -0.9994110, 0.9178200, -0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0343397, 1.0272900, -1.0000000, -0.0465162, -1.0268110, 3.0565050, 0.0000000, 0.0000000, -1.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.9994110, -0.0343158, 0.0000000, 1.0000010, -0.0000000, -0.0291952, 0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0343158, -0.9994110, 0.0000000, -0.0000000, 1.0000010, -0.9173584, -0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.0023019, 0.9178200, -1.0000000, -0.0291952, -0.9173584, 1.8423998, 0.0000000, 
};

float step8_M_correct_data[] = {
1.4142139, 0.0000000, 0.0000000, -0.7070833, 0.0057868, 0.0237462, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.0000000, 
0.0000000, 1.4142139, 0.0000000, -0.0057868, -0.7070833, 0.7064688, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 1.4142139, 0.0000000, 0.0000000, -0.7071066, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0000000, 1.2247455, 0.0000000, -0.0170474, -0.8164831, 0.0046514, 0.0145958, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.2247455, -0.4077544, -0.0046514, -0.8164831, 0.8193375, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.3539237, -0.0116813, -0.2458375, -0.4916545, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.1546422, -0.0024871, -0.0117849, -0.8655599, 0.0296997, 0.0146583, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.1282256, -0.4035994, -0.0323031, -0.8857610, 0.8627446, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.3896079, -0.0167227, -0.2570096, -0.4689266, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.1178033, -0.0064447, -0.0157254, -0.8945509, 0.0104335, -0.0027161, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.0716583, -0.3070486, -0.0162624, -0.9330071, 0.9686227, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.3743018, -0.0138693, -0.2083345, -0.5112621, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.0951360, -0.0079708, -0.0025923, -0.9130902, -0.0083591, 0.0559998, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.0420763, -0.2310663, 0.0018006, -0.9596462, 0.9530505, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.3507528, -0.0014443, -0.1641778, -0.5771870, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.0799360, -0.0056872, -0.0033836, -0.0169318, 0.9258261, -1.0300146, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.0256820, -0.1686155, -0.9748921, -0.0126939, 0.0226795, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.3100369, -0.1255226, 0.0007574, -0.7630785, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.0166337, 0.0033402, -0.0982546, -0.9835696, 0.0116614, 0.0337778, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.0689596, -0.1477026, -0.0080172, -0.9354600, 0.9609130, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.2501872, -0.0782478, -0.1096029, -0.6836991, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.0131166, -0.0045465, -0.0583225, -0.9864718, 0.0338715, -0.0022721, -0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.0548680, -0.1929244, -0.0367826, -0.9472815, 0.8700706, 0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.2743407, -0.0507162, -0.1418602, -0.6531022, 0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.1514902, -0.0569329, -0.2149048, -0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.2795455, -0.7081636, -0.0000000, 
0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.3333912, -0.0000000, 
};

