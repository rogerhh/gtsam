#pragma once

const bool step8_is_reconstruct = true;

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

const int step8_factor9_height = 7;
const int step8_factor9_width = 3;
int step8_factor9_ridx[] = {9, 10, 11, 12, 13, 14, 27, };
float step8_factor9_data[] = {
1.0000000, 0.0000000, 0.0000000, -0.9999320, 0.0116626, -0.0030360, -0.0000000, 
0.0000000, 1.0000000, 0.0000000, -0.0116626, -0.9999320, 1.0380500, -0.0000000, 
0.0000000, 0.0000000, 1.0000000, -0.0000000, -0.0000000, -1.0000000, 0.0000000, 
};

const int step8_factor9_num_blks = 1;
int step8_factor9_A_blk_start[] = {0, };
int step8_factor9_B_blk_start[] = {9, };
int step8_factor9_blk_width[] = {6, };

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

const int step8_node0_num_factors = 18;
const bool step8_node0_marked = true;
const bool step8_node0_fixed = false;
int step8_node0_factor_height[] = {step8_factor0_height, step8_factor1_height, step8_factor2_height, step8_factor3_height, step8_factor4_height, step8_factor5_height, step8_factor6_height, step8_factor7_height, step8_factor8_height, step8_factor9_height, step8_factor10_height, step8_factor11_height, step8_factor12_height, step8_factor13_height, step8_factor14_height, step8_factor15_height, step8_factor16_height, step8_factor17_height, };
int step8_node0_factor_width[] = {step8_factor0_width, step8_factor1_width, step8_factor2_width, step8_factor3_width, step8_factor4_width, step8_factor5_width, step8_factor6_width, step8_factor7_width, step8_factor8_width, step8_factor9_width, step8_factor10_width, step8_factor11_width, step8_factor12_width, step8_factor13_width, step8_factor14_width, step8_factor15_width, step8_factor16_width, step8_factor17_width, };
int* step8_node0_factor_ridx[] = {step8_factor0_ridx, step8_factor1_ridx, step8_factor2_ridx, step8_factor3_ridx, step8_factor4_ridx, step8_factor5_ridx, step8_factor6_ridx, step8_factor7_ridx, step8_factor8_ridx, step8_factor9_ridx, step8_factor10_ridx, step8_factor11_ridx, step8_factor12_ridx, step8_factor13_ridx, step8_factor14_ridx, step8_factor15_ridx, step8_factor16_ridx, step8_factor17_ridx, };
float* step8_node0_factor_data[] = {step8_factor0_data, step8_factor1_data, step8_factor2_data, step8_factor3_data, step8_factor4_data, step8_factor5_data, step8_factor6_data, step8_factor7_data, step8_factor8_data, step8_factor9_data, step8_factor10_data, step8_factor11_data, step8_factor12_data, step8_factor13_data, step8_factor14_data, step8_factor15_data, step8_factor16_data, step8_factor17_data, };
int step8_node0_factor_num_blks[] = {step8_factor0_num_blks, step8_factor1_num_blks, step8_factor2_num_blks, step8_factor3_num_blks, step8_factor4_num_blks, step8_factor5_num_blks, step8_factor6_num_blks, step8_factor7_num_blks, step8_factor8_num_blks, step8_factor9_num_blks, step8_factor10_num_blks, step8_factor11_num_blks, step8_factor12_num_blks, step8_factor13_num_blks, step8_factor14_num_blks, step8_factor15_num_blks, step8_factor16_num_blks, step8_factor17_num_blks, };
int* step8_node0_factor_A_blk_start[] = {step8_factor0_A_blk_start, step8_factor1_A_blk_start, step8_factor2_A_blk_start, step8_factor3_A_blk_start, step8_factor4_A_blk_start, step8_factor5_A_blk_start, step8_factor6_A_blk_start, step8_factor7_A_blk_start, step8_factor8_A_blk_start, step8_factor9_A_blk_start, step8_factor10_A_blk_start, step8_factor11_A_blk_start, step8_factor12_A_blk_start, step8_factor13_A_blk_start, step8_factor14_A_blk_start, step8_factor15_A_blk_start, step8_factor16_A_blk_start, step8_factor17_A_blk_start, };
int* step8_node0_factor_B_blk_start[] = {step8_factor0_B_blk_start, step8_factor1_B_blk_start, step8_factor2_B_blk_start, step8_factor3_B_blk_start, step8_factor4_B_blk_start, step8_factor5_B_blk_start, step8_factor6_B_blk_start, step8_factor7_B_blk_start, step8_factor8_B_blk_start, step8_factor9_B_blk_start, step8_factor10_B_blk_start, step8_factor11_B_blk_start, step8_factor12_B_blk_start, step8_factor13_B_blk_start, step8_factor14_B_blk_start, step8_factor15_B_blk_start, step8_factor16_B_blk_start, step8_factor17_B_blk_start, };
int* step8_node0_factor_blk_width[] = {step8_factor0_blk_width, step8_factor1_blk_width, step8_factor2_blk_width, step8_factor3_blk_width, step8_factor4_blk_width, step8_factor5_blk_width, step8_factor6_blk_width, step8_factor7_blk_width, step8_factor8_blk_width, step8_factor9_blk_width, step8_factor10_blk_width, step8_factor11_blk_width, step8_factor12_blk_width, step8_factor13_blk_width, step8_factor14_blk_width, step8_factor15_blk_width, step8_factor16_blk_width, step8_factor17_blk_width, };
const int step8_node0_parent = 1;
const int step8_node0_height = 28;
const int step8_node0_width = 27;
float step8_node0_data[] = {
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
};
const int step8_node0_num_blks = 0;
int step8_node0_A_blk_start[] = {};
int step8_node0_B_blk_start[] = {};
int step8_node0_blk_width[] = {};
const float step8_node0_H_correct_cond = 3403.4586031886065;
float step8_node0_H_correct_data[] = {
1.9999899240999999, 0.0, 0.0, -0.99996384943, 0.0081836938596, 0.033582113501999994, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.11021990366e-16, 
0.0, 1.9999899240999999, 0.0, -0.0081836938596, -0.99996384943, 0.9990955244899999, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
0.0, 0.0, 1.9999899240999999, 0.0, 0.0, -0.99999779047, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.364393481e-19, 
-0.99996384943, -0.0081836938596, 0.0, 2.000012417980298, -2.638723878179629e-10, -0.04175750403504, -0.99998755425, 0.0056967899025, 0.01787620605, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.110261131682e-16, 
0.0081836938596, -0.99996384943, 0.0, -2.638723878179629e-10, 2.000012417980298, -0.9987915178630152, -0.005696789726588738, -0.999987554251002, 1.0034817659968553, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.030347585199814e-18, 
0.033582113501999994, 0.9990955244899999, -0.99999779047, -0.04175750403504, -0.9987915178630152, 2.9993159441872, -1.0523740000612605e-08, 1.2840361139644401e-06, -0.99999788938492, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -5.47364860083e-18, 
0.0, 0.0, 0.0, -0.99998755425, -0.005696789726588738, -1.0523740000612605e-08, 1.9999961070876222, 1.3696899999681651e-08, -0.02359242123294, -0.9994101984000001, 0.034292461608, 0.016925059512000004, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.440450847989736e-16, 
0.0, 0.0, 0.0, 0.0056967899025, -0.999987554251002, 1.2840361139644401e-06, 1.3696899999681651e-08, 2.000011073753859, -1.0033622240757458, -0.0342926095482, -0.9994159985598761, 0.973338334985236, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.203732932500001e-18, 
0.0, 0.0, 0.0, 0.01787620605, 1.0034817659968553, -0.99999788938492, -0.02359242123294, -1.0033622240757458, 3.0072961460676497, 5.7538999967844174e-09, -3.06455299668711e-08, -1.00000141432467, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.80101559214e-18, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.9994101984000001, -0.0342926095482, 5.7538999967844174e-09, 1.9999940925648998, -4.742859999747845e-08, -0.050293102734600005, -0.9999291077999999, 0.0116625663, -0.003036011868, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.880462940770579e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.034292461608, -0.9994159985598761, -3.06455299668711e-08, -4.742859999747845e-08, 2.0000049371401802, -0.9721817816841101, -0.011662650754300001, -0.99993352239745, 1.038046670071882, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.7471019827199995e-18, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.016925059512000004, 0.973338334985236, -1.00000141432467, -0.050293102734600005, -0.9721817816841101, 2.94766258462205, -6.929939999957639e-08, 4.461750999768985e-07, -0.999998249434076, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0376575840800005e-18, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.9999291077999999, -0.011662650754300001, -6.929939999957639e-08, 2.00000993433725, -1.6181899999890792e-08, -0.00907050352444, -0.9999613826, -0.009154428579599999, 0.061327620972, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.332235379185914e-15, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0116625663, -0.99993352239745, 4.461750999768985e-07, -1.6181899999890792e-08, 2.0000082349422277, -1.0380108409357391, 0.0091544009172, -0.9999612749812535, 0.9927090251141519, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.7344735349999947e-20, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.003036011868, 1.038046670071882, -0.999998249434076, -0.00907050352444, -1.0380108409357391, 3.0775418411057522, 1.579439999830137e-09, -2.0121460378662482e-07, -0.9999995433275339, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.9410669521999992e-18, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.9999613826, 0.0091544009172, 1.579439999830137e-09, 2.0000090800721813, -2.52721800003472e-08, -0.05223720194521, -0.018285328091999998, 0.99983653044, -1.1123597987999998, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -9.358100123903231e-17, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.009154428579599999, -0.9999612749812535, -2.0121460378662482e-07, -2.52721800003472e-08, 1.999996543092868, -0.99322390857258, -0.9998309316884039, -0.018285255495719997, 0.029119859904399997, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.3508105846796e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.061327620972, 0.9927090251141519, -0.9999995433275339, -0.05223720194521, -0.99322390857258, 2.9892354421190004, -1.3586152001886828e-07, -6.818740000454233e-09, -1.0000029413805, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0179947131419995e-17, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.018285328091999998, -0.9998309316884039, -1.3586152001886828e-07, 1.9999934268992399, 7.757001998833058e-09, -0.00877524057499999, -0.9999257524699999, 0.011855430745, 0.034339524813999996, 0.0, 0.0, 0.0, 8.776710652568e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.99983653044, -0.018285255495719997, -6.818740000454233e-09, 7.757001998833058e-09, 2.0000021299116657, -1.112701944200986, -0.01185543246425, -0.9999303692746249, 1.02729038677645, 0.0, 0.0, 0.0, 1.353579972152051e-15, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.1123597987999998, 0.029119859904399997, -1.0000029413805, -0.00877524057499999, -1.112701944200986, 3.238190279440499, -2.472608999951815e-07, -1.210417500080674e-07, -1.0000022064990999, 0.0, 0.0, 0.0, -6.330888939024999e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.9999257524699999, -0.01185543246425, -2.472608999951815e-07, 2.0000070902121303, -6.713599999983426e-08, -0.04651634431750001, -0.99941451264, 0.03431589408, -0.0023019403455999997, -7.863223145167301e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.011855430745, -0.9999303692746249, -1.210417500080674e-07, -6.713599999983426e-08, 2.0000056071716403, -1.026810891923148, -0.03431588604344, -0.9994133594373199, 0.9178210710936023, -6.898410060010201e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.034339524813999996, 1.02729038677645, -1.0000022064990999, -0.04651634431750001, -1.026810891923148, 3.0565026636606003, -2.1412800002286285e-08, 9.29963999871553e-08, -1.000000145625288, 8.073228126103001e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.99941451264, -0.03431588604344, -2.1412800002286285e-08, 1.0000013194892, 1.3358199999543015e-08, -0.029195100238040006, 1.0823329625610002e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.03431589408, -0.9994133594373199, 9.29963999871553e-08, 1.3358199999543015e-08, 1.0000015105368598, -0.917357916135795, -8.312660567825e-17, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0023019403455999997, 0.9178210710936023, -1.000000145625288, -0.029195100238040006, -0.917357916135795, 1.842399131325737, 7.7704978258215e-17, 
};
float step8_node0_M_correct_data[] = {
1.41421, 0.0, 0.0, -0.707083, 0.00578676, 0.0237462, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -7.85046e-17, 
0.0, 1.41421, 0.0, -0.00578676, -0.707083, 0.706469, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
0.0, 0.0, 1.41421, 0.0, 0.0, -0.707107, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0861e-19, 
0.0, 0.0, 0.0, 1.22475, -2.1545e-10, -0.0170474, -0.816483, 0.00465139, 0.0145958, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.35975e-16, 
0.0, 0.0, 0.0, 0.0, 1.22475, -0.407754, -0.00465139, -0.816483, 0.819336, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.29466e-18, 
0.0, 0.0, 0.0, 0.0, 0.0, 1.35392, -0.0116813, -0.245837, -0.491655, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.62227e-18, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.15464, -0.00248708, -0.0117849, -0.86556, 0.0296997, 0.0146583, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.88417e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.12823, -0.403599, -0.0323031, -0.885761, 0.862745, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.81467e-18, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.38961, -0.0167227, -0.257009, -0.468927, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.46205e-18, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.1178, -0.0064447, -0.0157254, -0.894551, 0.0104335, -0.00271606, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -5.7117e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.07166, -0.307048, -0.0162624, -0.933007, 0.968618, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -6.86801e-18, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.3743, -0.0138693, -0.208334, -0.511263, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -5.6013e-18, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.09514, -0.00797076, -0.00259233, -0.91309, -0.00835914, 0.0559998, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.49772e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.04208, -0.231066, 0.00180061, -0.959646, 0.953051, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.22046e-18, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.35075, -0.00144436, -0.164178, -0.577188, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.63998e-18, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.07994, -0.00568722, -0.0033836, -0.0169318, 0.925826, -1.03002, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.47274e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.02568, -0.168615, -0.974892, -0.0126939, 0.0226795, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.15839e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.31004, -0.125522, 0.000757409, -0.763079, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.2461e-17, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.01663, 0.00334025, -0.0982545, -0.983569, 0.0116615, 0.0337778, 0.0, 0.0, 0.0, 6.62678e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.06896, -0.147703, -0.0080172, -0.93546, 0.960913, 0.0, 0.0, 0.0, 7.87647e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.25019, -0.0782477, -0.109603, -0.683699, 0.0, 0.0, 0.0, 7.98429e-17, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.01312, -0.00454648, -0.0583224, -0.986472, 0.0338715, -0.00227213, -1.20391e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.05487, -0.192924, -0.0367826, -0.947282, 0.87007, 4.49791e-17, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.27434, -0.0507162, -0.14186, -0.653103, 6.61702e-17, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.15149, -0.0569331, -0.214905, -3.64295e-17, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.279545, -0.708164, -1.04198e-16, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.333389, -3.17135e-19, 
};


const int step8_node1_num_factors = 0;
const bool step8_node1_marked = true;
const bool step8_node1_fixed = false;
int step8_node1_factor_height[] = {};
int step8_node1_factor_width[] = {};
int* step8_node1_factor_ridx[] = {};
float* step8_node1_factor_data[] = {};
int step8_node1_factor_num_blks[] = {};
int* step8_node1_factor_A_blk_start[] = {};
int* step8_node1_factor_B_blk_start[] = {};
int* step8_node1_factor_blk_width[] = {};
const int step8_node1_parent = -1;
const int step8_node1_height = 1;
const int step8_node1_width = 1;
float step8_node1_data[] = {
0, 
};
const int step8_node1_num_blks = 0;
int step8_node1_A_blk_start[] = {};
int step8_node1_B_blk_start[] = {};
int step8_node1_blk_width[] = {};
const float step8_node1_H_correct_cond = 1.0;
float step8_node1_H_correct_data[] = {
5.963168641600001e-60, 
};
float step8_node1_M_correct_data[] = {
-2.44196e-30, 
};


const int step8_nnodes = 2;
bool step8_node_marked[] = {step8_node0_marked, step8_node1_marked, };
bool step8_node_fixed[] = {step8_node0_fixed, step8_node1_fixed, };
int step8_node_num_factors[] = {step8_node0_num_factors, step8_node1_num_factors, };
int* step8_node_factor_height[] = {step8_node0_factor_height, step8_node1_factor_height, };
int* step8_node_factor_width[] = {step8_node0_factor_width, step8_node1_factor_width, };
int** step8_node_factor_ridx[] = {step8_node0_factor_ridx, step8_node1_factor_ridx, };
float** step8_node_factor_data[] = {step8_node0_factor_data, step8_node1_factor_data, };
int* step8_node_factor_num_blks[] = {step8_node0_factor_num_blks, step8_node1_factor_num_blks, };
int** step8_node_factor_A_blk_start[] = {step8_node0_factor_A_blk_start, step8_node1_factor_A_blk_start, };
int** step8_node_factor_B_blk_start[] = {step8_node0_factor_B_blk_start, step8_node1_factor_B_blk_start, };
int** step8_node_factor_blk_width[] = {step8_node0_factor_blk_width, step8_node1_factor_blk_width, };
int step8_node_parent[] = {step8_node0_parent, step8_node1_parent, };
int step8_node_height[] = {step8_node0_height, step8_node1_height, };
int step8_node_width[] = {step8_node0_width, step8_node1_width, };
float* step8_node_data[] = {step8_node0_data, step8_node1_data, };
int step8_node_num_blks[] = {step8_node0_num_blks, step8_node1_num_blks, };
int* step8_node_A_blk_start[] = {step8_node0_A_blk_start, step8_node1_A_blk_start, };
int* step8_node_B_blk_start[] = {step8_node0_B_blk_start, step8_node1_B_blk_start, };
int* step8_node_blk_width[] = {step8_node0_blk_width, step8_node1_blk_width, };
float step8_node_H_correct_cond[] = {step8_node0_H_correct_cond, step8_node1_H_correct_cond, };
float* step8_node_H_correct_data[] = {step8_node0_H_correct_data, step8_node1_H_correct_data, };
float* step8_node_M_correct_data[] = {step8_node0_M_correct_data, step8_node1_M_correct_data, };


