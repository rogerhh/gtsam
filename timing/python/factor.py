import numpy as np

def group_block_indices(reduced_sorted_keys, A_block_indices, B_block_indices):
    in_blk = False
    A_blk_start = []
    B_blk_start = []
    blk_width = []
    for key in reduced_sorted_keys:
        if key == 0:
            continue

        _, A_row, A_width = A_block_indices[key]
        _, B_row, B_width = B_block_indices[key]

        if not in_blk or B_blk_start[-1] + blk_width[-1] != B_row:
            in_blk = True
            A_blk_start.append(A_row)
            B_blk_start.append(B_row)
            blk_width.append(A_width)
        else:
            blk_width[-1] += A_width

        if blk_width[-1] >= Factor.MAX_DENSE_BLOCK_SIZE:
            in_blk = False

    return A_blk_start, B_blk_start, blk_width

class Factor:
    no_values = False

    def __init__(self, fin, key_to_ordering, key_width):

        ########################
        # Init
        ########################

        self.relin = False

        ########################
        # End Init
        ########################


        line = fin.readline()
        arr = line.split()

        self.key_width = key_width
        self.key_to_ordering = key_to_ordering

        self.index = int(arr[0])
        self.height = int(arr[1])
        self.width = int(arr[2])
        self.num_keys = int(arr[3])

        self.keys = []

        line = fin.readline()
        arr = line.split()

        assert(len(arr) == self.num_keys)

        for k in range(self.num_keys):
            self.keys.append(int(arr[k]))

        self.matrix = np.zeros((self.height, self.width))

        if not Factor.no_values:
            for i in range(self.height):
                line = fin.readline()
                arr = line.split()
                print(self.width)
                print(arr)
                assert(len(arr) == self.width)
                for j in range(self.width):
                    self.matrix[i, j] = float(arr[j])

        # Sort keys and factor matrix
        self.sorted_keys = sorted(self.keys, key=lambda key : key_to_ordering[key])
        self.sorted_matrix = np.zeros((self.height, self.width))

        self.block_indices = {}
        src_col = 0
        for key in self.keys:
            width = self.key_width[key]
            self.block_indices[key] = [key, src_col, width]
            src_col += width

        self.sorted_block_indices = {}
        cur_col = 0
        for sorted_key in self.sorted_keys:
            _, src_col, width = self.block_indices[sorted_key]
            self.sorted_matrix[:, cur_col:cur_col+width] = self.matrix[:, src_col:src_col+width]
            self.sorted_block_indices[sorted_key] = [sorted_key, cur_col, width]
            cur_col += width


    def __str__(self):
        return f"Factor {self.index} {self.keys}"

    def print_factor(self, fout, inverse_block_indices, prefix, lowest_key=-1):
        # Print the transposed factor in a header file

        # Factor may be truncated because not all keys are needed
        reduced_col_start = self.width
        reduced_index_start = len(self.sorted_keys)
        for key in reversed(self.sorted_keys):
            reduced_col_start -= self.key_width[key]
            reduced_index_start -= 1
            if key == lowest_key:
                break

        reduced_width = self.width - reduced_col_start
        self.reduced_keys = [self.sorted_keys[i] for i in range(reduced_index_start, len(self.sorted_keys))]
        self.reduced_matrix = self.sorted_matrix[:, reduced_col_start:self.width]

        reduced_block_indices = {}
        reduced_cur_col = 0
        for key in self.reduced_keys:
            _, col, width = self.sorted_block_indices[key]
            reduced_block_indices[key] = [key, reduced_cur_col, width]
            reduced_cur_col += width

        # Write transpose of factor
        fout.write(f"const int {prefix}factor{self.index}_height = {reduced_width};\n")
        fout.write(f"const int {prefix}factor{self.index}_width = {self.height};\n")
        fout.write(f"int {prefix}factor{self.index}_ridx[] = {{")

        for key in self.reduced_keys:
            _, row, height = inverse_block_indices[key]
            for r in range(row, row+height):
                fout.write(f"{r}, ")

        fout.write("};\n")

        reduced = self.reduced_matrix.astype(np.float32)

        fout.write(f"float {prefix}factor{self.index}_data[{reduced_width * self.height}] = {{")
        if not Factor.no_values:
            fout.write("\n")
            for j in range(self.height):
                for i in range(reduced_width):
                    fout.write(f"{self.reduced_matrix[j, i]:.7f}, ")

                fout.write("\n")
        else:
            fout.write("0")


        fout.write("};\n")
        fout.write("\n")

        # Group block indices
        A_blk_start, B_blk_start, blk_width = group_block_indices(self.reduced_keys, reduced_block_indices, inverse_block_indices)

        fout.write(f"const int {prefix}factor{self.index}_num_blks = {len(blk_width)};\n")
        fout.write(f"int {prefix}factor{self.index}_A_blk_start[] = {{")
        for a in A_blk_start:
            fout.write(f"{a}, ")
        fout.write("};\n")
        fout.write(f"int {prefix}factor{self.index}_B_blk_start[] = {{")
        for a in B_blk_start:
            fout.write(f"{a}, ")
        fout.write("};\n")
        fout.write(f"int {prefix}factor{self.index}_blk_width[] = {{")
        for a in blk_width:
            fout.write(f"{a}, ")
        fout.write("};\n")
        fout.write("\n")

        self.num_blks = len(blk_width)

    @staticmethod
    def print_metadata(fout, factor_indices, prefix):
        fout.write(f"const int {prefix}num_factors = {len(factor_indices)};\n")
        fout.write(f"const int {prefix}factor_height[] = {{")
        for index in factor_indices:
            fout.write(f"{prefix}factor{index}_height, ")
        fout.write("};\n")
        fout.write(f"const int {prefix}factor_width[] = {{")
        for index in factor_indices:
            fout.write(f"{prefix}factor{index}_width, ")
        fout.write("};\n")
        fout.write(f"int* {prefix}factor_ridx[] = {{")
        for index in factor_indices:
            fout.write(f"{prefix}factor{index}_ridx, ")
        fout.write("};\n")
        fout.write(f"float* {prefix}factor_data[] = {{")
        for index in factor_indices:
            fout.write(f"{prefix}factor{index}_data, ")
        fout.write("};\n")
        fout.write(f"int {prefix}factor_num_blks[] = {{")
        for index in factor_indices:
            fout.write(f"{prefix}factor{index}_num_blks, ")
        fout.write("};\n")
        fout.write(f"int* {prefix}factor_A_blk_start[] = {{")
        for index in factor_indices:
            fout.write(f"{prefix}factor{index}_A_blk_start, ")
        fout.write("};\n")
        fout.write(f"int* {prefix}factor_B_blk_start[] = {{")
        for index in factor_indices:
            fout.write(f"{prefix}factor{index}_B_blk_start, ")
        fout.write("};\n")
        fout.write(f"int* {prefix}factor_blk_width[] = {{")
        for index in factor_indices:
            fout.write(f"{prefix}factor{index}_blk_width, ")
        fout.write("};\n")
        fout.write("\n")






