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

        if blk_width[-1] >= Clique.MAX_DENSE_BLOCK_SIZE:
            in_blk = False

    return A_blk_start, B_blk_start, blk_width

class Clique:

    no_values = False

    RELIN_COST_7_3 = 15000
    RELIN_COST_13_6 = 20000
    SYM_COST = 2000
    SYM_COST_REORDER = 4000

    def __init__(self, fin, index=0):
        # ======================================== #
        # Variable initialization
        # ======================================== #

        self.marked = False
        self.fixed = False
        self.parent = None
        self.inverse_block_indices = None
        self.active_factor_indices = []
        self.index = index
        self.num_relin_factors = 0

        self.workspace_needed = 0
        self.children = []

        # ======================================== #
        # END Variable initialization
        # ======================================== #

        line = fin.readline()
        arr = line.split()

        self.height = int(arr[0])
        self.width = int(arr[1])
        self.clique_size = int(arr[2])
        self.block_indices_size = int(arr[3])

        self.block_indices = []

        for i in range(self.block_indices_size):
            line = fin.readline()
            arr = line.split()

            self.block_indices.append([int(x) for x in arr])

        self.keys = []
        for i in range(self.clique_size):
            self.keys.append(self.block_indices[i][0])

        self.matrix = np.zeros((self.height, self.width))

        if not Clique.no_values:
            for i in range(self.height):
                line = fin.readline()
                arr = line.split()
                assert(len(arr) == self.width)
                for j in range(self.width):
                    self.matrix[i, j] = float(arr[j])

        self.ancestor_block_indices = {}
        for key, row, height in self.block_indices:
            self.ancestor_block_indices[key] = [key, row, height]

        self.ancestor_keys = []
        for i in range(self.clique_size, self.block_indices_size):
            key, row, height = self.block_indices[i]
            self.ancestor_keys.append(key)

        self.C_block_indices = {}
        cur_row = 0
        for key in self.ancestor_keys:
            height = self.ancestor_block_indices[key][2]
            self.C_block_indices[key] = [key, cur_row, height]
            cur_row += height

        self.workspace_size = self.height * self.height * 4 + 100 # Padding


    def __str__(self):
        if self.marked:
            s = "(marked) "
        elif self.fixed:
            s = "(fixed) "
        else:
            s = "(none) "
        for k in self.keys:
            s += str(k) + " "
        return s

    def mark_clique_and_ancestors(self):
        if self.marked:
            return
        
        self.marked = True

        if self.parent is not None:
            self.parent.mark_clique_and_ancestors()

    def build_inverse_block_indices(self):
        self.inverse_block_indices = {}
        for key, row, height in self.block_indices:
            self.inverse_block_indices[key] = [key, row, height]


    def print_clique(self, fout, step, num_threads, is3D, is_reorder):
        self.active_factor_indices = sorted(self.active_factor_indices)

        # print(f"step = {step}, clique = {self.index}, block_indices = {self.block_indices}")

        # Print all the factors related to this clique
        fout.write(f"const int step{step}_node{self.index}_num_factors = {len(self.active_factor_indices)};\n")
        relin_cost = self.num_relin_factors * (Clique.RELIN_COST_13_6 if is3D else Clique.RELIN_COST_7_3)
        fout.write(f"const int step{step}_node{self.index}_relin_cost = {relin_cost};\n")

        sym_cost = self.clique_size * (Clique.SYM_COST_REORDER if is_reorder else Clique.SYM_COST)
        fout.write(f"const int step{step}_node{self.index}_sym_cost = {sym_cost};\n")

        marked_str = "true" if self.marked else "false"
        fixed_str = "true" if self.fixed else "false"
        fout.write(f"const bool step{step}_node{self.index}_marked = {marked_str};\n")
        fout.write(f"const bool step{step}_node{self.index}_fixed = {fixed_str};\n")

        fout.write(f"int step{step}_node{self.index}_factor_height[] = {{")
        for factor_index in self.active_factor_indices:
            fout.write(f"step{step}_factor{factor_index}_height, ")
        fout.write("};\n")

        fout.write(f"int step{step}_node{self.index}_factor_width[] = {{")
        for factor_index in self.active_factor_indices:
            fout.write(f"step{step}_factor{factor_index}_width, ")
        fout.write("};\n")

        fout.write(f"int* step{step}_node{self.index}_factor_ridx[] = {{")
        for factor_index in self.active_factor_indices:
            fout.write(f"step{step}_factor{factor_index}_ridx, ")
        fout.write("};\n")

        fout.write(f"float* step{step}_node{self.index}_factor_data[] = {{")
        for factor_index in self.active_factor_indices:
            fout.write(f"step{step}_factor{factor_index}_data, ")
        fout.write("};\n")

        fout.write(f"int step{step}_node{self.index}_factor_num_blks[] = {{")
        for factor_index in self.active_factor_indices:
            fout.write(f"step{step}_factor{factor_index}_num_blks, ")
        fout.write("};\n")

        fout.write(f"int* step{step}_node{self.index}_factor_A_blk_start[] = {{")
        for factor_index in self.active_factor_indices:
            fout.write(f"step{step}_factor{factor_index}_A_blk_start, ")
        fout.write("};\n")

        fout.write(f"int* step{step}_node{self.index}_factor_B_blk_start[] = {{")
        for factor_index in self.active_factor_indices:
            fout.write(f"step{step}_factor{factor_index}_B_blk_start, ")
        fout.write("};\n")

        fout.write(f"int* step{step}_node{self.index}_factor_blk_width[] = {{")
        for factor_index in self.active_factor_indices:
            fout.write(f"step{step}_factor{factor_index}_blk_width, ")
        fout.write("};\n")

        if self.parent is not None:
            fout.write(f"const int step{step}_node{self.index}_parent = {self.parent.index};\n")
        else:
            fout.write(f"const int step{step}_node{self.index}_parent = -1;\n")

        fout.write(f"const int step{step}_node{self.index}_height = {self.height};\n")
        fout.write(f"const int step{step}_node{self.index}_width = {self.width};\n")


        # fout.write(f"float step{step}_node{self.index}_data[{self.height * self.width}] = {{")
        # if not Clique.no_values:
        #     fout.write("\n")
        #     if self.marked:
        #         for j in range(self.width):
        #             for i in range(self.height):
        #                 fout.write("0, ")
        #             fout.write("\n")
        #     elif self.fixed:
        #         for j in range(self.width):
        #             for i in range(self.height):
        #                 fout.write(f"{self.matrix[i, j]}, ")
        #             fout.write("\n")
        # else:
        #     fout.write("0")

        # fout.write("};\n")

        fout.write(f"float step{step}_node{self.index}_data[{self.height * self.height}] = {{")
        if not Clique.no_values:
            fout.write("\n")
            if self.marked:
                for j in range(self.height):
                    for i in range(self.height):
                        fout.write("0, ")
                    fout.write("\n")
            elif self.fixed:
                for j in range(self.width):
                    for i in range(self.height):
                        fout.write(f"{self.matrix[i, j]}, ")
                    fout.write("\n")
                for j in range(self.width, self.height):
                    for i in range(self.height):
                        fout.write(f"0, ")
                    fout.write("\n")
        else:
            fout.write("0")

        fout.write("};\n")

        if self.parent is not None:
            A_blk_start, B_blk_start, blk_width = group_block_indices(self.ancestor_keys, self.C_block_indices, self.parent.ancestor_block_indices)
        else:
            A_blk_start = []
            B_blk_start = []
            blk_width = []

        fout.write(f"const int step{step}_node{self.index}_num_blks = {len(blk_width)};\n")
        fout.write(f"int step{step}_node{self.index}_A_blk_start[] = {{")
        for a in A_blk_start:
            fout.write(f"{a}, ")
        fout.write("};\n")

        fout.write(f"int step{step}_node{self.index}_B_blk_start[] = {{")
        for a in B_blk_start:
            fout.write(f"{a}, ")
        fout.write("};\n")

        fout.write(f"int step{step}_node{self.index}_blk_width[] = {{")
        for a in blk_width:
            fout.write(f"{a}, ")
        fout.write("};\n")

        if not Clique.no_values:
            LA = self.matrix[:self.width, :self.width]
            LAB = self.matrix[:, :self.width]

            H_correct = LAB @ LA.T

            cond = np.linalg.cond(H_correct)

            fout.write(f"const float step{step}_node{self.index}_H_correct_cond = {cond};\n")

            fout.write(f"float step{step}_node{self.index}_H_correct_data[] = {{\n")
            for j in range(self.width):
                for i in range(self.height):
                    fout.write(f"{H_correct[i, j]}, ")
                fout.write("\n")
            fout.write("};\n")

            fout.write(f"float step{step}_node{self.index}_M_correct_data[] = {{\n")
            for j in range(self.width):
                for i in range(self.height):
                    fout.write(f"{self.matrix[i, j]}, ")
                fout.write("\n")
            fout.write("};\n")

        fout.write("\n\n")

    def print_clique_ridx(self, fout, step, block_indices):
        fout.write(f"int step{step}_node{self.index}_ridx[] = {{\n")
        for key, row, height in self.block_indices:
            _, total_row, _ = block_indices[key]
            for i in range(height):
                fout.write(f"{total_row + i}, ")
        fout.write("};\n")

    @staticmethod
    def print_clique_variable(fout, t, prefix, postfix, pred, default, max_clique):
        fout.write(f"{t} {prefix}_node_{postfix}[] = {{")
        for i in range(max_clique):
            if pred(i):
                fout.write(f"{prefix}_node{i}_{postfix}, ")
            else:
                fout.write(f"{default}, ")
        fout.write("};\n")


    @staticmethod
    def print_clique_metadata(fout, step, active_clique_indices, max_clique):
        fout.write(f"const int step{step}_nnodes = {max_clique};\n");

        pred = lambda clique : clique in active_clique_indices

        Clique.print_clique_variable(fout, t="bool", prefix=f"step{step}", pred=pred, max_clique=max_clique, default="false", postfix="marked")

        Clique.print_clique_variable(fout, t="bool", prefix=f"step{step}", pred=pred, max_clique=max_clique, default="false", postfix="fixed")

        Clique.print_clique_variable(fout, t="int", prefix=f"step{step}", pred=pred, max_clique=max_clique, default="0", postfix="num_factors")

        Clique.print_clique_variable(fout, t="int", prefix=f"step{step}", pred=pred, max_clique=max_clique, default="0", postfix="relin_cost")

        Clique.print_clique_variable(fout, t="int", prefix=f"step{step}", pred=pred, max_clique=max_clique, default="0", postfix="sym_cost")

        Clique.print_clique_variable(fout, t="int*", prefix=f"step{step}", pred=pred, max_clique=max_clique, default="0", postfix="factor_height")

        Clique.print_clique_variable(fout, t="int*", prefix=f"step{step}", pred=pred, max_clique=max_clique, default="0", postfix="factor_width")

        Clique.print_clique_variable(fout, t="int**", prefix=f"step{step}", pred=pred, max_clique=max_clique, default="0", postfix="factor_ridx")

        Clique.print_clique_variable(fout, t="float**", prefix=f"step{step}", pred=pred, max_clique=max_clique, default="0", postfix="factor_data")

        Clique.print_clique_variable(fout, t="int*", prefix=f"step{step}", pred=pred, max_clique=max_clique, default="0", postfix="factor_num_blks")

        Clique.print_clique_variable(fout, t="int**", prefix=f"step{step}", pred=pred, max_clique=max_clique, default="0", postfix="factor_A_blk_start")

        Clique.print_clique_variable(fout, t="int**", prefix=f"step{step}", pred=pred, max_clique=max_clique, default="0", postfix="factor_B_blk_start")

        Clique.print_clique_variable(fout, t="int**", prefix=f"step{step}", pred=pred, max_clique=max_clique, default="0", postfix="factor_blk_width")

        Clique.print_clique_variable(fout, t="int", prefix=f"step{step}", pred=pred, max_clique=max_clique, default="0", postfix="parent")

        Clique.print_clique_variable(fout, t="int", prefix=f"step{step}", pred=pred, max_clique=max_clique, default="0", postfix="height")

        Clique.print_clique_variable(fout, t="int", prefix=f"step{step}", pred=pred, max_clique=max_clique, default="0", postfix="width")

        Clique.print_clique_variable(fout, t="float*", prefix=f"step{step}", pred=pred, max_clique=max_clique, default="0", postfix="data")

        Clique.print_clique_variable(fout, t="int", prefix=f"step{step}", pred=pred, max_clique=max_clique, default="0", postfix="num_blks")

        Clique.print_clique_variable(fout, t="int*", prefix=f"step{step}", pred=pred, max_clique=max_clique, default="0", postfix="A_blk_start")

        Clique.print_clique_variable(fout, t="int*", prefix=f"step{step}", pred=pred, max_clique=max_clique, default="0", postfix="B_blk_start")

        Clique.print_clique_variable(fout, t="int*", prefix=f"step{step}", pred=pred, max_clique=max_clique, default="0", postfix="blk_width")

        if not Clique.no_values:
            Clique.print_clique_variable(fout, t="float", prefix=f"step{step}", pred=pred, max_clique=max_clique, default="0", postfix="H_correct_cond")

            Clique.print_clique_variable(fout, t="float*", prefix=f"step{step}", pred=pred, max_clique=max_clique, default="0", postfix="H_correct_data")

            Clique.print_clique_variable(fout, t="float*", prefix=f"step{step}", pred=pred, max_clique=max_clique, default="0", postfix="M_correct_data")

        Clique.print_clique_variable(fout, t="int*", prefix=f"step{step}", pred=lambda i : True, max_clique=max_clique, default="0", postfix="ridx")


        fout.write("\n\n")





