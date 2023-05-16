import sys
from optparse import OptionParser
from helper_functions import *
import copy

def addNodeData(n1, n2):
    if n1 is None:
        return copy.deepcopy(n2)

    print("Not copying node")
    n1 = copy.deepcopy(n1)
    for d1, d2 in zip([n1.data, n1.flat_data], [n2.data, n2.flat_data]):
        for k in d1.keys():
            d1[k] += d2[k]
    return n1


if __name__ == "__main__":  
    parser = OptionParser()
    parser.add_option("--input_file1", dest="input_file1",
                      default="", help="The filename that contains the first csv.")
    parser.add_option("--input_file2", dest="input_file2",
                      default="", help="The filename that contains the second csv.")
    parser.add_option("--output_file", dest="output_file",
                      default="", help="The filename that the output is printed to.")
    (options, args) = parser.parse_args()

    l1 = readCSV(options.input_file1)
    l2 = readCSV(options.input_file2)

    l3 = mergeLists(l1, l2)

    flattenData(l3, ["Instructions Retired: Self", "Microarchitecture Usage: Total", "CPU Time: Self", "Microarchitecture Usage: Self", "Instructions Retired: Total"])

    total_cpu_time = l3[0].data["CPU Time: Total"]

    for node in l3:
        CI(node.data)
        CI(node.flat_data)
        node.data["CPU Time Percentage"] = node.data["CPU Time: Total"] / total_cpu_time

    delim = ";"


    with open(options.output_file, 'w') as fout:

        header_names = ["CPU Time: Total", "CPU Time Percentage", "flops", "mems", "CI"]

        operation_names = ["make_shared<gtsam::HessianFactor", "HessianFactor::eliminateCholesky", "GaussianBayesTree::optimize", "optimizeWildfire", "NonlinearFactorGraph::linearize"]
        operation_found = [False for _ in operation_names]
        operation_nodes = [None for _ in operation_names]

        stack = [l3[0]]
        while len(stack) != 0:
            cur_node = stack.pop(0)

            found_flag = False
            for i in range(len(operation_names)):
                op = operation_names[i]
                if op in cur_node.data["Source Function Stack"]:
                    # assert operation_found[i] == False, "Operation {} found twice".format(op)
                    if operation_found[i]:
                        print("Operation {} found twice".format(op))
                    operation_found[i] = True
                    operation_nodes[i] = addNodeData(operation_nodes[i], cur_node)
                    found_flag = True

            if not found_flag:
                stack = [*stack, *cur_node.children]

        fout.write("{}".format("Source Function Stack"))
        for header in header_names:
            fout.write("{} {}".format(delim, header))
        fout.write("\n")

        for i, node in enumerate(operation_nodes):
            # assert node is not None, f"{operation_names[i]} not found"
            if node is None:
                print(f"{operation_names[i]} not found")
                continue
            fout.write("{}".format(node.data["Source Function Stack"]))
            for header in header_names:
                fout.write("{} {}".format(delim, node.data[header]))

            fout.write("\n")


                    











