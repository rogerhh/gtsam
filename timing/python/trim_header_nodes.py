import re
import argparse

def process_line(i, lines, new_lines, add):
    while True:
        if add:
            new_lines.append(lines[i])
        if ";" in lines[i]:
            i += 1
            return i
        i += 1

def keep_indices(file_path, output_path, indices):

    MAX_INDEX=300
    discard_indices = list(range(0, MAX_INDEX))
    for index in indices:
        print(index)
        discard_indices.remove(int(index))

    # Build patterns based on the list of indices
    metadata_pattern = re.compile(rf'node\d')
    patterns = [re.compile(rf'\b\w*node{index}\D\w*\b') for index in indices]
    discard_patterns = [re.compile(rf'\b\w*node{index}\D\w*\b') for index in discard_indices]
    print(discard_patterns)
    print(patterns)

    new_lines = []

    # Read the contents of the file
    with open(file_path, 'r') as file:
        lines = file.readlines()

        # First find nnodes
        for line in lines:
            if "nnodes" in line:
                nnodes = int(line.split()[-1].split(";")[0])

        i = 0;
        while i < len(lines):
            print(i)
            line = lines[i]
            if not line.strip():
                i += 1
                continue

            factor_data_pattern42 = re.compile(rf'step\d+_factor\d+_data\[42\]')
            if factor_data_pattern42.findall(line):
                var = line.split("[")[0].split()[-1]
                new_lines.append(f"#define {var} fake_factor_data42\n");
                i = process_line(i, lines, new_lines, add=False);
                continue
            factor_data_pattern78 = re.compile(rf'step\d+_factor\d+_data\[78\]')
            if factor_data_pattern78.findall(line):
                var = line.split("[")[0].split()[-1]
                new_lines.append(f"#define {var} fake_factor_data78\n");
                i = process_line(i, lines, new_lines, add=False);
                continue

            node_factor_data_pattern = re.compile(rf'step\d+_node\d+_factor_data\[\]')
            step_node_factor_data_pattern = re.compile(rf'step\d+_node_factor_data\[\]')
            if node_factor_data_pattern.findall(line):
                line = line.replace("float", "const float")
                lines[i] = lines[i].replace("float", "const float")
            elif step_node_factor_data_pattern.findall(line):
                line = line.replace("float", "const float")
                lines[i] = lines[i].replace("float", "const float")

            parent_pattern = re.compile(rf'node\d+_parent = ')
            if parent_pattern.findall(line):
                parent = int(line.split()[-1].split(";")[0])
                if parent not in indices:
                    line = line.split("=")[0] + f"= {nnodes - 1};\n"
                    lines[i] = line


            if not metadata_pattern.findall(line):
                while True:
                    new_lines.append(lines[i])
                    if ";" in lines[i]:
                        break
                    i += 1
                i += 1
                continue

            node_found = False
            for pattern in patterns:
                if pattern.findall(line):
                    node_found = True
                    break
            print("node_found = ", node_found)

            if not node_found or "M_correct" in line or "H_correct" in line:
                while True:
                    if ";" in lines[i]:
                        break
                    i += 1
                i += 1
                continue

            if node_found:
                discard_matched = False
                for pattern in discard_patterns:
                    discard_matches = pattern.findall(line)
                    if discard_matches:
                        discard_matched = True
                        break
                # This is a metadata line
                if discard_matched:
                    print("metadata")
                    while True:
                        line = lines[i]
                        for pattern in discard_patterns:
                            line = pattern.sub("0", line)
                        new_lines.append(line)
                        if ";" in line:
                            break
                        i += 1

                # This is a normal line
                else:
                    print("normal")
                    while True:
                        new_lines.append(lines[i])
                        if ";" in lines[i]:
                            break
                        i += 1
            i += 1



    # Write the updated content to the output file
    with open(output_path, 'w') as file:
        for i, line in enumerate(new_lines):
            file.write(line)
            if i == 1:
                file.write("const float fake_factor_data42[42] = {\n")
                file.write("0.0100000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, \n")
                file.write("0.0000000, 0.0100000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, \n")
                file.write("0.0000000, 0.0000000, 0.0100000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, \n")
                file.write("0.0000000, 0.0000000, 0.0000000, 0.0100000, 0.0000000, 0.0000000, 0.0000000, \n")
                file.write("0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0100000, 0.0000000, 0.0000000, \n")
                file.write("0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0100000, 0.0000000, \n")
                file.write("};\n\n")
                file.write("const float fake_factor_data78[78] = {\n")
                file.write("-6.7509500, 0.0086870, -2.1034900, 0.0000000, 0.0000000, 0.0000000, 7.0710700, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0052129, \n")
                file.write("0.0074352, -7.0708600, -0.0530642, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 7.0710700, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0496537, \n")
                file.write("2.1034900, 0.0528737, -6.7507400, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 7.0710700, 0.0000000, 0.0000000, 0.0000000, -0.0305234, \n")
                file.write("-0.0203507, -0.7914150, 0.0620453, -6.7509500, 0.0086870, -2.1034900, 0.0000000, 0.0000000, 0.0000000, 7.0710700, 0.0000000, 0.0000000, 0.0131322, \n")
                file.write("0.4387230, -0.0089242, 1.2506300, 0.0074352, -7.0708600, -0.0530642, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 7.0710700, 0.0000000, -0.0150843, \n")
                file.write("-0.0668644, -1.0634100, -0.0291635, 2.1034900, 0.0528737, -6.7507400, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 7.0710700, -0.2542740, \n")
                file.write("};\n\n")


    print(f"All occurrences of whole words containing specified indices have been replaced with 'false' and written to {output_path}")

def main():
    # Create the parser
    parser = argparse.ArgumentParser(description='Replace whole words containing specific patterns with "false" in a file.')
        
    # Add arguments
    parser.add_argument('--infile', required=True, help='Path to the input file')
    parser.add_argument('--outfile', required=True, help='Path to the output file')
    parser.add_argument('--indices', required=True, nargs='+', help='List of indices to replace (e.g., 0 1 2)')
        
    # Parse the arguments
    args = parser.parse_args()
        
    # Convert indices to a list of strings
    indices = args.indices

        
    # Call the function with the provided file paths and indices
    keep_indices(args.infile, args.outfile, indices)

if __name__ == '__main__':
        main()

