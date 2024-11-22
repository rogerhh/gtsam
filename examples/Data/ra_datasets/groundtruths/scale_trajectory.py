import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Scale a trajectory')
    parser.add_argument('trajectory', type=str, help='Trajectory file')
    parser.add_argument('scale', type=float, help='Scale factor')
    parser.add_argument('output', type=str, help='Output file')
    args = parser.parse_args()

    outlines = []

    # Read file in kiiti format
    with open(args.trajectory, 'r') as f:
        for line in f:
            arr = line.split()
            if len(arr) == 12:
                arr[3] = str(float(arr[3]) * args.scale)
                arr[7] = str(float(arr[7]) * args.scale)
                arr[11] = str(float(arr[11]) * args.scale)
                outlines.append(' '.join(arr) + '\n')
            else:
                outlines.append(line)

    # Write to file
    with open(args.output, 'w') as f:
        for line in outlines:
            f.write(line)
