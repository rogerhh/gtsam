import numpy as np
import matplotlib.pyplot as plt
import argparse

def plot_trajectory(infile, outfile):
    # Load the data from the text file
    data = []

    with open(infile, 'r') as file:
        for line in file:
            # Each line contains 12 values corresponding to the first 3 rows of a 4x4 transformation matrix
            values = list(map(float, line.split()))
            transformation_matrix = np.array(values).reshape(3, 4)
            
            # Extract the x, y coordinates (elements [0, 3] and [1, 3] in the matrix)
            x, y = transformation_matrix[0, 3], transformation_matrix[2, 3]
            print(x, y)
            data.append((x, y))

    # Separate the x and y coordinates
    x_coords = [point[0] for point in data]
    y_coords = [point[1] for point in data]

    # Plot the trajectory
    plt.figure(figsize=(10, 6))
    plt.plot(x_coords, y_coords, marker='o', linestyle='-', color='b')
    plt.xlabel('X coordinate')
    plt.ylabel('Y coordinate')
    plt.title('Trajectory Plot')
    plt.grid(True)

    plt.show()
    
    # Save the plot to the specified output file
    plt.savefig(outfile)
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot x,y coordinates of a trajectory from a KITTI pose file.")
    parser.add_argument('--infile', type=str, required=True, help="Path to the input text file containing KITTI poses.")
    parser.add_argument('--outfile', type=str, required=True, help="Path to save the output plot image.")
    
    args = parser.parse_args()

    plot_trajectory(args.infile, args.outfile)

