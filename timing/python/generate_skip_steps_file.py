import argparse

WITH_VALUE = 0
NO_NUMERIC = 1
NO_VALUE = 2
WITH_VALUE_GEMMINI = 3
NO_SETUP = 4

def generate_skips_steps_file(output_fname, max_num_steps, period):
    with open(output_fname, "w") as fout:
        fout.write(f"{max_num_steps + 1}\n")
        for i in range(0, max_num_steps + 1, 1):
            if i % period == 0:
                fout.write(f"{NO_NUMERIC}\n")
            else:
                fout.write(f"{NO_SETUP}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate skip steps file")
    parser.add_argument("output_fname", type=str, help="Output file name")
    parser.add_argument("max_num_steps", type=int, help="Maximum number of steps")
    parser.add_argument("period", type=int, help="Period")
    args = parser.parse_args()

    generate_skips_steps_file(args.output_fname, args.max_num_steps, args.period)
