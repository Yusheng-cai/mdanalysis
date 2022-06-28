import sys
import os
import argparse
import subprocess
import numpy as np

def ReadTabulatedData(file:str):
    """
    Function that reads tabulated data 

    Assumes data is of the following format 
    # data1     data2   data3 ...
    """
    with open(file, "r") as f:
        lines = f.readlines()
        lines = [l.rstrip("\n").lstrip().split() for l in lines]
        lines = [[float(f) for f in l] for l in lines]

    return np.array(lines)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="For testing mdanalysis")
    parser.add_argument("-p", "--program", type=str, help="Program name")
    parser.add_argument("-nt", "--num_threads", type=str, help="Number of threads to run the program, will override system ENV")
    parser.add_argument("-i", "--input_file", type=str, help="Input file name for the program")
    parser.add_argument("-a", "--absolute_path", type=str, help="The absolute path in which to read the files in the input file")

    # parse the arguments passed in
    args = parser.parse_args()

    # set the omp of threads
    os.environ["OMP_NUM_THREADS"] = str(args.num_threads)

    # run the process
    p = subprocess.run([args.program, args.input_file, "-apath", args.absolute_path])

    if p.returncode != 0:
        raise Exception("Invalid result ", p.returncode)