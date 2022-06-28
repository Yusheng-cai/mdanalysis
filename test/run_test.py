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
        lines = [l for l in lines if not l[0].startswith("#") and not l[0].startswith(";")]
        lines = [[float(f) for f in l] for l in lines]

    return np.array(lines)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="For testing mdanalysis")
    parser.add_argument("-p", "--program", type=str, help="Program name")
    parser.add_argument("-nt", "--num_threads", type=str, help="Number of threads to run the program, will override system ENV")
    parser.add_argument("-i", "--input_file", type=str, help="Input file name for the program")
    parser.add_argument("-a", "--absolute_path", type=str, help="The absolute path in which to read the files in the input file")
    parser.add_argument("-r", "--ref_folder", type=str, help="The reference folder which contains all the results to be compared with")

    # parse the arguments passed in
    args = parser.parse_args()

    # set the omp of threads
    os.environ["OMP_NUM_THREADS"] = str(args.num_threads)

    # run the process
    program_to_be_run = [args.program, args.input_file, "-apath", args.absolute_path]
    p = subprocess.run(program_to_be_run)

    # see if there's error, raise if there's error
    if p.returncode != 0:
        raise Exception("Invalid result ", p.returncode)

    # compare the results 
    files_to_be_checked = os.path.join(args.ref_folder, "file_to_be_checked.dat")
    with open(files_to_be_checked) as f:
        filesref = f.readlines()
        filesref = [f.rstrip("\n").lstrip() for f in filesref]
        files = [f.split("_")[0] + ".out" for f in filesref]
    
    # compare files to filesref
    for i in range(len(files)):
        ref = os.path.join(args.ref_folder, filesref[i])
        f   = files[i]
        
        data_ref = ReadTabulatedData(ref)
        data = ReadTabulatedData(f)

        isclose = np.allclose(data, data_ref, equal_nan=True)

        if not isclose:
            raise Exception("Invalid result")