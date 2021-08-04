#/bin/bash

num_threads=0
test_dir=""
# We can read in multiple instance of the same flag by doing it multiple times, -o a -o b -o c etc.
while getopts ":p:i:d:n:r:o:s:" opt; do
    case $opt in 
        p) program=${OPTARG};;
        i) input=${OPTARG};;
	    n) num_threads=${OPTARG};;
        d) test_dir=${OPTARG};;
        r) output_ref+=(${OPTARG});;
        o) output_test+=(${OPTARG});;
        s) stdout_comp=${OPTARG};;
    esac
done 

if [[ ${num_threads} -eq 0 ]]
then
    num_threads=1
fi

# check if the test_dir string is null
if [[ -z ${test_dir} ]]
then 
    test_dir="."
fi

# list number of failed tests
failed_test=0

# set number of threads
export OMP_NUM_THREADS=${num_threads}

# Find the length of the output file 
len=${#output_test[@]}

if [[ ${len} -lt 1 ]]
then
    exit 1
fi

# run the program with the input
${program} ${input} > stdout

for ((i=0; i<${len}; ++i)) do
    test_file=${output_test[$i]}
    ref_file=${test_dir}/${output_ref[$i]}
    # if file not found, then exit right away
    if [[ ! -f ${test_file} || ! -f ${ref_file} ]]
    then
        echo "File not found."
        exit 1
    fi

    # if diff fails, then we can directly exit 
    diff ${test_file} ${ref_file}
    if [[ $? -eq 1 ]]
    then 
        ((failed_test++))
        exit 1
    fi
    rm ${test_file}
done

if [[ -f ${test_dir}/stdout_ref ]]
then
    diff stdout ${test_dir}/stdout_ref
fi

rm stdout

if [[ ${failed_test} -gt 0 ]]
then
    failed_test=1
else
    failed_test=0
fi

exit ${failed_test}