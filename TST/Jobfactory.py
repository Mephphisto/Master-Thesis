import os

header = """#$-q idefix.q,obelix.q
#$ -l h_vmem=3G
#$ -l h_cpu=999:00:00
#$ -pe smp {0}
#$-M j.teuffel@icloud.com -m ase
#$ -cwd

path="run_Time_Eval_Benchmark_G{1}_Res{1}_Threads{0}"
export BOOST_ROOT="~/boost_1_72_0/"
echo $path
module load intel/2020
mkdir $path
cd $path
cmake cmake .. -DCMAKE_BUILD_TYPE=Release  -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++  -DTHREADS={0} -DMKL_PARALEL=TRUE -DTARGET_ARCHITECTURE=skylake -DGRID={1} -DT_RES={2} -DT_ROT=1.0 -DT_C=-2.0 -DW_START=-3 -DW_END=-3 -DW_RES=1 -DUSE_MAGMA=FALSE -DVERBOSE=FALSE  -DUSE_GPU=FALSE -DTIME_EVOLUTION=TRUE -DCMAKE_PREFIX_PATH="~/boost_1_72_0/"
make -j 8
echo "running {0} threads" 
echo "Running Program"
./TST
echo "finished"
"""

# {0} Threads
# {1} Grid
# {2} T_RES


try:
    os.mkdir("Jobs")
except:
    pass
for i in range(7):
    file = open("Jobs/job" + str(i) + ".sh", "w+")
    file.write(header.format(str(2 ** i), str(35), str(500)))
    file.close()
    os.system("qsub Jobs/job" + str(i) + ".sh")
