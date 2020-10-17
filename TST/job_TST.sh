#$-q idefix.q,obelix.q
#$ -l h_vmem=3G
#$ -l h_cpu=999:00:00
#$ -pe smp 85
#$-M j.teuffel@icloud.com -m ase
#$ -cwd
path="run_Time_Eval_Release_G85_Res500_Th2pi"
export BOOST_ROOT="~/boost_1_72_0/"
module load intel/2020
mkdir $path
cd $path
cmake cmake .. -DCMAKE_BUILD_TYPE=Release  -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++  -DTHREADS=85 -DMKL_PARALEL=TRUE -DGRID=35 -DT_RES=500 -DT_ROT=1.0 -DT_C=-2.0 -DW_START=-3 -DW_END=3 -DW_RES=4 -DUSE_MAGMA=FALSE -DVERBOSE=TRUE  -DUSE_GPU=FALSE -DTIME_EVOLUTION=TRUE -DCMAKE_PREFIX_PATH="~/boost_1_72_0/"
make -j 8
echo "Running Program"
./TST
