#$-q idefix.q,obelix.q
#$ -l h_vmem=3G
#$ -l h_cpu=999:00:00
#$ -pe smp 85
#$-M j.teuffel@icloud.com -m ase
#$ -cwd
<<<<<<< HEAD
path="rune_G35_Res500_NegMu"
=======
path="run_Time_Eval_Release_G85_Res500_Th2pi"
>>>>>>> 681774e4aadd278002bf4b37513dab9bf4d6039b
export BOOST_ROOT="~/boost_1_72_0/"
echo $path
module load intel/2020
mkdir $path
cd $path
<<<<<<< HEAD
cmake cmake .. -DCMAKE_BUILD_TYPE=Release  -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++  -DTHREADS=85 -DGRID=85 -DT_RES=85 -DUSE_MAGMA=FALSE -DVERBOSE=FALSE  -DUSE_GPU=FALSE -DTIME_EVOLUTION=FALSE -DCMAKE_PREFIX_PATH="~/boost_1_72_0/"
=======
cmake cmake .. -DCMAKE_BUILD_TYPE=Release  -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++  -DTHREADS=85 -DMKL_PARALEL=TRUE -DTARGET_ARCHITECTURE=skylake -DGRID=35 -DT_RES=500 -DT_ROT=1.0 -DT_C=-2.0 -DW_START=-3 -DW_END=3 -DW_RES=4 -DUSE_MAGMA=FALSE -DVERBOSE=TRUE  -DUSE_GPU=FALSE -DTIME_EVOLUTION=TRUE -DCMAKE_PREFIX_PATH="~/boost_1_72_0/"
>>>>>>> 681774e4aadd278002bf4b37513dab9bf4d6039b
make -j 8
echo $path
echo "Running Program"
./TST
<<<<<<< HEAD
echo "finished"
=======
>>>>>>> 681774e4aadd278002bf4b37513dab9bf4d6039b
