#$-q idefix.q,obelix.q
#$ -l h_vmem=3G
#$ -l h_cpu=999:00:00
#$ -pe smp 85
#$-M j.teuffel@icloud.com -m ase
#$ -cwd
path="rune_G35_Res500_NegMu"
export BOOST_ROOT="~/boost_1_72_0/"
echo $path
module load intel/2020
mkdir $path
cd $path
cmake cmake .. -DCMAKE_BUILD_TYPE=Release  -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++  -DTHREADS=85 -DGRID=85 -DT_RES=85 -DUSE_MAGMA=FALSE -DVERBOSE=FALSE  -DUSE_GPU=FALSE -DTIME_EVOLUTION=FALSE -DCMAKE_PREFIX_PATH="~/boost_1_72_0/"
make -j 8
echo "Running Program"
./TST
echo "finished"
