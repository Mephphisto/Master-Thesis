#$-q idefix.q,obelix.q
#$ -l h_vmem=1G
#$ -l h_cpu=672:00:00
#$ -pe smp 35
#$-M j.teuffel@icloud.com -m ase
#$ -cwd
path="run_Time_Eval_Release_G35_Res500"
export BOOST_ROOT="~/boost_1_72_0/"
module load intel/2020
mkdir $path
cd $path
cmake cmake .. -DCMAKE_BUILD_TYPE=Release  -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++  -DTHREADS=35 -DGRID=85 -DT_RES=500 -DUSE_MAGMA=FALSE -DVERBOSE=FALSE  -DUSE_GPU=FALSE -DTIME_EVOLUTION=TRUE -DCMAKE_PREFIX_PATH="~/boost_1_72_0/"
make -j 8
echo "Running Program"
./TST
git add .
git commit -am "Run finished $date"
git push
