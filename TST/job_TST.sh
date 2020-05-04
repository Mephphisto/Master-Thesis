#$-q idefix.q,obelix.q
#$ -l h_vmem=2G
#$ -l h_cpu=672:00:00
#$ -pe smp 48
#$-M j.teuffel@icloud.com -m ase
#$ -cwd
export BOOST_ROOT="~/boost_1_72_0/"
git commit -am " Prepare Run $date"
git pull
module load intel/2020
mkdir "run Time Eval"
cd "run Time Eval"
cmake cmake .. -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++  -DTHREADS=45 -DGRID=35 -DT_RES=100 -DUSE_MAGMA=FALSE -DVERBOSE=FALSE  -DUSE_GPU=FALSE -DTIME_EVOLUTION=TRUE -DCMAKE_PREFIX_PATH="~/boost_1_72_0/" -DCMAKE_IGNORE_PATH="/opt/sw/"
make -j 4
echo "Running Program"
./TST
git add .
git commit -am "Run finished $date"
git push
