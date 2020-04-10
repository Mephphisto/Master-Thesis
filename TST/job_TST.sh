#$-q idefix.q,obelix.q
#$ -l h_vmem=2G
#$ -l h_cpu=168:00:00
#$ -pe smp 48
#$-M j.teuffel@icloud.com -m ase
#$ -cwd
export BOOST_ROOT="~/boost_1_72_0/boost"
git commit -am " Prepare Run $date"
git pull
module load intel/2020
mkdir "run Time Eval"
cd "run Time Eval"
cmake cmake .. -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++  -DTHREADS=48 -DGRID=85 -DT_RES=800 -DUSE_MAGMA=FALSE -DVERBOSE="" -DUSE_GPU=FALSE -DTIME_EVOLUTION=TRUE
make -j 4
echo "Running Program"
./TST
git add .
git commit -am "Run finished $date"
git push
