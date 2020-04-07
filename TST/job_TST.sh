#$-q obelix.q
#$ -l h_vmem=6G
#$ -l h_cpu=168:00:00
#$ -pe smp 96
#$-M j.teuffel@icloud.com -m ase
#$ -cwd
git commit -am " Prepare Run $date"
git pull
module load intel/2020
mkdir "run $date"
cd run
cmake cmake .. -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++  -DTHREADS=96 -DGRID=100 -DT_RES=800 -DUSE_MAGMA=FALSE -DVERBOSE=""
make -j 16
echo "Running Program"
./TST
git add .
git commit -am "Run finished $date"
git push
