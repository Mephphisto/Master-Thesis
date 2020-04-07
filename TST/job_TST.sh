#$-q idefix.q,obelix.q
#$ -l h_vmem=1G
#$ -l h_cpu=168:00:00
#$ -pe smp 32
#$-M j.teuffel@icloud.com -m ase
git commit -am " Prepare Run $date"
git pull
module load intel/2020
mkdir "run $date"
cd run
cmake cmake .. -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++  -DTHREADS=1 -DGRID=10 -DT_RES=1200 -DUSE_MAGMA=FALSE -DVERBOSE=""
make -j 16
echo "Running Program"
./TST
git add .
git commit -am "Run finished $date"
git push
