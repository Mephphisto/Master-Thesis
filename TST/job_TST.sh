#$-q idefix.q,obelix.q
#$ -l h_vmem=2G
#$ -l h_cpu=48:00:00
#$ -pe smp 16
#$-M j.teuffel@icloud.com -m ase
git commit -am " Prepare Run $date"
git pull
module load intel
mkdir "run"
cd run
cmake cmake .. -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++  -DTHREADS=16 -DGRID=100 -DT_RES=800 -DUSE_MAGMA=FALSE
make -j 16
./TST
git add .
git commit -am "Run finished $date"
git push