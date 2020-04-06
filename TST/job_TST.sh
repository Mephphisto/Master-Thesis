#$ -q idefix.q
#$ -M j.teuffel@icloud.com and -m b e a
#$ -l h_vmem=4G
#$ -l h_cpu=48:00:00
#$ -pe smp 4
git commit -am " Prepare Run $date"
git pull
module load intel
mkdir "build"
cd build
cmake cmake .. -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++  -DTHREADS=4 -DGRID=10 -DT_RES=800 -DUSE_MAGMA=FALSE
make -j 16
./TST
git commit -am "Run finished $date"
git push