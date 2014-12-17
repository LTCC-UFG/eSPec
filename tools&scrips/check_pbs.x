#PBS -S /bin/sh         
#PBS -N eSPec tests          
#PBS -l cput=3:0:0      
#PBS -l nodes=opteron        
#PBS -m ae              
cd $PBS_O_WORKDIR      

time ./runtests.x


