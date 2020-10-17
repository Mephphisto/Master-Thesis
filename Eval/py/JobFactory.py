import os
header = """#$-q idefix.q,obelix.q
#$ -l h_vmem=1G
#$ -l h_cpu=672:00:00
#$ -pe smp 1
#$-M j.teuffel@icloud.com -m ase
#$ -cwd
"""

try:
    os.mkdir("Jobs")
except:
    pass
for i in range(10):
    file = open("Jobs/job" + str(i) + ".sh", "w+")
    file.write(header + "echo \" job" + str(i) + "\" ")
    file.close()
    os.system("qsub Jobs/job" +str( i) + ".sh")
