import os
header = """#$-q idefix.q,obelix.q
#$ -l h_vmem=3G
#$ -l h_cpu=672:00:00
#$ -pe smp 85
#$-M j.teuffel@icloud.com -m ase
#$ -cwd)
"""

os.mkdir("Jobs")
for i in range(10):
    file = open("Jobs/job" + i + ".sh", "w+")
    file.write(header + "echo \" job" + i + "\" ")
    file.close()
    os.system("qsub Jobs/job" + i + ".sh")
