PRO_seq project for Khalil lab measuring synthetic transcription factors. 

# submitting a job to the scc, the most basic version
# put what commands you want to run in a file, save it with .sh extension
# note: there are job options you can specify w/in your .sh file, there are BU SCC docs for this online
# to submit a job:

qsub -P PROJECTNAME file.sh

# PROJECTNAME = khalil, for example
# can specify full path to file.sh to run it from anywhere

# to check job status:

qstat -u yourUsername 

# this will display your submitted jobs