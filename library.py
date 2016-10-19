import os,time

def clusterController(scratchDir,maxJobs,waitingTime):

    '''
    this function interrogates the cluster and submit +1 job of total capacity (4+1 in our case)
    '''

    greenLight=False

    while greenLight == False:

        # f.1. interrogating cluster
        statusFile=scratchDir+'status.txt'
        os.system('qstat > %s'%statusFile)

        # f.2. retrieving number of jobs
        currentJobs=0
        f=open(statusFile,'r')
        lines=f.readlines()
        f.close()
        for line in lines:
            vector=line.split()
            if len(vector) > 3:
                identity=vector[3]
                if identity == 'alomana':
                    currentJobs=currentJobs+1

        # f.3. deciding what to do (with my life)
        print('%s jobs found'%(str(currentJobs)))
        if currentJobs >= maxJobs:
            time.sleep(waitingTime)
        else:
            greenLight=True
            
    return None
