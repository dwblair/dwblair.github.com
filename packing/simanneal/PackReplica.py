import numpy as np
import random
import time
import os
import sys


class PackReplica:
    def __init__(self,L,N,outFileNum,pressure,randomSeed,initialDensity):

        self.pressure=pressure
        self.randomSeed=randomSeed
        self.L=L
        self.N=N
        self.outFileNum=outFileNum

        # assume particles are of radius "1"
        self.r=np.sqrt(initialDensity*L**2/(N*3.14159))
        self.V=self.N*np.pi/initialDensity

        self.density=initialDensity

        self.t=0
        self.maxt=1000000

        self.coords=np.zeros((N,2))+L/2.

        # for MSD
        self.prevCoords=np.copy(self.coords)
        self.MSD=0.

        # for assessing equilibrium
        #self.prevDensity=self.density
        
        self.colors=np.zeros(N)
        self.cellNums=np.zeros(N)
        self.maxNumberParticlesPerCell=10  # assumes no more than 10 particles per cell! may need to modify if polydisperse ...
        #self.radii=np.zeros(N)+self.r
        #self.stepSize=self.r*.1
        #self.rScale=1.
        self.cellSize=10


        #variables for MSD calc
        self.prevCoords=np.zeros((N,2))

        #variables for assessing equilibrium
        self.equilibrated=False
        self.densitySTD=1. 

        # params for translation algorithm
        self.mu = 0.
        self.translationStep=self.r*2
        self.expansionStep=self.V*.1
        
        self.ratioTranslateExpansion=100

        self.transAcceptRatio=0.
        self.expanAcceptRatio=0.
        self.expanRatioCount=0.
        self.transRatioCount=0.
        self.PIDfactor=0.1
        self.idealAcceptRatio=.4
        self.ratioStats=np.zeros(2)
        self.requiredSweepsToAchieveOptimalAcceptRatio=0
        
        #density stats
        self.highestDensity=0.

        #random seed
        random.seed(self.randomSeed)


    def attemptTranslate(self,i): #random translate, particle i
        self.transRatioCount=self.transRatioCount+1 
        oldCoords=np.copy(self.coords[i])
        dr=np.random.normal(self.mu,self.translationStep,2) #array size 2 -- one random value for ea coord
        self.coords[i]=self.coords[i]+dr

        if (self.checkIntersections(i))>0:
            self.coords[i]=oldCoords
        else:
            self.transAcceptRatio=self.transAcceptRatio+1

    def attemptExpansion(self):
        #particles=np.arange(0,self.N)
        #result=self.checkIntersections(particles)

        #print "attempt expand: r=",self.r

        
        #dVScale=np.random.normal(self.mu,self.expansionStep)+1.
        
        #dVScale=(2*random.random()-1)*self.expansionStep+1.
        #dVScale=(2*random.random()-1)
        
        #print dVScale
        #dV=(dVScale-1.)*self.V*self.expansionStep
        
        #make sure expansionStep isn't too big:
        #while (self.V-self.expansionStep)<0:
            # then expansStep is too large
        #    print "#oops!", "self.V=",self.V, "self.expansionStep=",self.expansionStep
        #   self.expansionStep=self.expansionStep*.9
         #   
        #dV=(2*random.random()-1)*self.expansionStep

        
        #while (self.V-dV)<0:
        dV=np.random.normal(self.mu,self.expansionStep,1)[0]
        
        
        expansionFlag=False
        if (dV < 0. and dV > -self.V): # is this second conditional biasing the system??
            expansionFlag=True
        if (dV > 0.):
            #print dV, (1-np.exp(self.pressure*dV))
            #print dV, np.exp(-self.pressure*dV)
            p=random.random()
            boltz=min(1,np.exp(-self.pressure*dV))

            #print "pressure=",self.pressure,"; p=",p,"; boltz=",boltz 
            if p<boltz:
                expansionFlag=True
                #print "growing!"

        if expansionFlag==True:
            #print "Expansion!"
            self.expanRatioCount=self.expanRatioCount+1 
            oldVolume=self.V
            newVolume=self.V+dV
            
            #volumeRatio=newVolume/oldVolume
            self.V=newVolume
            oldR=self.r

            oldDensity=self.N*3.14159*self.r**2/self.L**2
            newDensity=self.N*3.14159/self.V #assumes particle radius 1
            #print "oldDensity=",oldDensity,"; newDensity=",newDensity

            if newDensity*self.L**2/(self.N*3.14159) <0:
                print "less than zero error!"
                print "V=",self.V,"dV=",dV

            self.r=np.sqrt(newDensity*self.L**2/(self.N*3.14159))
            #self.r=np.sqrt(self.L/self.V)
            #print "2=",np.sqrt(4.)
            #print "oldRatio=",np.sqrt(self.L/oldVolume)
            
            #print "oldR=",oldR,"; newR=",self.r
            
            #self.r=np.sqrt(oldR**2.)*oldVolume/newVolume
            
            intersections=0
            for i in range(0,self.N):
                intersections=intersections+self.checkIntersections(i)

            if intersections>0:
                self.r=oldR
                self.V=oldVolume
            else: # we accepted
                self.expanAcceptRatio=self.expanAcceptRatio+1
                #self.acceptRatioMeasureLength=1000.
                #self.acceptRatioCount=0
                
    
    def getDensity(self):
        self.density=self.N*3.14159*self.r**2/self.L**2
        if self.density>self.highestDensity:
            self.highestDensity=self.density
        return self.density
    
    def sweep(self,numSweeps):
        #translation attempts

        # assess density over range of numSweeps
        prevDensities=np.zeros(numSweeps)
        
        for s in range(0,numSweeps):
            for q in range(0,self.ratioTranslateExpansion): #makes sure we have a sweep of as many expansion attempts as we do the ratio between expansion attempts and translations
                for p in range(0,self.ratioTranslateExpansion):
                    i=random.randint(0,self.N-1)
                    self.attemptTranslate(i)
                #expansion attempts
                self.attemptExpansion()
            self.t=self.t+1 #however we define a sweep, we count t in 'sweeps'
            prevDensities[s]=self.getDensity()
        
        # update the MSD
        self.MSD=self.getMSD()
        # update the step sizes using the PID loop
        self.updatePID() #modify the translation and expansion steps to get closer to ideal acceptance ratio

        # assess Equilibrium
        #print prevDensities
        self.densitySTD=np.std(prevDensities)
        #print self.t, self.PIDFactor, densitySTD
        
        
        
    def getMSD(self):
        # calcuates MSD, normalized by R**2
        #coordDiff=np.zeros((N,2))
        coordDiff=self.coords-self.prevCoords
        #print "coordDiff=",coordDiff
        diffX=coordDiff[:,0]
        diffY=coordDiff[:,1]
        SD=diffX*diffX+diffY*diffY
        MSD=np.sum(SD)/float(self.N)/(self.r*self.r) #normalized by radius of particle
    
        #update prevCoords
        self.prevCoords=np.copy(self.coords)

        #return value
        return MSD
     
    def updatePID(self):
        
        self.transAcceptRatio=float(self.transAcceptRatio)/float(self.transRatioCount)
        self.expanAcceptRatio=float(self.expanAcceptRatio)/float(self.expanRatioCount)

        
        print self.t, self.PIDfactor, self.pressure, self.getDensity(), self.transAcceptRatio,self.expanAcceptRatio, self.translationStep, self.expansionStep, self.getMSD(), self.V

        #update the step sizes
        transDiff=self.transAcceptRatio-self.idealAcceptRatio
        expanDiff=self.expanAcceptRatio-self.idealAcceptRatio
        #self.PIDfactor=.1
        
        if transDiff>0: # accept ratio too high --> step too small
            self.translationStep=self.translationStep*(1.+self.PIDfactor)
            #print "accept too high; mod=",(1.+factor)
        if transDiff<0:  # accept ratio too low -->  step too big
            self.translationStep=self.translationStep*(1.-self.PIDfactor)

        if expanDiff>0: # accept ratio too high --> step too small
            self.expansionStep=self.expansionStep*(1.+self.PIDfactor)
        if  expanDiff<0: # accept ratio too low -->  step too big
            self.expansionStep=self.expansionStep*(1.-self.PIDfactor)
        #if self.expansionStep>.5:
        #    self.expansionStep=.5
            
        #reset values
        self.ratioStats[0]=self.transAcceptRatio
        self.ratioStats[1]=self.expanAcceptRatio
        self.transRatioCount=0
        self.expanRatioCount=0
        self.transAcceptRatio=0
        self.expanAcceptRatio=0
        
        
    def checkIntersections(self,i):
        neighborIndices=np.concatenate((np.arange(0,i,1),np.arange(i+1,self.N,1)),axis=0)
        neighbors=self.coords[neighborIndices]
        dx=neighbors[:,0]-self.coords[i][0]
        dy=neighbors[:,1]-self.coords[i][1]
        r2=dx*dx+dy*dy
        # want to change this next to incorporate the diameter = radii[i]+radii[j], not global "r"
        #diams=self.radii*self.rScale
        
        intersectionIndices=r2<((2*self.r)**2) #i.e. less than a diameter
        numIntersections=len(r2[intersectionIndices])

        x=self.coords[i][0]
        y=self.coords[i][1]

        # in a box
        if x<(0+self.r) or x>(self.L-self.r):
            numIntersections=numIntersections+1
        if y<(0+self.r) or y>(self.L-self.r):
            numIntersections=numIntersections+1
        return numIntersections


        
    def initialize(self):

        #write out the params.in file
        p=open('params.in','w')
        thisline="#N,L,r\n"
        p.write(thisline)
        thisline=str(self.N)+","+str(self.L)+","+str(self.r)
        p.write(thisline)
        
        # write out a params.in file
        
        # place particles in a grid
        factor=1.5
        x=factor*2*self.r
        y=factor*2*self.r

        for i in range(0,self.N):
            self.coords[i][0]=x;
            self.coords[i][1]=y;
            x=x+factor*2*self.r
            if x>(self.L-factor*self.r):
                x=factor*2*self.r
                y=y+factor*2*self.r
            #self.radii[i]=self.r

        
    def printOut(self):

        #the coords file
        coordsOutFileName="N_"+str(self.N)+"_coords_"+str("%03d" % self.outFileNum)+".txt"
        densityOutFileName="N_"+str(self.N)+"_stats_"+str("%03d" % self.outFileNum)+".txt"
        #densityOutFileName="density"+str("%03d" % self.outFileNum)+".txt"
        
        f=open(coordsOutFileName,'a')

        for i in range(0,self.N):
            thisline=str(self.coords[i][0])+","+str(self.coords[i][1])+","+str(self.r)+","+str(self.cellNums[i])+"\n"
            f.write(thisline)

        thisline="@ "+str(self.t)+" "+str(self.cellSize)+" "+str(self.pressure)+" "+str(self.getDensity())+"\n"
        f.write(thisline)
        f.close()

        f=open(densityOutFileName,'a')
        
        thisline=str(self.t)+" "+str(self.N)+" "+str(self.pressure)+" "+str(self.getDensity())+" "+str(self.densitySTD)+" "+str(self.MSD)+" "+str(self.randomSeed)+" "+str(self.PIDfactor)+" "+str(self.ratioStats[0])+" "+str(self.ratioStats[1])+" "+str(self.V)+" "+str(self.requiredSweepsToAchieveOptimalAcceptRatio)+" "+str(self.highestDensity)+"\n"

        f.write(thisline)
        f.close()


#res = os.system(sys.argv[1], sys.argv[2])
#print sys.argv[1]

#N=sys.argv[1]

# run as "python PackReplica.py numSweeps tmax numSims N"
numSweeps=int(sys.argv[1])
tmax=int(sys.argv[2])
numSims=int(sys.argv[3])
N=int(sys.argv[4])

#simNum=1
#pressure=.1
#pressureList=[10.,3.,2.,1.,.5,.1,0.1,.01,.001]
#NList=[5,6,7,8,10,11]
simNum=0
#numSims=100
#pressureList=[.01, .1, .2, .4, .8, 2.,4.,8.,16.]
invplist=np.linspace(40,.03,50)
plist=1./invplist
pressureList=plist[35:]

randomSeed=0
initialDensity=.05
#N=11
L=2**9
#tmax=2000
#numSweeps=100
#while (pressure<10000):
#pr.initialize()

#for N in NList:
while simNum<numSims:
    
    pressure=pressureList[0]
    pr = PackReplica(L=L,N=N,outFileNum=simNum,pressure=pressure,randomSeed=randomSeed,initialDensity=initialDensity)
    pr.initialize()

    
    for pressure in pressureList:
        
        pr.pressure=pressure

        #equilibrate at first
        print "equilibrating ..."
        acceptRatiosGood=False

        # take a snapshot of the simulation before we let things evolve
        tSnap=pr.t
        pressureSnap=pr.pressure
        rSnap=pr.r
        VSnap=pr.V
        coordsSnap=np.copy(pr.coords)
        prevCoordsSnap=np.copy(pr.prevCoords)

             

        thisT=0
        while(acceptRatiosGood==False):
            pr.sweep(numSweeps=numSweeps)
            #pr.printOut()

            tr=pr.ratioStats[0] #translation accept ratio
            ex=pr.ratioStats[1] #expansion accept ratio
            print tr,ex, pr.MSD
            if (tr>.3 and tr<.5 and ex>.3 and ex<.5):
                acceptRatiosGood=True
                print "finished equilibrating in ", thisT, " steps"
                pr.requiredSweepsToAchieveOptimalAcceptRatio=thisT
                #thisT=tmax
            thisT=thisT+numSweeps

        #now begin actual simulation at this pressure value
        print "beginning serious sim!"

        #first retrieve original configurations
        pr.t=tSnap
        pr.pressure=pressureSnap
        pr.r=rSnap
        pr.V=VSnap
        pr.coords=np.copy(coordsSnap)
        pr.prevCoords=np.copy(prevCoordsSnap)

        #now carry out the sim
        thisT=0
        while (thisT<tmax):
            pr.sweep(numSweeps=numSweeps)
            pr.printOut()
            thisT=thisT+numSweeps
        
        #print pr.t, pr.PIDfactor, pr.pressure, pr.getDensity(), pr.transAcceptRatio,pr.expanAcceptRatio, pr.translationStep, pr.expansionStep, pr.getMSD(), pr.V

    randomSeed=randomSeed+1
    simNum=simNum+1
