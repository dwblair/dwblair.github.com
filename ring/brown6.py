from numpy import *
import random

#### general brownian dynamics params #####

L=700 #system size, compatible with size of display in processing
N=10 #number of particles
r=10 #particle radius
t=0 #time
maxt=1000000 #max time
stepSize=r/10. #random step size
dt=.1 #timestep
gamma=5. #friction coeff
kbt=.1 # KbT
MAXFORCE=10000000. #in case something goes wrong
timeGapForPrintout=100 #number of timesteps to skip before printing out coordinates

# radial position
rad=L/4.
# spring constant
k=1.

##### lennard jones params #####
epsilon=1. #depth of lennard jones potential
s=r # 'radius' of potential well
s6=s**6 
s12=s6**2

#### gaussian random number parameters ######
sigma=1.
mu=0.

#coordinates array
coords=zeros((N,2))+L/2. #start all particles in the center

####### initially put coordinates in a grid #######
factor = 1.1
x=factor*2*r
y=factor*2*r

rad=L/4;
theta=0.

for i in range(0,N):
    x=L/2+rad*cos(theta);
    y=L/2+rad*sin(theta);
    coords[i][0]=x;
    coords[i][1]=y;
    theta=theta+(2*3.14159)/(N+1);


############ run the simulation #########
for t in range(0,maxt): #time loop

    for i in range(0,N): #loop over all particles

        x=coords[i][0]
        y=coords[i][1]

        LJx=0.
        LJy=0.

        for j in range(0,N): #loop over all particles except the ith particle
            if j!=i:
                
                ##### get the distance between particles
                dx=coords[j][0]-coords[i][0]
                dy=coords[j][1]-coords[i][1]
                
                dr=sqrt(dx**2+dy**2)
                dr7=dr**7
                dr13=dr**13
               
                ###### calculate the LJ force in x and y

                ## note -- the neighbors need to be changed to reflect periodic boundaries (if not on compact surface)
                LJ=-24*epsilon*(2*s12/dr13 - s6/dr7)

                
                
                if (LJ>MAXFORCE):
                    LJ=MAXFORCE

                LJx=LJx+(dx/dr)*LJ
                LJy=LJy+(dy/dr)*LJ

            # get the radial force

            dx=coords[i][0]-L/2.
            dy=coords[i][1]-L/2.
            r=sqrt(dx**2+dy**2)
            dr=r-rad

            k=.0001
            Spring=-k*dr
	    Springx=(-dx/r)*Spring
	    Springy=(-dy/r)*Spring

            #print LJx, Springx, LJy, Springy
            #print LJ
            #print  LJx/Springx, LJy/Springy
	    #Fx=LJx+Springx
	    #Fy=LJy+Springy

            kbt=0.1
            
            Fx=Springx
            Fy=Springy
            
            
        #### update the particle positions
        x=x+(dt/gamma)*Fx+sqrt(2*kbt*dt/gamma)*random.gauss(mu,sigma)
        y=y+(dt/gamma)*Fy+sqrt(2*kbt*dt/gamma)*random.gauss(mu,sigma)
        
        ##### enforce periodic boundaries (not on periodic surface)
        x=x%L
        y=y%L

        coords[i][0]=x
        coords[i][1]=y
        
        #output the particle positions
        if t%timeGapForPrintout==0:
            print coords[i][0],",",coords[i][1]

    #mark end of timestep
    if t%timeGapForPrintout==0:
        print "@", t
