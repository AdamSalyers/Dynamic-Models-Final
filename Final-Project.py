###########################################################################
# March 2019, Orit Peleg, orit.peleg@colorado.edu
# Code for HW3 CSCI 4314/5314 Dynamic Models in Biology
###########################################################################

import numpy as np
import math
import matplotlib.pyplot as plt

class flock():
    def init_wrm(self, ln_wrm, p):
        directions = np.array(['N','S','E','W'])
        #We want to add segments to the worm in random directions
        my_wrms = []#a matrix of worm locations
        for pos in p:
            pos = pos.tolist()
            #print("Pos is:", pos)
            visited = [pos] #keep track of position already in the worm and to will be the list of the final worm.
            last = pos #Set our starting point
            #print("visited is: ", visited)
            ln = 1
            while ln < ln_wrm: #loop through until we have reached the length of our worm
                np.random.shuffle(directions)
                for direction in directions:#randomly shuffle the directions
                    if direction == 'N':
                        if [last[0]+1, last[1]] not in visited:
                            ln += 1
                            visited.append([last[0]+1, last[1]])#If the point is not already visited we add it 
                            break
                    elif direction == 'S':
                        if [last[0]-1, last[1]] not in visited:
                            ln += 1
                            visited.append([last[0]-1, last[1]])
                            break
                    elif direction == 'E':
                        if [last[0], last[1]+1] not in visited:
                            ln += 1
                            visited.append([last[0], last[1]+1])
                            break
                    elif direction == 'W':
                        if [last[0], last[1]-1] not in visited:
                            ln += 1
                            visited.append((last[0], last[1]-1))
                            break
            my_wrms.append(visited)#append our worm to the list of worms
        return my_wrms
                
    def upd_wrm(self, wrm_pos, v, L):
        for i in wrm_pos:
            pos = i[0][0]#Get cordinates of head of worm
            pos += v#add v then make sure we stay in bounds
            if pos[0][0] > L/2: pos[0] -= L
            if pos[0][0] < -L/2: pos[0] += L
            if pos[0][1] > L/2: pos[1] -= L
            if pos[0][1] < -L/2: pos[1] += L
            wrm_pos.pop()
            wrm_pos.insert(0,pos)
        return wrm_pos
        
    def flocking_python(self,c1=0.00001,c2=0.01,c3=1,c4=0.01):
        N = 400 #No. of Worms
        frames = 10 #No. of frames
        limit = 100 #Axis Limits
        L  = limit*2
        ln_wrm = 5 #length of each worm
        mu = 0.0 #starting point of our cluster
        kappa = 10.0 #starting dispersion
        pfinal = []
        P = 10 #Spread of initial position (gaussian)
        V = 10 #Spread of initial velocity (gaussian)
        delta = 1 #Time Step
        #c1 = 0.00001 #Attraction Scaling factor
        #c2 = 0.01 #Repulsion scaling factor
        #c3 = 1 #Heading scaling factor
        #c4 = 0.01 #Randomness scaling factor
        vlimit = 1 #Maximum velocity

        #Initialize
        p = P*np.random.vonmises(mu,kappa, size=(N,2))#A random circular probabibility distribution
        wrm_pos = self.init_wrm(ln_wrm,p)#Get our initial Worm Positions
        v = V*np.random.randn(N,2)#Not really sure how we should set initial velocities
        #Initializing plot
        plt.ion()


        for i in range(0, frames):
            v1 = np.zeros((N,2))
            v2 = np.zeros((N,2))
            
            #YOUR CODE HERE
            #Calculate Average Velocity v3 
            trans = np.transpose(v)
            v3 = [sum(trans[0]),sum(trans[1])]*c3
            
            if (np.linalg.norm(v3) > vlimit): #limit maximum velocity
                v3 = v3*vlimit/np.linalg.norm(v3)

            for n in range(0, N):
                for m in range(0, N):
                    if m!=n:
                        #YOUR CODE HERE
                        #Compute vector r from one agent to the next
                        #print("p[m] is",p[m])
                        r = p[m]-p[n]
                        
                        if r[0] > L/2:
                            r[0] = r[0]-L
                        elif r[0] < -L/2:
                            r[0] = r[0]+L

                        if r[1] > L/2:
                            r[1] = r[1]-L
                        elif r[1] < -L/2:
                            r[1] = r[1]+L

                        #YOUR CODE HERE
                        #Compute distance between agents rmag
                        rmag = math.sqrt(r[0]**2+r[1]**2)
                        #Compute attraction v1
                        v1[n] += c1*r
                        #Compute Repulsion [non-linear scaling] v2
                        v2[n] -= (c2*r)/(rmag**2)
                   
                #YOUR CODE HERE
                #Compute random velocity component v4
                v4 = c4*np.random.randn(1,2)
                #Update velocity
                #print("Shapes: v1:",v1[n].shape, "v2:",v2[n].shape,"v3:",v3.shape,v)
                v[n] = v1[n]+v2[n]+v3+v4
                
            #YOUR CODE HERE
            #Update position
            v *= delta
            wrm_pos = self.upd_wrm(wrm_pos, v, L)#update worm position
            
            plt.xlim(-limit, limit)
            plt.ylim(-limit, limit)
            for worms in wrm_pos:
                plt.plot(worms)
            plt.show()

            pfinal.append(wrm_pos)
        return pfinal


flock_py = flock()

c1dist = []

pfinal = flock_py.flocking_python()
    

