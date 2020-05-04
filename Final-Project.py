###########################################################################
# March 2019, Orit Peleg, orit.peleg@colorado.edu
# Code for HW3 CSCI 4314/5314 Dynamic Models in Biology
###########################################################################

import numpy as np
import math
import matplotlib.pyplot as plt
import random


class flock():
    def init_wrm(self, ln_wrm, p):
        directions = np.array(['N','S','E','W'])
        #We want to add segments to the worm in random directions
        my_wrms = []#a matrix of worm locations
        for pos in p:
            #print(pos)
            
            pos = pos.tolist()
            #print("Pos is:", pos)
            visited = [pos] #keep track of position already in the worm and to will be the list of the final worm.
            #print("visited is: ", visited)
            ln = 1
            while ln < ln_wrm: #loop through until we have reached the length of our worm
                np.random.shuffle(directions)
                for direction in directions:#randomly shuffle the directions
                    if direction == 'N':
                        if [visited[-1][0]+1, visited[-1][1]] not in visited:
                            ln += 1
                            visited.append([visited[-1][0]+1, visited[-1][1]])#If the point is not already visited we add it 
                            break
                    elif direction == 'S':
                        if [visited[-1][0]-1, visited[-1][1]] not in visited:
                            ln += 1
                            visited.append([visited[-1][0]-1, visited[-1][1]])
                            break
                    elif direction == 'E':
                        if [visited[-1][0], visited[-1][1]+1] not in visited:
                            ln += 1
                            visited.append([visited[-1][0], visited[-1][1]+1])
                            break
                    elif direction == 'W':
                        if [visited[-1][0], visited[-1][1]-1] not in visited:
                            ln += 1
                            visited.append([visited[-1][0], visited[-1][1]-1])
                            break
            
            my_wrms.append(visited)#append our worm to the list of worms
        return my_wrms
                
    def upd_wrm(self, wrm_pos, v, L):
        print("Original max position: ", wrm_pos[0][0])
        
        for i in range(0,len(wrm_pos)):
            new_head = [0.0, 0.0]

            #used as a temporary variable to work on the current worms positioning
            curr_worm = wrm_pos[i]

            #Get cordinates of head of worm
            curr_worm_head = curr_worm[0]
            #print(curr_worm_head)
            
            #print("velocity change for x:", v[i][0])
            new_head[0] = curr_worm_head[0] + v[i][0]#add v then make sure we stay in bounds
            new_head[1] = curr_worm_head[1] + v[i][1]#add v then make sure we stay in bounds
            
            #make checks on the boundaries
            if new_head[0] > L/2: new_head[0] -= L
            if new_head[0] < -L/2: new_head[0] += L
            if new_head[1] > L/2: new_head[1] -= L
            if new_head[1] < -L/2: new_head[1] += L

            #remove the last position of the current worm, i.e. the place its moving from
            curr_worm.pop()
            #add the new head position to the front of the current worm
            curr_worm.insert(0,new_head)
            #update position in array that is returned
            wrm_pos[i] = curr_worm
        print("New head: ",wrm_pos[0][0])
        
        return wrm_pos
    

            
    def check_to_branch(self, wrm_pos, max_branches, branches, frame):
        x_avg = 0.0
        y_avg = 0.0
        all_x = []
        all_y = []
        worm_branch_index_counter = 0
        print("max branches: ", max_branches)
        print("frame: ", frame)
        if (len(branches) >=  max_branches or frame <= 10):
            return branches, (x_avg, y_avg)
        else:
            for i in range(0, len(wrm_pos)):
                #get all of worms x position
                wx = [wrm_pos[i][0][0],wrm_pos[i][1][0],wrm_pos[i][2][0],wrm_pos[i][3][0],wrm_pos[i][4][0]]
                #get all of worms y position
                wy = [wrm_pos[i][0][1],wrm_pos[i][1][1],wrm_pos[i][2][1],wrm_pos[i][3][1],wrm_pos[i][4][1]]
                #print("sum: ", wx)
                #print("avg: ", sum(wx)/5)
                #average the worms points to find center of mass of worm
                all_x.append(sum(wx)/5)
                all_y.append(sum(wy)/5)
            #average all worms to find center of mass of blob
            x_avg = sum(all_x)/len(wrm_pos)
            y_avg = sum(all_y)/len(wrm_pos)
            #find differnce in ranges between x and y, this is the rough radius
            #average the two to find a better approximation and divide in two for radius
            max_x = max(all_x)
            min_x = min(all_x)
            #print("x diff: ", max_x-min_x)
            max_y = max(all_y)
            min_y = min(all_y)
            #print("y diff: ", max_y-min_y)
            #drb = distance from center of mass of blob required to branch
            drb = (sum([max_x-min_x, max_y-min_y])/4)
            #print("drb: ", drb)
            #loop through all worms to see if they are able to branch
            for tempx,tempy in zip(all_x, all_y):
                incir = self.in_circle(x_avg, y_avg, drb, tempx, tempy)
                #branches that pass this if statement are able to branch and are appened to branches
                if incir == False and len(branches) < max_branches and worm_branch_index_counter not in branches:
                    branches.append(worm_branch_index_counter)
                worm_branch_index_counter += 1
        
            return branches, (x_avg, y_avg)
        
    def branching_velocity(self, wrm, center):
        #get all of worms x position
        wx = [wrm[0][0],wrm[1][0],wrm[2][0],wrm[3][0],wrm[4][0]]
        x_pos = sum(wx)/5
        #get all of worms y position
        wy = [wrm[0][1],wrm[1][1],wrm[2][1],wrm[3][1],wrm[4][1]]
        y_pos = sum(wy)/5
        print("center: ",center)
        print("c x: ",center[0])
        print("c y: ",center[1])
        #get a random number between 0.5 and 1.0, this is what the branching worms velocity will sum to
        total_vel = random.uniform(0.5,1.0)
        print("total velocity: ",total_vel)
        #find rise and run between worm center of mass and worm blob center of mass then sum
        pos_sum = ((x_pos-center[0])+(y_pos-center[1]))
        #divide total velocity by pos_sum
        ratio = total_vel/pos_sum
        #result is a velocity in the opposite direction of the worm blob center of mass
        x_vel = x_pos*ratio
        y_vel = y_pos*ratio
        print("returned array: ", np.array([[x_vel, y_vel]]))
        return np.array([[x_vel, y_vel]])
        
    def in_circle(self, center_x, center_y, radius, x, y):
        square_dist = (center_x - x) ** 2 + (center_y - y) ** 2
        return square_dist <= radius ** 2
        
                
    
    def flocking_python(self,c1=0.001,c2=0.01,c3=1,c4=.90):
        N = 200 #No. of Worms
        frames = 20 #No. of frames
        limit = 40 #Axis Limits
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
        
        max_branches = int(N*0.10)
        branches = []

        #Initialize
        p = P*np.random.vonmises(mu,kappa, size=(N,2))#A random circular probabibility distribution
        wrm_pos = self.init_wrm(ln_wrm,p)#Get our initial Worm Positions
        print(len(wrm_pos))
        print(len(wrm_pos[0]))
        print(len(wrm_pos[0][0]))
        v = V*np.random.randn(N,2)#Not really sure how we should set initial velocities
        #Initializing plot
        plt.ion()
        
            
        for i in range(0, frames):
            v1 = np.zeros((N,2))
            v2 = np.zeros((N,2))
            
            
            #Calculate Average Velocity v3 
            trans = np.transpose(v)
            v3 = [sum(trans[0]),sum(trans[1])]*c3
            
            if (np.linalg.norm(v3) > vlimit): #limit maximum velocity
               v3 = v3*vlimit/np.linalg.norm(v3)
            
            branches, blob_center = self.check_to_branch(wrm_pos, max_branches, branches, i)
            print("branches length: ", len(branches))
            print("branches: ", branches)
            for n in range(0, N):
                p = np.array(wrm_pos[n])
                pn = p.mean(axis=0)
                for m in range(0, N):
                    p = np.array(wrm_pos[m])
                    pm = p.mean(axis=0)
                    if m!=n:
                        
                        #Compute vector r from one agent to the next
                        #print("p[m] is",p[m])
                        r = pm-pn
                        
                        if r[0] > L/2:
                            r[0] = r[0]-L
                        elif r[0] < -L/2:
                            r[0] = r[0]+L

                        if r[1] > L/2:
                            r[1] = r[1]-L
                        elif r[1] < -L/2:
                            r[1] = r[1]+L
                            
                        #Compute distance between agents rmag
                        rmag = math.sqrt(r[0]**2+r[1]**2)
                        #Compute attraction v1
                        v1[n] += c1*r
                        #Compute Repulsion [non-linear scaling] v2
                        v2[n] -= (c2*r)/(rmag**2)
                   
                #Compute random velocity component v4
                v4 = c4*np.random.randn(1,2)
                
                #Update velocity
                #print("Shapes: v1:",v1[n].shape, "v2:",v2[n].shape,"v3:",v3.shape,v)
                if n in branches:
                    v[n] = self.branching_velocity(wrm_pos[n], blob_center)
                    print("branching v[n]: ", v[n])
                else:
                    v[n] = v1[n]+v2[n]+v3+v4
                    print("non-branching v[n]: ", v[n])
                
            #Update position
            v *= delta
            wrm_pos = self.upd_wrm(wrm_pos, v, L)#update worm position
            #print(len(wrm_pos[0]))
            
            plt.xlim(-limit, limit)
            plt.ylim(-limit, limit)
            for worms in wrm_pos:
                x = [worms[0][0],worms[1][0],worms[2][0],worms[3][0],worms[4][0]]
                y = [worms[0][1],worms[1][1],worms[2][1],worms[3][1],worms[4][1]]
                plt.plot(x, y)
            plt.show()
            
            pfinal.append(wrm_pos)
        return pfinal


flock_py = flock()

c1dist = []

pfinal = flock_py.flocking_python()
