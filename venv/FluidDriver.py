# Class for integrating fluid simulation
import numpy as np
import scipy as sp
import sys
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation


np.set_printoptions(threshold=sys.maxsize)


class Fluid:
    def __init__(self,Re,dt,T,Lx,Ly,Nx,Ny):
        self.Re=Re
        self.dt=dt
        self.T=T
        self.Lx=Lx
        self.Ly=Ly
        self.Nx=Nx
        self.Ny=Ny

        self.dx=Lx/(Nx-1)
        self.dy = Ly / (Ny - 1)

    def SetGrid(self):
        y=self.Ny
        x=self.Nx
        self.vort_n = np.zeros((y,x))
        self.stream_n= np.zeros((y,x))
        self.vort_n1 = np.zeros((y, x))
        self.stream_n1 = np.zeros((y, x))

    #def InitialConditions(self):
        #something

    def BoundaryConditions(self):

        for i in range(self.Nx):
            self.vort_n[0,i]=(self.stream_n[0,i] - self.stream_n[1,i]) * 2 / pow(self.dy,2)   # Top
            self.vort_n[self.Ny-1,i]=(self.stream_n[self.Ny-1,i] - self.stream_n[self.Ny-2,i]) * 2 / pow(self.dy,2)  - 0.01  # Bottom
        for j in range(self.Ny):
            self.vort_n[j, 0] = (self.stream_n[j,0] - self.stream_n[j,1]) * 2 / pow(self.dx,2)     # Left
            self.vort_n[j, self.Nx - 1] = (self.stream_n[j,self.Nx-1] - self.stream_n[j,self.Nx-2]) * 2 / pow(self.dx,2)   # Right

    def InteriorVorticityT(self):
        for i in range(1,self.Nx-1):
            for j in range(1,self.Ny-1):
                self.vort_n[j,i]=-((self.stream_n[j,i+1]-2*self.stream_n[j,i]+self.stream_n[j,i-1])/(pow(self.dx,2))+(self.stream_n[j+1,i]-2*self.stream_n[j,i]+self.stream_n[j-1,i])/(pow(self.dy,2)))

    def InteriorVorticityT1(self):
        for i in range(1,self.Nx-1):
            for j in range(1,self.Ny-1):
                self.vort_n1[j,i] = self.vort_n[j,i] + self.dt*(( (self.stream_n[j,i+1]-self.stream_n[j,i-1])*(self.vort_n[j+1,i]-self.vort_n[j-1,i])/(4*self.dx*self.dy) ) - ( (self.stream_n[j+1,i]-self.stream_n[j-1,i])*(self.vort_n[j,i+1]-self.vort_n[j,i-1])/(4*self.dx*self.dy) ) + 1/self.Re*( (self.vort_n[j,i+1]-2*self.vort_n[j,i]+self.vort_n[j,i-1]) / (pow(self.dx,2)) + (self.vort_n[j+1,i]-2*self.vort_n[j,i]+self.vort_n[j-1,i]) / (pow(self.dy,2)) ))

    def Poisson(self):


        # Matrix coefficients
        a=-1.0/(self.dx*self.dx)
        b=-1.0/(self.dy*self.dy)
        c=(2.0/(self.dx*self.dx)+2.0/(self.dy*self.dy))

        A=np.zeros((self.Ny*self.Nx,self.Ny*self.Nx),dtype='uint8')
        B=self.vort_n1.flatten()

        # Assembling A matrix
        for i in range(0,self.Ny*self.Nx):
            A[i, i] = c
            if i>0:
                A[i, i - 1] = a
            if i<self.Ny*self.Nx-1:
                A[i, i + 1] = a
            if i>self.Nx:
                A[i,i-self.Nx] = b
            if i<self.Nx*self.Ny-self.Nx:
                A[i,i+self.Nx] = b

        #print(A)
        stream=spsolve(A,B)

        # Have to check
        j=0
        count=0
        for i in stream:
            if count==self.Nx:
                j=j+1
                count=0
            self.stream_n1[j,count]=stream[count]
            count = count + 1

    def Integrate(self):

        # For plotting
        x = np.linspace(0, 1, self.Nx)
        y = np.linspace(0, 1, self.Ny)
        X, Y = np.meshgrid(x, y)

        tsteps=int(self.T/self.dt)
        t_real=0
        #fig = plt.subplots()

        self.BoundaryConditions()
        self.InteriorVorticityT()
        # -- Integrating over time
        for t in range(0,tsteps):
            self.InteriorVorticityT1()
            self.Poisson()
            # Insert back into n
            self.stream_n=self.stream_n1
            self.vort_n=self.vort_n1

            # Checking time
            t_real = t_real + self.dt
            t_real=round(t_real,2)
            plt.title("t=" + str(t_real) + " s")

            Z = self.vort_n
            plt.contourf(X, Y, Z, 20, cmap='RdGy')
            #plt.colorbar()
            plt.pause(0.1)




        plt.show()


            #FuncAnimation(fig,self.Contour,blit=True,repeat=True)

        #plt.draw()



    def PrintGrid(self):

        for i in self.vort_n:
            print(i)

    def Contour(self):
        plt.contour(self.vort_n)
        plt.show()

    def Test(self):
        x = np.linspace(0, 1, self.Nx)
        y = np.linspace(0, 1, self.Ny)
        X, Y = np.meshgrid(x, y)
        Z = self.vort_n
        plt.contourf(X, Y, Z, 20,cmap='RdGy')
        plt.colorbar()
        plt.show()


