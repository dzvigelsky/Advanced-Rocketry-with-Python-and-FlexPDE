# Dennis Zvigelsky
# 2020

import subprocess
import scipy as sp
import matplotlib.pyplot as plt

# FlexPDE code with changing launch angle, theta
FlexCode = """TITLE 'OptimizingLaunchAngle'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
    xd (threshold=0.1)
    yd (threshold=0.1)
    vx (threshold=0.1)
    vy (threshold=0.1)
SELECT         { method controls }
	ngrid=1

DEFINITIONS    { parameter definitions }
    g = 9.81
    vi = 20
    v=sqrt(vx^2+vy^2)
    r=sqrt(xd^2+yd^2)
    theta=%s*pi/180
    
    mi =  1000 !inital mass
    mfuel=0.8*mi !mass of fuel
    mrocket = mi-mfuel !mass of rocket
    q=50 !fuel rate 
    vfuel=2e3 !fuel velocity 
    tfuel=mfuel/q
    
	m = if t<tfuel then mi-q*t else mrocket
    
    
    Fg = m*g !gravitational force
    
    !drag force
    vwindx=-25
    vwindy=10
    vrelx=vx-vwindx
    vrely=vy-vwindy
    vrel=sqrt(vrelx^2+vrely^2)
    Cd=0.6
    A=5
    rho = 2
    Fd=0.5*rho*Cd*A*vrel^2
    
    !thrust force
    Ft=if t<tfuel then q*vfuel else 0
    
    Fx = Ft*vx/v - Fd*vrelx/vrel
    Fy = -Fg + Ft*vy/v - Fd*vrely/vrel
    
    ax=Fx/m
    ay=Fy/m    
    
    KE=0.5*mrocket*(vx^2 + vy^2)

INITIAL VALUES
xd=0

vx=vi*cos(theta)
vy=vi*sin(theta)

EQUATIONS        { PDE's, one for each variable }
xd: dt(xd)=vx
yd: dt(yd)=vy
vx: dt(vx)=ax
vy: dt(vy)=ay

BOUNDARIES       { The domain definition }
  REGION 1       { For each material region }
    START(0,0)   { Walk the domain boundary }
    LINE TO (1,0) TO (1,1) TO (0,1) TO CLOSE

TIME 0 TO 600 halt yd<0    { if time dependent }

PLOTS            { save result displays }
	for t = 0 by endtime/1e3 to endtime
        history(xd,KE,yd) at (0,0) PrintOnly Export Format '#t#b#1#b#2#b#3' file = 'trajectory_data2.txt' 
SUMMARY 

END


"""

FlexFileName = "C:\\Users\\denni\\OneDrive\\Documents\\McMaster\\Winter 2020\\Computational Multiphysics - 2CM4\\M1-Scripting_AdvancedDynamics\\2CM4_Assignment1_automated.pde"

max_d=0 # maximum distance
d_angle=0 # corresponding launch angle (distance)
max_KE=0 # maximum Kinetic Energy
KE_angle=0 # corresponding launch angle (Kinetic Energy)
AngleRange = sp.arange(5,51,5) # for angles in increments 5 from 5 to 50
print("Time         ", "      x-disp        ", "    KE        ", "   Angle")
for Angle in AngleRange:
    with open(FlexFileName, "w") as f:
        print(FlexCode%Angle, file=f)
        
    completed =subprocess.run(["FlexPDE6s","-S",FlexFileName],timeout=100, shell=True)#,shell=True)
    
    with open("C:\\Users\\denni\\OneDrive\\Documents\\McMaster\\Winter 2020\\Computational Multiphysics - 2CM4\\M1-Scripting_AdvancedDynamics\\trajectory_data2.txt","r") as f:
        data=sp.loadtxt(f, skiprows=7) # read the output text file
    t=data[:,0] # time
    xd=data[:,1] # x-displacement
    KE=data[:,2] # Kinetic Energy
    yd=data[:,3] # y-displacement
    
    space="     "
    # print the values for the launch angle
    print(t[-1],space,xd[-1],space,KE[-1],space,Angle)
    
    # finding the max x-displacement and corresponding angle
    if xd[-1]>max_d:
        max_d=xd[-1]
        d_angle=Angle
    if abs(KE[-1])>max_KE:
        max_KE=KE[-1]
        KE_angle=Angle
    
    # add data to plot
    plt.plot(xd,yd)

print("")
print("The Maximum xd was ", max_d, "at", d_angle,"degrees, ", "and Maximum KE was ", max_KE, "at", KE_angle)
plt.title('Trajectory for various launch angles')
plt.legend(AngleRange)
plt.show()

