# Generated with SMOP  0.41-beta
from libsmop import *
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m


@function
def spiralling113MESA15(M=None,m=None,a=None,e=None,__=None,Tend=None,*args,**kwargs):
    varargin = spiralling113MESA15.varargin
    nargin = spiralling113MESA15.nargin

    #Solves the equations of motion for a two-body system in a common envelope
#with polytropic index p, i.e. P ~ rho^(1+1/p).
    G=6.67e-08
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:4
    Rsol=dot(696342,100000.0)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:5
    Msol=dot(1.98855,1e+33)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:6
    # -------------------------------------------------------------------------

    #Dynamical friction law

    global rho
    global C
    global M_e
    # load profile15Msol.mat rho q c_s;
    rho=concat([0,0,0])
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:20
    q=concat([0,0,0])
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:21
    c_s=concat([0,0,0])
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:22
    Rho=copy(rho)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:23
    q=dot(q,Rsol)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:24
    c_s=dot(c_s,100)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:24
    R=max(q)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:25
    q=smooth(q)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:25
    #Smoothing

    rho=lambda r=None: multiply((r <= R),real(interp1(q,Rho,r,'linear','extrap'))) + multiply((r > R),0)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:30
    C=lambda r=None: abs(interp1(q,c_s,r,'linear','extrap'))
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:32
    M_e=lambda x=None: multiply((x <= R),integral(lambda r=None: dot(dot(dot(4,pi),r ** 2.0),rho(r)),0,x)) + multiply((x > R),integral(lambda r=None: dot(dot(dot(4,pi),r ** 2.0),rho(r)),0,R))
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:34
    # -------------------------------------------------------------------------

    #initial conditions
    T2=dot(dot(4,pi ** 2),a ** 3) / (dot(G,(M + m)))
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:40
    T=sqrt(T2)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:41
    L=sqrt(dot(dot(dot(G,(M + m + M_e(a))),a),(1 - e ** 2)))
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:42
    y0=concat([(1 + e),0,0,dot(L,T) / (dot(a ** 2,(1 + e)))])
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:43
    #ode113 stuff

    tspan=concat([0,Tend]) / T
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:47
    opts=odeset('RelTol',1e-10,'Stats','on')
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:48
    #Integrator
    t,y=ode113(odefun,tspan,y0,opts,nargout=2)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:52
    X=dot(a,y(arange(),1))
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:55
    Y=dot(a,y(arange(),3))
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:56
    v_x=dot(a,y(arange(),2)) / T
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:57
    v_y=dot(a,y(arange(),4)) / T
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:58
    t=dot(t,T)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:59
    return X,Y,v_x,v_y,t


@function
def odefun(t=None,y=None,*args,**kwargs):
    varargin = odefun.varargin
    nargin = odefun.nargin

    dydt=zeros(1,4).T
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:63
    #parameters
    G=6.67e-08
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:65
    c=29979245800
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:66
    Msol=dot(1.98855,1e+33)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:67
    Rsol=dot(696342,100000.0)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:68
    R=dot(111.3642,Rsol)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:69
    M=dot(15,Msol)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:70
    # M = 1.4*Msol;
    m=dot(1,Msol)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:72
    mu=dot(M,m) / (M + m)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:73
    a=dot(0.7,R)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:74
    T=sqrt(dot(dot(4,pi ** 2),a ** 3) / (dot(G,(M + m))))
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:75
    global rho
    global C
    global M_e
    r=sqrt(y(1) ** 2 + y(3) ** 2)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:82
    v=sqrt(y(2) ** 2 + y(4) ** 2)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:83
    b_90=max(concat([dot(G,M) / (dot(v,a) / T) ** 2,dot(0.1,Rsol)]))
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:85
    N=R / b_90
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:86
    c_s=C(dot(a,r))
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:88
    Mach=(dot(v,a) / T) / c_s
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:89
    # Note f is cut at 1 so as not to diverge
    h1=max(concat([1 / N,0.01]))
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:91
    h=(1 - dot((2 - h1),exp(- 2 + dot(2,h1))) / (dot(N ** 2,h1))) ** (- 1 / 2) - 1
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:93
    #Ostriker
    if (Mach < 1 - h1):
        I=dot(0.5,log((1 + Mach) / (1 - Mach))) - Mach
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:97
    else:
        if ((Mach >= 1 - h1) and (Mach < 1 + h)):
            I=dot(0.5,log((2 - h1) / h1)) - 1 + h1
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:100
            #         disp('Speed of sound reached');
        else:
            I=dot(0.5,log(1 - 1 / Mach ** 2)) + log(N)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:103
            if (Mach < 1):
                disp('Wrong condition')

    f=dot(dot(dot(dot(dot(I,4),pi),G ** 2),M ** 2),rho(dot(a,r))) / ((dot(a,v) / T) ** 3)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:109
    if (I < 0):
        disp('boo')

    rdot=(dot(y(1),y(2)) + dot(y(3),y(4))) / r
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:116
    rdot=dot(rdot,a) / T
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:117
    nu=mu / (M + m)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:118
    v=dot(a,v) / T
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:120
    r=dot(a,r)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:121
    A=dot(1 / c ** 2,(dot(dot(- 3,rdot ** 2),nu) / 2 + v ** 2 + dot(dot(3,nu),v ** 2) - dot(dot(G,(M + m)),(4 + dot(2,nu))) / r)) + dot(1 / c ** 4,(dot(dot(15,rdot ** 4),nu) / 8 - dot(dot(45,rdot ** 4),nu ** 2) / 8 - dot(dot(dot(9,rdot ** 2),nu),v ** 2) / 2 + dot(dot(dot(6,rdot ** 2),nu ** 2),v ** 2) + dot(dot(3,nu),v ** 4) - dot(dot(4,nu ** 2),v ** 4) + dot(dot(G,(M + m)) / r,(dot(- 2,rdot ** 2) - dot(dot(25,rdot ** 2),nu) - dot(dot(2,rdot ** 2),nu ** 2) - dot(dot(13,nu),v ** 2) / 2 + dot(dot(2,nu ** 2),v ** 2))) + dot((dot(G,(M + m)) / r) ** 2,(9 + dot(87,nu) / 4)))) + dot(1 / c ** 5,(dot(dot(dot(dot(dot(- 24,nu),rdot),v ** 2),G),(M + m)) / (dot(5,r)) - dot((dot(dot(136,rdot),nu) / 15),(dot(G,(M + m)) / r) ** 2)))
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:123
    B=dot(1 / c ** 2,(dot(- 4,rdot) + dot(dot(2,rdot),nu))) + dot(1 / c ** 4,(dot(dot(9,rdot ** 3),nu) / 2 + dot(dot(3,rdot ** 3),nu ** 2) - dot(dot(dot(15,rdot),nu),v ** 2) / 2 - dot(dot(dot(2,rdot),nu ** 2),v ** 2) + dot(dot(G,(M + m)) / r,(dot(2,rdot) + dot(dot(41,rdot),nu) / 2 + dot(dot(4,rdot),nu ** 2))))) + dot(1 / c ** 5,(dot(dot((dot(dot(8,nu),v ** 2) / 5),G),(M + m)) / r + dot(dot(24,nu) / 5,(dot(G,(M + m)) / r) ** 2)))
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:129
    r=r / a
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:133
    dydt[1]=y(2)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:135
    dydt[2]=dot(dot(- 4,pi ** 2),(dot((1 + A),y(1)) / r + dot(dot(B,y(2)),a) / T)) / r ** 2 - dot(dot(T,f),y(2)) / M - dot(dot(dot(4,pi ** 2),((M_e(dot(a,r)) - m) / (M + m))),y(1)) / r ** 3
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:136
    dydt[3]=y(4)
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:137
    dydt[4]=dot(dot(- 4,pi ** 2),(dot((1 + A),y(3)) / r + dot(dot(B,y(4)),a) / T)) / r ** 2 - dot(dot(T,f),y(4)) / M - dot(dot(dot(4,pi ** 2),((M_e(dot(a,r)) - m) / (M + m))),y(3)) / r ** 3
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:138
    Rt=dot(dot(0.1,Rsol),max(concat([(M / m) ** (1 / 3),(m / M) ** (1 / 3)])))
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:142
    if (dot(r,a) <= Rt):
        dydt[1]=0
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:146
        dydt[2]=0
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:147
        dydt[3]=0
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:148
        dydt[4]=0
# ../gassy/studies/matlab_scripts/spiralling113MESA15.m:149
        disp(dot(t,T))

    return dydt

if __name__ == '__main__':
    pass
