#ipython27 wing_minimizer6.py data dimofPSF dimofPSFext x1 y1 x2 y2 x3 y3 x4 y4
import numpy as np
from sympy import *
import pyfits, sys
config = sys.argv[1]  #input the file

array=[]
with open("config", "r") as f: 
    tmp3 = f.readlines()        
    for line in tmp3:           
        words=line.split()      
        array.append(words)


data = pyfits.getdata(array[0][1])

database = Matrix(data.shape[0],data.shape[1],lambda i,j:0)
data = database + data #change data type to sympy matrix
pixelsize=float(array[2][1]) ###############################################impo
   


###############################################important input
x1=float(array[6][1])
y1=float(array[7][1])
x2=float(array[8][1])
y2=float(array[9][1])


print x1,y1,x2,y2

dim=int(array[4][1])
dim_ext=int(array[5][1])   ###############################################important input
print "the size of the PSF is",dim
xc1=int(x1/pixelsize)
x1=x1/pixelsize-xc1
yc1=int(y1/pixelsize)
y1=y1/pixelsize-yc1
xc2=int(x2/pixelsize)
x2=x2/pixelsize-xc2
yc2=int(y2/pixelsize)
y2=y2/pixelsize-yc2


variables = Matrix(dim,dim,lambda i,j:Symbol('a_%d'%(j+i*dim+1)))


pos_1=(yc1-1-(variables.shape[0]-1)/2,xc1-1-(variables.shape[1]-1)/2) #the upper left position of psf
pos_2=(yc2-1-(variables.shape[0]-1)/2,xc2-1-(variables.shape[1]-1)/2)



c=float(array[10][1])/float(array[11][1])######################important input


print "The ratio of c is",c




difference=data


sigma=pyfits.getdata(array[1][1])
for i in range(sigma.shape[0]):
    for j in range(sigma.shape[1]):
        sigma[i,j]=sigma[i,j]**2
##if the input file is the wht map
#for wht
#wht=1/sigma**2
#wht = np.zeros(shape=(data.shape[0],data.shape[1]))
#for i in range(data.shape[0]):
#    for j in range(data.shape[1]):
#        wht[i,j] = wht[i,j]+ 1/(0.002038*data[i,j]+0.12988)




#######################################################################create c bilinear interpolation matrix

variables_c = variables*c
def preinterpolate_c(i,j):
    if i==0 and j==0:
        return ZeroMatrix(1,1)
    if i==0 and j==1:
        return ZeroMatrix(1,variables.shape[1])
    if i==0 and j==2:
        return ZeroMatrix(1,1)
    if i==1 and j==0:
        return ZeroMatrix(variables.shape[0],1)
    if i==1 and j==1:
        return variables_c
    if i==1 and j==2:
        return ZeroMatrix(variables.shape[0],1)
    if i==2 and j==0:
        return ZeroMatrix(1,1)
    if i==2 and j==1:
        return ZeroMatrix(1,variables.shape[1])
    if i==2 and j==2:
        return ZeroMatrix(1,1)
#preNew_var_c=BlockMatrix(3,3,preinterpolate_c)
preNew_var_c=BlockMatrix([[preinterpolate_c(0,0),preinterpolate_c(0,1),preinterpolate_c(0,2)],[preinterpolate_c(1,0),preinterpolate_c(1,1),preinterpolate_c(1,2)],[preinterpolate_c(2,0),preinterpolate_c(2,1),preinterpolate_c(2,2)]])
preNew_var_c=preNew_var_c.as_explicit()
New_var_c=ZeroMatrix(variables.shape[0]+1,variables.shape[0]+1).as_mutable()
for i in range(New_var_c.shape[0]):
    for j in range(New_var_c.shape[1]):
        New_var_c[i,j]=preNew_var_c[i,j]*x1*y1+preNew_var_c[i,j+1]*(1-x1)*y1+preNew_var_c[i+1,j]*x1*(1-y1)+preNew_var_c[i+1,j+1]*(1-x1)*(1-y1)

print "Finish the c matrix setup"
#######################################################################create M bilinear interpolation matrix

variables_M = variables*1
def preinterpolate_M(i,j):
    if i==0 and j==0:
        return ZeroMatrix(1,1)
    if i==0 and j==1:
        return ZeroMatrix(1,variables.shape[1])
    if i==0 and j==2:
        return ZeroMatrix(1,1)
    if i==1 and j==0:
        return ZeroMatrix(variables.shape[0],1)
    if i==1 and j==1:
        return variables_M
    if i==1 and j==2:
        return ZeroMatrix(variables.shape[0],1)
    if i==2 and j==0:
        return ZeroMatrix(1,1)
    if i==2 and j==1:
        return ZeroMatrix(1,variables.shape[1])
    if i==2 and j==2:
        return ZeroMatrix(1,1)
#preNew_var_M=BlockMatrix(3,3,preinterpolate_M)
preNew_var_M=BlockMatrix([[preinterpolate_M(0,0),preinterpolate_M(0,1),preinterpolate_M(0,2)],[preinterpolate_M(1,0),preinterpolate_M(1,1),preinterpolate_M(1,2)],[preinterpolate_M(2,0),preinterpolate_M(2,1),preinterpolate_M(2,2)]])
preNew_var_M=preNew_var_M.as_explicit()
New_var_M=ZeroMatrix(variables.shape[0]+1,variables.shape[0]+1).as_mutable()
for i in range(New_var_M.shape[0]):
    for j in range(New_var_M.shape[1]):
        New_var_M[i,j]=preNew_var_M[i,j]*x2*y2+preNew_var_M[i,j+1]*(1-x2)*y2+preNew_var_M[i+1,j]*x2*(1-y2)+preNew_var_M[i+1,j+1]*(1-x2)*(1-y2)

print "Finish the M matrix setup"
#######################################################################create d bilinear interpolation matrix

for i in range(New_var_c.shape[0]):
    for j in range(New_var_c.shape[1]):
        difference[pos_1[0]+i,pos_1[1]+j]=difference[pos_1[0]+i,pos_1[1]+j]-New_var_c[i,j]
        difference[pos_2[0]+i,pos_2[1]+j]=difference[pos_2[0]+i,pos_2[1]+j]-New_var_M[i,j]

total=difference
print "Finish all matrix setup"

Mat = np.zeros(shape=(variables.shape[0]*variables.shape[1],variables.shape[0]*variables.shape[1]))
Mat = np.asmatrix(Mat)



b= np.zeros(shape=(variables.shape[0]*variables.shape[1],1))
b= np.asmatrix(b)
a= ZeroMatrix(variables.shape[0]*variables.shape[1],1).as_mutable()


for i in range(variables.shape[0]*variables.shape[1]):
    for k in range(pos_1[0]+i/dim,pos_1[0]+i/dim+2):
        for l in range(pos_1[1]+i%dim,pos_1[1]+i%dim+2):
            a[i]=a[i]+diff(total[k,l]**2/sigma[k,l],variables[i])

    for k in range(pos_2[0]+i/dim,pos_2[0]+i/dim+2):
        for l in range(pos_2[1]+i%dim,pos_2[1]+i%dim+2):
            a[i]=a[i]+diff(total[k,l]**2/sigma[k,l],variables[i])

                
    for k in range(3):
        for m in range(3):
            if 0 <= i/dim+k-1 < variables.shape[0] and 0 <= i%dim+m-1 < variables.shape[1]:
                var1=i/dim+k-1
                var2=i%dim+m-1
                Mat[i,var2+var1*variables.shape[1]]=Mat[i,var2+var1*variables.shape[1]]+diff(a[i],variables[var1,var2])
                #M[i,y+p*n for vairiables]= M[i,y+p*n for variables]+diff(a[i],variables[p,y])
                a[i]=a[i]-diff(a[i],variables[var1,var2])*variables[var1,var2]
        
            if 0 <= pos_1[0]+i/dim-pos_2[0]+k-1 < variables.shape[0] and 0 <= pos_1[1]+i%dim-pos_2[1]+m-1 < variables.shape[1]:
                var1=pos_1[0]+i/dim-pos_2[0]+k-1
                var2=pos_1[1]+i%dim-pos_2[1]+m-1
                Mat[i,var2+var1*variables.shape[1]]=Mat[i,var2+var1*variables.shape[1]]+diff(a[i],variables[var1,var2])
                a[i]=a[i]-diff(a[i],variables[var1,var2])*variables[var1,var2]

            if 0 <= pos_2[0]+i/dim-pos_1[0]+k-1 < variables.shape[0] and 0 <= pos_2[1]+i%dim-pos_1[1]+m-1 < variables.shape[1]:
                var1=pos_2[0]+i/dim-pos_1[0]+k-1
                var2=pos_2[1]+i%dim-pos_1[1]+m-1
                Mat[i,var2+var1*variables.shape[1]]=Mat[i,var2+var1*variables.shape[1]]+diff(a[i],variables[var1,var2])
                a[i]=a[i]-diff(a[i],variables[var1,var2])*variables[var1,var2]
        
    b[i]=-a[i]
    print i/float(variables.shape[0]*variables.shape[1])*100,"%"


############################################################################################## Normailized character
############################################################################################# creating the regularization Matrix


Mat_s = np.zeros(shape=(variables.shape[0]*variables.shape[1],variables.shape[0]*variables.shape[1]))
Mat_s = np.asmatrix(Mat_s)


a_s= ZeroMatrix(variables.shape[0]*variables.shape[1],1).as_mutable()
h= np.zeros(shape=(variables.shape[0]*variables.shape[1],1))
h= np.asmatrix(h)

variables_s = variables*1
def extend_s(i,j):
    if i==0 and j==0:
        return ZeroMatrix(2,2)
    if i==0 and j==1:
        return ZeroMatrix(2,variables.shape[1])
    if i==0 and j==2:
        return ZeroMatrix(2,2)
    if i==1 and j==0:
        return ZeroMatrix(variables.shape[0],2)
    if i==1 and j==1:
        return variables_s
    if i==1 and j==2:
        return ZeroMatrix(variables.shape[0],2)
    if i==2 and j==0:
        return ZeroMatrix(2,2)
    if i==2 and j==1:
        return ZeroMatrix(2,variables.shape[1])
    if i==2 and j==2:
        return ZeroMatrix(2,2)
#New_var_s=BlockMatrix(3,3,extend_s)
New_var_s=BlockMatrix([[extend_s(0,0),extend_s(0,1),extend_s(0,2)],[extend_s(1,0),extend_s(1,1),extend_s(1,2)],[extend_s(2,0),extend_s(2,1),extend_s(2,2)]])
New_var_s=New_var_s.as_explicit()

for i in range(2,New_var_s.shape[0]-2):
    for j in range(2,New_var_s.shape[1]-2):
        a_s[j-2+(i-2)*dim]=a_s[j-2+(i-2)*dim]+0.5*((New_var_s[i,j-2]-2*New_var_s[i,j-1]+New_var_s[i,j])**2)
        a_s[j-2+(i-2)*dim]=a_s[j-2+(i-2)*dim]+0.5*((New_var_s[i,j-1]-2*New_var_s[i,j]+New_var_s[i,j+1])**2)
        a_s[j-2+(i-2)*dim]=a_s[j-2+(i-2)*dim]+0.5*((New_var_s[i,j]-2*New_var_s[i,j+1]+New_var_s[i,j+2])**2)
        a_s[j-2+(i-2)*dim]=a_s[j-2+(i-2)*dim]+0.5*((New_var_s[i-2,j]-2*New_var_s[i-1,j]+New_var_s[i,j])**2)
        a_s[j-2+(i-2)*dim]=a_s[j-2+(i-2)*dim]+0.5*((New_var_s[i-1,j]-2*New_var_s[i,j]+New_var_s[i+1,j])**2)
        a_s[j-2+(i-2)*dim]=a_s[j-2+(i-2)*dim]+0.5*((New_var_s[i,j]-2*New_var_s[i+1,j]+New_var_s[i+2,j])**2)
        a_s[j-2+(i-2)*dim]=diff(a_s[j-2+(i-2)*dim],New_var_s[i,j])
        if New_var_s[i,j] != 0:
            Mat_s[j-2+(i-2)*dim,j-2+(i-2)*dim]=Mat_s[j-2+(i-2)*dim,j-2+(i-2)*dim]+diff(a_s[j-2+(i-2)*dim],New_var_s[i,j])
            a_s[j-2+(i-2)*dim]=a_s[j-2+(i-2)*dim]-diff(a_s[j-2+(i-2)*dim],New_var_s[i,j])*New_var_s[i,j]
                
        if New_var_s[i-1,j] != 0:
            Mat_s[j-2+(i-2)*dim,j-2+(i-1-2)*dim]=Mat_s[j-2+(i-2)*dim,j-2+(i-1-2)*dim]+diff(a_s[j-2+(i-2)*dim],New_var_s[i-1,j])
            a_s[j-2+(i-2)*dim]=a_s[j-2+(i-2)*dim]-diff(a_s[j-2+(i-2)*dim],New_var_s[i-1,j])*New_var_s[i-1,j]
        
        if New_var_s[i-2,j] != 0:
            Mat_s[j-2+(i-2)*dim,j-2+(i-2-2)*dim]=Mat_s[j-2+(i-2)*dim,j-2+(i-2-2)*dim]+diff(a_s[j-2+(i-2)*dim],New_var_s[i-2,j])
            a_s[j-2+(i-2)*dim]=a_s[j-2+(i-2)*dim]-diff(a_s[j-2+(i-2)*dim],New_var_s[i-2,j])*New_var_s[i-2,j]
        
        if New_var_s[i+1,j] != 0:
            Mat_s[j-2+(i-2)*dim,j-2+(i+1-2)*dim]=Mat_s[j-2+(i-2)*dim,j-2+(i+1-2)*dim]+diff(a_s[j-2+(i-2)*dim],New_var_s[i+1,j])
            a_s[j-2+(i-2)*dim]=a_s[j-2+(i-2)*dim]-diff(a_s[j-2+(i-2)*dim],New_var_s[i+1,j])*New_var_s[i+1,j]
        
        if New_var_s[i+2,j] != 0:
            Mat_s[j-2+(i-2)*dim,j-2+(i+2-2)*dim]=Mat_s[j-2+(i-2)*dim,j-2+(i+2-2)*dim]+diff(a_s[j-2+(i-2)*dim],New_var_s[i+2,j])
            a_s[j-2+(i-2)*dim]=a_s[j-2+(i-2)*dim]-diff(a_s[j-2+(i-2)*dim],New_var_s[i+2,j])*New_var_s[i+2,j]
        
        if New_var_s[i,j-2] != 0:
            Mat_s[j-2+(i-2)*dim,(j-2-2)+(i-2)*dim]=Mat_s[j-2+(i-2)*dim,(j-2-2)+(i-2)*dim]+diff(a_s[j-2+(i-2)*dim],New_var_s[i,j-2])
            a_s[j-2+(i-2)*dim]=a_s[j-2+(i-2)*dim]-diff(a_s[j-2+(i-2)*dim],New_var_s[i,j-2])*New_var_s[i,j-2]
        
        if New_var_s[i,j-1] != 0:
            Mat_s[j-2+(i-2)*dim,(j-1-2)+(i-2)*dim]=Mat_s[j-2+(i-2)*dim,(j-1-2)+(i-2)*dim]+diff(a_s[j-2+(i-2)*dim],New_var_s[i,j-1])
            a_s[j-2+(i-2)*dim]=a_s[j-2+(i-2)*dim]-diff(a_s[j-2+(i-2)*dim],New_var_s[i,j-1])*New_var_s[i,j-1]

        if New_var_s[i,j+1] != 0:
            Mat_s[j-2+(i-2)*dim,(j+1-2)+(i-2)*dim]=Mat_s[j-2+(i-2)*dim,(j+1-2)+(i-2)*dim]+diff(a_s[j-2+(i-2)*dim],New_var_s[i,j+1])
            a_s[j-2+(i-2)*dim]=a_s[j-2+(i-2)*dim]-diff(a_s[j-2+(i-2)*dim],New_var_s[i,j+1])*New_var_s[i,j+1]

        if New_var_s[i,j+2] != 0:
            Mat_s[j-2+(i-2)*dim,(j+2-2)+(i-2)*dim]=Mat_s[j-2+(i-2)*dim,(j+2-2)+(i-2)*dim]+diff(a_s[j-2+(i-2)*dim],New_var_s[i,j+2])
            a_s[j-2+(i-2)*dim]=a_s[j-2+(i-2)*dim]-diff(a_s[j-2+(i-2)*dim],New_var_s[i,j+2])*New_var_s[i,j+2]
        h[j-2+(i-2)*dim]=a_s[j-2+(i-2)*dim]
        a_s[j-2+(i-2)*dim]=None
    print (i-2)/float(dim)*100,"%"


############################################################################################# creating the regularization Matrix
################################################################################################################# regularization

from scipy.optimize import minimize_scalar
import scipy
from scipy import linalg

def optimizelambda(w):
    z=float(10**w)
    #factor=0.002
    print "The regularized factor is",z
    Mat_total=Mat+z*Mat_s
    cholesky = np.linalg.cholesky(Mat_total)
    Mat_total_inv=linalg.inv(cholesky).T*cholesky.I
    u=np.zeros(dim**2)
    for i in range(dim**2):
        u[i]=Mat_total_inv[i,i]
    u=u.reshape(variables.shape[0],variables.shape[1])
    #pyfits.writeto('covariance_matrix.fits',u,clobber=True)
    #Mat_total_inv=linalg.inv(Mat_total)
    correction=Mat_total_inv*(b)
    correction=correction.reshape(variables.shape[0],variables.shape[1])
    pyfits.writeto('correction.fits',correction,clobber=True)
    
    correction_Es = pyfits.getdata('correction.fits')
    correction_Es = np.append(correction_Es,np.zeros([2,correction_Es.shape[1]]),0)
    correction_Es = np.append(np.zeros([2,correction_Es.shape[1]]),correction_Es,0)
    correction_Es = np.append(correction_Es,np.zeros([correction_Es.shape[0],2]),1)
    correction_Es = np.append(np.zeros([correction_Es.shape[0],2]),correction_Es,1)
    Es=0
    for i in range(1,correction_Es.shape[1]-1):
        for j in range(1,correction_Es.shape[0]-1):
            Es=Es+(correction_Es[i,j-1]-2*correction_Es[i,j]+correction_Es[i,j+1])**2*0.5
            Es=Es+(correction_Es[i-1,j]-2*correction_Es[i,j]+correction_Es[i+1,j])**2*0.5
    Ns=dim**2
    n=np.trace(Mat_total_inv*Mat_s)
    diff=2*z*Es-Ns+z*n
    diff_2=diff**2
    #print "Ns is",Ns
    #print "trace is",n
    #print "Es is",Es
    print "diff is",diff
    return diff



def optimizelambda_abs(w):
    z=float(10**w)
    #factor=0.002
    print "The regularized factor is",z
    Mat_total=Mat+z*Mat_s
    cholesky = np.linalg.cholesky(Mat_total)
    Mat_total_inv=linalg.inv(cholesky).T*cholesky.I
    u=np.zeros(dim**2)
    for i in range(dim**2):
        u[i]=Mat_total_inv[i,i]
    u=u.reshape(variables.shape[0],variables.shape[1])
    #pyfits.writeto('covariance_matrix.fits',u,clobber=True)
    #Mat_total_inv=linalg.inv(Mat_total)
    correction=Mat_total_inv*(b)
    correction=correction.reshape(variables.shape[0],variables.shape[1])
    pyfits.writeto('correction.fits',correction,clobber=True)
    
    correction_Es = pyfits.getdata('correction.fits')
    correction_Es = np.append(correction_Es,np.zeros([2,correction_Es.shape[1]]),0)
    correction_Es = np.append(np.zeros([2,correction_Es.shape[1]]),correction_Es,0)
    correction_Es = np.append(correction_Es,np.zeros([correction_Es.shape[0],2]),1)
    correction_Es = np.append(np.zeros([correction_Es.shape[0],2]),correction_Es,1)
    Es=0
    for i in range(1,correction_Es.shape[1]-1):
        for j in range(1,correction_Es.shape[0]-1):
            Es=Es+(correction_Es[i,j-1]-2*correction_Es[i,j]+correction_Es[i,j+1])**2*0.5
            Es=Es+(correction_Es[i-1,j]-2*correction_Es[i,j]+correction_Es[i+1,j])**2*0.5
    Ns=dim**2
    n=np.trace(Mat_total_inv*Mat_s)
    diff=2*z*Es-Ns+z*n
    diff_2=diff**2
    #print "Ns is",Ns
    #print "trace is",n
    #print "Es is",Es
    print "diff is",diff
    return diff**2

tmp=-10
while (optimizelambda(tmp)*optimizelambda(tmp+3)>0 and tmp<=10):
    tmp +=3
    print 1
while (optimizelambda(tmp)*optimizelambda(tmp+0.4)>0 and tmp<=10):
    tmp +=0.4
    print 2

res=minimize_scalar(optimizelambda_abs,bounds=(tmp,tmp+0.4),method='bounded')
print res.x


#bnds=(0,None)
#factor=0.002
#xopt= opt.fmin_1_bfgs_b(optimizelambda, factor, method='SLSQP', bounds=bnds, options={'xtol':0.0001,'disp':True})
#print xopt


correction=pyfits.getdata('correction.fits')
old_psf=pyfits.getdata(array[3][1])


if correction.shape[0]>old_psf.shape[0]:
    tmp=(dim-old_psf.shape[0])/2
    old_psf = np.append(old_psf,np.zeros([tmp,old_psf.shape[1]]),0)
    old_psf = np.append(np.zeros([tmp,old_psf.shape[1]]),old_psf,0)
    old_psf = np.append(old_psf,np.zeros([old_psf.shape[0],tmp]),1)
    old_psf = np.append(np.zeros([old_psf.shape[0],tmp]),old_psf,1)
  
if correction.shape[0]<old_psf.shape[0]:
    tmp=(old_psf.shape[0]-correction.shape[0])/2
    correction = np.append(correction,np.zeros([tmp,correction.shape[1]]),0)
    correction = np.append(np.zeros([tmp,correction.shape[1]]),correction,0)
    correction = np.append(correction,np.zeros([correction.shape[0],tmp]),1)
    correction = np.append(np.zeros([correction.shape[0],tmp]),correction,1)
    


solution=correction+old_psf

for i in range(solution.shape[0]):
    for j in range(solution.shape[1]):
        if solution[i,j]<0:
            solution[i,j]=0
            print "solution is not all positive; correct it to 0"
pyfits.writeto('new.fits',solution,clobber=True)
solution_ext=solution[(solution.shape[0]+1)/2-((dim_ext+1)/2):(solution.shape[1]+1)/2+((dim_ext-1)/2),(solution.shape[0]+1)/2-((dim_ext+1)/2):(solution.shape[1]+1)/2+((dim_ext-1)/2)]
for i in range(solution_ext.shape[0]):
    for j in range(solution_ext.shape[1]):
        if solution_ext[i,j]<0:
            solution_ext[i,j]=0
pyfits.writeto('new_ext.fits',solution_ext,clobber=True)
