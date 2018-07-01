from math import log
import matplotlib.pyplot as plt
import numpy as np
V = [ -10, -13, -14,  -15, -16, -18, -19, -20, -22, -30, -40]
I = [ -9.05763773978, -4.59557622352, -3.55149121283
, -2.86789591171, -2.24838767008,  -1.43661825,-1.16157795307
,  -0.950624715614
, -0.614167653346, -0.104141471654, -0.0160217648699
]

I_sat = 39.0530518704 

logI=[]

for i in I:
	logI.append(log(i*-1))


print logI

grad , b = np.poly1d(np.polyfit(V, logI, 1))

grad = round(grad,4)
temp = 1/grad
tempe = temp
temp = round(temp,3)
#LABELS PLOTS, ADD LABEL TE =  SAVE AS PDF , EMAIL THIS FILE PLUS PDF TO ADD TO THESIS, DO SAME FOR MAXWELLIAN.  Plot IV curve too. Run a floating wall , show it sits perfectly where it should, proof floating works. Equatte to STANGEBY?
plt.plot(V, logI,'x' )
plt.plot(np.unique(V), np.poly1d(np.polyfit(V, logI, 1))(np.unique(V)), label = 'Gradient = %s  Temperature = %s eV '%(grad, temp ))
plt.xlabel('Probe Bias (V)')
plt.ylabel('ln($I_{electron}$)')
plt.legend()
plt.show() 

V = [-10, -13, -14,  -15, -16, -18, -19, -20, -21, -22, -30, -40]
I = [-60.2418359108, -32.5722479805, -26.7323146854, -21.9231149303, -18.344920776, -12.6438427765,-9.81066068867,  -8.1978030251, -6.46745241915, -5.37530211385, -0.988008833644, -0.133514707252
]

totI=[]
for i in I:
	totI.append(i+I_sat)


#From theory I = Isat*(1-exp(Vprobe-Vfloat)/Te)

theoryI = []


plt.plot(V, totI, 'x')
plt.show()
