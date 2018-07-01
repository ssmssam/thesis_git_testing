from math import log, exp
import matplotlib.pyplot as plt
import numpy as np
V = [-10, -13, -14, -15, -16, -18,  -19,-20, -21 ,  -22]
I = [-60.2418359108, -32.5722479805
, -26.7323146854,   -21.9231149303 ,-18.344920776
, -12.6438427765
, -9.81066068867, -8.1978030251,  -6.46745241915 , -5.37530211385]

float_I = -38.3294021579
float_logI = log(float_I*-1)
print "FLOAT" , float_logI
float_V = -12.14576056
I_sat = 39.0530518704 
print "LENGHTTT",len(V), len(I)
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
plt.plot(float_V, float_logI, '.', label = 'Floating')
plt.plot(np.unique(V), np.poly1d(np.polyfit(V, logI, 1))(np.unique(V)))#, label = 'Gradient = %s  Temperature = %s eV '%(grad, temp ))
plt.xlabel('Probe Bias (V)')
plt.ylabel('ln($I^-$)')
plt.legend()
plt.show() 

V = [-10, -13, -14,  -15, -16, -18, -19, -20, -21, -22, -30, -40]
I = [-60.2418359108, -32.5722479805, -26.7323146854, -21.9231149303, -18.344920776, -12.6438427765,-9.81066068867,  -8.1978030251, -6.46745241915, -5.37530211385, -0.988008833644, -0.133514707252
]

totI=[]
for i in I:
	totI.append(i+I_sat)


#From theory I = Isat*(1-exp(Vprobe-Vfloat)/Te)
V_float = -12.75317435
theoryI = []

for i in V:
	theoryI.append(I_sat*(1-exp((i-V_float)/tempe)))

plt.plot(V, totI, 'x')
plt.plot(V,theoryI)
plt.show()
