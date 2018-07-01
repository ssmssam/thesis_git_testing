import matplotlib.pyplot as plt
import numpy as np
from math import exp, log,sqrt, sin, tan
from scipy.optimize import curve_fit


v_float = -27.21946874
TE = 6


def func_a(x, a_fit, c):
	return x*a_fit + c 



probe_v = [-150, -70, -65, -60, -55, -61, -62, -63, -64, -66, -67, -68 , -69, -75, -80, -85, -90, -95, -100, -105, -110, -115, -120, -125 , -130, -135, -140, -145 ] #, -50 ,-45, -40, -35] 
#i_current = [32.6823464766, 29.89952737, 29.71293462, 29.44990157, 29.12218826, 29., , 29.6772624795,29.5925901414 , 29.7027425812 , 29.7203826516 , 29.7952549505 , 29.8830633011 , 29.9661676329  ]                                 #, 28.9304994982, 28.6133702321, 28.0720160709 ,27.737638736, ] 

i_current = [32.8260542503 , 29.89952737 , 29.71293462, 29.4499015718 , 29.0815769013, 29.4846329104 , 29.5094074094,  29.5736172657 , 29.6009789749 , 29.7212450551 , 29.7490771662 , 29.7977637606 , 29.8894137264, 30.1539363825 , 30.3787100798, 30.6041109797, 30.8000333619,  31.0341366965, 31.3478163488 , 31.4399367165, 31.6610255992 , 31.7185714289, 31.895677736, 32.0844656897 , 32.2810152743, 32.4465183351 , 32.6064549736 , 32.7826204769]



e_current = [0, -0.0468049868579, -0.0976083897, -0.222108086714 , -0.48169152304, -0.193335171845  ,-0.166679065427  ,  -0.132888530526  , -0.118933274813  , -0.0921203677689 
, -0.0760483036049  , -0.0666402660456 , -0.0488433949957,  -0.0214816857606 ,  -0.0114464456972  , -0.00509602034466 , -0.0029792118938  , -0.0014896059469 , -0.000470401877969 , -0.000313601251979 , -0.000392001564974, 0, 0 ,0 ,0 , 0, 0, 0]
 






probe_v = np.asarray(probe_v)
i_current = np.asarray(i_current)
e_current = np.asarray(e_current)








# order the arrays
inds = probe_v.argsort()

probe_v = probe_v[inds]
i_current = i_current[inds]
e_current = e_current[inds]

probe_v = np.asarray(probe_v)
i_current = np.asarray(i_current)


total_current = i_current +e_current

V = []
for i in probe_v:
	V.append((abs(v_float-i)/TE)**0.75)

V = np.asarray(V)

delete = [0, 10, 9, 8,27]

V = np.delete(V,delete)
i_current = np.delete(i_current,delete)



popt, pcov = curve_fit(func_a, V, i_current)

a = popt[0]/popt[1]
I0 = popt[1]

print "a", a , "I0", I0
a_fitted = func_a(V, *popt)


a = round(a,4)

I0 = round(I0,3)


plt.plot(V, i_current,'x')

plt.plot(V,a_fitted, '--', label = 'a = %s and I0 = %s'%(a, I0))

plt.legend()
plt.xlabel("V")
plt.ylabel("Ion Current (A)")
plt.show()






cot_theta = [4.4862463893, 3.4874144425 ,2.7474774194]
y = [2.6275912689 , 2.2048886947 , 1.8389081041] 







def func_a(x, c1,  c2):
	return x*c2 + c1


cot_theta = np.asarray(cot_theta)
y = np.asarray(y)

popt, pcov = curve_fit(func_a, cot_theta, y, (0.5,  0.6 ))

c2 = popt[0]
c1 = popt[1]

c1 = round(c1,1)

c2 = round(c2,2)

print "c1", c1 , "c2", c2
a_fitted = func_a(cot_theta, *popt)




plt.rc('text', usetex=True)
plt.rc('font', family='serif')


plt.plot(cot_theta,y, 'x')
plt.xlabel(r'{\cot \theta}')
plt.ylabel(r'{$a L {(\sin \theta )}^{1/2}$}')
plt.plot(cot_theta,a_fitted, linestyle = '--', label = r'$c_1$ = %g , $c_2$ = %g'%(c1,c2))
plt.legend()
plt.show()



	
	
	


