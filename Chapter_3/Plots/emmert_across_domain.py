import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import tables 

###At Start
h5=tables.openFile("./MAST_no_gap_electrons_0.h5" )
elecs=h5.root.electrons.read()
data = elecs[:,1]
temp1 = np.std(data)
#pos, data= np.loadtxt("felectron_final", dtype='d', unpack='true')  #for electrons
VETH = 937768.6678921346

print np.std(data)/VETH
#pos, data, time= np.loadtxt("fion", dtype='d', unpack='true')  #for ions

maxwell = stats.maxwell

data/=VETH
params = maxwell.fit(data, floc=0)
print(params)




# best fit of data
(mu, sigma) = norm.fit(data)

# the histogram of the data
#n, bins, patches = plt.hist(data, 1000, normed=1, facecolor='green', alpha=0.75)
n, bins = np.histogram(data, 1000, normed=1)
# add a 'best fit' line
y1 = mlab.normpdf( bins, mu, sigma)
y1/=np.max(y1)

temperature = temp1*temp1*9.11e-31/1.6e-19
temperature = round(temperature,1)
l = plt.plot(bins, y1, 'r--', linewidth=2, label = ' Initial: $T_e$ = 5.0eV')








##At end For Emmert Start of domain

dx =   0.0001

start= 200 *dx
end_1 = 350*dx

h5=tables.openFile("./MAST_no_gap_electrons_83.h5" )
elecs=h5.root.electrons.read()

a = np.where((elecs[:,0]<end_1) & (elecs[:,0]>start))
data_1 = []

for i in a:
	data_1.append(elecs[i,1])
print "AAAAA" , a

#data = elecs[:,1]
temp3 = np.std(data_1)
#pos, data= np.loadtxt("felectron_final", dtype='d', unpack='true')  #for electrons
VETH = 937099.713300 
print np.std(data_1)/VETH
#pos, data, time= np.loadtxt("fion", dtype='d', unpack='true')  #for ions

maxwell = stats.maxwell
i=0 
while i<len(data_1):
	data_1[i]/=VETH
	i+=1
params = maxwell.fit(data_1, floc=0)
print(params)



# best fit of data
(mu, sigma) = norm.fit(data_1)

# the histogram of the data
#n, bins, patches = plt.hist(data, 1000, normed=1, facecolor='green', alpha=0.75)
n, bins = np.histogram(data_1, 1000, normed=1)
# add a 'best fit' line
y3 = mlab.normpdf( bins, mu, sigma)
y3/=np.max(y3)

temperature = temp3*temp3*9.11e-31/1.6e-19
temperature = round(temperature,1)
l = plt.plot(bins, y3, 'b--', linewidth=2, label = 'Emmert: $T_e$ = %seV' %temperature)











#plot
plt.xlabel('Velocity ($V_{th}$)')
plt.ylabel('Probability')
#plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=%.3f,\ \sigma=%.3f$' %(mu, sigma))
plt.legend(loc=2)
plt.grid(True)

plt.show()
