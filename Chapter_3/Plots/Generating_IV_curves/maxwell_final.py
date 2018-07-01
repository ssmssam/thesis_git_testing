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
l = plt.plot(bins, y1, 'r--', linewidth=2, label = ' Initial: $T_e$ = %seV' %temperature)





##At end For Maxwell

h5=tables.openFile("./MAST_no_gap_electrons_84.h5" )
elecs=h5.root.electrons.read()
data = elecs[:,1]
temp2=np.std(data)
#pos, data= np.loadtxt("felectron_final", dtype='d', unpack='true')  #for electrons
VETH = 937099.713300 
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
y2 = mlab.normpdf( bins, mu, sigma)
y2/=np.max(y2)
temperature = temp2*temp2*9.11e-31/1.6e-19
temperature = round(temperature,1)
l = plt.plot(bins, y2, 'g--', linewidth=2, label = 'Maxwell: $T_e$ = %seV' %temperature)



temp2 = np.std(data)/VETH


##At end For Emmert

h5=tables.openFile("./MAST_no_gap_electrons_85.h5" )
elecs=h5.root.electrons.read()
data = elecs[:,1]
temp3 = np.std(data)
#pos, data= np.loadtxt("felectron_final", dtype='d', unpack='true')  #for electrons
VETH = 937099.713300 
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
