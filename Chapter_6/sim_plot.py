import matplotlib.pyplot as plt

depth_0 = [0, -1 , -2 ,-4]

potential_0 = [6.61892967, 6.82605274, 5.24332681, 3.6389412]



potential_250 = [2.79392426 , 22.1705316 , 24.04097461 ,23.80870974]


plt.plot(depth_0, potential_0, color = 'blue')
plt.plot(depth_0, potential_0, 'o', color = 'blue', label = "0 mT")
plt.plot(depth_0, potential_250, color = 'red')
plt.plot(depth_0, potential_250, 'o', color = 'red' , label = "250 mT")

plt.ylabel("Floating Potential (V)")
plt.xlabel("Position (mm)")
plt.legend()
plt.show() 


