import matplotlib.pyplot as plt
import numpy as np
import sys

# Load results
results = np.loadtxt(open(sys.argv[1]))

#Load in Global Average and Temperature
glob = results[:, 0]
temp = results[:, 1]

# Print Results

plt.plot(temp, abs(glob))

plt.title("Change in Absolute Magnetization with Temperature")
plt.xlabel("Temperature(KbT)")
plt.ylabel("Absolute Magnetization")


plt.show()


