import numpy as np


# read in binary file
numbers = np.fromfile("numbers.dat", dtype="longdouble")[1:]


result1 = np.longdouble(0)
for item in np.nditer(numbers[np.isfinite(numbers)], order="K"):
    result1 += item


print(result1)
