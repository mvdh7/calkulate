import calkulate as calk
from time import time
import numpy as np
from numba import jit

vec = np.full((100000,1),5.)
# vec = np.float_(5)

@jit
def test1(reps):
    for i in range(reps):
        x = vec**0.5

@jit
def test2(reps):
    for _ in range(reps):
        x = np.sqrt(vec)

def test3(reps):
    for _ in range(reps):
        x = np.sqrt(vec)

reps = int(1e5)

go = time()

test1(reps)

print(time() - go)

# print(x)

go = time()

test2(reps)

print(time() - go)

# print(x)

go = time()

test3(reps)

print(time() - go)
