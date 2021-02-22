import multiprocessing as mp, time
import numpy as np

# print(mp.cpu_count())


def myfunc(x):
    return x ** 2


xs = list(np.random.normal(size=1000000))


#%% Normal
def loop_normal(xs):
    return list(map(myfunc, xs))

#
def loop_pool(xs):
    return Parallel(n_jobs=8)(delayed(myfunc)(x) for x in xs)

#%%
from joblib import Parallel, delayed
import numpy as np

test = Parallel(n_jobs=8, prefer='threads')(delayed(np.sqrt)(i ** 2) for i in range(10))
