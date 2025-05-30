# %%
import calkulate as calk
from calkulate.classes import Titration


ds = calk.read_csv("tests/data/titration_table.csv")
tt = Titration(ds.loc[0])
