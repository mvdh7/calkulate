# %%
import calkulate as calk
from calkulate.classes import Titration


ds = calk.read_csv("tests/data/titration_table.csv")
tt = Titration(ds.loc[0])
print(tt)
tt.plot_emf()
tt.plot_pH()
tt.plot_alkalinity()
tt.plot_components()
tt.plot_alkalinity_gran()
tt.plot_emf0_gran()
