import pandas as pd


class Solution:
    def __init__(self, salinity=None, temperature=None):
        self.__dict__["salinity"] = salinity
        self.__dict__["temperature"] = temperature

    def __setattr__(self, name, value):
        if name in self.__dict__:
            super(Solution, self).__setattr__(name, value)
        else:
            raise AttributeError(f"{self.__class__.__name__} has no attribute {name}")

    def __repr__(self):
        return f"Solution(salinity={self.salinity}, temperature={self.temperature})"


test = Solution()


analyte = pd.Series(
    {
        "alkalinity": 2300,
        "dic": 2100,
        "salinity": 34.1,
        "temperature": 25.0,
    }
)
