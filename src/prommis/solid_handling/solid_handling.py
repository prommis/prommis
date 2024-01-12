#This module including power consumption for solid crushing; breakage probability and distribution.

__author__ = "Lingyan Deng"
__version__ = "1.0.0"

import math

class CrushPower():
    """
    This is the class for solid crushing power required.
    """
    def power(self,BWI,F80,P80,Massflow):
        self.solidcrushpower = 10*Massflow*BWI*(1/math.sqrt(P80)-1/math.sqrt(F80))
        """
        BWI: bond work index of rock. Its a measure of ore resistance to crushing and grinding. unit: kWh/tonne
        F80: the 80% passing size of the feed measured in micro meter.
        P80: the 80% passing size of the product measured in micro meter.
        Massflow: solid mass flow rate. Unit: tonne
        ref: https://help.syscad.net/Crusher2
        """

    def get_power(self):
        return self.solidcrushpower

"""
example:
my_instance=CrushPower()
my_instance.power(12,200,80,2)
print(my_instance.get_power())
results should be: 9.862252981520335
"""

class BreakageDistribution():
    """
    This is the class for accumulative fraction of solid breakage probability distribution smaller than size x
    """
    def rosinrammlerdistribution(self,x,x50,n):
        self.soliddistribution = 1 - math.exp(-(x / x50) ** n)
        """
        x: particle size. Unit: micro meter
        x50: the median particle size. unit: micro meter.
        n: a parameter related to the width of the size distribution.
        x50, n are determined experimentally. by collecting data and do linear regression to find n:
        ln(-ln(1-p(x))) = nln(x)-nln(x50). p(x) is the self.soliddistribution.
        draw a linear line: Ln(x) as x-axis, ln(-ln(1-p(x))) as y-axis: 
        slop is n, nln(x50) is the intersection. Validate, reevaluate, get final n and x_50.
        reference: Rosin, P., & Rammler, E. "The laws governing the fineness of powdered coal." Journal of the Institute of Fuel, 7, 29-36 (1933).
        """
    def get_rosinrammlerdistribution(self):
        return self.soliddistribution
    
"""
Example:
my_instance1=BreakageDistribution()
my_instance1.rosinrammlerdistribution(60,80,1.5)
print(my_instance1.get_rosinrammlerdistribution())
"""

