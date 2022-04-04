# Here we do some encapsulation. Also, this meets the requirement of standard input.
# Associating library is TideCal_Lib.py.
from TideCal_Lib import *

PowerSys=PowerSystem()
path1='input-HW-branch.txt'#支路数据
path2='input-HW-generator.txt'#发电机数据
path3='input-HW-tide.txt'#系统潮流数据
PowerSys.init_fromsile(path1,path2,path3)
PowerSys.layout()
PowerSys.initialize()
PowerSys.iterate()
PowerSys.PowerCal()
PowerSys.ultraV()
PowerSys.result()
PowerSys.layout2()

