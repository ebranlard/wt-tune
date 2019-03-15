all:
	python Main_WT_Tune_Parametric.py
# 	#ipython Main_AirfoilTuning_GA_TorqueFlap.py
#	ipython PowerCurve.py
# 	$u,ccipython Main_AirfoilTuning_GA.py
m:
	ipython Main_MIT.py

i:
	ipython Main_AirfoilTuning_GA_TorqueFlap.py
d:
	ipython DEBUG.py

p:
	_FastView.bat ../hubData.csv
test:
	python -m unittest discover

# 	_FastView.bat ../OpenFAST_V27_Power_PowerCurve/Turbine_ws10.0.out
# 	_FastView.bat ../wtgData.csv
