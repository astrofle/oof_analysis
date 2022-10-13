Module to mine AutoOOF results and M&C sampler data [1].
This is part of the High Frequency Task Force efforts to characterize the surface error during different observing conditions, as well as improve the Gravity-Zernike model.

Sub-modules:
   * mineOOFResults.py : Searches for projects that used AutoOOF.
   * SamplerData.py : Talks with the M&C sampler data.
   * ScanLog.py : Extracts times from the scan logs.

There's an exmaple in notebooks/mine_OOF_n_temperature.ipynb
and a python environment psalas_venv with the required dependencies.
To use do (in bash):
   source psalas_venv/bin/activate


[1]: Paul knows the proper name of the weather stations and others from the M&C perspective
