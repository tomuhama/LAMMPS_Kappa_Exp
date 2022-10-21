# LAMMPS_Kappa_Exp
LAMMPSの熱伝導率計算の例

- in.mp
  - Muller-Plathe法でNEMD(rNEMD)
- in.heatflux
  - Green-Kuboで計算するためのEMD
- cal_thermal_cond.py
  - 熱流束のlog(thermo)のデータから熱伝導率計算
- Si_NEMD_Fitting.png
  - ある長さにおけるNEMDの温度勾配Fitting
- Si_NEMD_Fitting_length.png
  - 複数の温度で計算したNEMDの結果をFittingしてバルクの熱伝導率を計算(各長さについて3回計算)
- Si_EMD.png
  - 複数のEMD計算により求められた $\kappa$ の相関時間依存性(30個)
