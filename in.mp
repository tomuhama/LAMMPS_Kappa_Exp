
# 後々変えたり何度も使いたい変数を定義
# lc:格子定数, x,y,z:repeatの個数, nbin:層数，nw:rNEMDのパラメータで何ステップごとに加熱(エネルギーの交換)をするか，t_eq:平均温度, thi,tlo:高温と低温，tdamp:加熱の頻度，kb:ボルツマン定数
variable   lc equal 5.431
variable   x equal 5
variable   y equal 5
variable   z equal 600
variable   nbin equal ${z}/2
variable   nw equal 100
variable   dt equal 0.001
variable   t_eq equal 300
variable   thi equal ${t_eq}+5
variable   tlo equal ${t_eq}-5
variable   tdamp equal ${dt}*100
variable   kb equal 8.617333262e-5 

# metal単位系
units       metal
atom_style  atomic
boundary     p p p

# 構造の定義
# lattice 結晶の種類と格子定数(慣用単位胞)
# region:boxという名前の領域を定義, create_box:何種類の元素でどの領域をシミュレーションセルにするか
# create_atoms:どこに何種類の原子を作るか, mass:何番目の元素が質量どれくらいか
lattice         diamond ${lc}
region          box block 0 ${x} 0 ${y} 0 ${z}
create_box      1 box
create_atoms    1 box
mass            1 28.085

# 初期速度と時間ステップ
velocity        all create ${t_eq} 87287
timestep        ${dt}

# potentialの情報，ここではOptimized SW
pair_style      sw
pair_coeff      * * Opt_Si.sw Si

# NPTで緩和計算 2 ns
fix             npt1 all npt temp ${t_eq} ${t_eq} ${tdamp} iso 0.0 0.0 ${tdamp} fixedpoint 0.0 0.0 0.0
thermo_style    custom step temp pe ke etotal press vol lx ly lz
thermo          100000

run             2000000
unfix  npt1

reset_timestep  0

compute         ke all ke/atom
variable        temp atom c_ke/1.5/${kb}

# NVE+MPで定常状態になるまで緩和計算 5 ns(大体公式Exampleの通り)
fix             1 all nve
compute         layers all chunk/atom bin/1d z lower 0.05 units reduced
fix             2 all ave/chunk 10 100 1000 layers v_temp file profile.mp
fix             3 all thermal/conductivity ${nw} z ${nbin}

variable        tdiff equal f_2[11][3]-f_2[1][3]
thermo_style    custom step temp epair etotal f_3 v_tdiff

thermo          50000
run             5000000

# 熱伝導率計算 5ns

fix             3 all thermal/conductivity ${nw} z ${nbin}

variable        start_time equal time
variable        lengthunit equal 1e-10
variable        timeunit equal 1e-12
variable        eneunit equal 1.602176634e-19
variable        scale equal ${eneunit}/${timeunit}/${lengthunit}
variable kappa equal (f_3/(time-${start_time})/(lx*ly)/2.0)*(lz/2.0)/f_ave*${scale}

fix             ave all ave/time 1 1 1000 v_tdiff ave running
thermo_style    custom step temp epair etotal f_3 v_tdiff f_ave

run             5000000
print           "Running average thermal conductivity: $(v_kappa:%.2f)"

