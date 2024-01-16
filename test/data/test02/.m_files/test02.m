function mpc = test02_ac()

%%Files are in Flexplan Github
%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

mpc.time_elapsed = 1.0;

%% AC bus
% type 1 = PQ, type 2 = PV, type 3 = reference, type 4 = isolated
% https://matpower.org/docs/ref/matpower5.0/caseformat.html
% bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
mpc.bus = [
      1    2 0   0  0  0    1 1   0    400    1  1.1  0.9;
      2    1 10  0  0  0    1 1   0    400    1  1.1  0.9;
];

%% generator
% bus Pg Qg Qmax Qmin Vg mBase status Pmax Pmin
mpc.gen = [
        1 10  0    10   -10  1   1      1   10    0;
];

%% non-dispatchable generators
%column_names% gen_bus   pref   qmax   qmin gen_status cost_gen cost_curt
mpc.ndgen = [
    1 0  0  -0   0    0    0;
];

%% generator cost
% 1 startup shutdown n x1 y1 ... xn yn
% 2 startup shutdown n c(n-1) ... c0
mpc.gencost = [
      2       0        0 2 50   0;
];

%% AC branch
% fbus tbus    r    x    b rateA rateB rateC ratio angle status angmin angmax
mpc.branch = [
     1    2 0.001 0.0001 0   20   20   20     0     0      1    -60     60;
];

%column_names% c_rating_a
%mpc.branch_currents = [
%                      1500;
%];

%% DC bus
%column_names% busdc_i grid Pdc Vdc basekVdc Vdcmax Vdcmin Cdc
%mpc.busdc = [];

%% converter
%column_names% busdc_i busac_i type_dc type_ac P_g Q_g islcc Vtar  rtf xtf transformer tm   bf filter   rc   xc reactor basekVac Vmmax Vmmin Imax status LossA LossB LossCrec LossCinv  droop   Pdcset Vdcset dVdcset Pacmax Pacmin Qacmax Qacmin conv_confi connect_at ground_type ground_z
%mpc.convdc = [];

%% DC branch
%column_names% fbusdc tbusdc     r l c rateA rateB rateC status line_confi return_type return_z connect_at
%mpc.branchdc = [];

%% load additional data
%column_names% load_id pf_angle pshift_up_rel_max pshift_down_rel_max tshift_up tshift_down eshift_rel_max pred_rel_max ered_rel_max cost_shift cost_red cost_curt  cost_inv flex co2_cost lifetime
%                    #        #          pu(load)            pu(load)         h           h            MWh          MWh          MWh      €/MWh    €/MWh     €/MWh     €/MWh  0/1    €/MWh     year
mpc.load_extra = [
                     1      0.0               0.2                 0.2         1           1            3.0          2.0          3.0        0.0      0.1      10.0       0.0    1     0.0        10;
];