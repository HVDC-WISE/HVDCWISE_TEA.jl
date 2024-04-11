function mpc = test01_ac()

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
      1    2 0 0  0  0    1 1  0     400    1  1.1  0.9;
      2    1 10 0  0  0    1 1  0     400    1  1.1  0.9;  % Pd : dimensioning power of the load
];

%% generator
% bus Pg Qg Qmax Qmin Vg mBase status Pmax Pmin
mpc.gen = [
    1 10  0    10   -10  1   1      1   10    0;
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
mpc.branch_currents = [
                      15.00;
];
