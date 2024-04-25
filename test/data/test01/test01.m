function mpc = test01()

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

mpc.time_elapsed = 1.0;

% AC buses
%column_names% bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
mpc.bus = [
1 2 0 0 0 0 1 1 0 400 1 1.1 0.9;
2 1 10 0 0 0 1 1 0 400 1 1.1 0.9;
];

% AC branches
%column_names% fbus tbus  r x b rateA rateB rateC ratio angle status angmin angmax
mpc.branch = [
1 2 0.001 0.0001 0 20 20 20 0 0 1 -60 60;
];

% AC branch rating
%column_names% c_rating_a
mpc.branch_currents = [
15;
];

% generators
%column_names% bus Pg Qg Qmax Qmin Vg mBase status Pmax Pmin
mpc.gen = [
1 10 0 10 -10 1 1 1 10 0;
];

% dispatchable generators operating cost model
%column_names% cost_model startup shutdown n c(1) c(0)
mpc.gencost = [
2 0 0 2 50 0;
];
end
