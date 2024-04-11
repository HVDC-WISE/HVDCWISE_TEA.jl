function mpc = test07()

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

mpc.time_elapsed = 1.0;

% AC buses
% bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
mpc.bus = [
1 1 100 0 0 0 1 1 0 400 1 1.1 0.9;
2 2 0 0 0 0 1 1 0 400 1 1.1 0.9;
3 1 0 0 0 0 1 1 0 400 1 1.1 0.9;
];

% AC branches
% fbus tbus  r x b rateA rateB rateC ratio angle status angmin angmax
mpc.branch = [
2 1 0.001 0.0001 0 12 12 12 0 0 1 -60 60;
2 3 0.001 0.0001 0 1 1 1 0 0 1 -60 60;
];

% AC branch rating
% c_rating_a
mpc.branch_currents = [
15;
15;
];

% generators
% bus Pg Qg Qmax Qmin Vg mBase status Pmax Pmin
mpc.gen = [
2 0 0 10 -10 1 100 1.0 10 0;
];

% dispatchable generators operating cost model
% 2={startup shutdown n c(n-1) ... C0} startup shutdown n c(1) c(0)
mpc.gencost = [
2 0 0 2 50 0;
];

% storage data
% storage_bus ps qs energy energy_rating charge_rating discharge_rating charge_efficiency discharge_efficiency thermal_rating  qmin qmax r x p_loss q_loss status
mpc.storage = [
3 0 0 0 1.0 1.0 1.0 0.9 0.9 1.0 -1 1 0 0 0.0 0.0 1;
];

% storage additional data
%column_names% max_energy_absorption stationary_energy_inflow stationary_energy_outflow self_discharge_rate
mpc.storage_extra = [
0 0.0 0.0 0.001;
];
end