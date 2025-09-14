%data for MaFan 3 Area ('Implementing Consensus Based Distributed Control in Power System Toolbox') case
disp('MaFan 3 Area data')

% bus data format
% bus:
% col1 number
% col2 voltage magnitude(pu)
% col3 voltage angle(degree)
% col4 p_gen(pu)
% col5 q_gen(pu),
% col6 p_load(pu)
% col7 q_load(pu)
% col8 G shunt(pu)
% col9 B shunt(pu)
% col10 bus_type
%       bus_type - 1, swing bus
%               - 2, generator bus (PV bus)
%               - 3, load bus (PQ bus)
% col11 q_gen_max(pu)
% col12 q_gen_min(pu)
% col13 v_rated (kV)
% col14 v_max  pu
% col15 v_min  pu
bus = [
    1  1.040 0.0  0.716  0.270  1.25  0.5  0.0  0.0  1  0.0  0.0  0.0  0.0  0.0;
    2  1.025 9.3  1.630  0.067  0.90  0.3  0.0  0.0  2  0.0  0.0  0.0  0.0  0.0;
    3  1.025 4.7  0.850 -0.109  1.00  0.35 0.0  0.0  2  0.0  0.0  0.0  0.0  0.0];

% line data format
% line: from bus, to bus, resistance(pu), reactance(pu),
%       line charging(pu), tap ratio, tap phase, tapmax, tapmin, tapsize
line = [
    1  2  0.0402  0.3661  0.5445 0.0000 0   0    0    0;
    2  3  0.0204  0.2939  0.358  0.0000 0   0    0    0;
    3  1  0.056   0.3782  0.516  0.0000 0   0    0    0];

% Machine data format
% Machine data format
%       1. machine number,
%       2. bus number,
%       3. base mva,
%       4. leakage reactance x_l(pu),
%       5. resistance r_a(pu),
%       6. d-axis sychronous reactance x_d(pu),
%       7. d-axis transient reactance x'_d(pu),
%       8. d-axis subtransient reactance x"_d(pu),
%       9. d-axis open-circuit time constant T'_do(sec),
%      10. d-axis open-circuit subtransient time constant
%                T"_do(sec),
%      11. q-axis sychronous reactance x_q(pu),
%      12. q-axis transient reactance x'_q(pu),
%      13. q-axis subtransient reactance x"_q(pu),
%      14. q-axis open-circuit time constant T'_qo(sec),
%      15. q-axis open circuit subtransient time constant
%                T"_qo(sec),
%      16. inertia constant H(sec),
%      17. damping coefficient d_o(pu),
%      18. dampling coefficient d_1(pu),
%      19. bus number
mac_con = [
    1 1 100 1/10*0.200  1/10*0.0025  1/10*1.8  1/10*0.30  1/10*0.25 8.00  0.03...
    1/10*1.7  1/10*0.55  1/10*0.25 0.4   0.05...
    6.5  0  0  1;
    2 2 100 1/10*0.200  1/10*0.0025  1/10*1.8  1/10*0.30  1/10*0.25 8.00  0.03...
    1.7  1/10*0.55  1/10*0.25 0.4   0.05...
    6.5  0  0  2;
    3 3 100 1/10*0.200  1/10*0.0025  1/10*1.8  1/10*0.30  1/10*0.25 8.00  0.03...
    1/10*1.7  1/10*0.55  1/10*0.25 0.4   0.05...
    6.5  0  0  3];


% governor model
% tg_con matrix format
%column	       data			unit
%  1	turbine model number (=1)
%  2	machine number
%  3	speed set point   wf		pu
%  4	steady state gain 1/R		pu
%  5	maximum power order  Tmax	pu on generator base
%  6	servo time constant   Ts	sec
%  7	governor time constant  Tc	sec
%  8	transient gain time constant T3	sec
%  9	HP section time constant   T4	sec
% 10	reheater time constant    T5	sec
tg_con = [
    1  1  1  25.1  3.0  0.45  0.1  0.0  1.25  5.0;
    1  2  1  25.5  3.0  0.5   0.5  0.0  1.25  5.0;
    1  3  1  25.4  3.0  0.24  0.18 0.0  1.25  5.0];

%active load modulation enabled
% col1 load number
% col2 bus number
% col3 MVA rating
% col4 maximum conductance pu
% col5 minimum conductance pu
% col6 regulator gain K
% col7  time constant T_R

lmod_con = [
    1  1 100  1  -1  1  0.01;
    2  2 100  1  -1  1  0.01;
    3  3 100  1  -1  1  0.01];

% non-conforming load
% col 1           bus number
% col 2           fraction const active power load
% col 3           fraction const reactive power load
% col 4           fraction const active current load
% col 5           fraction const reactive current load

load_con = [...
    1   0  0   .5  0.1;
    2   0  0   .5  0.1;
    3   0  0   .5  0.1];


%Switching file defines the simulation control
% row 1 col1  simulation start time (s) (cols 2 to 6 zeros)
%       col7  initial time step (s)
% row 2 col1  fault application time (s)
%       col2  bus number at which fault is applied
%       col3  bus number defining far end of faulted line
%       col4  zero sequence impedance in pu on system base
%       col5  negative sequence impedance in pu on system base
%       col6  type of fault  - 0 three phase
%                            - 1 line to ground
%                            - 2 line-to-line to ground
%                            - 3 line-to-line
%                            - 4 loss of line with no fault
%                            - 5 loss of load at bus
%       col7  time step for fault period (s)
% row 3 col1  near end fault clearing time (s) (cols 2 to 6 zeros)
%       col7  time step for second part of fault (s)
% row 4 col1  far end fault clearing time (s) (cols 2 to 6 zeros)
%       col7  time step for fault cleared simulation (s)
% row 5 col1  time to change step length (s)
%       col7  time step (s)
%
%
%
% row n col1 finishing time (s)  (n indicates that intermediate rows may be inserted)
sw_con = [...
    0    0    0    0    0    0    0.01; %sets intitial time step
    0.1  3    2    0    0    6    0.01; %apply three phase fault at bus 3 line 3-4
    0.15 0    0    0    0    0    0.01; %clear fault at bus 3
    0.18 0    0    0    0    0    0.01; %clear remote end
    1.0  0    0    0    0    0    0.01; % increase time step
    25  0    0    0    0    0    0]; % end simulation


