function f = mtg_sig(t,k)
addpath(genpath('/Users/kiguli/scots/mfiles'))
addpath(genpath('/Users/kiguli/Documents/My_papers/39-bus\ verification/NE39bus/Toy_System'))
% Syntax: f = mtg_sig(t,k)
% 12:37 PM 7/0/98
% defines modulation signal for turbine power reference
global tg_sig
%% no control
%global n_tg
%tg_sig(:,k) = zeros(n_tg,1);

%% SCOTS Controller
% global mac_spd
%
% BDD = SymbolicSet('toy_controller.bdd');
%
% if rem(t,2) == 0
%     for i = 1:3
%         if mac_spd(i,k)-1 < 0.002 && mac_spd(i,k)-1 > -0.002
%             tg_sig(i,k) = 0;
%         else
%             u = BDD.getInputs(round(mac_spd(:,k)-1,3));
%             tg_sig(i,k) = u(i);
%         end
%     end
% else
%     tg_sig(:,k) = tg_sig(:,k-1);
% end
%% Basic MPC
% global mac_spd xc MPCobj
% if t == 0
%     Ts = 0.02;%time step
%     load('sysFullMPC.mat')
%     CSTR = sysFull;
%     CSTR = c2d(CSTR,Ts);
%     CSTR.InputGroup.MV = [1:3];%inputs
%     CSTR.InputGroup.UD = [4:6];%unmeasured disturbances
%     MV = struct('Min',zeros(1,3));
%     old_status = mpcverbosity('off');
%     MPCobj = mpc(CSTR,Ts,[],[],[],MV);
%     MPCobj.PredictionHorizon = 2;
%     xc = mpcstate(MPCobj);
% end
% if k == 1
%    aa = [0;0;0];
% else
%    aa = tg_sig(:,k-1);
% end
% u = mpcmove(MPCobj,xc,mac_spd(:,k)-1,[0;0;0],aa);
% tg_sig(:,k) = u;
%% One Step MPC
% global mac_spd xmpc MPCobj
% Ts=2;
% if t == 0
%     tg_sig(:,k) = [0;0;0];
%     load('3StatesToyMPCObj.mat','MPCobj')
%     xmpc = mpcstate(MPCobj);
% end
% if rem(t,Ts)==0
%     u=mpcmove(MPCobj,xmpc,mac_spd(:,k)-1,r);
%     tg_sig(:,k) = u;
% else
%     tg_sig(:,k) = tg_sig(:,k-1);
% end
%% Secondary Frequency Response
global n_tg mac_spd intf timestep sw_con tg_con
K = [2;2;2;2;2;2;2;2;2];
if n_tg~=0
    if t == 0
        intf = zeros(size(tg_con,1),1);
        timestep = sw_con(1,7);
    end
    for j=1:size(tg_con,1)
        intf(j,k+1) = intf(j,k) + timestep*(mac_spd(j,k)-1);
        tg_sig(j,k) = -K(j)*intf(j,k);
    end
end

%% Consensus Subgradient Update Method
% global tg_sig n_tg lambda tg_pot lmod_sig P mac_con
% A = [1/3 1/3 1/3;
%     1/3 1/3 1/3;
%     1/3 1/3 1/3];
% 
% weights =  [1 2 0;
%     2 2 0;
%     3 2 0];
% alpha = 0.5;
% step = 2;
% 
% %ref = tg_pot(:,5);
% base_load = [1.25;0.9;1];
% 
% if n_tg~=0
%     if t == 0
%         tg_sig(:,k) = zeros(size(mac_con,1),1);
%         lambda(:,k) = zeros(3,1);
%         P(:,k) = tg_pot(:,5);
%         for j=1:3
%             lambda(j,k+1) = A(j,:)*lambda(:,k) + alpha*(sum(base_load)+sum(lmod_sig(:,k))-sum(P(:,k)));
%         end
%     else
%         if t < 2 || rem(t,step)==0
%             for j=1:3
%                 P(j,k) = (lambda(j,k) - weights(j,2))/(2*weights(j,1));
%             end
%             lambda(:,k+1) = A*lambda(:,k) + alpha*(sum(base_load)+sum(lmod_sig(:,k))-sum(P(:,k)));
%             tg_sig(:,k) = P(:,k)-tg_pot(:,5);
%         else
%             tg_sig(:,k) =tg_sig(:,k-1);
%             lambda(:,k+1)=lambda(:,k);
%         end
%     end
%     
% end
%%
return
