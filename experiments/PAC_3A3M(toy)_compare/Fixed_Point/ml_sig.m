function f = ml_sig(t,k)
load('SysFull');
load('disturbances');
% Syntax: f = ml_sig(t,k)
%4:40 PM 15/08/97
% defines modulation signal for lmod control
global mac_ang mac_spd eqprime edprime tg1 tg2 tg3 lmod_st lmod_sig n_lmod
global x_init check T Ti sysb psikd psikq noise
noise = disturbances;
f=0; %dummy variable

% you modify the following to do what you want with the load
% lmod_con must be specified in the data file
% and the load bus must be in the nonconforming load list.
if n_lmod~=0
    %% disturbances
    if t<1
        lmod_sig(:,k)= zeros(n_lmod,1);
    elseif rem(t,0.5)==0
        lmod_sig(2,k) = 0.3; %noise(1,k); %set to 3 with no noise
        lmod_sig(3,k) = 0.2; %noise(2,k); % set to 2 with no noise
    else
        lmod_sig(:,k)=lmod_sig(:,k-1);
    end
    %% balanced model reduction
    x = [mac_ang(1,k);
        mac_spd(1,k);
        eqprime(1,k);
        psikd(1,k);
        edprime(1,k);
        psikq(1,k);
        tg1(1,k);
        tg2(1,k);
        tg3(1,k);
        
        mac_ang(2,k);
        mac_spd(2,k);
        eqprime(2,k);
        psikd(2,k);
        edprime(2,k);
        psikq(2,k);
        tg1(2,k);
        tg2(2,k);
        tg3(2,k);
        
        mac_ang(3,k);
        mac_spd(3,k);
        eqprime(3,k);
        psikd(3,k);
        edprime(3,k);
        psikq(3,k);
        tg1(3,k);
        tg2(3,k);
        tg3(3,k);
        
        
        lmod_st(1,k);
        lmod_st(2,k);
        lmod_st(3,k)];
    
    if t == 0
        [sysb,g,T,Ti] = balreal(sysFull); % balance system
        x_init = x;
        check = [];
    end
    x = x-x_init;
    xb = T*x;
    xb = xb(1:3);
    %% SCOTS Controller
    BDD = SymbolicSet('04sec_t008a015.bdd');
    if rem(t,0.4) == 0
        if (-0.0115*xb(1) -0.2296*xb(2) + 0.0412*xb(3)) > -0.008
            lmod_sig(1,k) = 0;
        else
            disp('SCOTS control')
            u = -min(BDD.getInputs(xb));
            u = 0;
            lmod_sig(1,k) = u;
        end
    else
        lmod_sig(:,k) = lmod_sig(:,k-1);
    end
end
end
