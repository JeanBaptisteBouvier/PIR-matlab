function output = Targeting(choice, LVLH, SYN, times_lin, times_cont, settings, dro_T, phi0, n_transf)
% 
% Return the structur sol, it does not contain the name of the targeting algo
% 
% sol.err1   initial error due to the algo
% sol.fg     first guess
% sol.Lam    Lambert solution
% sol.Lam_t  Times of the Lambert solution
% sol.it     number of iterations
% 
% if choice == [1 1 1 0]
% sol.err_CW    CW error
% sol.err_SL    SL error
% sol.err_LR    LR error
% 
% The code begins with if/end/if/end/... and no ifelse because in the case
% choice = [1 1 1 0], the first three cases have to be done.
% 
% 
% 
% 
% 

%% Variables initialisation 

r_LVLH = LVLH.r;
v_LVLH = LVLH.v;
r_syn = SYN.r;
v_syn = SYN.v;
r_next = SYN.r_next;
T = times_lin.T;
tt = times_lin.tt;
tt_2 = times_lin.tt_2;
TOF_cont = times_cont.TOF;
int_cont = times_cont.int;

intervals = settings.Luq_int;
toll = settings.toll;
it_max = settings.it_max;
options = settings.options;
cr3bp = settings.cr3bp;   

choice_help = 0;

%% Computation of the first errors

if choice(1)        %% Clohessy Wiltshire (LVLH)
    [x_CW_LVLH, ~, ~] = Docking_CW(r_LVLH, v_LVLH, tt, tt_2, intervals, dro_T, cr3bp);
    x_CW_syn = rlvlh2rsyn(tt, x_CW_LVLH', state_time(tt,dro_T)', cr3bp.mu)';
    x_syn = x_CW_syn;
    err_CW = Lambert_cycle_errors(x_CW_syn, r_next, state_time(tt, dro_T), phi0, T, cr3bp, options);
    sol.err1 = err_CW;   % Initial error 
end

if  choice(2)  %% Straight Line
    [x_SL_LVLH, ~, ~] = Docking_SL(r_LVLH, v_LVLH, T);
    x_SL_syn  = rlvlh2rsyn(tt, x_SL_LVLH', state_time(tt,dro_T)', cr3bp.mu)';
    x_syn = x_SL_syn;
    err_SL = Lambert_cycle_errors(x_SL_syn, r_next, state_time(tt,dro_T), phi0, T, cr3bp, options);
    sol.err1 = err_SL;   % Initial error
end

if choice(3) %% Linearized Relative
    [x_LR_syn, ~, ~]  = Docking_LR(r_syn, v_syn, tt, tt_2, dro_T, intervals, cr3bp);
    x_syn = x_LR_syn;
    err_LR = Lambert_cycle_errors(x_LR_syn, r_next, state_time(tt,dro_T), phi0, T, cr3bp, options);
    sol.err1 = err_LR;   % Initial error
end



%% Complete Computation

if  sum(choice) == 1
    [fg, Transfer, Lam_t, it, err] = Lambert_cycle(x_syn, r_next, state_time(tt, dro_T), phi0, T, toll, it_max, cr3bp, options, tt);
    
elseif sum(choice) == 3
    sol.err1 = min([err_CW err_SL err_LR]);
    sol.err_CW = err_CW;
    sol.err_SL = err_SL;
    sol.err_LR = err_LR;
   
    if sol.err1 == err_CW
        [fg, Transfer, Lam_t, it, err] = Lambert_cycle(x_CW_syn, r_next, state_time(tt, dro_T), phi0, T, toll, it_max, cr3bp, options, tt);
        fprintf('C-W smallest error for transfer %d\n', n_transf) 
    elseif sol.err1 == err_SL
        [fg, Transfer, Lam_t, it, err] = Lambert_cycle(x_SL_syn, r_next, state_time(tt,dro_T), phi0, T, toll, it_max, cr3bp, options, tt);
        fprintf('SL smallest error for transfer %d\n', n_transf)
    else
        [fg, Transfer, Lam_t, it, err] = Lambert_cycle(x_LR_syn, r_next, state_time(tt,dro_T), phi0, T, toll, it_max, cr3bp, options, tt);
        fprintf('LR smallest error for transfer %d\n', n_transf)
    end
    
end

if it < it_max            
    if symmetric_sol(Transfer) 
        % Output generation
        sol.fg = fg;
        sol.Lam = Transfer;
        sol.Lam_t = Lam_t;
        sol.it = it;
        sol.err = err;
    else
        disp('Error: Symmetric Solution')
        choice_help = 1;
    end           
else
    disp('MAX iteration reached')
    choice_help = 1;
end
    

%% Continuation in case of problem

if choice(4) || choice_help 
    if choice_help
        disp('Continuation Algorithm is here to help')
    end
    
   [fg, Transfer, Lam_t, it, err] = Continuation(T, TOF_cont, int_cont, r_syn, v_syn, r_next, tt, dro_T, intervals, phi0, toll, options, it_max, cr3bp); 

   % Output generation 
   sol.fg = fg;
   sol.Lam = Transfer;
   sol.Lam_t = Lam_t;
   sol.it = it;
   sol.err = err;
end

%% Final output
sol.help = choice_help;
output.sol = sol;

end

