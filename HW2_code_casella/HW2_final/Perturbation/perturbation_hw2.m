%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Save empty dates and dseries objects in memory.
dates('initialize');
dseries('initialize');
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'perturbation_hw2';
M_.dynare_version = '4.5.6';
oo_.dynare_version = '4.5.6';
options_.dynare_version = '4.5.6';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('perturbation_hw2.log');
M_.exo_names = 'ez';
M_.exo_names_tex = 'ez';
M_.exo_names_long = 'ez';
M_.exo_names = char(M_.exo_names, 'ealpha');
M_.exo_names_tex = char(M_.exo_names_tex, 'ealpha');
M_.exo_names_long = char(M_.exo_names_long, 'ealpha');
M_.endo_names = 'c';
M_.endo_names_tex = 'c';
M_.endo_names_long = 'c';
M_.endo_names = char(M_.endo_names, 'l');
M_.endo_names_tex = char(M_.endo_names_tex, 'l');
M_.endo_names_long = char(M_.endo_names_long, 'l');
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names_long = char(M_.endo_names_long, 'k');
M_.endo_names = char(M_.endo_names, 'r');
M_.endo_names_tex = char(M_.endo_names_tex, 'r');
M_.endo_names_long = char(M_.endo_names_long, 'r');
M_.endo_names = char(M_.endo_names, 'u');
M_.endo_names_tex = char(M_.endo_names_tex, 'u');
M_.endo_names_long = char(M_.endo_names_long, 'u');
M_.endo_names = char(M_.endo_names, 'v');
M_.endo_names_tex = char(M_.endo_names_tex, 'v');
M_.endo_names_long = char(M_.endo_names_long, 'v');
M_.endo_names = char(M_.endo_names, 'ev');
M_.endo_names_tex = char(M_.endo_names_tex, 'ev');
M_.endo_names_long = char(M_.endo_names_long, 'ev');
M_.endo_names = char(M_.endo_names, 'm');
M_.endo_names_tex = char(M_.endo_names_tex, 'm');
M_.endo_names_long = char(M_.endo_names_long, 'm');
M_.endo_names = char(M_.endo_names, 'z');
M_.endo_names_tex = char(M_.endo_names_tex, 'z');
M_.endo_names_long = char(M_.endo_names_long, 'z');
M_.endo_names = char(M_.endo_names, 'alpha');
M_.endo_names_tex = char(M_.endo_names_tex, 'alpha');
M_.endo_names_long = char(M_.endo_names_long, 'alpha');
M_.endo_partitions = struct();
M_.param_names = 'beta';
M_.param_names_tex = 'beta';
M_.param_names_long = 'beta';
M_.param_names = char(M_.param_names, 'theta');
M_.param_names_tex = char(M_.param_names_tex, 'theta');
M_.param_names_long = char(M_.param_names_long, 'theta');
M_.param_names = char(M_.param_names, 'gamma');
M_.param_names_tex = char(M_.param_names_tex, 'gamma');
M_.param_names_long = char(M_.param_names_long, 'gamma');
M_.param_names = char(M_.param_names, 'eta');
M_.param_names_tex = char(M_.param_names_tex, 'eta');
M_.param_names_long = char(M_.param_names_long, 'eta');
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, 'delta');
M_.param_names_long = char(M_.param_names_long, 'delta');
M_.param_names = char(M_.param_names, 'rhoz');
M_.param_names_tex = char(M_.param_names_tex, 'rhoz');
M_.param_names_long = char(M_.param_names_long, 'rhoz');
M_.param_names = char(M_.param_names, 'sigmaz');
M_.param_names_tex = char(M_.param_names_tex, 'sigmaz');
M_.param_names_long = char(M_.param_names_long, 'sigmaz');
M_.param_names = char(M_.param_names, 'mualpha');
M_.param_names_tex = char(M_.param_names_tex, 'mualpha');
M_.param_names_long = char(M_.param_names_long, 'mualpha');
M_.param_names = char(M_.param_names, 'rhoalpha');
M_.param_names_tex = char(M_.param_names_tex, 'rhoalpha');
M_.param_names_long = char(M_.param_names_long, 'rhoalpha');
M_.param_names = char(M_.param_names, 'sigmaalpha');
M_.param_names_tex = char(M_.param_names_tex, 'sigmaalpha');
M_.param_names_long = char(M_.param_names_long, 'sigmaalpha');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 10;
M_.param_nbr = 10;
M_.orig_endo_nbr = 10;
M_.aux_vars = [];
M_.Sigma_e = zeros(2, 2);
M_.Correlation_matrix = eye(2, 2);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
M_.det_shocks = [];
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
M_.hessian_eq_zero = 0;
erase_compiled_function('perturbation_hw2_static');
erase_compiled_function('perturbation_hw2_dynamic');
M_.orig_eq_nbr = 10;
M_.eq_nbr = 10;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./' M_.fname '_set_auxiliary_variables.m'], 'file') == 2;
M_.lead_lag_incidence = [
 1 7 0;
 0 8 0;
 2 9 0;
 0 10 17;
 3 11 0;
 0 12 18;
 4 13 0;
 0 14 19;
 5 15 0;
 6 16 0;]';
M_.nstatic = 1;
M_.nfwrd   = 3;
M_.npred   = 6;
M_.nboth   = 0;
M_.nsfwrd   = 3;
M_.nspred   = 6;
M_.ndynamic   = 9;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(10, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(10, 1);
M_.NNZDerivatives = [39; 92; 402];
close all
M_.params( 1 ) = 0.99;
beta = M_.params( 1 );
M_.params( 2 ) = 0.5;
theta = M_.params( 2 );
M_.params( 3 ) = 10;
gamma = M_.params( 3 );
M_.params( 5 ) = 0.1;
delta = M_.params( 5 );
M_.params( 6 ) = 0.95;
rhoz = M_.params( 6 );
M_.params( 7 ) = 0.005;
sigmaz = M_.params( 7 );
M_.params( 8 ) = 0.03;
mualpha = M_.params( 8 );
M_.params( 9 ) = 0.9;
rhoalpha = M_.params( 9 );
M_.params( 10 ) = 0.01;
sigmaalpha = M_.params( 10 );
l_ss     = 100;
alpha_ss = 0.3;
k_ss     = (beta*alpha_ss*l_ss^(1-alpha_ss)/(1+beta*(delta-1)))^(1/(1-alpha_ss));
c_ss     = -delta*k_ss + k_ss^alpha_ss*l_ss^(1-alpha_ss);
M_.params( 4 ) = (1-alpha_ss)*k_ss^alpha_ss*l_ss^((-1)-alpha_ss)/c_ss;
eta = M_.params( 4 );
u_ss     = (log(c_ss)-eta*l_ss^2/2);
v_ss     = (log(c_ss)-eta*l_ss^2/2);
ev_ss    = v_ss^(1-gamma);
m_ss     = beta;
r_ss     = (alpha_ss*k_ss^(alpha_ss-1)*l_ss^(1-alpha_ss)+1-delta); 
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 1 ) = c_ss;
oo_.steady_state( 2 ) = l_ss;
oo_.steady_state( 3 ) = k_ss;
oo_.steady_state( 4 ) = r_ss;
oo_.steady_state( 6 ) = v_ss;
oo_.steady_state( 5 ) = u_ss;
oo_.steady_state( 7 ) = ev_ss;
oo_.steady_state( 8 ) = m_ss;
oo_.steady_state( 10 ) = alpha_ss;
oo_.exo_steady_state( 2 ) = 0;
oo_.steady_state( 9 ) = 0;
oo_.exo_steady_state( 1 ) = 0;
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 1;
M_.Sigma_e(2, 2) = 1;
steady;
oo_.dr.eigval = check(M_,options_,oo_);
options_.k_order_solver = 1;
options_.irf = 100;
options_.order = 3;
var_list_ = char('k','l','c','z','alpha');
info = stoch_simul(var_list_);
save('perturbation_hw2_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('perturbation_hw2_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('perturbation_hw2_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('perturbation_hw2_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('perturbation_hw2_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('perturbation_hw2_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('perturbation_hw2_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
