%% ------------------------------------------------------------------------
% This code was jointly written by Dr Andrew Rawlinson and A/Prof Roslyn Hickson
% The code was most recently updated in Januaray 2023
%% ------------------------------------------------------------------------
function stigma
% Simple population-level model to explore the effects of stigma on 
% prevalence of anxiety/depression in a community

more off;

%% ------------------------------------------------------------------------
% flags to control what is run and/or plotted. Absence of any given .mat
% will override a flag being set to zero
flag_plot_all = 1; % plot all results -- overrides individual settings if 1

flag_calibrate_beta = 0; % if 1, will calibrate beta to specified target 
    %  prevalence, otherwise loads all params from file

flag_run_1D_sweeps = 0;  % if 1, explores effects of p,h,g,k
flag_plot_1D_sweeps = 0;      

flag_run_multiple_1D_sweeps =0; % if 1, param sweep for pairs [beta, nu]
flag_plot_multi_1D_sweeps = 0;

flag_run_2D_sweeps = 0;  % if 1, explore combined effects of p-h and p-g
flag_plot_2D_sweeps = 0;

flag_fresh_sensanal = 0;  % if 1, will run the prcc, otherwise loads file
flag_plot_sensanal = 0;

%% ------------------------------------------------------------------------
% set params needed for simulation
params = initialiseParameters;  % model params - grab defaults from csv
sim.fonts = 16;     % set fontsize
sim.fignum = 1;     % increment fignum so don't overwrite
sim.tspan = 0:1:365.25*50;  % 50 years
sim.options = [];   % in case want to set tols or similar for DE solver
sim.alpha_0 = 0.214;    % proportion starting in "acute" stage
sim.total_pop = 1;  % assume N=1     
num_repeats = 1e4;  % for the PRCC 
num_fit_repeats = 1e3;  % for the calibration of beta and nu

% set params to explore effects of, including a default value
svec = 0:0.01:1;  % proportion experiencing stigma
stig_default = 0.15;  % default value prop'n commmunity experiencing stigma 
stig_plotname = 'Proportion of the population experiencing stigma (p)';

hvec = 0:0.01:1;  % affects gamma
h_default = params.h;
h_plotname = 'Effect of stigma on treatment seeking behaviour (h)';

gvec = 0:0.01:1;  % affects sigma
g_default = params.g;
g_plotname = 'Effect of stigma on rate become managed (g)';

kvec = 0.5:0.01:5;  % affects psi
k_default = params.k;
k_plotname = 'Effect of stigma on defaulting from treatment (k)';

%% ------------------------------------------------------------------------
% calibration of model parameters
if flag_calibrate_beta ==1  % otherwise use values in `inputs.csv`
    disp("Calibrating beta and nu")

    
    beta = zeros(1,num_fit_repeats);
    nu = zeros(1,num_fit_repeats);
    prevalences = zeros(1,num_fit_repeats);
    for i=1:num_fit_repeats
        disp(i);
        tmp = calibrate_transmission(params, sim);
        beta(i) = tmp(1);
        nu(i) = tmp(2);

        prevalences(i) = calc_ss_prevalence(beta(i), nu(i), params, sim);
    end

    save('calibrate_beta_nu.mat','beta', 'nu', 'prevalences', 'params',...
        'sim');

    % update here for use in these runs
    [exp_beta, exp_nu] = clean_calibration();
    params.beta = exp_beta;
    params.nu = exp_nu;
    params.eta = 2*exp_nu;
    params = updateParameters(params);

end

%% ------------------------------------------------------------------------
% 1D sweeps for expected values of [beta,nu]
if flag_run_1D_sweeps==1 || exist('1D_sweeps.mat','file')==0
    disp("running 1D stigma parameter sweeps")
    % run through param variations
    [sim, s_prev, ~,~,~,~,~,~,~,~] = run_stigma(params, svec, 'p', ...
        stig_default, stig_plotname, sim, 1);
    params.p = stig_default;
    params = updateParameters(params);
    
    [sim, h_prev, ~,~,~,~,~,~,~,~] = run_stigma(params, hvec, 'h', ...
        h_default, h_plotname, sim, 1);
    params.h = h_default;
    params = updateParameters(params);
    
    [sim, g_prev, ~,~,~,~,~,~,~,~] = run_stigma(params, gvec, 'g', ...
        g_default, g_plotname, sim, 1);
    params.g = g_default;
    params = updateParameters(params);
    
    [sim, k_prev, ~,~,~,~,~,~,~,~] = run_stigma(params, kvec, 'k', ...
        k_default, k_plotname, sim, 1);
    params.k = k_default;
    params = updateParameters(params);

    save('1D_sweeps.mat', 'svec', 's_prev', 'hvec', 'h_prev', 'gvec', ...
        'g_prev', 'kvec', 'k_prev','stig_plotname', 'h_plotname',...
        'g_plotname','k_plotname', 'params', 'sim');  
end

if flag_plot_all==1 || flag_plot_1D_sweeps==1                  
    sim = plot_1D_sweeps(sim);  % done in function to avoid overwriting inputs for 
        % subsequent simulations in this file
    
end

%% ------------------------------------------------------------------------
% 1D sweeps for pairs of refined [beta,nu] values
if flag_run_multiple_1D_sweeps ==1 || exist('repeat_betas_prev.mat', 'file')==0
    disp("running 1D stigma parameter sweeps for pairs of [beta,nu]")

    % storing to reset after this set of runs
    old_beta = params.beta;
    old_nu = params.nu;

    load("calibrate_beta_nu.mat", 'select_betas','select_nus');
    num_runs = length(select_betas);
    p_mult = zeros(num_runs, length(svec));
    h_mult = zeros(num_runs, length(hvec));
    g_mult = zeros(num_runs, length(gvec));

    num_ts = length(svec);
    p_Un = zeros(num_runs, num_ts);
    p_An = zeros(num_runs, num_ts);
    p_Tn = zeros(num_runs, num_ts);
    p_Mn = zeros(num_runs, num_ts);
    p_Us = zeros(num_runs, num_ts);
    p_As = zeros(num_runs, num_ts);
    p_Ts = zeros(num_runs, num_ts);
    p_Ms = zeros(num_runs, num_ts);

    for ind=1:num_runs

        params.beta = select_betas(ind);
        params.nu = select_nus(ind);
        params.eta = 2 * select_nus(ind);  % not strictly needed
        params = updateParameters(params);
    
        % run through param variations
        [sim, s_prev, Un, An, Tn, Mn, Us, As, Ts, Ms] = ...
            run_stigma(params, svec, 'p', stig_default, stig_plotname, ...
            sim, 0);
        p_mult(ind,:) = s_prev;
        p_Un(ind,:) = Un; 
        p_An(ind,:) = An;
        p_Tn(ind,:) = Tn;
        p_Mn(ind,:) = Mn;
        p_Us(ind,:) = Us;
        p_As(ind,:) = As;
        p_Ts(ind,:) = Ts;
        p_Ms(ind,:) = Ms;

        params.p = stig_default;
        params = updateParameters(params);

        [sim, h_prev, ~,~,~,~,~,~,~,~] = run_stigma(params, hvec, 'h', ...
            h_default, h_plotname, sim, 0);
        h_mult(ind,:) = h_prev;
        params.h = h_default;
        params = updateParameters(params);

        [sim, g_prev, ~,~,~,~,~,~,~,~] = run_stigma(params, gvec, 'g', ...
            g_default, g_plotname, sim, 0);
        g_mult(ind,:) = g_prev;
        params.g = g_default;
        params = updateParameters(params);

    end
    save('repeat_betas_prev.mat','p_mult','h_mult','g_mult', 'p_Un', ...
        'p_An','p_Tn','p_Mn','p_Us','p_As','p_Ts','p_Ms', ...
        'svec','hvec','gvec','stig_default',...
        'select_betas','select_nus', 'params', 'sim'); 

    % resetting to expected values for subsequent simulations
    params.beta = old_beta;
    params.nu = old_nu;
    params.eta = 2 * select_nus(ind);  % not strictly needed
    params = updateParameters(params);

end

if flag_plot_all ==1 || flag_plot_multi_1D_sweeps==1

    sim = plot_multi_1D_sweeps(sim);

end

%% ------------------------------------------------------------------------
% 2D sweeps of p,g,h,k for expected values [beta,nu]
if flag_run_2D_sweeps==1 || exist('2D_sweeps.mat', 'file')==0
    disp("running 2D stigma parameter sweeps")
    % combined effects of p-h and p-g on ss prevalence
    sim = run_2D_stigma(params, svec, hvec, gvec, ...
        kvec, [stig_default, h_default, g_default, k_default], ...
        {stig_plotname, h_plotname, g_plotname, k_plotname}, sim);

end

if flag_plot_all==1 || flag_plot_2D_sweeps==1
    sim = plot_2D_sweeps(sim);

end


%% ------------------------------------------------------------------------
% Sensitivity analysis
if flag_fresh_sensanal == 1 || exist('prcc_results.mat','file')==0
    disp("running sensitivity analysis")
    params.p = stig_default; 
    ipops = initialConditions(params.p, sim.alpha_0, sim.total_pop);
    lhsonly(num_repeats);
    runlhs(sim.tspan,sim.options,ipops);  
    [delta, varnamesGreek] = prcconly;
end

if flag_plot_all ==1 || flag_plot_sensanal==1
    load('prcc_results.mat') % loads the Greek variable names                                   
    tornadoplot(delta, varnamesGreek); 
end

end

% =========================================================================

function sim = plot_multi_1D_sweeps(sim)
% Plots the results of the 1D sweeps for each of the values in the 
% [beta,nu] chain, here and not elsewhere to avoid overwriting key params
%
% Inputs:
%    `repeat_betas_prev.mat` = saved inputs and outputs required to create 
%     the plots
%
% Outputs:
%    saved figures:
%         prevalences vs all params with distributions
%         proportion population vs p with distributions for no stigma
%         proportion population vs p with distributions with             stigma


    load('repeat_betas_prev.mat', 'p_mult','h_mult','g_mult', 'p_Un', ...
        'p_An','p_Tn','p_Mn','p_Us','p_As','p_Ts','p_Ms', ...
        'svec','hvec','gvec','stig_default');
    
    
    sim.fignum = sim.fignum + 1;
    plot_mult_tx(p_mult, svec, '-b', sim.fignum);
    plot_mult_tx(h_mult, hvec, '--r', sim.fignum);
    plot_mult_tx(g_mult, gvec, '-.k', sim.fignum);
    legend('','Proportion with stigma (p)', '',...
        'Effect of stigma on treatment seeking (h)', '',...
        'Effect of stigma on becoming managed (g)', 'location','best')
    set(gca,'FontSize',sim.fonts);
    axis([0, 1, 0, 0.8]);
    box on;
    savefig('./figures/all_params_prev_dist.fig');

    % time series for varying p -- no stigma strata
    sim.fignum = sim.fignum + 1;
    plot_mult_tx(p_Un, svec, '-b', sim.fignum);
    plot_mult_tx(p_An, svec, '--r', sim.fignum);
    plot_mult_tx(p_Tn, svec, '-.k', sim.fignum);
    plot_mult_tx(p_Mn, svec, ':g', sim.fignum);
    hold on;
    plot([stig_default, stig_default], [0, 1], ':k');
    xlabel('Proportion of the population experiencing stigma (p)');
    ylabel('Proportion of the population');
    def_text = ['Default value (' num2str(stig_default) ')'];
    legend('','U_n','','A_n','','T_n','','M_n', ...
        def_text, 'Location','best');
    set(gca,'FontSize',sim.fonts);
    box on;
    savefig('./figures/p_dist_n.fig');

    % time series for varying p -- with stigma strata
    sim.fignum = sim.fignum + 1;
    plot_mult_tx(p_Us, svec, '-b', sim.fignum);
    plot_mult_tx(p_As, svec, '--r', sim.fignum);
    plot_mult_tx(p_Ts, svec, '-.k', sim.fignum);
    plot_mult_tx(p_Ms, svec, ':g', sim.fignum);
    hold on;
    plot([stig_default, stig_default], [0, 1], ':k');
    xlabel('Proportion of the population experiencing stigma (p)');
    ylabel('Proportion of the population');
    legend('','U_s','','A_s','', 'T_s','','M_s', ...
        def_text, 'Location','best');
    set(gca,'FontSize',sim.fonts);
    box on;
    savefig('./figures/p_dist_s.fig');
    
end

% =========================================================================

function sim = plot_1D_sweeps(sim)
% Plots the results of the 1D sweeps, here and not elsewhere to avoid
% overwriting key params
%
% Inputs:
%    `1D_sweeps.mat` = saved inputs and outputs required to create the
%    plots
%
% Outputs:
%    saved figure of prevalence vs all params

    load('1D_sweeps.mat','svec', 's_prev', 'hvec', 'h_prev', 'gvec', ...
        'g_prev', 'kvec', 'k_prev','stig_plotname', 'h_plotname',...
        'g_plotname','k_plotname');
    
    sim.fignum = sim.fignum + 1;
    figure(sim.fignum), clf, box on, hold on;
    plot(svec,s_prev,'-b',hvec,h_prev,'--r',gvec,g_prev,'-.g',kvec, ...
        k_prev,':k','LineWidth',2);
    xlabel('Parameter value');
    ylabel('Prevalence of anxiety and/or depression');
    legend(stig_plotname,h_plotname,g_plotname,k_plotname, ...
                    'Location','best');
    axis([0,1.1,0,0.6]);
    set(gca,'FontSize',sim.fonts);
    savefig('./figures/all_params_prev.fig')
end

% =========================================================================

function plot_mult_tx(mult,xvec,colour,fignum)
% Plots of the runs of 1D parameter sweeps for the chain of [beta,nu]
% values calibrated to target prevalence (and refined in
% `clean_calibration`).
%
% Inputs:
%   mult = the steady state prevalences from multiple [beta,nu] combos
%   xvec = vector of the parameter value being considered
%   colour = chosen colour for this parameter
%   fignum = specified figure number, so can hold on from previous combos
%
% Output:
%   not yet saved figure, in `fignum`

    % determine boundaries
    mult_min = min(mult);
    mult_max = max(mult);
    mult_exp = median(mult);

    % plot shaded regions
    figure(fignum), hold on, box on; 
    patch([xvec, fliplr(xvec)],[mult_min, fliplr(mult_max)],colour, ...
        'FaceAlpha',.2,'EdgeColor','none');
    plot(xvec,mult_exp,colour,'LineWidth',2);
    xlabel('Parameter value (see legend)');
    ylabel('Prevalence of anxiety and/or depression');
    
end

% =========================================================================

function [exp_beta, exp_nu] = clean_calibration()
% the optimiser results in a wide range of steady state prevalences. Here
% we reduce the pairs to those within a target range, perform linear
% regression, and randomly pick one of the final pairs to become the
% expected values.
%
% Inputs:
%   `calibrate_beta_nu.mat` = this is where all key values have been
%       stored from the multiple runs through `calibrate_transmission` for
%       different starting values
% 
% Outputs:
%   exp_beta = randomly selected expected value of beta from subset
%       resulting in refined range around target prevalence
%   exp_nu = randomly selected expected value of nu from subset
%       resulting in refined range around target prevalence
%   `calibrate_beta_nu.mat` = where chain of cleaned params etc stored
%   linear regressions = printed to command window

    load('calibrate_beta_nu.mat');

%     % plot things
%     figure; plot(prevalences);
%     figure; plot(beta, 'xb'); hold on; 
%     plot(nu, 'or');
%     figure; histogram(beta); title('beta')
%     figure; histogram(nu); title('nu')
%     figure; figure; plot(beta, nu, 'x');    

    % linear regression of nu= m* beta + b
    p1 = polyfit(beta, nu, 1);
    SSresid1 = sum((nu - polyval(p1, beta)).^2);
    SStotal1 = (length(nu)-1) * var(nu);
    rsq1 = 1 - SSresid1 / SStotal1;
    fprintf(['linear regression between raw params:\n ' ...
        'nu = %4.8f * beta + %4.8f,\n    with an R^2 of %4.8f'], p1(1), ...
        p1(2), rsq1);

    % clean up for saving
    cleaner_indices = find(prevalences<(sim.alpha_0+0.001) & ...
        prevalences>(sim.alpha_0-0.1)); 
    cleaner_betas = beta(cleaner_indices);
    cleaner_nus = nu(cleaner_indices);
    cleaner_prevalences = prevalences(cleaner_indices);

    % linear regression of nu= m* beta + b
    p = polyfit(cleaner_betas, cleaner_nus, 1);
    SSresid = sum((cleaner_nus - polyval(p, cleaner_betas)).^2);
    SStotal = (length(cleaner_nus)-1) * var(cleaner_nus);
    rsq = 1 - SSresid / SStotal;
    fprintf(['linear regression between cleaned params:\n' ...
        ' nu = %4.8f * beta + %4.8f,\n' ...
        ' with an R^2 of %4.8f\n'], p(1), p(2), rsq);

    % find distribution estimate of cleaner_nus
    pd = fitdist(cleaner_nus', 'Gamma');
    playind = find(cleaner_nus<(median(pd)+std(pd)/2) & ...
        cleaner_nus>(median(pd)-std(pd)/2));
    select_nus = cleaner_nus(playind);
    select_betas = cleaner_betas(playind);
    select_prevalences = cleaner_prevalences(playind);

    % print details wanted for inputs.csv
    rand_ind = randi(size(select_betas));
    min_beta = min(beta);
    max_beta = max(beta);
    exp_beta = select_betas(rand_ind);

    min_nu = min(nu);
    max_nu = max(nu);
    exp_nu = select_nus(rand_ind);

    % overwrite inputs.csv
    readindata = readtable('inputs.csv', 'ReadVariableNames', true);

    readindata(1,'beta') = {min_beta};
    readindata(2,'beta') = {max_beta};
    readindata(3,'beta') = {exp_beta};
    readindata(1,'nu') = {min_nu};
    readindata(2,'nu') = {max_nu};
    readindata(3,'nu') = {exp_nu};
    readindata(1,'eta') = {2 * min_nu};
    readindata(2,'eta') = {2 * max_nu};
    readindata(3,'eta') = {2 * exp_nu};

    writetable(readindata, 'inputs.csv');


    save('calibrate_beta_nu.mat','beta', 'nu', 'prevalences', ...
        'select_betas', 'select_nus', 'select_prevalences',  ...
        'params', 'sim');
end

% ====================

function calibrated_params = calibrate_transmission(params, sim)
% calibrate the transmission rate based on target population prevalence
% specified as one of the `sim` parameters
%
% Inputs:
%   params = model parameter values
%   sim = various inputs needed to run the simulations, such as tspan
%
% Outputs:
%   calibrated_params = estimated values for [beta, nu]

% reduce parmeter dependency
newfun = @(calvals) calc_sse(calvals, params, sim);
lb = [0,0]; % lower bounds of possible parameter ranges to use as 
        % conditions in the optimiser
ub = [1,1]; % upper bounds of possible parameter ranges to use as 
        % conditions in the optimiser
x0 = - lb + (ub-lb).*rand(1,2);  % randomise starting point based on ranges

calibrated_params = fmincon(newfun,x0,[],[],[],[],lb,ub);  % estimate

end

% ====================

function sse = calc_sse(calvals, params, sim)
% calculate the sum of squared error between the steady state prevalence 
% for this combination of parameter values for the system of ODEs and the
% target prevalence specified as one of the `sim` parameters
%
% Inputs:
%   calvals = vector for updated values of [beta, nu]
%   params = model parameter values
%   sim = various inputs needed to run the simulations, such as tspan
%   
% Output:
%   sum of squared error (NB in percentages not proportion)

    beta = calvals(1);
    nu = calvals(2);

    new_prevalence = calc_ss_prevalence(beta, nu, params, sim);

    sse = (100*sim.alpha_0 - 100*new_prevalence)^2;

end

% =========================================================================

function prevalence = calc_ss_prevalence(beta, nu, params, sim)
% calculate the steady state prevalence of the model for a given set of
% parameter values
%
% Inputs:
%   beta = value of "transmission" rate to explore, to override `params`
%   nu = value of spontaneous rate, to override `params`
%   params = model parameter values
%   sim = various inputs needed to run the simulations, such as tspan
%
% Output:
%   prevalence = steady state prevalence for these parameter values

    params.beta = beta;
    params.nu = nu;
    params = updateParameters(params);
    ipops = initialConditions(params.p, sim.alpha_0, sim.total_pop);
    [~,pops] = ode15s(@stigmaODEmodel,sim.tspan,ipops,sim.options,params);
%     sim.fignum = sim.fignum + 1;
%     figure(sim.fignum), clf, box on, hold on;
%     plot(t,pops,'LineWidth',2);
%     xlabel('time (days)');
%     ylabel('population');
%     legend('U_n','A_n','T_n','M_n','U_s','A_s','T_s','M_s',...
%            'Location','best');
%     axis([0,20000,-0.2,1]);
%     set(gca,'FontSize',sim.fonts);
%     drawnow('update');
    prevalence = 1 - pops(end,1) - pops(end,5); % 1 - Un - Us
end

% =========================================================================

function sim = run_2D_stigma(params, pvec, hvec, ...
    gvec, kvec, default_values, plot_name, sim)
% Solves the system of ODEs for variations in two sets of model parameters.
% Later plotted in `plot_2D_sweeps` function.
%
% Inputs:
%   params = model parameter values
%   pvec = vector of `p` values for varying
%   hvec = vector of `h` values for varying
%   gvec = vector of `g` values for varying
%   kvec = vector of `k` values for varying
%   default_values = vector of default values of the parameters being 
%       varied, used for resetting 
%   plot_name = cells of plot names for saving
%   sim = various inputs needed to run the simulations, such as tspan
%
% Outputs:
%   sim = updated values of inputs to run the simulations, such as
%       incremented figure number
%   `2D_sweeps.mat` where all the key outputs are saved, such as:
%       p_h = matrix of prevalences for combos of `p` and `h`
%       p_g = matrix of prevalences for combos of `p` and `h`
%       p_k = matrix of prevalences for combos of `p` and `h`
%       g_h = matrix of prevalences for combos of `g` and `h`


num_pvec = length(pvec);
% combined effects of p-h and p-g
p_h = zeros(num_pvec, length(hvec));
p_g = zeros(num_pvec, length(gvec));
p_k = zeros(num_pvec, length(kvec));
g_h = zeros(length(hvec), length(gvec));

for pind=1:num_pvec
    
    params.p = pvec(pind);
    params = updateParameters(params);

    [sim, h_prev, ~,~,~,~,~,~,~,~] = run_stigma(params, hvec, 'h', ...
        default_values(2), plot_name(2), sim, 0);
    p_h(pind,:) = h_prev;
    params.h = default_values(2);
    params = updateParameters(params);
    
    [sim, g_prev, ~,~,~,~,~,~,~,~] = run_stigma(params, gvec, 'g', ...
        default_values(3), plot_name(3), sim, 0);
    p_g(pind,:) = g_prev;
    params.g = default_values(3);
    params = updateParameters(params);
    
    [sim, k_prev, ~,~,~,~,~,~,~,~] = run_stigma(params, kvec, 'k', ...
        default_values(4), plot_name(4), sim, 0);
    p_k(pind,:) = k_prev;
    params.k = default_values(4);
    params = updateParameters(params);
end

params.p = default_values(1);  % reset p to default value

for hind = 1:length(hvec)
    params.h = hvec(hind);
    params = updateParameters(params);
    [sim, g_prev, ~,~,~,~,~,~,~,~] = run_stigma(params, gvec, 'g', ...
        default_values(3), plot_name(3), sim, 0);
    g_h(hind,:) = g_prev;
    params.g = default_values(3);
    params = updateParameters(params);

end

save('2D_sweeps.mat','p_h', 'p_g', 'p_k', 'g_h', 'sim', 'pvec', 'hvec', ...
    'gvec', 'kvec', 'plot_name');

end

% ====================

function sim = plot_2D_sweeps(sim)
% creates the combined contour and heatmap plots of the 2D parameter sweeps
% showing anxiety/depression prevalence in the steady state for the
% expected values
%
% No explicit inputs or outputs as loads outputs from file (see 
% `run_2D_stigma` function) and saves figures.

load('2D_sweeps.mat','p_h', 'p_g', 'p_k', 'g_h', 'pvec', 'hvec', ...
    'gvec', 'kvec', 'plot_name');

sim.fignum = sim.fignum + 1;
figure(sim.fignum), clf, box on, hold on;
imagesc([pvec(1),pvec(end)], [hvec(1),hvec(end)],p_h');
set(gca,'FontSize',sim.fonts);
cRange = caxis;
[C,h] = contour(pvec, hvec, p_h', '-k', 'ShowText','on'); 
caxis(cRange);
clabel(C,h,'FontSize',sim.fonts)
set(gca,'YDir','normal');
xlabel(plot_name{1});
ylabel(plot_name{2});
set(gca,'FontSize',sim.fonts);
clim(gca,[0 1]);
a = colorbar();
a.Label.String = 'Prevalence of anxiety and/or depression';
set(gca,'FontSize',sim.fonts);
savefig('./figures/prevalences_for_p_and_h.fig')

sim.fignum = sim.fignum + 1;
figure(sim.fignum), clf, box on, hold on;
imagesc([pvec(1),pvec(end)], [gvec(1),gvec(end)],p_g');
cRange = caxis;
[C,h] = contour(pvec, gvec, p_g', '-k', 'ShowText','on'); 
clabel(C,h,'FontSize',sim.fonts)
caxis(cRange);
set(gca,'YDir','normal');
clim(gca,[0 1]);
xlabel(plot_name{1});
ylabel(plot_name{3});
set(gca,'FontSize',sim.fonts);
a = colorbar();
a.Label.String = 'Prevalence of anxiety and/or depression';
set(gca,'FontSize',sim.fonts);
savefig('./figures/prevalences_for_p_and_g.fig')

sim.fignum = sim.fignum + 1;
figure(sim.fignum), clf, box on, hold on;
imagesc([pvec(1),pvec(end)], [kvec(1),kvec(end)],p_k');
cRange = caxis;
[C,h] = contour(pvec, kvec, p_k', '-k', 'ShowText','on'); 
clabel(C,h,'FontSize',sim.fonts)
caxis(cRange);
set(gca,'YDir','normal');
clim(gca,[0 1]);
xlabel(plot_name{1});
ylabel(plot_name{4});
set(gca,'FontSize',sim.fonts);
a = colorbar();
a.Label.String = 'Prevalence of anxiety and/or depression';
set(gca,'FontSize',sim.fonts);
savefig('./figures/prevalences_for_p_and_k.fig')

sim.fignum = sim.fignum + 1;
figure(sim.fignum), clf, box on, hold on;
imagesc([hvec(1),hvec(end)], [gvec(1),gvec(end)],g_h');
cRange = caxis;
[C,h] = contour(hvec, gvec, g_h', '-k', 'ShowText','on'); 
clabel(C,h,'FontSize',sim.fonts)
caxis(cRange);
set(gca,'YDir','normal');
clim(gca,[0 1]);
xlabel(plot_name{2});
ylabel(plot_name{3});
set(gca,'FontSize',sim.fonts);
a = colorbar();
a.Label.String = 'Prevalence of anxiety and/or depression';
set(gca,'FontSize',sim.fonts);
savefig('./figures/prevalences_for_h_and_g.fig')


end

% =========================================================================

function [sim, prevalence, Un, An, Tn, Mn, Us, As, Ts, Ms] = ...
    run_stigma(params, param_varying, param_name, default_value, ...
    plot_name, sim, plotflag)
% runs the solving of the system of ODEs for sweeps of parameter values
%
% Inputs:
%   params = model parameter values
%   param_varying = vector of the parameter being varied
%   param_name = string for name of parameter being varied
%   default_value = default value of the parameter being varied for
%       plotting and to enable resets
%   plot_name = string for saving plot
%   sim = various inputs needed to run the simulations, such as tspan
%   plotflag = whether or not to plot the outputs
%       
% Outputs:
%   sim = updated values of inputs to run the simulations, such as
%       incremented figure number
%   prevalence = vector of steady state anxiety/depression prevalences for
%       each value in the vector of `param_varying`
%   Un, An, Tn, Mn, Us, As, Ts, Ms = vectors of the steady state prevalence
%        from the ODEs for each value of `param_varying`

num_param_vals = length(param_varying);

% iterate over different levels of stigma prevalence in community
Un = zeros(1,num_param_vals);
An = Un;
Tn = Un;
Mn = Un;
Us = Un;
As = Un;
Ts = Un;
Ms = Un;
for i = 1:num_param_vals
    params.(param_name) = param_varying(i);
    params = updateParameters(params);
    ipops = initialConditions(params.p, sim.alpha_0, sim.total_pop);
    [t,pops] = ode15s(@stigmaODEmodel,sim.tspan,ipops,sim.options,params);
    
    if param_varying(i) == default_value && plotflag==1
        figure(sim.fignum),clf,box on;
        plot(t,pops,'LineWidth',2);
        xlabel('time (days)');
        ylabel('population');
        legend('U_n','A_n','T_n','M_n','U_s','A_s','T_s','M_s', ...
            'Location','best');
        axis([0,20000,-0.2,1]);
        title([plot_name ' = ',num2str(param_varying(i))]);
        set(gca,'FontSize',sim.fonts);
        drawnow('update');
        savefig(['./figures/time_series_' param_name '=' ...
            num2str(default_value) '.fig'])
    end
    
    Un(i) = pops(end,1);
    An(i) = pops(end,2);
    Tn(i) = pops(end,3);
    Mn(i) = pops(end,4);
    Us(i) = pops(end,5);
    As(i) = pops(end,6);
    Ts(i) = pops(end,7);
    Ms(i) = pops(end,8);
end

if plotflag==1
    sim.fignum = sim.fignum + 1;
    figure(sim.fignum), clf, box on, hold on;
    plot(param_varying,Un,'-b',param_varying,An,'--r', ...
        param_varying,Tn,'-.g',param_varying,Mn,':k','LineWidth',2);
    xlabel(plot_name);
    ylabel('Proportion of population');
    plot([default_value, default_value],[0,1],':k')
    def_text = ['Default value (' num2str(default_value) ')'];
    legend('U_n','A_n','T_n','M_n',def_text,'Location','best');
    axis([0,max(param_varying),0,1]);
    set(gca,'FontSize',sim.fonts);
    savefig(['./figures/' param_name '_n.fig'])
    
    
    sim.fignum = sim.fignum + 1;
    figure(sim.fignum), clf, box on, hold on;
    plot(param_varying,Us,'-b',param_varying,As,'--r', ...
        param_varying,Ts,'-.g',param_varying,Ms,':k','LineWidth',2);
    xlabel(plot_name);
    ylabel('Proportion of population');
    plot([default_value, default_value],[0,1],':k')
    legend('U_s','A_s','T_s','M_s', def_text,'Location','best');
    axis([0,max(param_varying),0,1]);
    set(gca,'FontSize',sim.fonts);
    savefig(['./figures/' param_name '_s.fig'])
end

prevalence = 1 - Us - Un;

end

% =========================================================================

function runlhs(tSpan,options,ipops)
% determine values used in the prcc sensitivity analysis using a latin
% hypercube sampling approach
%
% Inputs:
%    tSpan = time to run the DE solver for
%    options = DE solver options
%    ipops = initial conditions for the DE solver
%   `trilhspoints.txt` = generated in `lhsonly.m`
%   `fixedParams.txt` = fixed parameter values from `lhsonly.m`
% 
% Outputs:
%    saved in `lhsoutput.mat`                               


% Variable names are as defined by the first row in this txt file. 
readindata = readtable('trilhspoints.txt', 'ReadVariableNames', true);
readinParams = readtable('fixedParams.txt', 'ReadVariableNames', true);

isize = size(readindata(:,1));
output = zeros(isize);

for i = 1:isize
    a = sprintf( '%f %s', [i, 'in master loop'] ) ;
    disp(a);
    paramss = [readindata(i,:) readinParams];
    paramss.sensAnalFlag = 1;
    % set up gammas, psis, sigmas
    paramss = updateParameters(paramss);

    [~,pops] = ode15s(@stigmaODEmodel,tSpan,ipops,options,paramss); % solve 
                            % the ODEs   NB not guaranteed steady state
    
    output(i) = max(0, 1 - pops(end,1) - pops(end,5));  % prevalence
end %master for loop.

save('lhsoutput.mat', 'output');

end % function runlhs

% =========================================================================

function pops = stigmaODEmodel(t,prev_pops,in)
% The system of ODEs to explore the effect of stigma on the population
% level prevalence of anxiety/depression
%
% Inputs:
%    t = current time. Note this is a DE solver requirement
%    prev_pops = population at previous time point
%    in = various model input parameters
%
% Compartment definitions:
%   U = susceptible/unaffected
%   A = affected/(anxious/depressed)
%   T = undergoing treatment (less infectious)
%   M = managed anxiety/depresssion
%
% Subscripts:
%   n = neutral
%   s = stigma

Un = prev_pops(1);
An = prev_pops(2);
Tn = prev_pops(3);
Mn = prev_pops(4);
Us = prev_pops(5);
As = prev_pops(6);
Ts = prev_pops(7);
Ms = prev_pops(8);

N = Un+An+Tn+Mn + Us+As+Ts+Ms;
if in.sensAnalFlag == 1
    eta = in.eta;
else
    eta = in.etaMultiplier * in.nu; 
end
% "forces of infection" 
lambdan = (in.beta*An + in.beta*As + in.a*in.beta*Tn + in.a*in.beta*Ts)/N;
lambdas = (in.beta*An + in.beta*As + in.a*in.beta*Tn + in.a*in.beta*Ts)/N;

% system of ODEs, not dropped the `/dt` on LHS for simplicity
dUn = (1-in.p)*in.mu*N -lambdan*Un + in.omega*Mn - in.nu*Un - in.mu * Un;
dAn = lambdan*Un - in.gamman*An + in.psin*Tn + in.c*lambdan*Mn + ...
    eta*Mn + in.nu*Un - in.mu*An;
dTn = in.gamman*An - in.psin*Tn - in.sigman*Tn - in.mu*Tn;
dMn = in.sigman*Tn - in.omega*Mn - in.c*lambdan*Mn - eta*Mn - in.mu*Mn;

dUs = in.p*in.mu*N-lambdas*Us + in.omega*Ms - in.nu*Us - in.mu*Us;
dAs = lambdas*Us - in.gammas*As + in.psis*Ts + in.c*lambdas*Ms + ...
    eta*Ms + in.nu*Us - in.mu*As;
dTs = in.gammas*As - in.psis*Ts - in.sigmas*Ts - in.mu*Ts;
dMs = in.sigmas*Ts - in.omega*Ms - in.c*lambdas*Ms - eta*Ms - in.mu*Ms;

pops = [dUn; dAn; dTn; dMn; dUs; dAs; dTs; dMs]; % needed by de solver

end

% =========================================================================

function ipops = initialConditions(prop_stig, alpha_0, total_pop)
% calculate the initial conditions based on key inputs 
%
% Inputs:
%    prop_stig = proportion of population undergoing stigma
%    alpha_0 = target steady state prevalence, used here as initial
%         proportion with "active" anxiety/depression for simplicity
%    total_pop = total population size. For MODSIM paper = 1
%
% Output:
%    ipops = initial condition for the system of ODEs, the "initial
%        populations"

% initial population neutral
An0 = alpha_0 * (1-prop_stig) * total_pop; 
Tn0 = 0;
Mn0 = 0;
Un0 = (1-prop_stig) * total_pop-An0-Tn0-Mn0;

% initial population stigmatised
As0 = alpha_0 * prop_stig * total_pop;
Ts0 = 0;
Ms0 = 0;
Us0 = prop_stig* total_pop-As0-Ts0-Ms0;

ipops = [Un0,An0,Tn0,Mn0,Us0,As0,Ts0,Ms0]; 
end

% =========================================================================

function out = initialiseParameters
% Read in parameter values from file (`inputs.csv`) where:
% row 1 = header
% row 2 = min values
% row 3 = max values
% row 4 = expected values
% row 5 = flag for whether or not to include in sensitivity analysis
% NB since use default header true, indices are reduced by one below

% use expected values from `inputs.csv`
readindata = readtable('inputs.csv', 'ReadVariableNames', true);
out = readindata(3,:);  % expected values are row 3

out.sensAnalFlag = 0;  % for whether to use eta=2*nu or LHS values

% set up gammas, psis, sigmas
out = updateParameters(out);

end

% =========================================================================

function params = updateParameters(params)
% updating relative effect of stigma on treatment parameters 
                    
params.gammas = params.h * params.gamman;
params.sigmas = params.g * params.sigman;
params.psis = params.k * params.psin;

end
