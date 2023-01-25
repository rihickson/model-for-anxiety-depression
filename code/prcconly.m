function [delta, varnamesGreek] = prcconly
% Modified from Roslyn Hickson by Tim Lynar in 2015 to read a comma seperated file and 
% based on that file perform PRCC.
% Modified by Andrew Rawlinson circa 2017

%params set in `stigma.m` using latin hypercube sampling
%Using dlmread to preserve compatability.
%note assumes no header.

load('lhsoutput.mat') % loads the variable named "output"

dcol = 2; %columns to start at, if it is greater than 1 then we are discarding columns.
drow = 1; %rows to start at, if it is greater than 1 then we are discarding rows.

%end params set in `stigma.m` using latin hypercube sampling

%This is the matrix with all the data we want and will perform a prcc on. 
numberofsamples = length(output);  % for the prcc test
tspan = 0:numberofsamples;           % useful to specify the outputs of the ode caller


%------------------------------
% The Partial Rank Correlation Coefficient (PRCC) sensitivity analysis
%------------------------------
% initialise for efficient memory use etc

Cumi =output; % out of interest. Must have monotonic relationship with params

%params needs to be everything but the output of interest
params = readtable('trilhspoints.txt', 'ReadVariableNames', true);
tparams = table2array(params);
varnames = params.Properties.VariableNames;

% Model parameter values
%Variable names are as defined by the first row in this txt file. 
readindata = readtable('inputs.csv', 'ReadVariableNames', true);
%readindata = readtable1('inputs.txt', 'ReadVariableNames');
inputs = table2array(readindata);
flag = inputs(4,:); % 0 for triangular, 1 for uniform, 2 for remove

% hack to get greek names working for the plotting
varnamesGreek={'$\beta$','$a$','$\sigma_n$','$\gamma_n$','$\psi_n$',...
    '$\omega$','$c$','$\eta$','$\nu$', '$h$', '$g$', '$k$', '$\mu$',...
    '$p$', '$\eta_{\nu}$'};

% remove fixed params from vars
varnamesGreek(flag==2) = [];

% and the prcc:
% B =[input1 input2 ... inputK output]
% where output is the variable of interest
B=[tparams Cumi];
[delta,ttest]=prccgnm(B);   % NB: delta IS used
disp('prcc delta:')
disp(delta)
disp('prcc ttest:')
disp(ttest)

% save results to file. 
string_headers = [strjoin(string(varnames),','), '\n'];
fid = fopen('prcc_results.txt', 'w');
fprintf(fid, strjoin(string_headers,','), '\n');
fclose(fid);
dlmwrite('prcc_results.txt', delta,'-append', 'delimiter', ',', 'precision', 4);
dlmwrite('prcc_results.txt', ttest,'-append', 'delimiter', ',', 'precision', 4);

save('prcc_results.mat','varnamesGreek','delta')

return  


%==========================================================================
%==========================================================================

function x=trirandu(a,c,b,u)
%====================================
% u is a random uniform variable on [0 1]
%====================================

fc=(c-a)./(b-a);
fc(isnan(fc))=0;  % this is necessary for zero diagonals in matrices

if u<fc
    x=a+sqrt(u.*(b-a).*(c-a));
else
    x=b-sqrt((1-u).*(b-a).*(b-c));
end

return

%==========================================================================
%==========================================================================

function  [gamma,t]=prccgnm(D)
%====================================
% Geoff Mercer 1/4/10
% Implementation of the Partial Rank Correlation Coefficient as given in
% S.M. Blower and H. Dowlatabadi,
% Sensitivity and uncertainty analysis of complex models of disease
% transmission: An HIV model as an example
% International Stat. Rev. 62(2) 229-243, 1994
%
% B=[input1 input2 .... inputK output]
% only does one output variable at present
%
% gamma=[PRCCinput1 PRCCinput2 ... PRCCinputK]
% t is the corresponding student T test statistic, I think they had a typo
% in their formula, the denominator should be 1-gamma^2
%====================================

[norows,nocols]=size(D);

K=nocols-1;

% find the ranks of each column, store in R
% doesnt handle equal values well as it just ranks them
% not a problem if data is randomly generated
R=zeros(norows,nocols);
for i=1:nocols
    A=D(:,i);
    [~,order]=sort(A);
    temp=1:norows;
    myrank(order)=temp;
    R(:,i)=myrank;
end

% average rank is mu
mu=(1+norows)/2;
C = zeros(nocols);
tmp = (R-mu).^2;
tmp2 = R-mu;
sumR = sum(tmp);
for i=1:K+1
    newdenom1 = sumR(i);
    for j=1:K+1
        numer=0;
        newdenom2=sumR(j);
        for t=1:norows
            numer=numer+tmp2(t,i)*tmp2(t,j);
        end

        C(i,j)=numer/sqrt(newdenom1*newdenom2);

    end
end

B=inv(C);

gamma = zeros(1,K);
t = zeros(1,K);
for i=1:K
    gamma(i)=-B(i,K+1)/sqrt(B(i,i)*B(K+1,K+1));
    t(i)=gamma(i)*sqrt((norows-2)/(1-gamma(i)^2));
end

return

%==========================================================================
%==========================================================================
