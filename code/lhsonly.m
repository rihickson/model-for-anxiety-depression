function lhsonly(varargin)
%==========================================================================
% lhsonly is  Latin Hypercube Sampling. 
% 
% functions within this file to note:
%   * lhsu.m        = Latin Hypercube Sampling function, by Geoff N. Mercer
%
% inputs: optionally can specify the number of samples, otherwise defaults
% to 1e3
% 
% Also you should check out "Modeling Infectious Disease in humans and
% animals" by Keeling & Rohani. Especially the associated website with
% MATLAB code.
%
%  
% TO DO AN LHS CHANGE:
% numberofsamples - the number of samples to make.
% 
% Need 
% specified parameters in `inputs.csv`
% xmin - the minimum value for each variable
% Xmax - the maximum value for each variable
% Xexpected - the expected value for each variable
%
% Author: Roslyn Hickson
% Date: 07 Feb 2011
%
% modified by Tim Lynar : 19 Nov 2015.
% modified by Andrew Rawlinson: circa 2017
% modified by Roslyn Hickson: Januaray 2023
%==========================================================================


%------------------------------
% Setting up the sensitivity analysis: the Latin Hypercube Sampling
%------------------------------
if nargin == 0
    numberofsamples = 1e3;          % for the prcc test (500 samples = 5e2)
else
    numberofsamples = varargin{1};
end


% Model parameter values
%Variable names are as defined by the first row in this txt file. 
readindata = readtable('inputs.csv', 'ReadVariableNames', true);
%readindata = readtable1('inputs.txt', 'ReadVariableNames');
varnames = readindata.Properties.VariableNames;
inputs = table2array(readindata);
xmin = inputs(1,:);
xmax = inputs(2,:);
xexpected = inputs(3,:);
flag = inputs(4,:); % 0 for triangular, 1 for uniform, 2 for remove

% store fixed params for runlhs 
fixedParamNames = varnames(flag==2);
fixedVarheader = [strjoin(fixedParamNames,',') '\n'];
fixedParamVals = xexpected(flag==2);
writeToFile('fixedParams.txt', fixedVarheader, fixedParamVals);

% remove fixed params from vars
xmin(flag==2) = [];
xmax(flag==2) = [];
xexpected(flag==2) = [];
varnames(flag==2) = [];
flag(flag==2) = [];

% create header
varheader = [strjoin(varnames,',') '\n'];

% ruv is random uniform variables on [0 1] in a latin cube
ruv = lhsu(zeros(size(xmin)),ones(size(xmax)),numberofsamples);

vfromdisp = zeros(size(ruv)); % sample variable from specified (by flag) distribution

%E.g. Tri distro is ruv(sample, variable)
for sample = 1:numberofsamples
    for var = 1:length(xmin)
        if flag(var) ==0  % triangular
            vfromdisp(sample,var) = trirandu(xmin(var), xexpected(var), xmax(var),ruv(sample,var));
        elseif flag(var)==1  % uniform
            vfromdisp(sample,var) = randu(xmin(var), xmax(var),ruv(sample,var));
        end
    end
end

writeToFile('trilhspoints.txt', varheader, vfromdisp)

disp('All done output written to trilhspoints.txt')


return  % this is not strictly necessary, but it feels wrong not to

function writeToFile(fileName, header, values)

%Write the variable names as the header of the file. 
fid = fopen(fileName, 'w');
fprintf(fid, header, '\n');
fclose(fid);

%Now write the data to the file. 
dlmwrite(fileName, values,'-append', 'delimiter', ',', 'precision', 4);

return




function s=lhsu(xmin,xmax,nsample)
%====================================
% s=lhsu(xmin,xmax,nsample)
% LHS from uniform distribution
% Input:
%   xmin    : min of data (1,nvar)
%   xmax    : max of data (1,nvar)
%   nsample : no. of samples
% Output:
%   s       : random sample (nsample,nvar)
%   Budiman (2003)
%====================================

nvar=length(xmin);
ran=rand(nsample,nvar);
s=zeros(nsample,nvar);
for j=1: nvar
    idx=randperm(nsample);
    P =(idx'-ran(:,j))/nsample;
    s(:,j) = xmin(j) + P.* (xmax(j)-xmin(j));
end

return

%==========================================================================
%==========================================================================

%==========================================================================
%==========================================================================

function x=trirandu(a,c,b,u)
%====================================
%MIGHT BE FROM http://au.mathworks.com/matlabcentral/answers/251154-how-to-generate-a-traingular-distrbution-sawtooth-model
% u is a random uniform variable on [0 1]
%====================================
%a=lower bound, c=expected value, b=upper bound, u=random selection within

fc=(c-a)./(b-a);
fc(isnan(fc))=0;  % this is necessary for zero diagonals in matrices

if u<fc
    x=a+sqrt(u.*(b-a).*(c-a));
else
    x=b-sqrt((1-u).*(b-a).*(b-c));
end

return

function x=randu(a,b,u)
% u random uniform variable on [0 1], needs to be on [a,b]
    x=u.*(b-a)+a;
return
%==========================================================================
%==========================================================================

