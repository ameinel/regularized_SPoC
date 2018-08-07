function [W, A, lambda_values, p_values_lambda, Cxx, Cxxz, Cxxe] = reg_spoc(X, z, varargin)
% Regularized version of Source Power Co-modulation Analysis (SPoC).
%
% Optimizes spatial filters such that the power of the filtered
% signal maximally covaries with the univariate target function z, as
% described in Dahne et al., 2014a.
% 
% [W, A, lambda_values, p_values_lambda, Cxx, Cxxz, Cxxe] = reg_spoc(X, z, <'keyword1', value1, 'keyword2', value2, ...>)
%
%
% Input:
% X     - epoched data, with size(X) = [n_samples_per_epoch, n_channels, n_epochs]
%           X should already have been band-pass filtered!
% z     - univariate target function, with length(z) == n_epochs
%
% Optional keyword input:
% 'n_bootstrapping_iterations' - The number of repetitions to bootstrap the
%                           distribution of lambda values with shuffled z.
%                           Default value is 0, i.e. no bootstrapping and
%                           therefore no p-values. 
% 'pca_X_var_explained' - Optional dimensionality reduction with PCA, using
%                       the number of PCA components is determined by how
%                       much cumulative variance they should explain. This
%                       number must be between 0 and 1. Default is 1, i.e.
%                       no dimensionality reduction. 
% 'verbose' - If this value is larger than 0, some info will be printed
%               during bootstrapping. Default value is 1.
% 'alpha' - Regularization parameter, if alpha > 0, then the Tikhonov
%               regularization is applied
% 'autoshrink_Cxxe' - Flag, that applies automatic covariance shrinkage 
%               (Ledoit & Wolf, 2014) to the trial-wise covariance matrices 
%
% Output:
%   W       - filter matrix, each column is a spatial filter, corresponding
%               the lambda values.
%   A       - pattern matrix, each column is a spatial pattern, corresponding
%               the lambda values.
%   lambda  - Lambda values indicate the covariance between the
%   Cxx     - data covariance matrix
%   Cxxz    - data covariance matrix weighted by z
%   Cxxe    - trialwise covariance matrix computed from X
%   lambda_values  - vector with lambda(i) = covariance between the power of the
%                   data after projection on W(:,i) and the target function
%   p_values_lambda - vector with p_values for each lambda. Only meaningfull if
%                   the number of bootstrapping iterations has been set.
%
%
% Example function calls:
%
% Calling the function without arguments will evoke a small example
% >> spoc()
% 
% The most simple way:
% >> W = spoc(X,z); 
%
% Bootstrap the lambda distribution to get some p-values:
% >> [W,A,lambdas,p_values] = spoc(X,z, 'n_bootstrapping_iterations',500); 
%
%
% Reference:
% A. Meinel, S. Castano, B. Blankertz, F. Lotte, M. Tangermann, "Characterizing 
% Regularization Techniques for Spatial Filter Optimization in Oscillatory 
% EEG Regression Problems",
% T.B.D: Journal + Vol etc.
%
% S. Daehne, F. C. Meinecke, S. Haufe, J. Hoehne, M. Tangermann, K. R. Mueller, 
% V. V. Nikulin, "SPoC: a novel framework for relating the amplitude of neuronal
% oscillations to behaviorally relevant parameters",
% NeuroImage, 86(0):111-122, 2014
%
%
% original function by: sven.daehne@tu-berlin.de, 2014
% 12/2017: included regularization techniques, andreas.meinel@blbt.uni-freiburg.de

%% check input

if not(exist('X','var')) && not(exist('z','var'))
    spoc_example;
    return
end

opt= propertylist2struct(varargin{:});
opt= set_defaults(opt, ...
                  'n_bootstrapping_iterations',0,...
                  'pca_X_var_explained', 1, ...
                  'verbose',1, ...
                  'Cxx', [], ...
                  'Cxxz', [], ...
                  'Cxxe', [], ...
                  'autoshrink_Cxxe', 0,...
                  'alpha', 0);
              
if not(isempty(X))
    [N_sam, N_c , N_e] = size(X);
else
    [N_sam, N_c ,N_e] = size(opt.Cxxe);
end
if not(length(z) == N_e)
    error('X and z must have the same number of epochs!!!');
end

%% some preprocessing

% normalize the target funtion to have zero mean and unit-variance and
% average over all epochs
z = (z-mean(z(:)))./std(z(:));
   
if isempty(opt.Cxxe) 
    % create mean-free trial-wise covariance matrices
    [Cxxe, Cxx] = create_Cxxe(X,opt);
else
    Cxxe = opt.Cxxe;
end

if isempty(opt.Cxxz)
    % create the z-weighted covariance matrix
    Cxxz = create_Cxxz(Cxxe, z);
else
    Cxxz = opt.Cxxz;
end

% add penalty term to the denominator
if opt.alpha > 0 && ~opt.autoshrink_Cxxe
    Cxx_reg = eye(N_c,N_c);
    Cxx = (1-opt.alpha)*Cxx + opt.alpha*Cxx_reg;
end

% whiten the data, and possibly reduce dimensionality
[foo, M] = whiten_data(X, opt.pca_X_var_explained, Cxx);
Cxxz_white = M * Cxxz * M';

%% SPoC
% compute SPoC in whitened space. Here the covariance matrix is the
% identity and thus the generalized eigenvalue problem is reduced to an
% ordinary eigenvalue problem
[W, D] = eig(Cxxz_white);
[lambda_values, sorted_idx] = sort(diag(D), 'descend');
W = W(:, sorted_idx);
W = M'*W; % project back to original (un-whitened) channel space

%% some postprocessing

% Normalize the filters such that the extracted components have unit
% variance. This is necessary because the scaling of the eigenvectors is
% arbitrary. Every multiple of an eigenvector fullfils the eigenvector
% equation with the same eigenvalue.
% Cxz*w = lambda*Cxx*w   <=>  Cxz*(c*w) = lambda*Cxx*(c*w), with c in R 
% Thus, the eigenvalues are unique but the scaling of the eigenvectors is not.
for k=1:size(W,2)
    W(:,k) = W(:,k) / sqrt(squeeze(W(:,k)'*Cxx*W(:,k)));
end

A = Cxx * W / (W'* Cxx * W); % compute patterns


%% bootstrap the eigenvalue distribution
n_bootstrapping_iterations = opt.n_bootstrapping_iterations;
z_amps = [];
n_ev = length(lambda_values);
p_values_lambda = inf(1,n_ev);
p_values_r = inf(1,n_ev);
if n_bootstrapping_iterations > 0
    lambda_samples = zeros(1, n_bootstrapping_iterations);
    r_samples = zeros(1, n_bootstrapping_iterations);
    for k=1:n_bootstrapping_iterations
        if opt.verbose && mod(k, ceil(n_bootstrapping_iterations/25)) == 0
            fprintf('bootstrapping iteration %d/%d\n', k, n_bootstrapping_iterations);
        end
        % shuffle the target function
%         z_shuffled = z(randperm(length(z)));
        [z_shuffled, z_amps] = random_phase_surrogate(z, 'z_amps', z_amps);
        % re-compute SPoC
        Cxxz_s = create_Cxxz(Cxxe, z_shuffled);
        Cxxz_s_white = M * Cxxz_s * M';
        [W_s, D_s] = eig(Cxxz_s_white);
        [lambda_values_s, sorted_idx] = sort(diag(D_s), 'descend');
        W_s = W_s(:, sorted_idx);
        W_s = M'*W_s; 
        
        % compute and store lambda and correlation values
        lambda_samples(k) = max(abs(lambda_values_s));
        fv = get_var_features(W_s, Cxxe);
        R = corrcoef([z_shuffled', fv']);
        r_samples(k) = max(abs(R(1,2:end)));
    end
    
    % compute bootstrapped p-values
    for n=1:n_ev
        p_values_lambda(n) = sum(abs(lambda_samples(:))>=abs(lambda_values(n)))/n_bootstrapping_iterations;
    end
end