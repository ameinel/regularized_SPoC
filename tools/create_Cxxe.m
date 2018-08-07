function [Cxxe, Cxx] = create_Cxxe(X,varargin)
% create trial-wise and averaged covariance matrix from epoched and
% bandpass-filtered data X
%
% [Cxxe, Cxx] = create_Cxxe(X,'keyword1',value1)
%
% Input:
%   X          - epoched data of size [samples per epoch, channels, epochs]
%
% Optional keyword input:
% 'normalize_Cxxe'  - boolean flag for trace normalization of single-trial 
%                     covariance matrix
% 'normalize_Cxx'   - boolean flag for trace normalization of averaged 
%                     covariance matrix
% 'autoshrink_Cxxe' - boolean flag for automatic shrinkage of covariance
%                     matrix using the Ledoit & Wolf estimator
%
% Output:
%   Cxxe      - mean free trial-wise covariance matrix
%   Cxx       - averaged covariance matrix
%
% andreas.meinel@blbt.uni-freiburg.de, 12/2017

%% loading default options

opt = propertylist2struct(varargin{:});
opt = set_defaults(opt,...
    'normalize_Cxxe',0,...
    'normalize_Cxx',0,...
    'autoshrink_Cxxe',0);

%% Compute covariance matrices

[N_s, N_c,N_e] = size(X);
Cxx = zeros(N_c, N_c);
Cxxe = zeros(N_c, N_c, N_e);

for e=1:N_e
    X_e = squeeze(X(:,:,e));
    if opt.autoshrink_Cxxe 
        X_e=X_e./sqrt(trace(cov(X_e)));
        [C_tmp, ~, ~] = clsutil_shrinkage(X_e','Target','B');
        Cxxe(:,:,e) = C_tmp;
    else
        C_tmp = cov(X_e);
        if opt.normalize_Cxxe
             Cxxe(:,:,e) = C_tmp/trace(C_tmp);
        else
            Cxxe(:,:,e) = C_tmp;
        end       
    end
    Cxx = Cxx + C_tmp;
end
Cxx = Cxx/N_e;

if opt.normalize_Cxx 
   Cxx = Cxx/trace(Cxx);
end
