%% Test script: Regularized SPoC vs. standard SPoC on simulation data
% 
% In this script, a regularized variant (you can choose between 'NTik-SPoC',
% 'Tik-SPoC' or 'AS-SPoC') and standard SPoC is both trained on a simulation 
% dataset. Their performances (correlation, angle between spatial filters, 
% spatial pattern) are then reported. 
%
% NOTE: for 'Tik-SPoC' and 'NTik-SPoC' a selection of the regularization 
% parameter 'opts.alpha' is required which can be done e.g. via cross-
% validation. For simplicity, a fixed value is chosen here.
%
% andreas.meinel@blbt.uni-freiburg.de
% December 2017

clear; close all; clc
addpath('tools')
addpath(fullfile('external','bbci_public-master'))
addpath(genpath(fullfile('external','Post-HocLabeling-master')))

startup_bbci_toolbox()

%% load eeg data and generate set of possible target sources (post-hoc labeling)
% here: we use an ICA decomposition for the generating a continous variable
subject_str = 'S2';
options = struct(...
    'ford', 5, ...
    'cFreq', [10] , ...
    'wFreq', [4], ...
    'windowLength', 1000, ...
    'N_compICA', 10,...
    'type', 'ica', ...
    'select_sources', 'all');

eeg_data = load(fullfile('data',subject_str));
[epo_sources, Ax] = load_simulation(eeg_data, options);

%% bandpass-filter EEG data to target frequency bands

[zfilt,pfilt,kfilt] = butter(5,[8,12]/(eeg_data.cnt.fs/2));
cnt_filt = proc_filt(eeg_data.cnt,zfilt,pfilt,kfilt);
epo = proc_segmentation(cnt_filt, eeg_data.vmrk, [0,1000]);
epo = proc_selectEpochs(epo, 'not', eeg_data.iart);

%% Generate target variable, train/test split

close all
rng(111)
% select random source as target and extract labels
ix_targetIndex = randi(size(epo_sources.x,2));
epo_target = proc_selectChannels(epo_sources,ix_targetIndex);
z = squeeze(mean(epo_target.x,1));
epo.y = z';

% split train/test set
Ne = size(epo.x,3);
[ix_train,ix_val,~] = divideblock(Ne,0.1,0.9);
epo_tr = proc_selectEpochs(epo, ix_train, 'RemoveVoidClasses', 0);
epo_val = proc_selectEpochs(epo, ix_val,  'RemoveVoidClasses', 0);

fprintf('# of training data: %i \n',length(ix_train))

%% Run regularized SPoC vs. standard SPoC

% Choose between 'NTik-SPoC','Tik-SPoC' and 'AS-SPoC'
reg_method = 'NTik-SPoC'; 
opts = select_RegularizationMethod(reg_method);

if isfield(opts,'alpha')
    fprintf('Fixed regularization parameter for %s: %.2e \n',reg_method,...
        opts.alpha)
end

% train rregularized/baseline SPoC
[W_reg, A_reg] = reg_spoc(epo_tr.x,epo_tr.y,opts);
w_reg = W_reg(:,1);
a_reg = A_reg(:,1);

[W, A] = spoc(epo_tr.x,epo_tr.y,opts);
w_bl = W(:,1);
a_bl = A(:,1);

% apply SPoC
epo_targetPred_reg = proc_linearDerivation(epo_val, w_reg, 'prependix','spoc');
epo_targetPred_bl = proc_linearDerivation(epo_val, w_bl, 'prependix','spoc');

%% Evaluation of decoding performance

% correlation rho
z_val = z(ix_val);
z_est_reg = squeeze(var(epo_targetPred_reg.x,[],1));
z_est_bl = squeeze(var(epo_targetPred_bl.x,[],1));
rho_reg = corr(z_val,z_est_reg);
rho_bl = corr(z_val,z_est_bl);

fprintf('Correlation for SPoC: %.3f | for %s: %.3f \n',rho_bl,...
    reg_method,rho_reg);

% angle between spatial filters 
[~,ix_order_channels] = ismember(epo.clab,Ax.clab);
val_ch = find(ix_order_channels ~= 0);
ix_order_channels = ix_order_channels(val_ch);

w_true = Ax.W_ica(ix_order_channels,ix_targetIndex);

theta_reg = acos(dot(w_reg(val_ch),w_true)/(norm(w_reg)*norm(w_true)));
if (theta_reg > pi/2)
    theta_reg = pi- theta_reg;
end
theta_bl = acos(dot(w_bl(val_ch),w_true)/(norm(w_bl)*norm(w_true)));
if (theta_bl > pi/2)
    theta_bl = pi- theta_bl;
end

fprintf('Angle theta bwetween filters for SPoC: %.3f | for %s: %.3f \n',...
    theta_bl,reg_method,theta_reg);

% plot original and estimated spatial patterns
mnt = mnt_setElectrodePositions(Ax.clab(ix_order_channels));
a_true = Ax.Ax_all(ix_order_channels,ix_targetIndex);

figure('Position',[0,0,800,300]);
subplot(1,3,1); 
plot_scalp(mnt, a_true, 'ScalePos', 'none');
title('true pattern')

subplot(1,3,2); 
plot_scalp(mnt, a_bl(val_ch), 'ScalePos', 'none');
title('estimated pattern: SPoC')

subplot(1,3,3); 
plot_scalp(mnt, a_reg(val_ch), 'ScalePos', 'none');
title(['estimated pattern: ',reg_method])

