function [opts] = select_RegularizationMethod(method)
% sets the configurations for the different SPoC variants (regularized or
% non-regularized versions)
%
% [opts] = select_RegularizationMethod(method)
% 
% Input:
% method - string describing the selected method, the following choices are
%          possible: 'Tik-SPoC','NTik-SPoC', 'AS-SPoC', 'SPoC'
%
% Output: 
% opts   - struct containing all optional arguments for the "reg_spoc()" 
%          function. Please note that 'Tik-SPoC' and 'NTik-SPoC' require a
%          selection of the regularization parameter alpha, e.g. via 
%          cross-validation.

opts=struct();

try
    switch method
        
        case 'Tik-SPoC'
            % Tikhonov Regularization without any trace normalization
            opts.normalize_Cxx = 0;
            opts.normalize_Cxxe = 0;
            opts.normalize_Cxxz = 0;
            opts.alpha = 5e-3;
            
        case 'NTik-SPoC'
            % Tikhonov Regularization with trace normalization (to Cxx and
            % Cxxe)
            opts.normalize_Cxx = 1;
            opts.normalize_Cxxe = 1;
            opts.alpha = 1e-4;
            
        case 'AS-SPoC'
            % automatic shrinkage for Cxxe (analytic solution for optimal
            % regularization parameter)
            opts.normalize_Cxx = 0;
            opts.normalize_Cxxe = 0;
            opts.autoshrink_Cxxe = 1;
            
        case 'SPoC'
            % standard SPoC algorithm without regularization/normalization
            opts.normalize_Cxx = 0;
            opts.normalize_Cxxe = 0;
            opts.autoshrink_Cxxe = 0;
            opts.alpha = 0;
         
    end
    
catch
    error('Input struct contains an invalid method')
end

end