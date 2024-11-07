% 2024-08-09, Test
% 2024-08-15, Add function of memoryPolynomialMatrixGen
% 2024-08-19, Change input of Degree and MemoryDepth to setDegree_cell and setMemoryDepth_cell
% 2024-08-28, Initial k1=1, kEnd=order, m1=0, mEnd=depth
% 2024-08-29, Add cross term settings
% 2024-09-03, Add isSignalGenerator: 'coefficientFinderInverse'
% 2024-09-11, Modify RMS Err.
%% 2024-09-20, add isEVM

function varargout = rfSim_pa_memoryModel(signal_coefs_cell,isSignalGenerator,setDegree_cell,setMemoryDepth_cell,setCrossDelay_cell,setCrossLag_cell,...
    setInputPowerRangeDbm,ReferenceImpedance,fnum)
%
% If isSignalGenerator == 'coefficientFinder',
% Generate polynomial memory pa model using pa's input and output
% measurement data
% Return varargout - cell array containing:
% coef: memory polynomial coefficients
% X_pm_matrix_cell: memory polynomial pa's input matrix
% y_memFit_cell: memory polynomial fitting pa's output
% rmsError: rms error between measurement and fitting data
%
% Elseif isSignalGenerator == 'signalGenerator',
% Generate polynomial memory fitting pa's output using pa input and
% polynominal fitting coefficients,
% Noticed the signal_coefs_cell is CELL type
% Return varargout = cell array containing:
% y_memFit_cell: memory polynomial fitting pa's output
% X_pm_matrix_cell: memory polynomial pa's input matrix
%
% 2024-08-19, setDegree_cell and setMemoryDepth_cell
% - If the type is cell, setDegree_cell = {k1, kStep, kMax_vec} and
%   setMemoryDepth_cell = {m1, mStep, mMax_vec}.
% - If the type is vector, setDegree_cell = [pMax_vec], with k1=1 and kStep=2 as defaults, and
%   setMemoryDepth_cell = [mMax_vec], with m1=0 and mStep=1 as defaults.
%
if iscell(setDegree_cell)
    k1 = setDegree_cell{1};
    kStep = setDegree_cell{2};
    kEnd_vec = setDegree_cell{3};
elseif isvector(setDegree_cell)
    k1 = 1; % 1st order initialization
    kStep = 2;
    kEnd_vec = setDegree_cell;
else
    error('setDegree_cell type ?')
end
if iscell(setMemoryDepth_cell)
    m1 = setMemoryDepth_cell{1};
    mStep = setMemoryDepth_cell{2};
    mEnd_vec = setMemoryDepth_cell{3};
elseif isvector(setMemoryDepth_cell)
    m1 = 0; % zero memory initialization
    mStep = 1;
    mEnd_vec = setMemoryDepth_cell;
else
    error('setMemoryDepth_cell type ?')
end
if ~exist('setCrossDelay_cell','var') || any(isempty(setCrossDelay_cell)) || any(cellfun(@isempty, setCrossDelay_cell))
    setCrossDelay_cell = [];
end
if ~exist('setCrossLag_cell','var') || any(isempty(setCrossLag_cell)) || any(cellfun(@isempty, setCrossLag_cell))
    setCrossLag_cell = [];
end

isFnum = exist('fnum', 'var') && ~isempty(fnum);
if exist('ReferenceImpedance','var') && ~isempty(ReferenceImpedance)
    Rohm = ReferenceImpedance;
else
    Rohm = 1;
end
if exist('setInputPowerRangeDbm','var') && ~isempty(setInputPowerRangeDbm)
    setInputPowerRangeDbm;
    isSetInputPowerRange = 1;
else
    setInputPowerRangeDbm = [];
    isSetInputPowerRange = 0;
end

switch isSignalGenerator
    case {'coefficientFinder', 'coefficientFinderInverse'}
        x = signal_coefs_cell(:,1); % pa input signal
        y = signal_coefs_cell(:,2); % pa output signal
    case 'signalGenerator'
        x = signal_coefs_cell{1};
        pmCoefficients = signal_coefs_cell{2};
        y = nan;
    case 'polynominalMemoryMatrix'
        x = signal_coefs_cell;
        y = nan;
end
signalIn = x; % Get input signal
Nsamps = length(x);
signalOut = y;

switch isSignalGenerator
    case {'coefficientFinder', 'coefficientFinderInverse'}
        if isSetInputPowerRange && ~isempty(setInputPowerRangeDbm)
            piDbm = 10*log10((abs(x)*sqrt(1000)./sqrt(Rohm)).^2);
            try
                piDbm_min = setInputPowerRangeDbm(1);
            catch
                piDbm_min = min(piDbm-eps);
            end
            try
                piDbm_max = setInputPowerRangeDbm(2);
            catch
                piDbm_max = max(piDbm+eps);
            end
            idx_pin = find(piDbm>=piDbm_min & piDbm<=piDbm_max);
            x = x(idx_pin);
            y = y(idx_pin);
        end
end
N = length(x); % updated xLen
Nm = numel(mEnd_vec);
Nk = numel(kEnd_vec);
Nkm = Nm * Nk;

PolynominalMatrixVsRmsError = zeros(Nkm, 4+~isempty(setCrossDelay_cell)+~isempty(setCrossLag_cell)); % [count, P, M, rmsError]
rmsError = zeros(Nm,Nk);
coef = cell(Nkm, 1);
signalIn_pm_matrix_cell = coef;
signalOut_pmFit_cell = coef;

count = 0;
% Loop through polynomial degrees and memory depths
for k=1:numel(kEnd_vec)
    for n=1:numel(mEnd_vec)
        count = count + 1;

        K = kEnd_vec(k);
        M = mEnd_vec(n);
        poly_vec = [k1:kStep:K];
        memory_vec = [m1:mStep:M];

        % Construct polynomial memory matrix
        switch isSignalGenerator
            case 'coefficientFinderInverse'
                [Y_pm_matrix, maxDelay, count_p_m]  = memoryPolynomialMatrixGen(y, poly_vec, memory_vec, setCrossDelay_cell, setCrossLag_cell);
                Y_pm_matrix_cell{count} = Y_pm_matrix;
            otherwise
                [X_pm_matrix, maxDelay, count_p_m]  = memoryPolynomialMatrixGen(x, poly_vec, memory_vec, setCrossDelay_cell, setCrossLag_cell);
                X_pm_matrix_cell{count} = X_pm_matrix;
        end

        if 0
            X_pm_matrix_tmp = X_pm_matrix(M:end,:);
            xrow = reshape((M:-1:1)' + (0:Nsamps:Nsamps*(K-1)),1,[]);
            xVec = (0:Nsamps-M)' + xrow;
            xPow = x.*(abs(x).^(0:K-1));
            xVec = xPow(xVec);
            find(xVec - X_pm_matrix_tmp~=0)
        end

        % Get coefficients using least square
        switch isSignalGenerator
            case 'coefficientFinder'
                isMemoryEffect = 'memory';
                if strcmpi(isMemoryEffect, 'memoryLess')
                    coef = zeros(K*M,1);
                    coef([1,2,3],:) = [c1 0 c3*(3/4)]';
                elseif  strcmpi(isMemoryEffect, 'memory')
                    if 1
                        xVec = X_pm_matrix(1+maxDelay:N,:);
                        yVec = y(1+maxDelay:N);
                        try
                            coef{count} = xVec \ yVec;
                        catch
                            coef{count} = ( xVec'*xVec + 0*(eps)*eye(size(xVec'*xVec)) ) \ (xVec'*yVec);
                        end

                        if 0
                            xVec_tmp = X_pm_matrix;
                            yVec_tmp = y;
                            memDepth = 1+maxDelay;
                            coef_tmp = reshape(xVec_tmp\(yVec_tmp(memDepth:end)),memDepth,[]);
                            coef_tmp2 = reshape(xVec\(yVec(1:end)),1,[]);
                        end
                    else
                        X_tmp = X_pm_matrix(M:Nsamps,:);
                        y_tmp = y(M:Nsamps);
                        X_tmp_H = X_tmp';
                        coef{count} = inv(X_tmp'*X_tmp + eps*eye(35)) * X_tmp'*y_tmp;
                    end
                else
                    error('isMemoryEffect ?')
                end
            case 'coefficientFinderInverse' % get coefficients of inverse between pa input/output (as comm.DPD toolbox)
                xVec = x(1+maxDelay:N);
                yVec = Y_pm_matrix(1+maxDelay:N,:);
                try
                    coef_inverse = yVec\ xVec;
                catch
                    coef_inverse = ( yVec'*yVec + 0*(eps)*eye(size(yVec'*yVec)) ) \ (yVec'*xVec);
                end

                coef{count} = coef_inverse;
            case 'signalGenerator'
                coef{1} = pmCoefficients;
            case 'polynominalMemoryMatrix'
                coef{1} = ones(K*M,1);
        end

        if 0
            xVec_tmp    = X_pm_matrix(M:Nsamps,:);
            Coef_tmp    = reshape(xVec_tmp \ y(M:Nsamps),M,[]);
            Coef        = reshape(xVec\(y(M:Nsamps)),M,[]);
        end

        % if strcmpi(isSignalGenerator, 'coefficientFinder')
        switch isSignalGenerator
            case {'coefficientFinder','signalGenerator','coefficientFinderInverse'}
                % Construct polynomial memory matrix
                signalIn_pm_matrix = memoryPolynomialMatrixGen(signalIn, poly_vec, memory_vec, setCrossDelay_cell, setCrossLag_cell);

                % Generate memory model fitting output
                signalOut_pmFit = zeros(size(signalIn));
                signalOut_pmFit(1+maxDelay:Nsamps,:) = signalIn_pm_matrix(1+maxDelay:Nsamps,:) * coef{count};

                if strcmpi(isSignalGenerator,'coefficientFinder')
                    % Calculate rms error
                    isRemoveDelay = 1;
                    try
                        isEVM;
                    catch
                        isEVM = 'evm'; % 2024-09-20, add isEVM
                    end
                    if 1
                        % 2024-09-11, Modify RMS Err.
                        switch isEVM
                            case 'rmsErr'
                                signalOut_norm = signalOut(isRemoveDelay*maxDelay+1:end)/rms(signalOut(isRemoveDelay*maxDelay+1:end));
                                signalOut_pmFit_norm = signalOut_pmFit(isRemoveDelay*maxDelay+1:end)/rms(signalOut_pmFit(isRemoveDelay*maxDelay+1:end));
                                errNorm = (signalOut_norm - signalOut_pmFit_norm);
                                rmsError(n,k) = sqrt(sum(abs(errNorm).^2)) / length(errNorm) * 100;
                            case 'evm'
                                signalOut_nonNorm = signalOut(isRemoveDelay*maxDelay+1:end)/1;
                                signalOut_pmFit_nonNorm = signalOut_pmFit(isRemoveDelay*maxDelay+1:end)/1;
                                rmsError(n,k) = evm(signalOut_nonNorm, signalOut_pmFit_nonNorm, '%');
                        end
                        if 0 % debug
                            figure(0911)
                            plot(real(signalOut_norm)), hold on
                            plot(real(signalOut_pmFit_norm), 'o','MarkerSize',5)
                        end
                        % % else
                        % %     switch isEVM
                        % %         case 'rmsErr'
                        % %             rmsError(n,k) = rms(err)/rms(signalOut(isRemoveDelay*maxDelay+1:end))*100;
                        % %         case 'evm'
                        % %             rmsError(n,k) = evm(signalOut_norm, signalOut_pmFit_norm, '%');
                        % %     end

                    end

                    if size(PolynominalMatrixVsRmsError,2) == 6
                        PolynominalMatrixVsRmsError(count,:) = [count, K, M, max(setCrossDelay_cell{end}), max(setCrossLag_cell{3}), rmsError(n,k)];
                    else
                        PolynominalMatrixVsRmsError(count,:) = [count, K, M, rmsError(n,k)];
                    end
                end

                % Return
                signalIn_pm_matrix_cell{count} = signalIn_pm_matrix;
                signalOut_pmFit_cell{count} = signalOut_pmFit;
            otherwise
                % Return
                signalIn_pm_matrix_cell{count} = X_pm_matrix_cell;
        end

    end % n
end % k

% Output results
switch isSignalGenerator
    case 'coefficientFinder'
        if isscalar(coef)
            varargout{1} = coef{1};
            varargout{2} = signalIn_pm_matrix_cell{1};
            varargout{3} = signalOut_pmFit_cell{1};
            varargout{4} = rmsError;
            varargout{5} = count_p_m;
        else
            varargout{1} = coef;
            varargout{2} = signalIn_pm_matrix_cell;
            varargout{3} = signalOut_pmFit_cell;
            varargout{4} = rmsError;
            varargout{5} = count_p_m;
        end
    case 'coefficientFinderInverse'
        varargout{1} = coef{1};
        varargout{2} = signalOut_pmFit_cell{1}; % pa input multiple with inverse coefficients for PA input dpd
        fnum = [];
    case 'signalGenerator'
        varargout{1} = signalOut_pmFit;
        varargout{2} = signalIn_pm_matrix_cell{1};
        varargout{3} = count_p_m;
    case 'polynominalMemoryMatrix'
        varargout{1} = signalIn_pm_matrix_cell{1};
        varargout{2} = count_p_m;
end

if isFnum
    figure(fnum)
    if all(size(rmsError)>1)
        contourf(rmsError,'ShowText','on'), colorbar, grid minor
        xlabel(sprintf('Polynomial Degree (Fundamental is 1)'))
        ylabel(sprintf('Memory Depth (Without memory is 0)'))
        title([upper(isEVM), ' [%]'])
        set(gca, 'XTick', 1:numel(kEnd_vec), 'XTickLabel', kEnd_vec);
        set(gca, 'YTick', 1:numel(mEnd_vec), 'YTickLabel', mEnd_vec);
    elseif size(rmsError,2)>1
        scatter(1:numel(kEnd_vec), rmsError); grid minor
        xlabel(sprintf('Polynomial Degree (Fundamental is 1)'))
        ylabel([upper(isEVM), ' [%]'])
        title([upper(isEVM), ' [%]'])
        set(gca, 'XTick', 1:numel(kEnd_vec), 'XTickLabel', kEnd_vec);
    else
        scatter(1:numel(kEnd_vec), rmsError); grid minor
        xlabel(sprintf('Memory Depth (Without memory is 0)'))
        ylabel([upper(isEVM), ' [%]'])
        title([upper(isEVM), ' [%]'])
        set(gca, 'XTick', 1:numel(mEnd_vec), 'XTickLabel', mEnd_vec);
    end
end

end

function [X_km_matrix, maxDelay, count_p_m] = memoryPolynomialMatrixGen(x, degree_vec, memoryDepth_vec, setCrossDelay_cell, setCrossLag_cell)
% [1]. Behavioral Modelling and Predistortion of Wideband Wireless Transmitters
% [2]. https://ieeexplore.ieee.org/abstract/document/1703853

% Ma, Ka - memory depth and polynomial order for first term
% Mb, Kb, P - parameters for second term
% Mc, Kc, Q - parameters for third term

% Initialize
x = x(:); % row
if 1
    Ka = degree_vec;
    Ma = memoryDepth_vec;

    if Ka(1)==0
        error('Fundermental start from polynomial order ONE')
    end
    if Ma(1)~=0
        error('Fundermental start from memory depth ZERO')
    end
end

% 2024-08-29, Add cross-delay term settings
if exist("setCrossDelay_cell","var")&&~isempty(setCrossDelay_cell)
    Kb = setCrossDelay_cell{1}; % polynomial order vector
    Mb = setCrossDelay_cell{2}; % memory depth vector
    Lb = setCrossDelay_cell{3}; % corss delay vector
    if Kb(1)<2
        error('Cross delay term start from polynomial order TWO')
    end
    if Mb(1)<0
        error('Cross delay term start from memory depth ZERO')
    end
    if Lb(1)<1
        error('Cross delay term start from corss-delay ONE')
    end
else
    Kb = nan;
    Mb = nan;
    Lb = nan;
end
% 2024-08-29, Add cross-lag term settings
if exist("setCrossLag_cell","var")&&~isempty(setCrossLag_cell)
    Kc = setCrossLag_cell{1}; % polynomial order vector
    Mc = setCrossLag_cell{2}; % memory depth vector
    Lc = setCrossLag_cell{3}; % corss lag vector
    if Kc(1)<2
        error('Cross lag term start from polynomial order TWO')
    end
    if Mc(1)<0
        error('Cross lag term start from memory depth ZERO')
    end
    if Lc(1)<1
        error('Cross lag term start from corss-delay ONE')
    end
else
    Kc = nan;
    Mc = nan;
    Lc = nan;
end

% Construct polynomial memory matrix
N = length(x);
NKa = numel(Ka);
NMa = numel(Ma);
NKb = ~isnan(Kb(1))*numel(Kb);
NMb = ~isnan(Mb(1))*numel(Mb);
NLb = ~isnan(Lb(1))*numel(Lb);
NKc = ~isnan(Kc(1))*numel(Kc);
NMc = ~isnan(Kc(1))*numel(Mc);
NLc = ~isnan(Kc(1))*numel(Lc);

X_km_matrix = zeros(N,NKa*NMa+NKb*NMb*NLb+NKc*NMc*NLc);
count_pm = 0;
count_p_m = [];

if 1
    % tic
    p = 0; % cross term delay
    q = 0; % cross term lag
    for k=degree_vec
        for m=memoryDepth_vec
            if m<N
                count_pm = count_pm + 1;
                delayed_x_m = [zeros(m,1); x(1:N-m)]; % delay the x by m
                poly_x_m = abs(delayed_x_m).^(k-1); % polynomial term by m
                X_km_matrix(:,count_pm) = delayed_x_m.*poly_x_m; % gmp matrix
                count_p_m = [count_p_m; [count_pm k m p q nan]];
            end
        end
    end
    % toc
else
    tic
    p = 0; % cross term delay
    q = 0; % cross term lag

    % memoryDepth_vec
    valid_m = memoryDepth_vec(memoryDepth_vec < N);

    % (k, m)combination
    [K_grid, M_grid] = ndgrid(degree_vec, valid_m);
    K_list = K_grid(:);
    M_list = M_grid(:);
    total_terms = length(K_list);
    delayed_x_m_cell = arrayfun(@(m) [zeros(m,1); x(1:N-m)], M_list, 'UniformOutput', false);
    poly_x_m_cell = arrayfun(@(idx) abs(delayed_x_m_cell{idx}).^(K_list(idx)-1), 1:total_terms, 'UniformOutput', false);
    X_terms = cellfun(@(dx, px) dx .* px, delayed_x_m_cell, poly_x_m_cell(:), 'UniformOutput', false);
    X_km_matrix(:, 1:total_terms) = cell2mat(X_terms');
    toc
end

if ~isnan(Lb(1))
    q = nan;
    % cross term delays
    for m=Mb
        for k=Kb
            for p=Lb
                if m<N
                    delayed_x_m = [zeros(m,1); x(1:N-m)]; % delay the x by m
                end
                if (m+p)<N
                    count_pm = count_pm + 1;
                    delayed_x_mp = [zeros(m+p,1); x(1:N-(m+p))]; % delay the x by m+p
                    poly_x_mp = abs(delayed_x_mp).^(k-1); % polynomial term by m and p
                    X_km_matrix(:,count_pm) = delayed_x_m.*poly_x_mp; % gmp matrix
                    count_p_m = [count_p_m; [count_pm k m p q (m+p)]];
                    % if count_pm==132 % debug
                    %     debug = 1;
                    % end
                end
            end
        end
    end
end

if ~isnan(Lc(1))
    p = nan;
    % cross term lags
    for m=Mc
        for k=Kc
            for q=Lc
                if m<N
                    delayed_x_m = [zeros(m,1); x(1:N-m)]; % delay the x by m
                end
                if (m-q)>=0 && (m-q)<N
                    count_pm = count_pm + 1;
                    delayed_x_mq = [zeros(m-q,1); x(1:N-(m-q))]; % delay the x by m+q
                    poly_x_mq = abs(delayed_x_mq).^(k-1); % polynomial term by m and p
                    X_km_matrix(:,count_pm) = delayed_x_m.*poly_x_mq; % gmp matrix
                    count_p_m = [count_p_m; [count_pm k m p q (m-q)]];
                end
            end
        end
    end
end

% remove empty term and return
X_km_matrix = X_km_matrix(:,1:count_pm);
maxDelay = max([count_p_m(:,3);count_p_m(:,end)]);
end

function [X_km_matrix, count_p_m] = memoryPolynomialMatrixGen_bk(x, setDegree, setMemoryDepth, setCrossDelay_cell, setCrossLag_cell)
% [1]. Behavioral Modelling and Predistortion of Wideband Wireless Transmitters
% [2]. https://ieeexplore.ieee.org/abstract/document/1703853

% Ma, Ka - memory depth and polynomial order for first term
% Mb, Kb, P - parameters for second term
% Mc, Kc, Q - parameters for third term

% Initialize
x = x(:); % row
try
    k1 = setDegree(1);
    kStep = setDegree(2);
    pEnd = setDegree(3);
catch
    k1 = 0;
    kStep = 1;
    pEnd = setDegree;
end
try
    m1 = setMemoryDepth(1);
    mStep = setMemoryDepth(2);
    mEnd = setMemoryDepth(3);
catch
    m1 = 0;
    mStep = 1;
    mEnd = setMemoryDepth;
end
% 2024-08-29, Add cross term settings
if exist("setCrossDelay_cell","var")&&~isempty(setCrossDelay_cell)
    switch class(setCrossDelay_cell)
        case 'cell'
            Kb = setCrossDelay_cell{1}; % polynomial order vector
            Mb = setCrossDelay_cell{2}; % memory depth vector
            Lb = setCrossDelay_cell{3}; % corss delay vector
        case 'double'
            Kb = k1:kStep:pEnd;
            Mb = m1:mStep:mEnd;
            if isscalar(setCrossDelay_cell)
                Lb = 2:1:setCrossLag_cell; % setCrossDelay_cell is max. of cross delay
            else
                Lb = setCrossLag_cell; % setCrossDelay_cell is vector of cross delay
            end
    end
else
    Lb = 0;
end
if exist("setCrossLag_cell","var")&&~isempty(setCrossLag_cell)
    switch class(setCrossLag_cell)
        case 'cell'
            Kc = setCrossLag_cell{1}; % polynomial order vector
            Mc = setCrossLag_cell{2}; % memory depth vector
            Lc = setCrossLag_cell{3}; % corss lagging vector
        case 'double'
            Kc = k1:kStep:pEnd;
            Mc = m1:mStep:mEnd;
            if isscalar(setCrossLag_cell)
                Lc = 2:1:setCrossLag_cell; % setCrossLag_cell is max. of cross lag
            else
                Lc = setCrossLag_cell; % setCrossLag_cell is vector of cross lag
            end
    end
else
    Lc = 0;
end

% Construct polynomial memory matrix
N = length(x);
X_km_matrix = zeros(N,numel(k1:kStep:pEnd-1)*numel(m1:mStep:mEnd-1));
count_pm = 0;
count_p_m = [];
if k1~=1 % 2024-08-28, Initial p1=1, pEnd=order, m1=0, mEnd=depth
    disp('polynominal order should start from 1 and end to maximum order-1')
end
if m1~=0
    disp('memory depth should start from 0 and end to maximum depth')
end

if 1
    p = 0; % cross term delay
    q = 0; % cross term lag
    for k=k1:kStep:pEnd
        for m=m1:mStep:mEnd
            count_pm = count_pm + 1;
            x_m = zeros(N,1);
            if m<N
                x_m(m+1:end) = x(1:N-m);
            end
            x_k = sqrt(x_m.*conj(x_m)).^(k-1);
            X_km_matrix(:,count_pm) = x_m.*x_k;
            count_p_m = [count_p_m; [count_pm k m p q]];
        end
    end
end

if Lb~=0
    % cross term delays
    for m=Mb
        for k=Kb
            for p=Lb
                count_pm = count_pm + 1;
                if m<N
                    delayed_x_m = [zeros(1,m); x(1:N-m)]; % delay the x by m
                end
                if (m+p)<N
                    delayed_x_mp = [zeros(1, m+p); x(1:N-m-p)]; % delay the x by m+p
                    poly_x_mp = abs(delayed_x_mp).^(k-1); % polynomial term by m and p
                end
                X_km_matrix(:,count_pm) = delayed_x_m.*poly_x_mp;
            end
        end
    end
end

if Lc~=0
    % cross term lags
    for m=Mc
        for k=Kc
            for q=Lc
                count_pm = count_pm + 1;
                if m<N
                    delayed_x_m = [zeros(1,m); x(1:N-m)]; % delay the x by m
                end
                if (m-q)>=0 && (m-q)<N
                    delayed_x_mq = [zeros(1, m-p); x(1:N-m-q)]; % delay the x by m+q
                    poly_x_mq = abs(delayed_x_mq).^(k-1); % polynomial term by m and p
                end
                X_km_matrix(:,count_pm) = delayed_x_m.*poly_x_mq;
            end
        end
    end
end

% return
X_km_matrix;
end