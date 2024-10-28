%% PA Modeling and DPD Implementation 
% This example shows to implement the PA model and apply DPD methology for correction.
% History
%% 
% * 2024-09-19, Draft
% * 2024-09-26, Release V1.0, Apply PA model fitting coefs. >> inverse Matrix 
% DPD >> CFR
% * 2024-10-18, Release V2.0, Apply PA modele flatness by isPaFlatness
% * 2024-10-28, Release V3.0
% Set Test Waveform

warning off

try
    Rohm;
    RbwKHz;
catch
    Rohm = 1 % set resistor in ohm
    RbwKHz = 10 % set resolution bw for aclr plot
end
testSignal = "WIFI"; % select signal

switch upper(testSignal) % open example
    case "WIFI"
        if 0
            open('Example_WIFI_Transmit_240919_v1.mlx')
            paInSig = paInSignal;
            paOutSig = paSignal;
            paOutSig_ideal = paSignal;
            wvParams = passParams;
            aclrVal_wv;
            pltComm = PLOT_COMM('aclr', aclrVal_wv, [], RbwKHz, Rohm); % set plt object
            dlSlots;
            fsMHz;
            signalFormat;
            MCS_index
        else
            clear paInSig paOutSig paOutSig_ideal
            load('WIFI_EHT320_MCS13_15.80dBm_NumPkts4_1280MHz_0527SIMCW_0dB.mat') % load waveform
            wvParams;
            paInSig = wvParams.PaInSignal;
            try
                paOutSig = wvParams.PaOutSignal;
            catch
                paOutSig = wvParams.PaOutSignal_ideal;
            end
            try
                paOutSig_ideal = wvParams.PaOutSignal_ideal;
            catch
                paOutSig_ideal = paOutSig;
            end
            try
                PaModelFlatnessPerCBW = wvParams.PaModelFlatnessPerCBW;
            catch
                PaModelFlatnessPerCBW = 0;
            end
            try
                isPaModelFlatCompensateLoss = wvParams.PaModelFlatCompensateLoss;
            catch
                isPaModelFlatCompensateLoss = 0;
            end
            if PaModelFlatnessPerCBW~=0
                isPaFlatness = 1;
                b_coefs_flat = wvParams.PaModelFlatFirCoefs;
            else
                isPaFlatness = 0;
                b_coefs_flat = 1;
            end
            aclrVal_wv = wvParams.PlotACLR;
            pltComm = wvParams.PlotComm;
            dlSlots = wvParams.usefulWvfm;
            fsMHz = wvParams.SampleRateMHz;
            signalFormat = wvParams.SignalFormat;
            MCS_index = wvParams.MCS;
        end

    case "CW"

end

% check power
paOutPwrDbm = powerDbm(paOutSig,'rms',[],'dBm',dlSlots,Rohm)
paInPwrDbm = powerDbm(paInSig,'rms',[],'dBm',dlSlots,Rohm)

% Plot Test Waveform

% plot
aclrVal_paIn = pltComm.ACLR(paInSig, {091901, 'paInput'}); % plot aclr
aclrVal_paOut = pltComm.ACLR(paOutSig, {091901, 'paOutput'}); % plot aclr
pltComm.TIME(paOutSig, {092001, 'paOutout'}, [], 'usec'); % plot time
pltComm.AMAM([paInSig, paOutSig], {092002, 'paOutput'}, 0.01); % plot amam

if 0
    wlanDemodulation(paOutSig, wvParams, 'tx')
end
% 
% PA Characterization - MemoryLess Polynomial Model
% Create coefficients of MemoryLess Polynomial Model using PA input and ouput 
% data.

inputPowerBackOffDb = 40; % set the minimum power range using peak power back-off
[~, paMemLess_obj, paOutput_fitMemLess] = rfSim_pa_memoryLessLUTModel([paInSig,paOutSig],...
    inputPowerBackOffDb,Rohm,070403);
% PA Characterization - Memory Polynomial Model
% Create coefficients of Memory Polynomial Model using PA input and ouput data.
% Sweep Polynomial Order and Memory Depth
% To search the minimum rms error in %

setInputPowerRangeDbm = -70; % set the mimimun power of pa input instantaneous signal
isSweepMemoryPolynomialFit = 0;
nSamps = length(paInSig); % signal length

% set signal captured length to compute the fitting coefficients
if 0
    nCapturedSamps = nSamps/2;
else
    nCapturedSamps = nSamps;
end

% set fitting model using the polynminal order and memory depth
if 1 % set sweeping range
    MemoryDepth_max = 7;
    Degree_max = 11;
    Degree_step = 1;
    Degree_start = 1;
end

MemDepth_sweep = [0:MemoryDepth_max]; % start from 0
Degree_sweep = {Degree_start,Degree_step,[Degree_start:2:Degree_max]}; % start from 0

if isSweepMemoryPolynomialFit
    % sweep memory polynomial parameters
    tic
    [~,~,~,evmList] = rfSim_pa_memoryModel([paInSig(1:nCapturedSamps),paOutSig_ideal(1:nCapturedSamps)],...
        'coefficientFinder',Degree_sweep,MemDepth_sweep,[],[],setInputPowerRangeDbm,Rohm); hold on
    toc
    % plot contourf - evm in %
    plt_contourf(evmList, [], {Degree_sweep{end}, MemDepth_sweep},...
        {[070404,1,1,1], 'Polynomial Degree', 'Memory Depth', 'EVM. Sweep [%%] for PA output fit signal'}, 'min')
end
% Sweep Polynomial Order, Memory Depth, Cross Delay and Cross Lag
% To search the minimum rms error in %

isSweepCrossTerm = 0;
% Create Memory Polynomial Fitting PA Output with/without Cross Terms
% To evaluate the rms error of fitting results

if isSweepCrossTerm % set memory polynomial cross term parameters
    Degree_fit_set={1,2,11};
    MemoryDepth_fit_set={0,1,4};
    CrossDelay_set{3}=1;
    CrossLag_set{3}=1;
else % set memory polynomial without cross term parameters
    Degree_fit_set={1,Degree_step,7}; % set max. degree based on min. rmsError
    MemoryDepth_fit_set={0,1,2}; % set max. memory based on min. rmsError
    CrossDelay_set=[];
    CrossLag_set=[];
end

% get max. delay
try
    maxDelay = MemoryDepth_fit_set{end}
catch
    maxDelay = MemoryDepth_fit_set;
end

isSelectDlSlots = 0;
if isSelectDlSlots
    idxFeed2PaModelFit = dlSlots;
else
    idxFeed2PaModelFit = 1:nCapturedSamps;
end
% captured pa signal for pa model fit
paInSig_feed2PaModelFit = paInSig(idxFeed2PaModelFit);
paOutSig_feed2PaModelFit = paOutSig_ideal(idxFeed2PaModelFit);

if 1 % memory less

    % calculate fit coefficients - memoryLess
    [coefs_memLessFit,~,~,evm_memLessFit] = rfSim_pa_memoryModel([paInSig_feed2PaModelFit,paOutSig_feed2PaModelFit],...
        'coefficientFinder',Degree_fit_set,[0],CrossDelay_set,CrossLag_set,setInputPowerRangeDbm,Rohm);

    % create pa output signal - memoryLess fit
    paOutSig_memLessFit = rfSim_pa_memoryModel({paInSig(1:nSamps),coefs_memLessFit},...
        'signalGenerator',Degree_fit_set,[0],CrossDelay_set,CrossLag_set);

    % apply flatness fir
    if isPaFlatness
        paOutSig_memLessFit_flat = conv(paOutSig_memLessFit, b_coefs_flat, 'same');
    end

    % apply flatness fir
    if isPaFlatness
        % apply flatness fir
        paOutSig_memLessFit_flat = conv(paOutSig_memLessFit_flat, b_coefs_flat, 'same');

        % adjust power
        paOutSig_memLessFit_flat = powerDbm(paOutSig_memLessFit_flat, 'set', paOutPwrDbm, 'dBm', dlSlots);
    else
        paOutSig_memLessFit_flat = paOutSig_memLessFit;
    end

    pltComm.AMAM([paInSig, paOutSig_memLessFit_flat], {092002, 'paOutSig_memLessFit_flat'}, 0.01, setInputPowerRangeDbm); % plot amam
end

% calculate fit coefficients - memory
[coefs_memoryFit,~,~,evm_memoryFit] = rfSim_pa_memoryModel([paInSig_feed2PaModelFit,paOutSig_feed2PaModelFit],...
    'coefficientFinder',Degree_fit_set,MemoryDepth_fit_set,CrossDelay_set,CrossLag_set,setInputPowerRangeDbm,Rohm);

% create pa output signal - memory fit
paOutSig_memoryFit = rfSim_pa_memoryModel({paInSig(1:nSamps),coefs_memoryFit},...
    'signalGenerator',Degree_fit_set,MemoryDepth_fit_set,CrossDelay_set,CrossLag_set);

% apply flatness fir
if isPaFlatness
    % apply flatness fir
    paOutSig_memoryFit_flat = conv(paOutSig_memoryFit, b_coefs_flat, 'same');

    % adjust power
    paOutSig_memoryFit_flat = powerDbm(paOutSig_memoryFit_flat, 'set', paOutPwrDbm, 'dBm', dlSlots);
else
    paOutSig_memoryFit_flat = paOutSig_memoryFit;
end

if 1
    % get dlSlots for memory fit
    dlSlots_fit = dlSlots;
    dlSlots_fit(1:maxDelay) = logical(0);
else
    error('!')
end

% check power
paOutPwrDbm = powerDbm(paOutSig(dlSlots));
paOutPwrDbm_memoryFit_flat = powerDbm(paOutSig_memoryFit_flat(dlSlots_fit));
paOutPwrDbm_memoryFit = powerDbm(paOutSig_memoryFit(dlSlots_fit));
deltaPowerDb_memoryFit = paOutPwrDbm - paOutPwrDbm_memoryFit_flat;
if abs(deltaPowerDb_memoryFit) > 0.5
    error('delta power is larger than 0.5dB')
end

if 1 % plot
    aclrVal_paOutFit = pltComm.ACLR(paOutSig_memoryFit(1+maxDelay:end), {091901, 'paOutSig_memoryFit'}) % plot aclr
    pltComm.TIME(paOutSig_memoryFit(1+0:end), {092001, 'paOutSig_memoryFit'}, [], 'usec'); % plot time
    pltComm.AMAM([paInSig(1+maxDelay:end), paOutSig_memoryFit(1+maxDelay:end)], {092002, 'paOutSig_memoryFit'}, 0.01, setInputPowerRangeDbm); % plot amam
end

if 1 % demodulation
    demod_paOutFit = wlanDemodulation(paOutSig_memoryFit, wvParams, 'tx')
    demod_paOutFit_flat = wlanDemodulation(paOutSig_memoryFit_flat, wvParams, 'tx')
end

%% DPD - comm.DPD (Inverse approach)

if 0
    dispConditions_dpdFit = "dpdMemLessFit";
else
    dispConditions_dpdFit = "dpdMemoryFit";
end

if 0 % set signal captured range
    error('!?'), idxFeed2dpd = dlSlots;
elseif 1*0
    idxFeed2dpd = 1+maxDelay:21976;
else
    idxFeed2dpd = 1+maxDelay:nCapturedSamps;
end
if idxFeed2dpd(1) ~= 1+maxDelay
    error('!')
end
paInSig_feed2Dpd = paInSig(idxFeed2dpd);
dlSlots_feed2Dpd = dlSlots(idxFeed2dpd);

if 1
    isDpdCoefUsePaSignal = 'memPolyFit'; % 1: memory polynomail fit, 2: memoryLess polynomila fit, 3: pa output signal
    paOutSig_feed2Dpd = paOutSig_memoryFit_flat(idxFeed2dpd);
    paOutSig_feed2Dpd_memLess = paOutSig_memLessFit_flat(idxFeed2dpd);
elseif 1*1
    isDpdCoefUsePaSignal = 'memLessPolyFit';
    paOutSig_feed2Dpd = paOutSig_memLessFit_flat(idxFeed2dpd);
    paOutSig_feed2Dpd_memLess = paOutSig_memLessFit_flat(idxFeed2dpd);
else
    isDpdCoefUsePaSignal = 'paOutput'; % 1: memory polynomail fit, 2: memoryLess polynomila fit, 3: pa output signal
    paOutSig_feed2Dpd = paOutSig_ideal(idxFeed2dpd);
end

% chck power
paInPwrDbm_feed2Dpd = powerDbm(paInSig_feed2Dpd(dlSlots_feed2Dpd));
paOutPwrDbm_feed2Dpd = powerDbm(paOutSig_feed2Dpd(dlSlots_feed2Dpd));

clear dpm % set DPD parameters
dpm.DesiredAmplitudeGaindB = (paOutPwrDbm_feed2Dpd - paInPwrDbm_feed2Dpd) / 2;
dpm.PolynomialType = 'Memory polynomial';
dpm.Algorithm = 'Least squares';

if 1
    isDpdCalculateMethod = "coefficientFinderInverse";
else
    isDpdCalculateMethod = "comm.DPD";
end

isSweepDPD = 1;
% set memory polynomial parameters for dpd sweep
if 1
    Degree_dpm_sweep = [7];
    MemoryDepth_dpm_sweep = [0];
else
    Degree_dpm_sweep = [5 7];
    MemoryDepth_dpm_sweep = [0 1];
end

if 1
    CrossDelay_dpm_set = []; % disable cross delay
    CrossLag_dpm_set = []; % disable cross lag
end

for K=Degree_dpm_sweep
    for M=MemoryDepth_dpm_sweep

        isDpdCoefsReused = 0; % set to reuse the dpd model/coeficients
        isCfrAfterDpd = 1; % set turn cfr on/off
        cfr_nIteration = 20; % set cfr_nIteration
        cfr_paprDb_vec = [8 8.5 9 9.5 10 10.5 11 11.5 12 12.5 13]; % set cfr vectors
        cfr_paprDb_vec = [8 9 10 11 12 13 14 15 16]; % set cfr vectors
        cfr_paprDb_vec = [12]; % set cfr vectors

        for idxCfr=1:numel(cfr_paprDb_vec)
            close all

            deltaPowerDb_dpd = 0; % init.
            SIMULATING_tuningPower = 1; % init. tuning the dpm.DesiredAmplitudeGaindB to meeting pa output power w/ dpd
            while SIMULATING_tuningPower

                if isDpdCoefsReused
                elseif deltaPowerDb_dpd > 1
                    dpm.DesiredAmplitudeGaindB = dpm.DesiredAmplitudeGaindB + 1/2;
                elseif deltaPowerDb_dpd < -1
                    dpm.DesiredAmplitudeGaindB = dpm.DesiredAmplitudeGaindB - 1/2;
                elseif deltaPowerDb_dpd > 0.5
                    dpm.DesiredAmplitudeGaindB = dpm.DesiredAmplitudeGaindB + 0.5/2;
                elseif deltaPowerDb_dpd < -0.5
                    dpm.DesiredAmplitudeGaindB = dpm.DesiredAmplitudeGaindB - 0.5/2;
                elseif deltaPowerDb_dpd > 0.1
                    dpm.DesiredAmplitudeGaindB = dpm.DesiredAmplitudeGaindB + 0.1/2;
                elseif deltaPowerDb_dpd < -0.1
                    dpm.DesiredAmplitudeGaindB = dpm.DesiredAmplitudeGaindB - 0.1/2;
                elseif deltaPowerDb_dpd > 0.05
                    dpm.DesiredAmplitudeGaindB = dpm.DesiredAmplitudeGaindB + 0.05/2;
                elseif deltaPowerDb_dpd < -0.05
                    dpm.DesiredAmplitudeGaindB = dpm.DesiredAmplitudeGaindB - 0.05/2;
                elseif deltaPowerDb_dpd > 0
                    dpm.DesiredAmplitudeGaindB = dpm.DesiredAmplitudeGaindB + 0.005;
                else
                    dpm.DesiredAmplitudeGaindB = dpm.DesiredAmplitudeGaindB - 0.005;
                end

                % set dpd parameters (max.)
                if isSweepDPD
                    % get initial
                    K1 = Degree_dpm_sweep(1);
                    M1 = MemoryDepth_dpm_sweep(1);
                    % set degree and depth from sweeping parameters
                    if 1
                        Degree_dpm_set = {1,1,K};
                        MemoryDepth_dpm_set = {0,1,M};
                    else
                        Degree_dpm_set = {1,2,K};
                        MemoryDepth_dpm_set = {0,1,M};
                    end
                else
                    % set initial
                    K1 = 1;
                    M1 = 0;
                    % set degree and depth from sweeping parameters
                    Degree_dpm_set = Degree_fit_set;
                    MemoryDepth_dpm_set = MemoryDepth_fit_set;
                end

                if 0 % check
                    isequal(Degree_fit_set, Degree_dpm_set)
                    isequal(MemoryDepth_fit_set, MemoryDepth_dpm_set)
                end

                % get maximum delay
                maxDelay_dpd = MemoryDepth_dpm_set{end};

                % update dpm
                dpm.Degree = Degree_dpm_set{end};
                dpm.MemoryDepth = MemoryDepth_dpm_set{end};

                switch isDpdCalculateMethod % Estimate the dpd memory polynomial coefficients
                    case "comm.DPD"

                        % Generate estimator model
                        estimator = comm.DPDCoefficientEstimator( ...
                            'DesiredAmplitudeGaindB',dpm.DesiredAmplitudeGaindB, ...
                            'PolynomialType',dpm.PolynomialType, ...
                            'Degree',dpm.Degree,'MemoryDepth',(dpm.MemoryDepth+1), ...
                            'Algorithm',dpm.Algorithm); % set 'MemoryDepth' to dpm.MemoryDepth+1 for length alignment

                        if ~isDpdCoefsReused
                            % Get coefficients
                            switch dispConditions_dpdFit
                                case "dpdMemLessFit"
                                    coefs_dpd = estimator(paInSig_feed2Dpd,paOutSig_feed2Dpd_memLess);
                                case "dpdMemoryFit"
                                    coefs_dpd = estimator(paInSig_feed2Dpd,paOutSig_feed2Dpd);
                            end
                            coefs_dpd_org = coefs_dpd; % store original
                        else
                            if deltaPowerDb_dpd == 0
                                coefs_dpd = coefs_dpd_org;
                            elseif deltaPowerDb_dpd > 1
                                coefs_dpd = coefs_dpd * 10^(0.5/20);
                            elseif deltaPowerDb_dpd < -1
                                coefs_dpd = coefs_dpd / 10^(0.5/20);
                            elseif deltaPowerDb_dpd > 0.5
                                coefs_dpd = coefs_dpd * 10^(0.01/20);
                            elseif deltaPowerDb_dpd < -0.5
                                coefs_dpd = coefs_dpd / 10^(0.01/20);
                            elseif deltaPowerDb_dpd > 0.1
                                coefs_dpd = coefs_dpd * 10^(0.05/20);
                            elseif deltaPowerDb_dpd < -0.1
                                coefs_dpd = coefs_dpd / 10^(0.05/20);
                            elseif deltaPowerDb_dpd > 0
                                coefs_dpd = coefs_dpd * 10^(0.002/20);
                            else
                                coefs_dpd = coefs_dpd / 10^(0.002/20);
                            end
                        end

                        % Create dpd model using comm.DPD object
                        dpd_model = comm.DPD('PolynomialType',dpm.PolynomialType,'Coefficients',coefs_dpd);

                        % Generate pa input dpd
                        paInSig_wiDpd = dpd_model(paInSig);

                        if 0
                            pltComm.AMAM([paInSig, paInSig_wiDpd], {1014, 'paInSig_wiDpd '+dispConditions_dpdFit}, 0.01); % plot amam
                        end

                        % maxDelay_dpd = 0;

                        coefs_dpd(:);

                    case "coefficientFinderInverse"
                        % Normalize pa output dpd feed
                        LinearGain = 10^(dpm.DesiredAmplitudeGaindB/10);
                        paOutSig_feed2Dpd_norm = paOutSig_feed2Dpd / LinearGain; % normalize

                        % Get the inverse coefs. between pa input/output for dpd used
                        [coefs_paInverse, ~] = rfSim_pa_memoryModel([paInSig_feed2Dpd,paOutSig_feed2Dpd_norm],...
                            'coefficientFinderInverse',Degree_dpm_set,MemoryDepth_dpm_set,CrossDelay_dpm_set,CrossLag_dpm_set);

                        % Create paInputWiDpd using paInputWiDpd and coefs_paInverse
                        paInSig_wiDpd = rfSim_pa_memoryModel({paInSig,coefs_paInverse},...
                            'signalGenerator',Degree_dpm_set,MemoryDepth_dpm_set,CrossDelay_dpm_set,CrossLag_dpm_set) * 1;

                        % return coefs.
                        coefs_dpd = coefs_paInverse;
                end

                % set idxMem
                idxMemDpd = 1+maxDelay_dpd:length(paInSig_wiDpd);
                dlSlots_dpd = dlSlots;
                dlSlots_dpd(1:maxDelay_dpd) = logical(0);

                if 1 % check power
                    paInPwrDbm
                    paInPwrDbm_wiDpd = powerDbm(paInSig_wiDpd(dlSlots_dpd));
                end

                if 0 % check dpd expansion
                    pltComm.TIME(paInSig_wiDpd(idxFeed2dpd), {092401, 'paInSig_wiDpd(idxFeed2dpd)'}, [], 'samples', [], 1); % plot time
                    pltComm.TIME(paInSig_feed2Dpd(:), {092401, 'paInSig_feed2Dpd'}, [], 'samples', [], 1); % plot time

                    [ccdfVal_paInWiDpd,~,idxPaprDbMax] = pltComm.CCDF(paInSig_wiDpd(idxMemDpd), {[092402], 'paInSig_wiDpd'}, [], 'max');
                    [ccdfVal_paInFeed2Dpd,~,idxPaprDbMax] = pltComm.CCDF(paInSig_feed2Dpd, {[092402], 'paInSig_feed2Dpd'}, [], 'max');
                end

                if isCfrAfterDpd || ~isnan(cfr_paprDb_vec(idxCfr))% apply CFR

                    if 1
                        cfr_paprDb = cfr_paprDb_vec(idxCfr);
                    else
                        cfr_paprDb = 10.0 + 0.5*(-1); % set cfr_paprDb for HT20
                        cfr_paprDb = 10.0 + 5.0*(1); % set cfr_paprDb
                    end

                    paInSig_wiDpd_org = paInSig_wiDpd;
                    sprintf('%% set transmit index, remove abnormal points (peak voltage out of dpd expansion point)');
                    % limit the dpd expansion in dB by CFR
                    % paInSig_wiDpdCfr = dfe_cfr_nco_2(paInSig_wiDpd_org, cfr_paprDb, [], [], cfr_nIteration, []);
                    paInSig_wiDpdCfr = dfe_cfr_nco_2(paInSig_wiDpd, cfr_paprDb, [], [], cfr_nIteration, []);

                    % reutrn cfr signal
                    paInSig_wiDpd = paInSig_wiDpdCfr;
                    % plot
                    if 0
                        pltComm.TIME(paInSig_wiDpdCfr(idxFeed2dpd), {092401, 'paInputWiDpdCfr'}, [], 'samples', [], 1); % plot time
                        [ccdfVal_paInSig_wiDpdCfr,~,idxPaprDbOverflow] = pltComm.CCDF(paInSig_wiDpdCfr(idxMemDpd), {092402, 'paInputWiDpdCfr'}, [], 14);
                    end
                else
                    cfr_paprDb = nan; % assign nan
                end

                % Apply paInputDPD to pa memoryless model
                switch isDpdCoefUsePaSignal
                    case 'memLessPolyFit'
                        paOutSig_wiDpd = paMemLess_obj(paInSig_wiDpd);
                    otherwise
                        % Calculate coefficients of pa model
                        if 1
                            [coefs_pa,~,~,evm_memoryFit] = rfSim_pa_memoryModel([paInSig_feed2Dpd,paOutSig_feed2Dpd],...
                                'coefficientFinder',Degree_dpm_set,MemoryDepth_dpm_set,CrossDelay_dpm_set,CrossLag_dpm_set,setInputPowerRangeDbm,Rohm);
                        else
                            [coefs_pa,~,~,evm_memoryFit] = rfSim_pa_memoryModel([paInSig_feed2PaModelFit,paOutSig_feed2PaModelFit],...
                                'coefficientFinder',Degree_dpm_set,MemoryDepth_dpm_set,CrossDelay_set,CrossLag_set,setInputPowerRangeDbm,Rohm);
                        end

                        % Create paOutputWiDpd using paInputWiDpd and coefs_pa
                        paOutSig_wiDpd = rfSim_pa_memoryModel({paInSig_wiDpd,coefs_pa},...
                            'signalGenerator',Degree_dpm_set,MemoryDepth_dpm_set,CrossDelay_dpm_set,CrossLag_dpm_set) * 1;

                        if 0 % debug
                            demod_paOutSig_wiDpdWoFlat = wlanDemodulation(paOutSig_wiDpd(:),{wvParams,['paOutSig_wiDpd']},'tx');
                        end
                end

                if 1 % check delta power
                    paOutPwrDbm_wiDpd = powerDbm(paOutSig_wiDpd(idxMemDpd));
                    deltaPowerDb_dpd = paOutPwrDbm_feed2Dpd - paOutPwrDbm_wiDpd % check deltaPower_dpd_dB
                    sprintf('delta power between pa output wi/wo DPD %.4f dB\n', deltaPowerDb_dpd);
                    deltaPowerDb_list_dpd(M-M1+1,K-K1+1) =deltaPowerDb_dpd;
                end

                if abs(deltaPowerDb_dpd) < 0.005
                    SIMULATING_tuningPower = 0;
                end
            end % SIMULATING_tuningPower

            if 1 % check evm
                evmInband_dpd_percent = evm(paInSig(idxMemDpd),paOutSig_wiDpd(idxMemDpd),'%','inband',...
                    [aclrVal_paIn.SampleRateMHz*1e6, -aclrVal_paIn.CarrierBwMHz*1e6/2, aclrVal_paIn.CarrierBwMHz*1e6/2]); % sampleRate_freqL_freqU
            end

            % Plot
            if idxCfr==1
                aclrVal_paInFeed2Dpd= pltComm.ACLR(paInSig_feed2Dpd(:), {091901, 'paInSig_feed2Dpd'}); % plot aclr
                aclrVal_paOutFeed2Dpd = pltComm.ACLR(paOutSig_feed2Dpd(:), {091901, 'paOutSig_feed2Dpd'}); % plot aclr
            end
            aclrVal_paOutWiDPD = pltComm.ACLR(paOutSig_wiDpd(idxMemDpd), {091901, 'paOutSig_wiDpd'}); % plot aclr

            if 0
                isNorm = 0;
                linearGain = 10^((paOutPwrDbm_feed2Dpd - paInPwrDbm_feed2Dpd)/20);
                pltComm.TIME(paInSig_feed2Dpd(:), {[092711], 'paInSig_feed2Dpd'}, [], 'samples', [], isNorm); % plot time
                pltComm.TIME(paOutSig_wiDpd(idxFeed2dpd)/linearGain, {[092711], 'paInSig_feed2Dpd/gain'}, [], 'samples', [], isNorm); % plot time
                pltComm.TIME(paOutSig_feed2Dpd(:)/linearGain, {[092711], 'paOutSig_feed2Dpd/gain'}, [], 'samples', [], isNorm); % plot time
            end
            if 0
                isNormGain = 20*log10(LinearGain);
                pltComm.AMAM([paInSig_feed2Dpd(:), paOutSig_feed2Dpd(:)], {09203, 'paOutSig_feed2Dpd'}, 0.01,setInputPowerRangeDbm,isNormGain); % plot amam
                pltComm.AMAM([paInSig(idxMemDpd), paInSig_wiDpd(idxMemDpd)], {09203, 'paInSig_wiDpd'}, 0.01,setInputPowerRangeDbm); % plot amam
                pltComm.AMAM([paInSig(idxMemDpd), paOutSig_wiDpd(idxMemDpd)], {09203, 'paOutSig_wiDpd'}, 0.01,setInputPowerRangeDbm,isNormGain); % plot amam
            end
            ccdfVal_paInFeed2Dpd = pltComm.CCDF(paInSig_feed2Dpd(:), {092003, 'paInSig_feed2Dpd'}, 0.01); % plot amam
            ccdfVal_paInWiDpd = pltComm.CCDF(paInSig_wiDpd(idxMemDpd), {092003, 'paInSig_wiDpd'}, 0.01); % plot amam
            ccdfVal_paOutFeed2Dpd = pltComm.CCDF(paOutSig_feed2Dpd(:), {092003, 'paOutSig_feed2Dpd'}, 0.01); % plot amam
            ccdfVal_paOutWiDpd = pltComm.CCDF(paOutSig_wiDpd(idxMemDpd), {092003, 'paOutSig_wiDpd'}, 0.01); % plot amam
            ccdfVal_paInWiDpd_org = pltComm.CCDF(paInSig_wiDpd_org(idxMemDpd), {092003, 'paInSig_wiDpd_org'}, 0.01); % plot amam

            % Demodulation
            if 1 % set condtions to dispLgnd
                dispConditions = sprintf(".CFR %ddB .Flat %ddB",cfr_paprDb, PaModelFlatnessPerCBW);
            end
            if idxCfr==1
                demod_paOut = wlanDemodulation(paOutSig(:),{wvParams,'paOutSig'},'tx');
            end
            demod_paOutWiDPD = wlanDemodulation(paOutSig_wiDpd(:),{wvParams,['paOutSig_wiDpd'+dispConditions]},'tx');

            if ~isSweepDPD
                dispLgnd_dpd = sprintf('DPD output,K=%d ,M=%d\nRMSErr. %.2f%%, \x0394Pout %.2fdB', K1, M1, evmInband_dpd_percent, deltaPowerDb_dpd);
            end

            if 1 % print results
                if 1 % conditions
                    result.MCS = MCS_index;
                    result.SignalFormat = signalFormat;
                    result.InputSignalPaprDb = round(ccdfVal_paInWiDpd,2);
                    % result.OutputSignalPaprDb = round(ccdfVal_paOutWiDpd,2);

                    if isDpdCoefsReused
                        result.DPD = "Coefs_Reused";
                    else
                        result.DPD = "Coefs_New";
                    end

                    result.CfrDb = round(cfr_paprDb,2);
                    result.SampleRateMHz = fsMHz;
                    result.CoefsModel = isDpdCoefUsePaSignal;
                    result.Degree = dpm.Degree;
                    result.Memory = dpm.MemoryDepth;
                    result.DpdMethod = isDpdCalculateMethod;

                    if isSweepDPD==0
                        result.DPD = "DPD_Off";
                        result.DpdMethod = "NaN";
                    end
                end
                if isPaModelFlatCompensateLoss
                    result.PaModelFlatnessPerCBW = "Cp" + string(PaModelFlatnessPerCBW)+"dB";
                else
                    result.PaModelFlatnessPerCBW = string(PaModelFlatnessPerCBW)+"dB";
                end
                % result.PaOutChPwrDbm = round(aclrVal_paOutWiDPD.ChannelPowerDbm,2);
                result.PaInPwrDbm = powerDbm(paInSig_wiDpd_org(dlSlots_dpd));
                result.PaOutEvmRmsDb = round(max(demod_paOutWiDPD.evmRMS),2);
                result.PaOutACLR1Db = round(aclrVal_paOutWiDPD.ACLR1,2);
                result.PaOutACLR2Db = round(aclrVal_paOutWiDPD.ACLR2,2);

                if idxCfr==1
                    result_paOut = result;
                    result_paOut.CfrDb = nan;
                    result_paOut.DPD = "DPD Off";
                    result_paOut.Degree = nan;
                    result_paOut.Memory = nan;
                    result_paOut.DpdMethod = "Nan";
                    % result_paOut.PaOutChPwrDbm = round(aclrVal_paOut.ChannelPowerDbm,2);
                    result.PaInPwrDbm = powerDbm(paInSig_feed2Dpd(dlSlots_feed2Dpd));
                    result_paOut.PaOutEvmRmsDb = round(max(demod_paOut.evmRMS),2);
                    result_paOut.PaOutACLR1Db = round(aclrVal_paOutFeed2Dpd.ACLR1,2);
                    result_paOut.PaOutACLR2Db = round(aclrVal_paOutFeed2Dpd.ACLR2,2);
                    try
                        result_tab = vertcat(result_tab,struct2table(result_paOut));
                    catch
                        result_tab = struct2table(result_paOut);
                    end
                end

                % vertcat table
                result_tab = vertcat(result_tab,struct2table(result));
            end

        end %idxCfr

    end %M
end %K

c = nan(10, 5);
if 0 % summary plot
    c(:,1) = [1:10];
    c(2:10,2) = [8:1:16];

    dispLgnd_format = "EHT320 15dBm w/"
    color = [0 0.4470 0.7410];

    dispLgnd_format = "HE160 18dBm w/"
    color = [0.8500 0.3250 0.0980];

    dispLgnd_format = "VHT80 23dBm w/"
    color = [0.9290 0.6940 0.1250];

    dispLgnd_format = "HT40 24dBm w/"
    color = [0.4940 0.1840 0.5560];

    dispLgnd_format = "HT20 27dBm w/"
    color = [0.4660 0.6740 0.1880];

    switch paModelFlatnessPerCBW
        case {4}
            paModelFlatnessPerCBW = 4
            mk = '+'
        case {2}
            paModelFlatnessPerCBW = 2
            mk = 'v'
        case {0}
            paModelFlatnessPerCBW = 0
            mk = 'o'
        case '2cp'
            paModelFlatnessPerCBW = 2
            mk = 'square'
            isCP = 'cp'
        case '4cp'
            paModelFlatnessPerCBW = 4
            mk = 'o'
            isCP = 'cp'
    end
    figure(102301)
    if paModelFlatnessPerCBW~=0
        plt = scatter(c(:,1), c(:,4), mk, 'CData', color, 'LineWidth', 2, 'DisplayName', sprintf(dispLgnd_format + " Flatness %ddB %s", paModelFlatnessPerCBW, isCP));
    else
        plt2 = plot(c(:,1), c(:,4), '--', 'LineWidth', 1, 'Color', color,...
            'DisplayName', sprintf(dispLgnd_format + " Flatness %ddB %s", paModelFlatnessPerCBW));
    end
    hold on, grid on
    xticks(c(:,1));
    xtklab = num2cell(c(:,2)); xtklab{1} = "DPD off"; % if 1, xtklab{end} = "CFR off"; end
    xticklabels(xtklab), xlabel('CFR [dB]'), ylabel('EVM RMS [dB]'), title('PA Output wi/wo DPD EVM vs CFR')
    legend, legend('Location','eastoutside')
    figure(102302)
    if paModelFlatnessPerCBW~=0
        plt = scatter(c(:,1), c(:,5), mk, 'CData', color, 'LineWidth', 2, 'DisplayName', sprintf(dispLgnd_format + " Flatness %ddB %s", paModelFlatnessPerCBW, isCP));
    else
        plt2 = plot(c(:,1), c(:,5), '--', 'LineWidth', 1, 'Color', color,...
            'DisplayName', sprintf(dispLgnd_format + " Flatness %ddB", paModelFlatnessPerCBW));
    end
    hold on, grid on
    xticklabels(xtklab), xlabel('CFR [dB]'), ylabel('EVM RMS [dB]'), title('PA Output wi/wo DPD ACLR vs CFR')
    legend, legend('Location','eastoutside')


end

if 0 % comparsion dpd coefs. reused
    b = nan(8, 8);
    b(:,1) = [1:8];
    b(:,2) = [7:1:14];
    figure(1008)

    dispLgnd_format = "EHT160 w/"
    color = [0 0.4470 0.7410];

    dispLgnd_format = "HE160 w/"
    color = [0.8500 0.3250 0.0980];

    dispLgnd_format = "VHT80 w/"
    color = [0.9290 0.6940 0.1250];

    dispLgnd_format = "HT40 w/"
    color = [0.4940 0.1840 0.5560];

    dispLgnd_format = "HT20 w/"
    color = [0.4660 0.6740 0.1880];

    plt = plot(b(:,1), b(:,3), 'Color', color, 'LineWidth', 2, 'DisplayName', dispLgnd_format + " DPD coefs. new"); hold on, grid on
    plt = scatter(b(:,1), b(:,3), 'filled', 'CData', color, 'LineWidth', 2, 'DisplayName', dispLgnd_format + " DPD coefs. new"); hold on, grid on
    plt = scatter(b(:,1), b(:,4),  'diamond', 'CData', color, 'LineWidth', 2, 'DisplayName', dispLgnd_format + " DPD coefs. reused from  HT20 27dBm"); hold on, grid on
    plt = scatter(b(:,1), b(:,5), '+', 'CData', color, 'LineWidth', 2, 'DisplayName', dispLgnd_format + " DPD coefs. reused from  HT40 24dBm"); hold on, grid on
    plt = scatter(b(:,1), b(:,6), '*', 'CData', color, 'LineWidth', 2, 'DisplayName', dispLgnd_format + " DPD coefs. reused from  VHT80 23dBm"); hold on, grid on
    plt = scatter(b(:,1), b(:,7), 'square', 'CData', color, 'LineWidth', 2, 'DisplayName', dispLgnd_format + " DPD coefs. reused from  HE160 18dBm"); hold on, grid on
    plt = scatter(b(:,1), b(:,8), 'v', 'CData', color, 'LineWidth', 2, 'DisplayName', dispLgnd_format + " DPD coefs. reused from  EHT160 15dBm"); hold on, grid on
    xticks(b(:,1));
    xtklab = num2cell(b(:,2));
    xticklabels(xtklab), xlabel('CFR [dB]'), ylabel('EVM RMS [dB]'), title(sprintf('PA Output wi/wo Resued DPD Coefficients\n EVM vs CFR'))
    lgd = legend('show', 'Location','eastoutside'),

    lgd.NumColumns = 5

    find_obj_sca = findobj(gca, 'Type', 'scatter')
    find_obj_plt = findobj(gca, 'Type', 'line')

    if 0, delete(find_obj_sca(7:12)), end
    if 0, delete(find_obj_plt(1)), end

end

a = nan(14, 5);
if 0 % summary plot
    a(:,1) = [1:14];
    a(2:13,2) = [7.5:0.5:13];
    figure(1001)
    plt = scatter(a(:,1), a(:,4), 'filled', 'DisplayName', signalFormat); hold on, grid on
    xticks(a(:,1));
    xtklab = num2cell(a(:,2)); xtklab{1} = "DPD off"; if 1, xtklab{end} = "CFR off"; end
    xticklabels(xtklab), xlabel('CFR [dB]'), ylabel('EVM RMS [dB]'), title('PA Output wi/wo DPD EVM vs CFR')
    plt2 = plot(a(:,1), specEVM*ones(size(a,1),1), '--', 'LineWidth', 2, 'Color', get(plt, 'CData'),...
        'DisplayName', sprintf("%s @%.2fdBm Req.",signalFormat, 27.68));
    legend, legend('Location','eastoutside')
end

if 0 % summary plot - aclr
    figure(1004)
    plt3 = scatter(a(:,1), a(:,5), 'filled', 'DisplayName', sprintf("%s %.2fdBm", signalFormat, 27.68)); hold on, grid on
    xticks(a(:,1));
    xticklabels(xtklab), xlabel('CFR [dB]'), ylabel('ACLR [dBc]'), title('PA Output wi/wo DPD ACLR vs CFR')
    legend, lgend('Location','eastoutside')
end

%% DPD - Plot Results

if isSweepDPD * 0
    % plot contourf
    plt_contourf(evm_memoryFit, [], {Degree_dpm_sweep, MemoryDepth_dpm_sweep},...
        {[091101,1,3,1], 'Polynomial Degree', 'Memory Depth', 'RMS Err. Sweep [%%] for PA output fit signal'}, 'min')

    plt_contourf(evmList_dpd, [], {Degree_dpm_sweep, MemoryDepth_dpm_sweep},...
        {[091101,1,3,2], 'Polynomial Degree', 'Memory Depth', 'EVM Sweep [%%] for PA output w/ DPD signal'}, 'min')

    plt_contourf(deltaPowerDb_list_dpd, [], {Degree_dpm_sweep, MemoryDepth_dpm_sweep},...
        {[091101,1,3,3], 'Polynomial Degree', 'Memory Depth', '\x0394Pout Sweep [dB] for PA output w/ DPD signal'}, 'absMin')
end
%%
function plt_contourf(inputMatrix, contourfLevels, pltTickLabels_cell ,fnum_dispLabels_dispTitle_cell, isMarker)
fnum = fnum_dispLabels_dispTitle_cell{1};
dispXLabel = fnum_dispLabels_dispTitle_cell{2};
dispYLabel = fnum_dispLabels_dispTitle_cell{3};
dispTitle = fnum_dispLabels_dispTitle_cell{4};
setXTicklabel = pltTickLabels_cell{1};
setYTicklabel = pltTickLabels_cell{2};
if isempty(contourfLevels)
    % Set contourfLevels
    minVal = min(inputMatrix,[],'all');
    maxVal = max(inputMatrix,[],'all');
    stepVal = (maxVal - minVal) / 8;
    contourfLevels = minVal:stepVal:maxVal;
end

% plot contourf
figure(fnum(1))
if numel(fnum)==4
    subplot(fnum(2),fnum(3),fnum(4))
end
contourf(inputMatrix, contourfLevels, 'ShowText','on'), grid minor, hold on
c = colorbar;
clim([min(contourfLevels) max(contourfLevels)]);% set colorbar range
c.Ticks = contourfLevels;% set colorbar ticks
xlabel(sprintf(dispXLabel))
ylabel(sprintf(dispYLabel))
title(sprintf(dispTitle))
set(gca, 'XTick', 1:size(inputMatrix,2), 'XTickLabel', setXTicklabel); % xtick
set(gca, 'YTick', 1:size(inputMatrix,1), 'YTickLabel', setYTicklabel); % ytick

if exist('isMarker','var')&&~isempty(isMarker)
    % marker value
    switch isMarker
        case 'max'
            % max
            [inputMatrix_mk, idx] = max(inputMatrix(:)');
        case 'absMax'
            [inputMatrix_mk, idx] = max(abs(inputMatrix(:)'));
        case 'abxMin'
            [inputMatrix_mk, idx] = min(abs(inputMatrix(:)'));
        otherwise
            % min
            [inputMatrix_mk, idx] = min(inputMatrix(:)');
    end
    [xId, yId] = ind2sub(size(inputMatrix), idx);
    scatter(yId, xId, 'red', 'square', 'LineWidth', 2, 'SizeData', 100)
    txt = text(yId+0.1, xId+0.1, [num2str(round(inputMatrix_mk,4)), ' ' ,isMarker]);
    set(txt, 'FontSize', 14, 'Color', 'red', 'FontName', 'Arial', 'FontWeight', 'bold');
end

end