%% 2024-09-20, remove nan


function varargout = rfSim_pa_memoryLessLUTModel(signals,inputPowerBackOffDb,ReferenceImpedance,fnum)
% Generate look-up table (LUT) for input of comm.MemorylessNonlinearity
% Return varargout - cell array containing:
% inOutTable: LUT
% pa_memoryLessNonlinearity_class: comm.MemorylessNonlinearity object
% paOutputFitMemless: pa output memory less fitting data

paInput = signals(:,1);
paOutput = signals(:,2);

isFnum = exist('fnum', 'var') && ~isempty(fnum);
if exist('ReferenceImpedance','var') && ~isempty(ReferenceImpedance)
    Rohm = ReferenceImpedance;
else
    Rohm = 1;
end
if exist('inputPowerBackOffDb','var') && ~isempty(inputPowerBackOffDb)
    isInputPowerBackOffDb = 1;
else
    isInputPowerBackOffDb = 0;
end

% Calculate pa psd
paInputDbm = powerDbm(paInput,'psd',[],'dBm',[],Rohm);
paOutputDbm = powerDbm(paOutput,'psd',[],'dBm',[],Rohm);

% Partitialize input power values into bins
[N,edges,idx] = histcounts(paInputDbm, 'BinWidth', 0.5);

% Set the input power minimum
if isInputPowerBackOffDb
    minInPowerDbm = max(paInputDbm) - inputPowerBackOffDb;
else
    minInPowerDbm = min(paInputDbm);
end
minBinIdx = find(edges < minInPowerDbm, 1, 'last');

% Calculate the midpoint of each bin and set to inOutTable
inOutTable = zeros(length(edges)-minBinIdx-1,3);
for p = minBinIdx+1:length(edges)-1
    inOutTable(p-minBinIdx,1) = mean(paInputDbm(idx == p)); % Average input power of p_th bin
    inOutTable(p-minBinIdx,2) = mean(paOutputDbm(idx == p)); % Average output power of p_th bin
    inOutTable(p-minBinIdx,3) = mean(angle(paOutput(idx == p)./paInput(idx == p)))*180/pi; % Average phase shift of p_th bin
end

% 2024-09-20, remove nan
inOutTable = rmmissing(inOutTable, 1);

% Create comm.MemorylessNonlinearity object
pa_memoryLessNonlinearity_class = comm.MemorylessNonlinearity('Method','Lookup table','Table',inOutTable,'ReferenceImpedance',Rohm);
paOutputFitMemless = pa_memoryLessNonlinearity_class(paInput);

% Return to varargout
varargout{1} = inOutTable;
varargout{2} = pa_memoryLessNonlinearity_class;
varargout{3} = paOutputFitMemless;

if isFnum
    figure(fnum)
    plt = PLOT_COMM([],Rohm);
    plt.amam(plt, [paInput,paOutput], {fnum, 'Measurement'});
    plt.amam(plt, [paInput,paOutputFitMemless], {fnum, 'MemoryLess LUT fit'});
end

end