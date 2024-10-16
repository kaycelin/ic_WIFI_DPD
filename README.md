# ic_WLAN_DPD
Simulation of WLAN coporate with DPD
## Main and History: 2024-10-15
- Example_WIFI_Transmit_240919_v1.mlx
- Example_PAModel_and_DPD_240926_v2.mlx
## Simulation flow
<img src="https://github.com/user-attachments/assets/8906efe5-eb29-4099-8f39-68f42f8f9cae" width="60%">

### Set Test Waveform
```js
Rohm = 1 // set resistor in ohm
RbwKHz = 10 // set resolution bw for aclr plot
testSignal = "WIFI"; // select signal
```
### PA Characterization
- Memory/MemoryLess Polynomianl Model with/without Cross Term 
```js
setInputPowerRangeDbm = -70; // set the mimimun power of pa input instantaneous signal
```
```js
//set fitting model using the polynminal order and memory depth - memory polynominal parameters
Degree_fit_set={1,Degree_step,7}; // set max. degree based on min. rmsError
MemoryDepth_fit_set={0,1,2}; // set max. memory based on min. rmsError
CrossDelay_set=[];
CrossLag_set=[];
```
- Compare the memory and memoryless models
<img src="https://github.com/user-attachments/assets/d9fd8cff-63a6-4173-9d86-875bd49bcc50" width="60%">

### DPD Processing
- Set signal captured range
- Set DPD parameters: LinearGain, Model and Algorithm (for comm.DPD tool)
- Set CFR parameter
```js
idxCaptured4Dpd = 1:nSamps;
dpm.DesiredAmplitudeGaindB = (powerDbm(paOutput_feed4Dpd) - powerDbm(paInput_feed4Dpd))/2 + 0.01*(1); // tuning the linear gain
dpm.Degree = 7;
dpm.MemoryDepth = 1;
cfr_paprDb_vec = [9]; // set cfr vectors
```
- Generate or reuse the dpd coefficients, `coefs_dpd`
- Create pa input signal with DPD, `paInputWiDpd`
- Limit the dpd expansion in dB by CFR, `paInputWiDpdCfr`
- Cenerate PA output signal by PA model, `paOutputWiDpd`
- Demodulated siganl and evaluation

<img src="https://github.com/user-attachments/assets/d818458e-ef3f-4296-b59c-fce14501da82" width="60%">
<img src="https://github.com/user-attachments/assets/1eabfcce-179c-476b-8f4b-4aa0ede54c9d" width="60%">
<img src="https://github.com/user-attachments/assets/10ce0449-df2c-4906-a667-684e48fa4991" width="60%">
<img src="https://github.com/user-attachments/assets/5c62e5d9-f3d0-4e9d-8255-491dcb23c5ef" width="60%">
<img src="https://github.com/user-attachments/assets/80190075-1a33-4286-8248-23d613065357" width="60%">




