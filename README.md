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

|**AMAM & AMPM**  |**ACLR**  |
|:--:|:--:|
| ![Image 1](https://github.com/user-attachments/assets/ef385d1f-33b6-419d-880c-b7ff5dc14e64) | ![Image 2](https://github.com/user-attachments/assets/45ffee83-900b-4a13-9422-5acbfda894d8) |        
|**Constellation**  |**EVM Subcarriers**  |
|:--:|:--:|
| ![Image 1](https://github.com/user-attachments/assets/31fadb37-8434-4b13-83e9-0de59bbe2cc1) | ![Image 2](https://github.com/user-attachments/assets/29ca5412-34d0-47fb-854f-8d402931c76d) |     
|**Flatness**       |**Spectral Mask**    |
| ![Image 3](https://github.com/user-attachments/assets/20321d3a-8c05-4984-bde5-eb89485b9323) | ![Image 4](https://github.com/user-attachments/assets/bea3c6fa-b334-41d8-a74d-6ea9a8f30524) |    

## Simulation Results
### Simulate PA wi/wo DPD vs CFR
|**EVM**  |**ACLR**  |
|:--:|:--:|
| ![Image 1](https://github.com/user-attachments/assets/b38d8074-34cd-4f3d-a205-c8f97daf3fd7) | ![Image 2](https://github.com/user-attachments/assets/969925b3-fb08-4651-953c-ca68e20e568d) |    

### Simulate the conditions for reusing DPD coefficients
<img src="https://github.com/user-attachments/assets/6595bc4f-396a-4f90-9f3a-9be36617cee4" width="60%">

### Simulate PA wi/wo DPD vs CFR, conscidering flatness    
|**EVM**  |**ACLR**  |
|:--:|:--:|
| ![Image 1](https://github.com/user-attachments/assets/b6258799-b35f-426b-9ed7-082b4d1f5b51) | ![Image 2](https://github.com/user-attachments/assets/1bcef2b3-ac60-4d6a-8fea-3dbd1053a84d) | 
    
- Demodulation of HT20 at 27dBm with flatness 2dB    

|**Constellation**  |**EVM Subcarriers**  |
|:--:|:--:|
| ![Image 1](https://github.com/user-attachments/assets/6407092e-d0f0-4628-a3fa-f8222be40194) | ![Image 2](https://github.com/user-attachments/assets/e9b0cfd6-1f9e-41e3-a18e-5c3440b96df9) |     
|**Flatness**       |**Spectral Mask**    |
| ![Image 3](https://github.com/user-attachments/assets/a8f98161-c3d7-42c0-95b4-2bdddfa703ab) | ![Image 4](https://github.com/user-attachments/assets/c0deadef-f35d-4720-89bb-40067c2ce708) |    

- Demodulation of EHT320 at 15dBm with flatness 2dB        
    
|**Constellation**  |**EVM Subcarriers**  |
|:--:|:--:|
| ![Image 1](https://github.com/user-attachments/assets/afdfb21d-e06d-47a1-8eca-f9f1c96abeff) | ![Image 2](https://github.com/user-attachments/assets/0944c7e9-e9f4-4c35-b5ea-8f05dd1ea0b3) |     
|**Flatness**       |**Spectral Mask**    |
| ![Image 3](https://github.com/user-attachments/assets/31b82c5a-3da2-4161-8fe4-285b45523e80) | ![Image 4](https://github.com/user-attachments/assets/bffccd54-3b86-464e-9161-9d6f6269fd03) |    

