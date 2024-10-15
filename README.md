# ic_WLAN_DPD
Simulation of WLAN coporate with DPD
## Main and History: 2024-10-15
- Example_WIFI_Transmit_240919_v1.mlx
- Example_PAModel_and_DPD_240926_v2.mlx
## Simulation flow
<img src="https://github.com/user-attachments/assets/09470ede-4011-461d-853e-14df9ced3d20" width="60%">

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
- 



