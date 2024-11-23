# **EEG in Sleep**

## **Objective**

The goal of this project is to analyze and process an EEG recording from a sleeping patient to detect and highlight different sleep stages, including:
- **REM sleep**
- **NREM sleep**
- **Awake state**

The main code is `BSP_Team02Topic01_Code.m`

---

## **Data Overview**

The project is based on a one-night, single-lead EEG recording stored in the file `data_eeg.mat`:
- **Sampling frequency**: 100 Hz
- **EEG placement**: Forehead
- The data contains approximately **6 hours** of EEG signals recorded during sleep.

---

## **Methodology**

### **1. Visualization of Raw EEG Data**
The raw EEG signal is plotted to understand its structure and identify potential artifacts.

### **2. Preprocessing**
- **Artifact removal**:
  - Application of stop-band filters to remove artifacts at 0.99 Hz and 2 Hz (`coefficients_SOS_G.mat`& `coefficients_SOS2_G2.mat`).
- **Band-pass filtering [0.5 - 30 Hz]**:
  - Isolation of relevant EEG frequency bands to focus on sleep-related activities without noise like baseline drift or eye blink. 
  - This range excludes gamma waves (>30 Hz), as they do not contribute to sleep analysis and helps reduce potential noise at higher frequencies.

### **3. Spectral Analysis**
- Computation of the Fast Fourier Transform (FFT) to analyze frequency components.
- Calculation of power distribution over 30-second epochs to study changes in frequency bands.

### **4. Compressed Spectral Array (CSA)**
- The **Compressed Spectral Array (CSA)** is generated to visualize the evolution of spectral power over time.  
- This representation displays:
  - Spectral power densities for each frequency band.
  - Temporal changes in power levels, helping identify transitions between sleep stages.

### **5. Sleep Stage Detection**
- Frequency bands analyzed:
  - **Delta (δ)**: 1–4 Hz (Deep sleep).
  - **Theta (θ)**: 4–7 Hz (Light sleep).
  - **Alpha (α)**: 7–10 Hz (Relaxation).
  - **Beta (β)**: 10–30 Hz (Wakefulness).
  - **Gamma (γ)**: >30 Hz (Cognitive processing and high-level brain activity). _-not taken into account for sleep stages detection_
- Power levels in these bands are compared to predefined thresholds to classify sleep stages (REM, NREM, awake).

### **6. Hypnogram Generation**
- An **hypnogram** is created to give an overview of sleep architecture, highlighting transitions between stages over time.

### **7. Data Normalization**
- A manual **Box-Cox transformation** is applied to normalize the data and ensure consistency across epochs.

---

## **Tools, Language, and Libraries**

- **Language**: MATLAB
- **Libraries and Functions Used**:
  - **Signal Processing Toolbox**:
    - `filterDesigner`: For designing stop-band and band-pass filters.
    - `filter`: For applying the filters to the EEG signal.
  - **MATLAB Math Functions**:
    - `fft`: For performing the Fast Fourier Transform.
    - `mean`, `std`: For statistical analysis of signal epochs.
  - **Visualization**:
    - `plot`: For plotting raw and filtered signals.
    - `imagesc`: For generating the Compressed Spectral Array (CSA).
---

## **Outputs**

### **1. Plots**
- **Raw EEG signal**:
  - Visualization of the unprocessed EEG data to observe its structure and identify artifacts.
- **Filtered EEG signal**:
  - Results after artifact removal (stop-band filters) and band-pass filtering.
- **Spectral Power Distribution**:
  - Power distribution across different frequency bands computed over 30-second epochs.

### **2. Sleep Stage Detection**
- **Classification**:
  - REM, NREM, and awake stages are identified based on frequency eeg-band power.
- **Representation**:
  - Sleep stages are detected for each epoch of the EEG signal, providing a timeline of the sleep phases.
- **Hypnogram**:
  - A graphical representation showing the progression of sleep stages throughout the night.
  - The x-axis represents time (epochs), and the y-axis indicates the detected sleep stage (e.g., REM, NREM, awake).

---

## **Results**

The script provides:
- **Visualization**:
  - A clear view of raw and filtered EEG signals.
- **Spectral analysis**:
  - Insights into power distribution across frequency bands over time.
- **Sleep stage classification**:
  - Identification of REM, NREM, and awake states based on spectral power thresholds.
- **Representation of the sleep stages during the night**
  - A hypnogram is generated, showing the progression of sleep stages throughout the night.

---


