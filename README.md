# **Data Processing for the <Exp 1 : Metro>. Analysis and Machine Learning**  

## **Overview**  
This repository contains Python scripts for data preprocessing, resampling, and analysis using various machine learning techniques, including XGBoost and deep learning models. The scripts handle raw data processing, time and space resampling, and feature extraction for model training and evaluation.  

## **Files and Functionalities**  

- **`DataPreprocVGG.py`**, **`DataPreprocXGBoost.py`**, **`DataPreprocDeepLearrn.py`** – Preprocess raw data for training machine learning models (VGG, XGBoost, Deep Learning). They include filtering, resampling, and feature extraction.  
- **`demo01_xgboost.py`** – Example script for training an XGBoost classifier on the breast cancer dataset.  
- **`interceptionPointP1P1.py`** – Calculates intersection points between different trajectory lines for data analysis.  
- **`SpaceTimeResampeMSdata.py`** – Resamples spatial and temporal data to a standard frame rate while applying oscillation filtering.  
- **`angleAnalysis.py`** – Computes and visualizes angle-related metrics from movement trajectories.  
- **`dataCleanupPolyReg.py`**, **`dataConvertAngleAnalysis.py`**, **`dataConvertLinRegression.py`** – Data transformation and regression-based cleanup scripts.  


# **Data Processing and Motion Analysis for < Exp2: Entropy analysis >**  

## **Overview**  
This repository contains Python scripts for preprocessing, resampling, and analyzing motion data. The scripts focus on **time and space resampling, filtering, feature extraction, and visualization**. These methods are particularly useful for **motion tracking experiments, biomechanics analysis, and machine learning applications**.  

---

## **Files and Functionalities**  

### **Preprocessing & Filtering Scripts**  
- **`dataLowPassFilter.py`** – Applies a **low-pass filter** (Butterworth filter) to smooth noisy motion data by removing high-frequency oscillations.  
- **`dataTimeResampling.py`** – Resamples **time-series motion data** to ensure a consistent frame rate, which is useful for synchronizing multiple datasets.  

### **Motion Trajectory Processing & Feature Extraction**  
- **`motionAnalysisVisExp2.py`** – Loads, processes, and visualizes motion trajectories from an experiment. Extracts movement patterns and generates plots to **analyze path characteristics**.  
- **`dataConvert.py`** – Converts raw **motion tracking data (CSV format)** into a structured format for further analysis, including coordinate transformations and basic filtering.  
- **`dataConvertExp2.py`** – A specialized version of `dataConvert.py` for an **experiment on motion randomness**, focusing on trajectory deviations and movement variability.  

### **Spatiotemporal Resampling & Analysis**  
- **`SpaceTimeResampeMSdata.py`** – Performs both **time resampling** (for synchronization) and **spatial resampling** (to normalize trajectory length), making motion data consistent across samples.  
- **`TimeSpaceResampeMSdata_visOnlyEpx2_v2.py`** – A visualization-optimized version of the above, designed for quickly displaying results with **trajectory overlays and movement plots**.  
- **`TimeSpaceResampeMSdata_visOnlyEpx2-7-2.py`** – Main data processing script for **Experiments conducted on June 13-15, 2023**, focusing on:  
  - **Time & space resampling**  
  - **Filtering out oscillations**  
  - **Extracting key movement metrics** (turning points, area covered, trajectory length, and time to complete tasks)  




## **Setup & Dependencies**  
Ensure you have the required Python libraries installed:  
```sh
pip install numpy pandas matplotlib scipy sklearn xgboost pwlf
```

## **Usage**  
Run any script as follows:  
```sh
python script_name.py
```
Replace `script_name.py` with the desired file.

