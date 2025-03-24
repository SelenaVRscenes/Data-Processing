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
