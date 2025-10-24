<p align="center">
  <a href="LICENSE">
    <img src="https://img.shields.io/badge/License-MIT-green.svg" alt="License: MIT">
  </a>
  <img src="https://img.shields.io/badge/R-4.2%2B-blue" alt="R Version">
</p>

# Deciphering Gut Microbial Impact in Coronary Artery Disease through Multimodal Computational Approaches

This paper presents an integrated computational framework for **gut microbiome analysis** in **Coronary Artery Disease (CAD)**.  
This repository contains reproducible **R code** implementing the full pipeline described in our manuscript:  

> *Deciphering Gut Microbial Impact in Coronary Artery Disease through a Multimodal Computational Study*  
> Alaeddini, Nemati, Gholizadeh (2025)

---

## 🧠 Overview

This framework combines **microbiome bioinformatics, stochastic modelling, and machine learning** to reveal CAD-associated microbial alterations.  

Key components:
- **Microbiome preprocessing & diversity analysis** (DADA2, phyloseq, vegan, DESeq2)  
- **Differential abundance testing**  
- **Stochastic logistic growth simulations** of microbial dynamics  
- **Supervised machine learning models** (Random Forest, SVM, XGBoost) for biomarker discovery  

<p align="center">
 <img width="975" height="658" alt="image" src="https://github.com/user-attachments/assets/d4f3900d-fe2a-409c-9e7a-e263afc9cfe0" />
</p>

---

## ⚙️ Key Features

✅ End-to-end pipeline fully implemented in **R**  
✅ **Taxonomic profiling** and **diversity metrics** (alpha & beta)  
✅ **Stochastic differential equation (SDE)** modelling of microbial dynamics  
✅ **Machine learning classifiers** for CAD risk prediction  
✅ Validation using **independent external dataset**  

---

## 📊 Machine Learning Models

We implemented three supervised models to classify **CAD vs. healthy controls** and top microbial genera identified as important features for distinguishing CAD patients from healthy individuals using the RF model:  
- 🌳 **Random Forest (RF)** – best performance (AUC ≈ 0.95)  
- ⚡ **XGBoost** – robust gradient boosting classifier  
- 📈 **Support Vector Machine (SVM, RBF kernel)**  

Evaluation metrics include:  
- Accuracy  
- Precision, Recall, F1-score  
- AUC (ROC analysis)  
- Cohen’s Kappa  

---

## 📦 Prerequisites

Please make sure you have **R >= 4.2** installed with the following R packages:

```r
install.packages(c("phyloseq", "vegan", "DESeq2", "igraph", "ggplot2", "e1071"))
install.packages(c("randomForest", "xgboost", "caret", "pROC"))
