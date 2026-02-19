# HIV Prevalence Modeling in Kenya – MATLAB Capstone Project

This repository contains the MATLAB code for a mathematical modeling project that simulates HIV prevalence in Kenya. The project uses a system of differential equations and optimization methods (fminsearch and GA) to fit the model to historical data.  

> **Disclaimer:** This project was initially developed in MATLAB and is **showcased here for reference only**. The code may not run outside MATLAB or without the appropriate toolboxes. This repository is intended for **demonstration and educational purposes**, not for execution.

---

## **Project Files and Descriptions**

| File | Description |
|------|-------------|
| `model4_function.m` | Defines the main system of differential equations for the epidemic model. Calculates the dynamics of susceptible, infected, and other compartments. |
| `model4_prevalence_only.m` | Solves the ODE system over time to simulate HIV prevalence based on model parameters. |
| `cost_function0.m` | Implements the cost function used for parameter fitting with fminsearch and GA. Computes the error between simulated prevalence and real data. |
| `main_model4_compare.m` | Main script that runs the model and fits parameters using `fminsearch`. Generates plots comparing simulated prevalence with observed data. |
| `main_model4_compare_GA.m` | Main script that runs the model and fits parameters using Genetic Algorithm (GA). Generates plots comparing simulated prevalence with observed data. |
| `gaEarlyStop.m` | Helper function used by the GA to implement early stopping criteria. |
| `fminsearch_calibration.png`  | Example output plot showing observed vs simulated HIV prevalence. |
| `GA_calibration.png`  | Example output plot showing observed vs simulated HIV prevalence. |
---

## **How to Use**

- These files are provided **as-is for showcase purposes**.  
- You can browse the MATLAB code to understand the model structure, differential equations, and fitting methodology.  
- Plots and results can be seen from saved images or by opening MATLAB if you have access.  

---

## **About the Project**

- Developed as a **capstone project** in Applied Mathematics.  
- Combines **ODE-based epidemic modeling** and **parameter fitting**.  
- Focused on HIV prevalence trends in Kenya from 1950–2022.  
- Demonstrates skills in **mathematical modeling, numerical simulation, and optimization**.

---