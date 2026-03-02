# Dimuon Invariant Mass Analysis

## Overview
This project contains analysis codes to study the dimuon mass distribution change of $€omega$ and $€phi$ mesons in the forward rapidity region ($-3.6 < y < -2.5$). 
The analysis is based on data from high-energy physics experiments and uses ROOT.

We separate the heavy data processing (C++) from the visual plotting (Python/PyROOT) to make the workflow faster and more flexible.

## Directory Structure
* `include/` : Header files (e.g., `DimuonAnalyzer.h`).
* `src/` : C++ source codes for data processing and fitting.
* `scripts/` : Python scripts using PyROOT to draw and save plots.

## Prerequisites
* ROOT (v6.24 or later recommended)
* C++ compiler supporting C++11 or later (g++, clang, etc.)
* Python 3 with PyROOT enabled

## How to Compile
Use `g++` with `root-config` to compile the C++ core logic. Run the following command in the `myAnalysis` directory:

```bash
g++ src/run_analysis.cxx src/DimuonAnalyzer.cxx -o run_analysis $(root-config --cflags --libs) -I./include/
