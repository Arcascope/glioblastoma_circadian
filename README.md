# Circadian rhythms in Glioblastoma

## Introduction

This is code to accompany the paper "Tracking Daily Variations in Rest-Wake to Guide Personalized Timing of Temozolomide for High-Grade Glioma Patients" for the Journal of Biological Rhythms by Gonzalez-Aponte et al.

The story of that paper is that even when people are highly consistent with taking their pills at a prescribed wall clock time, their biological time can still be (and often is) highly variable. This suggests that tracking personal biological time and not just wall clock time could be valuable for chronomedicine. 

## How to use this
- `run_preprocessing.py` prepares the data for analysis (e.g., sorting out time zones)
- `run_analysis.py` runs the circadian model on the data
- `make_all_figures.py` calls all the functions needed to make the figures in the paper

Questions? Reach out to the senior authors, Erik Herzog and Olivia Walch.