# N2O2NOeON
8-Component High Temp Air Equilibrium Properties

This code was originally developed for the MAE 250f - Hypersonic Flow graduate course at UCLA. See the report vanlew_report.pdf in the  main directory for a complete description of the formulation (and limitations) of this code.

## Introduction
We start from the partition functions to compute equilibrium thermodynamic properties for air as a function of pressure and temperature. We will assume that the air is a mixture of N 2 , and O 2 at low temperature (using the approximate mole ratio of 79% to 21%). At high temperature, we can consider the air is a mixture of the following species: N 2 , O 2 , NO, N, O, N + , O + , and e − . The objective is to write a computer program to compute mixture’s specific enthalpy, h, internal energy, e, and entropy, s, as well as the mass fractions of the species, c i , compressibility factor, Z, and total density, ρ. In addition we shall find the specific heat ratio, γ, and speed of sound, a, of the gas mixture as functions of pressure, p, and temperature, T. First the code will be validated through direct comparison of computed results with established results from literature. Following validation, more results will be given to prompt discussion of the thermodynamic properties at high temperatures.
