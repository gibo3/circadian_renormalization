# Waveform distortion for temperature compensation and synchronization in circadian rhythms: An approach based on the renormalization group method

This repository contains MATLAB code for simulations and data analyses from the manuscript 
"Waveform distortion for temperature compensation and synchronization in circadian rhythms: An approach based on the renormalization group method", 
by Shingo Gibo, Teiji Kunihiro, Tetsuo Hatsuda, Gen Kurosawa.

- **Goodwin_Hill.m** MATLAB code for simulation of Goodwin model. The transcriptional-translational function is Hill function $f(x_{3})=r/(1+(x_{3}/K)^{n})$.
- **Goodwin_Hill_with_ExternalForce.m** MATLAB code for simulation of Goodwin model with external force. The transcriptional-translational function is Hill function $f(x_{3})=r/(1+(x_{3}/K)^{n})$.
- **Goodwin.m** MATLAB code for simulation of Goodwin model. The transcriptional-translational function is $f(x_{3})=r/x_{3}^{n})$.
- **LotkaVolterramodel.m** MATLAB code for simulation of Lotka-Volterra model.
- **vanderPolmodel.m** MATLAB code for simulation of van der Pol model.
- **Relogio2011model.m** MATLAB code for simulation of Relogio (2011) model.
- **GHA.m** MATLAB code for generalized harmonic analyis (GHA).
- **NS_experimental_data.m** MATLAB code for estimating $NS$ value of experimental time-series data. Running this code quantifies the $NS$ of time series created by adding noise to 'timeseries_example.csv' as test data, and displays the distribution of $NS$.
- **timeseries_example.csv** Test data for NS_experimental_data.m.
