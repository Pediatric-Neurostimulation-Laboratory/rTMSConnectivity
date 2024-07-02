# rTMSConnectivity

This repository contains MATLAB scripts designed to analyze the impact of 1 Hz repetitive transcranial magnetic stimulation (rTMS) on TMS-evoked potentials (TEP) and weighted phase lag index (wPLI) connectivity.

Author: Xiwei She, Ph.D.

Usage: The scripts are named and indexed to be run sequentially. For instance, to obtain all wPLI results comparing pre- and post-rTMS under both real and sham conditions, execute the four scripts named rTMSEffects_1_*_wPLIConnectivity_*.m in order.

Dependencies:
Some scripts call functions located in the toolbox folder. Due to GitHub's file upload limitations, you need to download the eeglab and fieldtrip toolboxes separately from their official websites:
EEGLAB: https://sccn.ucsd.edu/eeglab/index.php
FieldTrip: https://www.fieldtriptoolbox.org/
After downloading, add them to your MATLAB path. The scripts include the necessary lines to facilitate this.

An example dataset is provided in the exampleData folder. This dataset contains four .set files corresponding to pre- and post-real-rTMS, and pre- and post-sham-rTMS conditions. Use this dataset to test the wPLI and TEP analysis methodologies described in the associated paper. While fully de-identified, this dataset is intended for testing purposes only. Please do not distribute it without the authors' permission.

Please contact xiweishe@stanford.edu or fbaumer@stanford.edu for any questions
