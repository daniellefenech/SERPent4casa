# SERPent4casa
A version of the <a href="https://github.com/daniellefenech/SERPent"> SERPent</a> RFI mitigation software for use with measurement set data.


This was written for calibrating e-MERLIN observations, however is suitable for use with any radio interferometric array data.

With the increase in bandwidth for e-MERLIN observations, the detection rate of Radio Frequency Interference (RFI) has dramatically increased, particularly at L-band (1.4-1.8 GHz) where the spectrum is actively used for other purposes.

In order to make the maximum use of the available bandwidth, the application of RFI mitigation techniques is necessary to remove the affected data.

The SERPent4casa software is written in python making use of the ms toolkit for CASA to interact directly with measurement set data. It utilises the SumThreshold flagging algorithms described in the articles "Post-correlation radio frequency interference classification methods'" (Offringa et al. 2010, MNRAS, 405, 155-167) and "A morphological algorithm for improving radio-frequency interference detection" (Offringa et al. 2012, A&A, 539, A95). The SumThreshold algorithm has been adapted to prevent over-flagging with the use of a sigma threshold.

In addition to the SumThreshold flagging, this software will optionally run Lovell-dropout and Zero-Level-dropout flagging. The first deals with the slow slew-time of the Lovell telescope and therefore the non-useable data recorded for every other phase calibration scan for Lovell baselines. The zero-level-dropout removes sections of data in time where the amplitudes drop to zero. This can occur for a number of reasons including when data is being recorded for a telescope before it is pointing at the target source.


## Requirements

CASA, python2.4+, numpy

## Usage

1. Download both SERPent4casa.py and SERPent4casa_input.py
2. Make the relevant changes to the input file 
3. Run with `casa -c SERPent4casa.py`
