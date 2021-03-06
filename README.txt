# shunt

## Purpose
This script calculates effective shunt for a given set of inputs from an arterial blood gas result (a physiological measurement often acquired by clinicians caring for sick patients in a hospital setting).

## Effective shunt
Effective shunt (ES) is an estimate of the proportion of cardiac output that would have to bypass the lungs entirely to obtain these blood gas results, at steady state. Obviously this is a rare event in clinical practice, but the measure provides an intuitively-interpretable result. This is analogous to an estimated glomerular filtration rate (eGFR) for the oxygenation function of the lungs.

## Quickstart

Open a terminal window and run this command:

python shunt.py -fio2 0.9 -pao2 8.3 -paco2 6.7 -pH 7.45 -gasunit kPa

## Inputs
Required inputs are:
- fio2
- pao2
- paco2
- pH

- Optional inputs are:
- gasunit=kPa
- acidunit=pH
- thisHb=8
- thisTemp=309.65
- thisVO2=0.25
- thisQ=5
- thismaxOER=0.8
- thisRER=0.8
- thisDPG=0.00465
