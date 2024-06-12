# gbtds_optimizer
Tool for optimizing yields or metrics for the Roman Galactic Bulge
Time Domain Survey

The current script (optimizeSlew.py) just computes best path around fields and a rough
scaling of microlensing planet detection rates. Caution should be
taken when increasing cadence beyond 15 minutes.

To run:

`python optimizeSlew.py <fields> <slew-times (short axis)> {<slew-times (diagonal)> <slew-times (long-axis)>`

The most up-to-date slew time file is
`slew_times_withResetReference_McEnery05232024.txt` - this file is the
`slew_times_McEnery05232024.txt` file provided by Julie McEnery with
6.12 seconds added to every slew to account for the reset read cycle
(3.08 s) and the first reference read that is subtracted from the ramp
(3.04 s). If you are using sample up the ramp signal to noise
estimates the reference read is already accounted for in the exposure
time, so you can remove it from the overheads. 

Fields files should have 3 columns with:
<Field_name> <l(deg)> <b(deg)>
