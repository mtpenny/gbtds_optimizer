# gbtds_optimizer
Tool for optimizing yields or metrics for the Roman Galactic Bulge
Time Domain Survey

The current script (optimizeSlew.py) just computes best path around fields and a rough
scaling of microlensing planet detection rates. Caution should be
taken when increasing cadence beyond 15 minutes.

To run:

`python optimizeSlew.py <fields> <slew-times (short axis)> {<slew-times (diagonal)> <slew-times (long-axis)>`

The most up-to-date slew time file is `slew_times_McEnery05232024.txt`

Fields files should have 3 columns with:
<Field_name> <l(deg)> <b(deg)>
