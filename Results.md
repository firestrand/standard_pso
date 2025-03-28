# A few results. For more explanation, see ReadMe.txt and the source code

# SPSO 2007

# KISS, seed=1294404794

# 100 runs

|Function|Dim.|Comment|S|% Mean best|% Mean best| | |
|---|---|---|---|---|---|---|---|
|4 Tripod|2| |12|56|5.03E-01|63|3.10E-01|
|11 Network|42|Partly binary|22|0|1.35E+02|0|1.06E+02|
|15 Step|10|Biased|16|99|1.00E-02|3|4.53E+00|
|17 Lennard-Jones|18|6 atoms|18|4|4.26E-01|3|6.40E-01|
|18 Gear train|4|Discrete|14|9|1.55E-09|16|2.47E-10|
|20 Perm|5|Discrete|14|16|5.10E+02|46|2.92E+02|
|21 Compression spring|3|Partly discrete|13|31|3.96E-02|72|1.91E-03|
|100 Sphere|30|Shifted|20|100|9.39E-07|100|9.00E-07|
|102 Rosenbrock|10|Shifted|16|68|5.75E+00|9|1.81E+00|
|103 Rastrigin|30|Shifted|20|0|5.39E+01|0|3.89E+01|
|104 Schwefel|10|Shifted|16|100|9.09E-05|100|8.57E-05|
|105 Griewank|10|Shifted|16|5|5.26E-02|18|3.05E-02|
|106 Ackley|30|Shifted|20|30|1.12E+00|98|1.87E-02|
| | |Total|518|7.07E+02|528|4.44E+02| |

# SPSO 2011, swarm size 40

# Uniform radius, BW=(0,0,0,0), Confinement

# KISS, seed=1294404794

# Mersenne, seed=1294404794

# Mersenne, seed=1234567890

# 100 runs

|Function|100 runs|100 runs|1000 runs|100 runs| | | | |
|---|---|---|---|---|---|---|---|---|
|4 Tripod|79|1.46E-01|72|1.54E-01|75|1.39E-01|74|1.45E-01|
|11 Network|0|1.09E+02|0|1.12E+02|0|1.11E+02|0|1.10E+02|
|15 Step|99|1.00E-02|99|1.00E-02|99|1.10E-02|98|2.00E-02|
|17 Lennard-Jones|0|9.18E-01|0|9.93E-01|0|9.34E-01|0|8.42E-01|
|18 Gear train|58|1.90E-11|46|2.61E-11|49|2.57E-11|64|3.19E-11|
|20 Perm|36|3.09E+02|29|3.43E+02|33|3.03E+02|32|3.38E+02|
|21 Compression spring|81|3.26E-03|79|3.55E-03|77|3.25E-03|78|4.51E-03|
|100 Sphere|100|0.00E+00|100|0.00E+00|100|0.00E+00|100|0.00E+00|
|102 Rosenbrock|50|5.77E+01|46|5.95E+01|45|6.62E+01|44|6.02E+01|
|103 Rastrigin|1|5.39E+00|0|5.26E+00|0|5.23E+00|0|5.24E+00|
|104 Schwefel|100|0.00E+00|100|0.00E+00|100|0.00E+00|100|0.00E+00|
|105 Griewank|9|2.15E-02|21|1.65E-02|12|2.07E-02|15|2.13E-02|
|106 Ackley|100|0.00E+00|100|0.00E+00|100|1.16E-03|100|0.00E+00|
|Total|713|4.82E+02|692|5.21E+02|690|4.86E+02|705|5.15E+02|

100 runs is not always enough for a good estimation of the success rate (and of the mean best). The result may depend on the RNG. Actually, even for a given RNG, it may depend on the seed..