This started out as a quick perl script to play with Covid numbers
The perl script is left here for interest - but (for reasons of speed) its been migrated to c++17.
(NB, much of the conversion was done simply by reusing the same code, so it's fairly hideous code to read!)

The simulation is current setup to play with the French and UK numbers from e.g. http://www.data.gouv.fr/
Numbers for France and UK are provided.

The script also takes into account mobility data as supplied by Google https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv

The core of the script algorithm runs accross each day using an adjusted gaussian distribution to model reported cases.
Each case generates new cases at the R0 value, but distributed with a mean and variance accross the days. 

These mean and variance numbers have a HUGE impact on the results. 

10% of Cases are then finally recorded (which seems to be what has been reported), some days after infection.

Cases are under-reported during weekends.

Hence the factors that are 'adjutable' are

- Base R0 for virus,
- Mobility multiply,
  - By how much to adjust the moblity data to effect R0.
- Daily Imported cases,
  - (starts the epidemic)
- Infectious day mean,
  - Day on which a patient is most infections
- Infectious day variance,
- Case reporting delay mean,
  - Day on which a case is most likely reported
- Case reporting delay variance,
- Social Distancing effect,
  - Effect on R0 of social distancing etc.
- Social Distance Inforduction day
  - The mean day on which social distancing is introduced
- Social Distancing Introduction day variance
- Start day

A modified simulated annealing approach is then used to find the best match using these variables.

The result shows a relatively close match to the actual data.

France:
<img src="./output-fr.svg">

UK:
<img src="./output-uk.svg">

To run the simulation first compile with your preffered compiler, e.g.
```bash
  > c++ -O3 -std=c++17 covid.c++ -o covid
```

The program can be run in a few ways:
``bash
 > ./covid country [-c][-o file.csv][-n runs][-q] [factors]
```
To find the best fit, simply run and wait ( a long time ) selecting the country you wish. This will use the supplied data for that country in mobility-<country>.csv and <country>.csv
```bash
  > ./covid france
```
  
To generate a csv file suitable for import into a spreadsheet for a specific set of numbers, run it thus:
```bash
  > ./covid france -o numbers-fr.csv 2.4207148 0.805779693 1.02442175 3.76362184 5 5 12.3217681 1.50674125 0.695714995 71.3405724 21.2957611 
```

To restart the annealing process you can add the '-c' flag, e.g.:
```bash
  > ./covid france -c 2.4207148 0.805779693 1.02442175 3.76362184 5 5 12.3217681 1.50674125 0.695714995 71.3405724 21.2957611
```
This can be helpful to start the process from a reasonabl position, given new data. It will also not only output the best guess, but also the 10 most promising alternatives.

To set the number of runs to try, use -n <number> (default 50000)

If you hit ^C during a run, then the program will finish the current runs and output the current best result.


