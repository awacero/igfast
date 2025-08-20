# GFAST_influxDB
A Python version of GFAST that utilizes influxDB for the data handling, based upon the database installed at IG-EPN in Ecuador.

Real-time Mode:

To run an event that is in the real-time database, you first need to check that gfast.props has the correct database set. From here, you will use the command line to enter the values specific to the earthquake. The order of variables is latitude, longitude, depth, origin time, earthquake name, seconds to model, and model style. The origin time, in UTC, needs to be in the following format:

"YYYY-MM-DD HH:MN:SC"

For model style, you will use 0 for the real-time database. A value of 1 is used for running off of the archive database. If using the archive database, an additional variable is needed, the channel file, since the coordinate information is not stored in the database.

Below is an example call to run data from the real-time database:

python3 GFAST_run.py 1.084 -79.532 18.0 '2025-04-25 11:44:52' esmeraldas_20250425 180 0

For the archive database, we currently only have the Pedernales earthquake available. To run that, here is an example call:

python3 GFAST_run.py -0.29 -80.18 20.6 '2016-04-16 23:58:36' pedernales_cent_est 120 1 'Ecuador2016_disp_pgd_v2.chan'

