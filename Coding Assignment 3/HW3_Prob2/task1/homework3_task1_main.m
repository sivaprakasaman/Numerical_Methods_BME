%%% ===========Instruction (running under windows, Mac or Linux shall be easier, not tested)================%%%%
%%% 1. Need latest python version (3.9/3.8)(from official website)
%%%  Check with this command line code (windows): py --version
%%% 2. Install the following python latest libraries: pandas, numpy.
%%% 3. use the python.m file that worked for you in homework 2.
%%% 4. Do not change the .py files and .csv files.
%%% 5. You are good to use the code now!
%%% =========== Brief description ===========
%%% The main simulator file is called "main.py".
%%% Inputs(2 types): type 1. l1, l2, l3: real parameters to be varied with
%%%                         domain (-infinity, infinity)
%%% Output: cumulative number of deaths
%%% Example:(first 3 parameters to be varied)
%%% 11 sec to run
outputTemp = python('main.py','1.5','-2.2','1'); %e.g. l1:'1.5', l2:'-2.2', l3:'1' 
output = sscanf(outputTemp,'%f'); %
%%% output two-line result (1.cumulative death counts 2.cumulative ICU patient counts)
%%% a set of initial states, the number is equivalent to the observations
%%% of the rolling sums.
%%% provide a set of 100 observations from raw data 14-day sum of H and ICU()

