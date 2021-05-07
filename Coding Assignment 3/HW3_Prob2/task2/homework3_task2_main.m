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
%%%                  type 2. seed: you do not need to do anything with this
%%%                  argument, just use the following code template to
%%%                  execute the function
%%% Output: cumulative(3-day)number of hospital inpatients and ICU
%%%         30*2 result: 
%%%                 (column 1.cumulative (3-day sum) H(inpatient) counts
%%%                  column 2.cumulative(3-day sum) ICU patient counts) 
%%% Note:(Only first 3 parameters to be varied after 'main.py' below)
%%% 30 replications, 28 seconds

seed = table2array(readtable('seeds100.csv','ReadVariableNames',false));
output = zeros(30,2);

parfor i = 1:30
    outputTemp = python('main.py','3.15','-0.4','8',num2str(seed(i))); %e.g. l1:'3.15', l2:'-0.2', l3:'8' 
    output(i,:) = transpose(sscanf(outputTemp,'%f')); %
end










