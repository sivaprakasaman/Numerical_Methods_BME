%%% ===========Instruction (running under windows, Mac or Linux shall be easier, not tested)================%%%%
%%% 1. Need latest python version (3.9)(from official website)
%%%  Check with this command line code (windows): py --version
%%% 2. Install the following python libraries: pandas, numpy
%%% 3. You are good to use the code now!
%%% =========== Brief description ===========
%%% The main simulator file is called "main.py".
%%% Inputs: 1. BH: number of beds allocated to Hospital inpatient (suggestion: >=10)
%%%         2. BU: number of beds allocated to ICU (suggestion: >=10)
%%%         3. seed (each run needs a different integer number, with a fixed seed, the simulator will give the same output on the same machine while fixing other inputs) 
%%% Output: cumulative number of deaths
%%% Example:
r = python('main.py','40','70','5');% BH:'40', BU:'70', seed: '5'
%%% output a number, takes 2-4 secs per run.