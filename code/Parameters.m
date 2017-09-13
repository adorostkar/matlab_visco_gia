function pars = Parameters()
  % Debug variable
  % 1 - print debug texts
  % 0 - no debug
  pars.debug = 1;

  % % Variables
  % 0 - No additional output(figure and text)
  % 1 - Only figures
  % 2 - everything
  pars.verbose = 0;

  % output directory
  pars.outdir = './out/';

  % years to output the results
  pars.t_pr = [1,2,3,4,5,6,7,8,9,10,50,100,200,300,400,500,600,1000,[1500:500:10000]];

  pars.beta = 0.1;
  pars.no_domains = 2;
  pars.Emagn = 1;

  pars.flag_incr = 1; % ==1, incremental form

  % pars.test_problem = 0; % pure elasticity in unit square, benchmark
  % pars.test_problem = 9; % visco-elasticity in unit square, benchmark
  % pars.test_problem = 8; % visco-elasticity in unit square, benchmark
  pars.test_problem = 10; % visco-elasticity
