%{
Runs the tests defined in test_multicolor19F
%}


fprintf('Running automated tests for multicolor 19F app\n')

% find tests path and set this to current directory
cd(fileparts(which('runtests_19F')))

% run all tests defined in test_multicolor19F
runtests('test_multicolor19F')

% run only one test like this: 
% runtests('test_multicolor19F/test_load_data_19F')