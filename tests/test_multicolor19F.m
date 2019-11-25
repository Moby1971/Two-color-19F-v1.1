classdef test_multicolor19F < matlab.uitest.TestCase
    % test function for multicolor 19F
    %{
    Goal of this script is to collect (unit) tests for the multicolor app
    This can be used for 

    1. debugging and development without manually clicking buttons
    2. Verify that a new installation works
    
    TO DO: 
    * automate loading of data (testdata is in tests/testdata)
    * add way more test cases (recons with different settings etc.) 
    
    Jasper 25/11/'19 
    %}
    
    methods (Test)
        function test_GPU_checkbox(testCase)
            fprintf('test_GPU_checkbox\n')
            app = multicolor19F;                    %create app
            testCase.addTeardown(@delete,app);      %delete app after test
            testCase.press(app.GPUCheckBox)
        end
         
        function test_load_data(testCase)
            fprintf('test CS data loading and recon\n')
            app = multicolor19F;
            testCase.addTeardown(@delete,app);
            testCase.press(app.Load19FDataButton);
            assert(app.valid19Ffile)
            testCase.press(app.PseudoInverseButton);
            testCase.press(app.DeconvolutionButton);
            pause(3)

        end
    end
    
end
