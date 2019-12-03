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
            testCase.press(app.GPUCheckBox);
        end
         
        function test_load_data_19F(testCase)
            fprintf('\nTest data loading\n')
            app = multicolor19F;
            testCase.addTeardown(@delete,app);
            testCase.press(app.UseTestDataCheckBox);
            testCase.press(app.Load19FDataButton);
            assert(app.valid19Ffile)
            pause(1)
        end
        
        function test_load_data_1H(testCase)
            fprintf('\nTest data loading\n')
            app = multicolor19F;
            testCase.addTeardown(@delete,app);
            testCase.press(app.UseTestDataCheckBox);
            testCase.press(app.Load1HDataButton);
            assert(app.valid1Hfile)
            pause(1)
        end
        
        function test_CS_recon(testCase)
            fprintf('\nTest CS data loading and recon\n')
            app = multicolor19F;
            testCase.addTeardown(@delete,app);
            testCase.press(app.UseTestDataCheckBox);
            app.data_import_path_test_19F = 'tests\testdata\MPCS\'
            testCase.press(app.Load19FDataButton);
            assert(app.valid19Ffile)
            testCase.press(app.TransCorrCheckBox); %turn off registration
            testCase.press(app.PseudoInverseButton);
            app.IterationsEditField.Value = 2; 
            pause(2)
            testCase.press(app.DeconvolutionButton);           

            
        end
        
        
    end
    
end
