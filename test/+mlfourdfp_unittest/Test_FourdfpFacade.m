classdef Test_FourdfpFacade < matlab.unittest.TestCase
	%% TEST_FOURDFPFACADE 

	%  Usage:  >> results = run(mlfourdfp_unittest.Test_FourdfpFacade)
 	%          >> result  = run(mlfourdfp_unittest.Test_FourdfpFacade, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 10-Mar-2016 21:23:54
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/test/+mlfourdfp_unittest.
 	%% It was developed on Matlab 9.0.0.307022 (R2016a) Prerelease for MACI64.
 	

	properties
 		registry
 		testObj
        view = true
 	end

	methods (Test)
        function test_t4ResolveSubject(this)
            studyd = mlpipeline.StudyDataSingletons.instance('raichle');
            sessd = mlraichle.SessionData( ...
                'studyData', studyd, 'sessionPath', fullfile(studyd.subjectsDir, 'NP995_09', ''));            
            t4b = mlfourdfp.T4ResolveBuilder('sessionData', sessd);
            ff  = mlfourdfp.FourdfpFacade( 'sessionData', sessd, 't4ResolveBuilder', t4b);
            fprintf('Test_FourdfpFacade.test_t4ResolveSubject:  running t4ResolveSubject which may requires hours of processing time..........');
            ff  = ff.t4ResolveSubject;
            this.verifyTrue(~isempty(ff.product));
            
            if (this.view)
                product.fdg.view;
            end
        end
        function test_t4ResolveSubject2(this)
            studyd = mlpipeline.StudyDataSingletons.instance('test_arbelaez');
            sessd = mlarbelaez.SessionData( ...
                'studyData', studyd, 'sessionPath', fullfile(studyd.subjectsDir, 'p7991_JJL', ''));
            t4b = mlfourdfp.T4ResolveBuilder('sessionData', sessd);
            ff  = mlfourdfp.FourdfpFacade('sessionData', sessd, 't4ResolveBuilder', t4b);
            fprintf('Test_FourdfpFacade.test_t4ResolveSubject2:  running t4ResolveSubject which may requires hours of processing time..........');
            ff = ff.t4ResolveSubject;
            this.verifyTrue(~isempty(ff.product));
            
            if (this.view)
                product.gluc.view;
            end
        end
	end

 	methods (TestClassSetup)
		function setupFourdfpFacade(this)
 			import mlfourdfp.*;
 			this.testObj_ = FourdfpFacade;
 		end
	end

 	methods (TestMethodSetup)
		function setupFourdfpFacadeTest(this)
 			this.testObj = this.testObj_;
 			this.addTeardown(@this.cleanFiles);
 		end
	end

	properties (Access = private)
 		testObj_
 	end

	methods (Access = private)
		function cleanFiles(this)
 		end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

