classdef Test_O15UmapResolveBuilder < matlab.unittest.TestCase
	%% TEST_O15UMAPRESOLVEBUILDER 

	%  Usage:  >> results = run(mlfourdfp_unittest.Test_O15UmapResolveBuilder)
 	%          >> result  = run(mlfourdfp_unittest.Test_O15UmapResolveBuilder, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 24-Oct-2016 21:03:19
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/test/+mlfourdfp_unittest.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.
 	

	properties
 		registry
        studyd
        sessd
 		testObj
 	end

	methods (Test)
        function test_ctor(this)
            this.verifyEqual(pwd, fullfile(getenv('PPG'),'jjlee', 'HYGLY09', 'V1', ''));
        end
        function test_buildO15NAC(this)
            this.testObj.buildO15NAC('OC', 1);
            this.testObj.buildO15NAC('HO', 1);
            this.testObj.buildO15NAC('OO', 1);
        end
        function test_convertNativeSumtTo4dfp(this)
        end
        function test_alignUmapToSumtResolved(this)
        end
        function test_convertUmapToE7Format(this)
        end
	end

 	methods (TestClassSetup)
		function setupO15UmapResolveBuilder(this)
            this.studyd = mlraichle.StudyData;
            this.sessd = mlraichle.SessionData( ...
                'studyData', this.studyd, ...
                'sessionPath', fullfile(this.studyd.subjectsDir, 'HYGLY09',  ''), ...
                'tracer', 'HO'); 
 			this.testObj_ = mlfourdfp.O15UmapResolveBuilder('sessionData', this.sessd);
 		end
	end

 	methods (TestMethodSetup)
		function setupO15UmapResolveBuilderTest(this)
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

