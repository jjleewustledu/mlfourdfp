classdef Test_SessionResolvingData < matlab.unittest.TestCase
	%% TEST_SESSIONRESOLVINGDATA 

	%  Usage:  >> results = run(mlfourdfp_unittest.Test_SessionResolvingData)
 	%          >> result  = run(mlfourdfp_unittest.Test_SessionResolvingData, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 27-May-2018 15:05:39 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/test/+mlfourdfp_unittest.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties
 		registry
 		testObj
 	end

	methods (Test)
		function test_afun(this)
 			import mlfourdfp.*;
 			this.assumeEqual(1,1);
 			this.verifyEqual(1,1);
 			this.assertEqual(1,1);
        end
        function test_compositeT4ResolveBuilderBlurArg(this)
        end
        function test_maxLengthEpoch(this)
        end
        function test_supEpoch(this)
        end
        function test_t4ResolveBuilderBlurArg(this)
        end       
	end

 	methods (TestClassSetup)
		function setupSessionResolvingData(this)
 			import mlfourdfp.*;
 			this.testObj_ = SessionResolvingData;
 		end
	end

 	methods (TestMethodSetup)
		function setupSessionResolvingDataTest(this)
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

