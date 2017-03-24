classdef Test_CochainResolveBuilder < matlab.unittest.TestCase
	%% TEST_COCHAINRESOLVEBUILDER 

	%  Usage:  >> results = run(mlfourdfp_unittest.Test_CochainResolveBuilder)
 	%          >> result  = run(mlfourdfp_unittest.Test_CochainResolveBuilder, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 07-Feb-2017 22:21:19
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/test/+mlfourdfp_unittest.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.  Copyright 2017 John Joowon Lee.
 	

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
	end

 	methods (TestClassSetup)
		function setupCochainResolveBuilder(this)
 			import mlfourdfp.*;
 			this.testObj_ = CochainResolveBuilder;
 		end
	end

 	methods (TestMethodSetup)
		function setupCochainResolveBuilderTest(this)
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

