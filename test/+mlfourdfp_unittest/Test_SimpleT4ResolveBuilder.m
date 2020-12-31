classdef Test_SimpleT4ResolveBuilder < matlab.unittest.TestCase
	%% TEST_SIMPLET4RESOLVEBUILDER 

	%  Usage:  >> results = run(mlfourdfp_unittest.Test_SimpleT4ResolveBuilder)
 	%          >> result  = run(mlfourdfp_unittest.Test_SimpleT4ResolveBuilder, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 30-Dec-2020 13:11:39 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/test/+mlfourdfp_unittest.
 	%% It was developed on Matlab 9.9.0.1538559 (R2020b) Update 3 for MACI64.  Copyright 2020 John Joowon Lee.
 	
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
		function setupSimpleT4ResolveBuilder(this)
 			import mlfourdfp.*;
 			this.testObj_ = SimpleT4ResolveBuilder;
 		end
	end

 	methods (TestMethodSetup)
		function setupSimpleT4ResolveBuilderTest(this)
 			this.testObj = this.testObj_;
 			this.addTeardown(@this.cleanTestMethod);
 		end
	end

	properties (Access = private)
 		testObj_
 	end

	methods (Access = private)
		function cleanTestMethod(this)
 		end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

