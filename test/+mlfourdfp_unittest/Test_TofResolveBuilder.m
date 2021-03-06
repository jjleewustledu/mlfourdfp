classdef Test_TofResolveBuilder < matlab.unittest.TestCase
	%% TEST_TOFRESOLVEBUILDER 

	%  Usage:  >> results = run(mlfourdfp_unittest.Test_TofResolveBuilder)
 	%          >> result  = run(mlfourdfp_unittest.Test_TofResolveBuilder, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 30-Dec-2020 00:37:06 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/test/+mlfourdfp_unittest.
 	%% It was developed on Matlab 9.9.0.1538559 (R2020b) Update 3 for MACI64.  Copyright 2020 John Joowon Lee.
 	
	properties
 		registry
        sesd
 		testObj
 	end

	methods (Test)
		function test_afun(this)
 			import mlfourdfp.*;
 			this.assumeEqual(1,1);
 			this.verifyEqual(1,1);
 			this.assertEqual(1,1);
        end
        function test_ctor(this)
            disp(this.testObj_)
        end
        function test_resolve(this)
            t4rb = this.testObj_.resolve();
            disp(t4rb)
        end
        function test_t4img_on_atl(this)
            this.testObj_.t4img_on_atl();
        end
	end

 	methods (TestClassSetup)
		function setupTofResolveBuilder(this)
 			import mlfourdfp.*;
            this.sesd = mlraichle.SessionData.create( ...
                'CCIR_00559/ses-E03056/OC_DT20190523112618.000000-Converted-AC');
 			this.testObj_ = TofResolveBuilder(this.sesd);
 		end
	end

 	methods (TestMethodSetup)
		function setupTofResolveBuilderTest(this)
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

