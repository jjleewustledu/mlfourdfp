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
        function test_resolveOnTof(this)
            cd(fullfile(getenv('SINGULARITY_HOME'), 'CCIR_00559_00754', 'derivatives', 'resolve', 'sub-S58163', 'anat', ''))
            bids = mlraichle.Ccir559754Bids();
            tof_b6 = bids.tof_ic.blurred(6);
            tof_bin = tof_b6.thresh(30);
            tof_bin = tof_bin.binarized();
            t1_b6 = bids.T1w_ic.blurred(6);
            t1_bin = t1_b6.thresh(10);
            t1_bin = t1_bin.binarized();
            masks = {tof_bin, t1_bin};
            images = {bids.tof_ic, bids.T1w_ic};
            rb = mlfourdfp.SimpleT4ResolveBuilder('maskForImages', masks, 'theImages', images);
            rb.resolve
        end
	end

 	methods (TestClassSetup)
		function setupSimpleT4ResolveBuilder(this)
 			import mlfourdfp.*;
 			this.testObj_ = []; % SimpleT4ResolveBuilder;
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

