classdef Test_InnerFourdfp < matlab.unittest.TestCase
	%% TEST_INNERFOURDFP 

	%  Usage:  >> results = run(mlfourdfp_unittest.Test_InnerFourdfp)
 	%          >> result  = run(mlfourdfp_unittest.Test_InnerFourdfp, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 27-Jul-2018 23:02:22 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/test/+mlfourdfp_unittest.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties
        pwd0
        ref
 		registry
 		testObj
 	end
    
    properties (Dependent)
        dataroot
        TmpDir
    end

    methods
        
        %% GET/SET
        
        function g = get.dataroot(~)
            g = fullfile(getenv('HOME'), 'MATLAB-Drive', 'mlfourd', 'data', '');
        end
        function g = get.TmpDir(~)
            g = fullfile(getenv('HOME'), 'Tmp', '');
        end
    end

	methods (Test)
		function test_afun(this)
 			import mlfourdfp.*;
 			this.assumeEqual(1,1);
 			this.verifyEqual(1,1);
 			this.assertEqual(1,1);
        end
        function test_ctor(this)
        end
        function test_save(this)
        end
        function test_ifh(this)
        end
        function test_imgreg(this)
        end
	end

 	methods (TestClassSetup)
		function setupInnerFourdfp(this)
 			import mlfourdfp.*;
            this.ref = mlfourd.ReferenceMprage;
            this.ref.copyfiles(this.TmpDir);
 			this.testObj_ = InnerFourdfp(FourdfpInfo(this.ref.dicomAsFourdfp));
 		end
	end

 	methods (TestMethodSetup)
		function setupInnerFourdfpTest(this)
            this.pwd0 = pushd(this.TmpDir);
 			this.testObj = this.testObj_;
 			this.addTeardown(@this.cleanTestMethod);
 		end
	end

	properties (Access = private)
 		testObj_
 	end

	methods (Access = private)
		function cleanTestMethod(this)
            popd(this.pwd0);
 		end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

