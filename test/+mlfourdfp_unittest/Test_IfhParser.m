classdef Test_IfhParser < matlab.unittest.TestCase
	%% TEST_IFHPARSER 

	%  Usage:  >> results = run(mlfourdfp_unittest.Test_IfhParser)
 	%          >> result  = run(mlfourdfp_unittest.Test_IfhParser, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 02-Aug-2018 13:17:25 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/test/+mlfourdfp_unittest.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties
        noDelete = true
        pwd0
 		registry
 		testObj
 	end
    
    properties (Dependent)
        dataroot
        fileprefix
        TmpDir
    end

    methods
        
        %% GET/SET
        
        function g = get.dataroot(~)
            g = fullfile(getenv('HOME'), 'MATLAB-Drive', 'mlfourd', 'data', '');
        end
        function g = get.fileprefix(~)
            g = 'fdgv2ConvertedLaForest_dcm2niix';
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
        function test_constructDenovo(this)
            import mlfourdfp.*;
            nh = InnerFourdfp(FourdfpInfo([this.fileprefix '.4dfp.hdr']));
            tmpfp = tempFileprefix(this.fileprefix);
            deno = IfhParser.constructDenovo( ...
                nh.hdr, ...
                'fqfileprefix', tmpfp);
            s = struct(deno);
            this.verifyEqual(s.conversion_program, 'mlfourdfp.IfhParser.constructDenovo');
            this.verifyTrue(lstrfind(s.name_of_data_file, 'dgv2ConvertedLaForest_dcm2niix'));
            this.verifyEqual(s.matrix_size, [344 344 127 1]);
            this.verifyEqual(s.scaling_factor, [2.086260080337524 2.086260080337524 2.031250000000000], 'RelTol', 1e-4);
            
            deno.save(nh);
            this.verifyTrue(lexist(fullfile(this.TmpDir, [this.fileprefix '.4dfp.ifh']), 'file'));
            this.deleteExisting([tmpfp '.4dfp.ifh']);
        end
        function test_read(this)
        end
        function test_write(this)
        end
	end

 	methods (TestClassSetup)
		function setupIfhParser(this)
 			import mlfourdfp.*;
 			this.testObj_ = IfhParser;
 		end
	end

 	methods (TestMethodSetup)
		function setupIfhParserTest(this)
            this.pwd0 = pushd(this.TmpDir);
            fv = mlfourdfp.FourdfpVisitor;
            if (~lexist_4dfp(this.fileprefix))
                fv.copyfile_4dfp(fullfile(this.dataroot, this.fileprefix), this.fileprefix);
            end
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
        function deleteExisting(this, varargin)
            if (this.noDelete)
                return
            end
            deleteExisting(varargin{:});
        end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

