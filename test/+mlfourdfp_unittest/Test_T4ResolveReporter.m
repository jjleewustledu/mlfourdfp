classdef Test_T4ResolveReporter < matlab.unittest.TestCase
	%% TEST_T4RESOLVEREPORTER 

	%  Usage:  >> results = run(mlfourdfp_unittest.Test_T4ResolveReporter)
 	%          >> result  = run(mlfourdfp_unittest.Test_T4ResolveReporter, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 04-May-2016 23:14:05
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/test/+mlfourdfp_unittest.
 	%% It was developed on Matlab 9.0.0.341360 (R2016a) for MACI64.
 	

	properties
 		registry
 		testObj
 	end

	methods (Test)
        function test_reports(this)
            for r = 1:this.testObj.length
                report = this.testObj.reports{r};
                report.pcolor('z(etas)');
            end
        end
        function test_reports1(this)
            report = this.testObj.reports{1};
            report.pcolor('z(etas)');
        end
        function test_ctor(this)
            this.verifyTrue(isa(this.testObj, 'mlfourdfp.T4ResolveReporter'));
        end
	end

 	methods (TestClassSetup)
		function setupT4ResolveReporter(this)
 			import mlfourdfp.*;
            studyd = mlpipeline.StudyDataSingletons.instance('raichle');
            sessd  = mlraichle.SessionData( ...
                'studyData',   studyd, ...
                'sessionPath', fullfile(studyd.subjectsDir, '..', 'REG_TEST2', ''), ...
                'snumber',     1, ...
                'vnumber',     1);
 			this.testObj_ = T4ResolveReporter( ...
                'sessionData', sessd, ...
                'imagingPath', sessd.sessionPath, ...
                'loggingPath', sessd.sessionPath, ...
                'imagingFileprefix', 'fdgv1', ...
                'loggingFileprefix', 'frame_reg_fdgv1');
 		end
	end

 	methods (TestMethodSetup)
		function setupT4ResolveReporterTest(this)
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

