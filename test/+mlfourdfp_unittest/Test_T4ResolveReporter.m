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
                report.pcolor('etas');
            end
        end
        function test_reports1(this)
            report = this.testObj.reports{1};
            report.pcolor('etas');
        end
        function test_ctor(this)
            this.verifyTrue(isa(this.testObj, 'mlfourdfp.T4ResolveReporter'));
        end
	end

 	methods (TestClassSetup)
		function setupT4ResolveReporter(this)
 			import mlfourdfp.*;
            studyd = mlraichle.SynthStudyData;
            sessd  = mlraichle.SynthSessionData( ...
                'studyData',   studyd, ...
                'sessionPath', fullfile(mlraichle.RaichleRegistry.instance.subjectsDir, 'HYGLY09', ''), ...
                'snumber',     1);
 			this.testObj_ = T4ResolveReporter( ...
                'sessionData', sessd, ...
                'imagingPath', sessd.fdgNACLocation, ...
                'loggingPath', fullfile(sessd.fdgNACLocation, 'Log', ''), ...
                'imagingFileprefix', 'fdgv1r1_resolved_b55', ...
                'loggingFileprefix', 'fdgv1r1_FdgBuilder_redoT4ResolveAndPasteForAll_D20170103T074515481', ...
                'frameLength', 72);
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

