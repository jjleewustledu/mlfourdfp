classdef Test_T4ResolveParser < matlab.unittest.TestCase
	%% TEST_T4RESOLVEPARSER 

	%  Usage:  >> results = run(mlfourdfp_unittest.Test_T4ResolveParser)
 	%          >> result  = run(mlfourdfp_unittest.Test_T4ResolveParser, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 28-Feb-2016 12:21:00
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/test/+mlfourdfp_unittest.
 	%% It was developed on Matlab 9.0.0.307022 (R2016a) Prerelease for MACI64.
 	

	properties
 		registry
 		testObj
 	end

	methods (Test)
        function test_parseLog(this)
            this = this.setupLog;
            this.verifyEqual(this.testObj.imgregLogParser.length, 144259);
            this.verifyEqual(this.testObj.etas{4,5}, 0.43914);
            this.verifyEqual(this.testObj.curves{4,5}, [615. 589. 746. 370. 434. 302.]);
        end
		function test_report(this)
            this = this.setupLog;
            r = this.testObj.report;
            this.verifyInstanceOf(r, 'mlfourdfp.T4ResolveReport');
            r.pcolor('etas',      this.testObj);
            r.pcolor('curves',    this.testObj);
            r.pcolor('z(etas)',   this.testObj);
            r.pcolor('z(curves)', this.testObj);
        end
        function test_report2(this)   
            this = this.setupLog;         
            r = this.testObj.report;
            this.verifyInstanceOf(r, 'mlfourdfp.T4ResolveReport');
            t4r = mlfourdfp.T4ResolveParser('sessionData', this.testObj.sessionData, 'frameLength', 28);
            t4r = t4r.parseLog( ...
                fullfile(t4r.sessionData.sessionPath, 'V1', 'imgreg_4dfp_pid3054.log'));
            t4r = t4r.shiftFrames(3);
            r.pcolor('etas',      t4r);
            r.pcolor('curves',    t4r);
            r.pcolor('z(etas)',   t4r);
            r.pcolor('z(curves)', t4r);
        end
        function test_reportD(this)   
            this = this.setupLog;         
            r = this.testObj.report;
            this.verifyInstanceOf(r, 'mlfourdfp.T4ResolveReport');
            t4r = mlfourdfp.T4ResolveParser('sessionData', this.testObj.sessionData, 'frameLength', 28);
            t4r = t4r.parseLog( ...
                fullfile(t4r.sessionData.sessionPath, 'V1', 'imgreg_4dfp_pid3054.log'));
            t4r = t4r.shiftFrames(3);
            r.d('etas',      t4r, this.testObj);
            r.d('curves',    t4r, this.testObj);
            r.d('z(etas)',   t4r, this.testObj);
            r.d('z(curves)', t4r, this.testObj);
        end
	end

 	methods (TestClassSetup)
		function setupT4Resolve(this)
 			import mlfourdfp.*;
            studyd = mlpipeline.StudyDataSingletons.instance('test_raichle');
            sessd = mlraichle.SessionData( ...
                'studyData',   studyd, ...
                'sessionPath', fullfile(studyd.subjectsDir, 'NP995_09'), ...
                'snumber',     1);
 			this.testObj_ = T4ResolveParser('sessionData', sessd, 'frameLength', 31);
 		end
	end

 	methods (TestMethodSetup)
		function setupT4ResolveTest(this)
 			this.testObj = this.testObj_;
 			this.addTeardown(@this.cleanFiles);
 		end
	end

	properties (Access = private)
 		testObj_
 	end

	methods (Access = private)
        function this = setupLog(this)
            this.testObj = this.testObj.parseLog( ...
                fullfile(this.testObj_.sessionData.sessionPath, 'V1', 'imgreg_4dfp.log'));
        end
		function cleanFiles(this)
 		end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

