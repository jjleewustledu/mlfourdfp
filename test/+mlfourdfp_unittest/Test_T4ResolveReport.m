classdef Test_T4ResolveReport < matlab.unittest.TestCase
	%% TEST_T4RESOLVEREPORT 

	%  Usage:  >> results = run(mlfourdfp_unittest.Test_T4ResolveReport)
 	%          >> result  = run(mlfourdfp_unittest.Test_T4ResolveReport, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 28-Feb-2016 15:06:12
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/test/+mlfourdfp_unittest.
 	%% It was developed on Matlab 9.0.0.307022 (R2016a) Prerelease for MACI64.
 	

	properties
 		registry
 		testObj
        frame_reg_fdgv1_log   = 'frame_reg_fdgv1_20160302T194001.log'
        frame_reg_fdgv1r1_log = 'frame_reg_fdgv1r1_20160303T010324.log'
        frame_reg_fdgv1r2_log = 'frame_reg_.sh_fdgv1r2_20160310T20171457662635.log' % w/o frame 20
 	end

	methods (Test)
		function test_report(this)
            r = this.testObj.report;
            r.pcolor('z(etas)',   this.testObj);
            %r.pcolor('z(curves)', this.testObj);
        end
        function test_report1(this)          
            r = this.testObj.report;
            t4r = mlfourdfp.T4ResolveParser('sessionData', this.testObj.sessionData);
            t4r = t4r.parseLog( ...
                fullfile(t4r.sessionData.sessionPath, this.frame_reg_fdgv1r1_log), ...
                'frameLength', 28);
            t4r = t4r.shiftFrames(3);
            r.pcolor('z(etas)',   t4r);
            %r.pcolor('z(curves)', t4r);
        end
        function test_report2(this)          
            r = this.testObj.report;
            t4r = mlfourdfp.T4ResolveParser('sessionData', this.testObj.sessionData);
            t4r = t4r.parseLog( ...
                fullfile(t4r.sessionData.sessionPath, this.frame_reg_fdgv1r2_log), ...
                'frameLength', 27);
            t4r = t4r.shiftFrames(3);
            r.pcolor('z(etas)',   t4r);
            %r.pcolor('z(curves)', t4r);
        end
%         function test_reportD(this)  
%             r = this.testObj.report;
%             t4r = mlfourdfp.T4ResolveParser('sessionData', this.testObj.sessionData);
%             t4r = t4r.parseLog( ...
%                 fullfile(t4r.sessionData.sessionPath, 'V1', 'frame_reg_fdgv1r2_20160303T201053.log'), ...
%                 'frameLength', 28);
%             t4r = t4r.shiftFrames(3);
%             r.d('z(etas)',   t4r, this.testObj);
%             r.d('z(curves)', t4r, this.testObj);
%         end
	end

 	methods (TestClassSetup)
		function setupT4ResolveReport(this)
            studyd = mlpipeline.StudyDataSingletons.instance('raichle');
            sessd  = mlraichle.SessionData( ...
                'studyData',   studyd, ...
                'sessionPath', fullfile(studyd.subjectsDir, '..', 'REG_TEST2', ''), ...
                'snumber',     1, ...
                'vnumber',     1);
 			this.testObj_ = mlfourdfp.T4ResolveParser('sessionData', sessd);
            this.testObj_ = this.testObj_.parseLog( ...
                fullfile(this.testObj_.sessionData.sessionPath, this.frame_reg_fdgv1_log), ...
                'frameLength', 31);
 		end
	end

 	methods (TestMethodSetup)
		function setupT4ResolveReportTest(this)
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

