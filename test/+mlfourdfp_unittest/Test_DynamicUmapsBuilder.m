classdef Test_DynamicUmapsBuilder < matlab.unittest.TestCase
	%% TEST_DYNAMICUMAPSBUILDER 

	%  Usage:  >> results = run(mlfourdfp_unittest.Test_DynamicUmapsBuilder)
 	%          >> result  = run(mlfourdfp_unittest.Test_DynamicUmapsBuilder, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 14-Dec-2016 17:51:42
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/test/+mlfourdfp_unittest.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.
 	

	properties
        hyglyNN = 'HYGLY09'
        pwd0
 		registry
        sessd
        studyd
 		testObj
 	end

	methods (Test)		
		function test_ctor(this)
            this.verifyClass(this.testObj, 'mlfourdfp.DynamicUmapsBuilder');
 		end
		function test_buildUmaps(this)
            this.testObj.buildUmaps;
            this.testObj.viewUmaps;
            this.verifyTrue(this.testObj_.isfinished);
        end
	end

 	methods (TestClassSetup)
		function setupDynamicUmapsBuilder(this)
            this.studyd = mlraichle.StudyData;
            this.sessd = mlraichle.SessionData( ...
                'studyData', this.studyd, ...
                'sessionPath', fullfile(this.studyd.subjectsDir, this.hyglyNN, ''), ...
                'rnumber', 1, ...
                'tracer', 'FDG');
 			this.testObj_ = mlfourdfp.DynamicUmapsBuilder('sessionData', this.sessd);
 		end
	end

 	methods (TestMethodSetup)
		function setupDynamicUmapsBuilderTest(this)
 			this.testObj = this.testObj_;
            this.pwd0 = pushd(this.sessd.sessionLocation);
 			this.addTeardown(@this.cleanFiles);
 		end
	end

	properties (Access = private)
 		testObj_
 	end

	methods (Access = private)
		function cleanFiles(this)
            if (lexist(this.testObj_.finished.markerFilename, 'file'))
                delete(this.testObj_.finished.markerFilename); end
            popd(this.pwd0);
 		end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

