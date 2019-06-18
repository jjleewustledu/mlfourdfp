classdef Test_RawDataSorter < matlab.unittest.TestCase
	%% TEST_RAWDATASORTER 

	%  Usage:  >> results = run(mlfourdfp_unittest.Test_RawDataSorter)
 	%          >> result  = run(mlfourdfp_unittest.Test_RawDataSorter, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 06-Sep-2016 15:36:58
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/test/+mlfourdfp_unittest.
 	%% It was developed on Matlab 9.0.0.341360 (R2016a) for MACI64.
 	

	properties
 		studyData
 		testObj
 	end

	methods (Test)
        function test_setup(this)
            this.verifyClass(this.testObj, 'mlfourdfp.RawDataSorter');
            this.verifyClass(this.testObj.studyData, 'mlraichle.StudyDataSingleton');
            this.verifyClass(this.testObj.sessionData, 'mlpatterns.CellComposite');
        end
        function test_copyRawData(this)
            import mlraichle.* mlsystem.*;
            
            srcLoc   = fullfile(getenv('PPG'), 'rawdata/SYNTH25_V1_1234/RESOURCES/RawData', '');
            sessPath = fullfile(mlraichle.StudyRegistry.instance.subjectsDir, 'SYNTH25', '');
            vPath    = fullfile(sessPath, 'V1', '');
            this.testObj.sessionData = SessionData( ...
                'sessionPath', sessPath, ...
                'studyData', this.studyData);
            this.testObj = this.testObj.copyRawData(srcLoc);
            
            tracers = {'FDG' 'HO1' 'HO2' 'OC1' 'OC2' 'OO1' 'OO2'};
            for t = 1:length(tracers)
                this.verifyTrue(isdir(fullfile(vPath, [tracers{t} '_V1'], '')));
            end
        end
        function test_copyUTE(this)
            import mlraichle.* mlsystem.*;
            
            srcLoc = fullfile(getenv('PPG'), 'rawdata/HYGLY25_VISIT_1', 'scans', '');
            sessPath = fullfile(mlraichle.StudyRegistry.instance.subjectsDir, 'HYGLY25', '');
            vPath = fullfile(sessPath, 'V1', '');
            this.testObj.sessionData = SessionData( ...
                'sessionPath', sessPath, ...
                'studyData', this.studyData);
            this.testObj = this.testObj.copyUTE(srcLoc);
            
            tracers = {'FDG' 'HO1' 'HO2' 'OC1' 'OC2' 'OO1' 'OO2'};
            for t = 1:length(tracers)
                destLoc = fullfile(vPath, [tracers{t} '_V1'], 'umap', '');
                dt = DirTool(fullfile(destLoc, '*'));
                this.verifyTrue(isdir(destLoc));
                this.verifyTrue(~isempty(dt.fns));
            end
        end
        function test_rawDataMatch(this)
            import mlfourdfp.* mlraichle.*;
            RawDataPath = fullfile(getenv('PPG'), 'rawdata', 'SYNTH25_v1_1234', 'RESOURCES', 'RawData', '');
            sessPath    = fullfile(mlraichle.StudyRegistry.instance.subjectsDir, 'SYNTH25', '');
            vPath       = fullfile(sessPath, 'V1', '');
 			this.testObj = RawDataSorter( ...
                           'studyData',   this.studyData, ...
                           'sessionData', SessionData('sessionPath', sessPath));
            [destLocs,srcLocs] = this.testObj.rawDataMatch(RawDataPath, vPath);
            
            this.verifyEqual(destLocs, ...
                {fullfile(vPath, 'HO1_V1', '') ...
                 fullfile(vPath, 'HO2_V1', '') ...
                 fullfile(vPath, 'FDG_V1', '') ...
                 fullfile(vPath, 'OC1_V1', '') ...
                 fullfile(vPath, 'OC2_V1', '') ...
                 fullfile(vPath, 'OO1_V1', '') ...
                 fullfile(vPath, 'OO2_V1', '')});
            this.verifyEqual(srcLocs, ...
                {fullfile(RawDataPath, 'Head_HO1_HD_PET_RawData', '') ...
                 fullfile(RawDataPath, 'Head_HO2_HD_PET_RawData', '') ...
                 fullfile(RawDataPath, 'Head_MRAC_PET_60min_RawData', '') ...
                 fullfile(RawDataPath, 'Head_OC1_HD_MRAC_RawData', '') ...
                 fullfile(RawDataPath, 'Head_OC2_HD_MRAC_RawData', '') ...
                 fullfile(RawDataPath, 'Head_OO1_HD_PET_RawData', '') ...
                 fullfile(RawDataPath, 'Head_OO2_HD_PET_RawData', '')});
        end
        function test_UTEMatch(this)
            import mlfourdfp.* mlraichle.*;
            scansPath = fullfile(getenv('PPG'), 'rawdata', 'HYGLY25_VISIT_1', 'scans', '');
            sessPath  = fullfile(mlraichle.StudyRegistry.instance.subjectsDir, 'HYGLY25', '');
            vPath     = fullfile(sessPath, 'V1', '');
 			this.testObj = RawDataSorter( ...
                           'studyData',   this.studyData, ...
                           'sessionData', SessionData('sessionPath', sessPath, 'tag', 'UTE_AC_only_UMAP'));
            [destLocs,srcLocs] = this.testObj.UTEMatch(scansPath, vPath);
            
            this.verifyEqual(destLocs, ...
                {fullfile(vPath, 'FDG_V1', 'umap', '') ...
                 fullfile(vPath, 'HO1_V1', 'umap', '') ...
                 fullfile(vPath, 'HO2_V1', 'umap', '') ...
                 fullfile(vPath, 'OC1_V1', 'umap', '') ...
                 fullfile(vPath, 'OC2_V1', 'umap', '') ...
                 fullfile(vPath, 'OO1_V1', 'umap', '') ...
                 fullfile(vPath, 'OO2_V1', 'umap', '')});
            for s = 1:7
                this.verifyEqual(srcLocs{s}, ...
                    fullfile(scansPath, '45_Head_UTE_AC_only_UMAP', 'DICOM', ''));
            end
        end
	end

 	methods (TestClassSetup)
		function setupRawDataSorter(this)
 		end
	end

 	methods (TestMethodSetup)
		function setupRawDataSorterTest(this)
 			import mlfourdfp.* mlraichle.*;
            this.studyData = StudyDataSingleton.instance('initialize');
 			this.testObj   = RawDataSorter('studyData', this.studyData);
 			this.addTeardown(@this.cleanFiles);
            cd(this.studyData.subjectsDir);
            
            synth = fullfile(this.studyData.subjectsDir, 'SYNTH25', '');
            if (hostnameMatch('ophthalmic') && ~isdir(synth))
                mkdir(synth);
            end
            hygly = fullfile(this.studyData.subjectsDir, 'HYGLY25', '');
            if (hostnameMatch('ophthalmic') && ~isdir(hygly))
                mkdir(hygly);
            end
 		end
    end
    
	methods (Access = private)
        function ensureDestinations(~)
            hygly = fullfile(this.studyData.subjectsDir, 'HYGLY25');
            if (~isdir(hygly))
                mkdir( hygly);
            end
            synth = fullfile(this.studyData.subjectsDir, 'SYNTH25');
            if (~isdir(synth))
                mkdir( synth);
            end
        end
		function cleanFiles(this)
            if (hostnameMatch('ophthalmic'))
                hygly = fullfile(this.studyData.subjectsDir, 'HYGLY25', '');
                if (isdir(hygly))
                    rmdir(hygly, 's');
                end
                synth = fullfile(this.studyData.subjectsDir, 'SYNTH25', '');
                if (isdir(synth))
                    rmdir(synth, 's');
                end
            end
 		end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

