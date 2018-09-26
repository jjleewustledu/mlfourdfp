classdef Test_DicomSorter < matlab.unittest.TestCase
	%% TEST_DICOMSORTER 

	%  Usage:  >> results = run(mlfourdfp_unittest.Test_DicomSorter)
 	%          >> result  = run(mlfourdfp_unittest.Test_DicomSorter, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 18-Sep-2016 15:18:39
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/test/+mlfourdfp_unittest.
 	%% It was developed on Matlab 9.0.0.341360 (R2016a) for MACI64.
 	

	properties
        sessionData
        sessionPath
        studyData
 		testObj
        ctFqfp  = fullfile(mlraichle.RaichleRegistry.instance.subjectsDir, 'HYGLY28', 'AC_CT')
        mprFqfp = fullfile(mlraichle.RaichleRegistry.instance.subjectsDir, 'HYGLY28', 'V1', 't1_mprage_sag')
 	end

	methods (Test)
        function test_dcmInfos(this)
            cd(fullfile(this.studyData.rawdataDir, 'HYGLY28_VISIT_1', ''));
            [infos,fqdns] = this.testObj.dcmInfos;
            for idx = 1:5
                this.verifyEqual(fullfile(pwd, num2str(infos{idx}.SeriesNumber)), fqdns{idx});
                fprintf('%s\t%g\t%s\t%s\n', fqdns{idx}, infos{idx}.SeriesNumber, infos{idx}.SeriesDescription, infos{idx}.ProtocolName);
            end
        end
        function test_findDcmInfos(this)
            cd(fullfile(this.studyData.rawdataDir, 'HYGLY28_VISIT_1', ''));
            infos = this.testObj.findDcmInfos('UMAP');
            for idx = 1:length(infos)
                fprintf('%g\t%s\t%s\t%s\t%s\n', ...
                    infos{idx}.SeriesNumber, infos{idx}.SeriesDescription, infos{idx}.ProtocolName, ...
                    infos{idx}.ImageType, infos{idx}.SequenceName);
            end
        end
        function test_findDcmInfos2(this)
            cd(fullfile(this.studyData.rawdataDir, 'HYGLY28_VISIT_1', ''));
            infos = this.testObj.findDcmInfos({'UMAP' 't2_spc_sag'});
            for idx = 1:length(infos)
                fprintf('%g\t%s\t%s\t%s\t%s\n', ...
                    infos{idx}.SeriesNumber, infos{idx}.SeriesDescription, infos{idx}.ProtocolName, ...
                    infos{idx}.ImageType, infos{idx}.SequenceName);
            end
        end
        function test_destPath(this)
            cd(this.testObj.studyData.rawdataDir);
            dt = mlsystem.DirTools('HYGLY*', 'NP*', 'hygly*', 'TW*', 'DT*');
            for idns = 1:length(dt.dns)
                fprintf('dns-> %s\ndestPath -> %s\n\n', dt.dns{idns}, this.testObj.destPath(dt.dns{idns}));
            end
        end
        function test_session_to_4dfp(this)
            src  = fullfile(getenv('PPG'), 'rawdata', 'HYGLY28_VISIT_1', '');
            dest = fullfile(mlraichle.RaichleRegistry.instance.subjectsDir, 'HYGLY28', 'V1', '');
            this.testObj.session_to_4dfp( ...
                src, dest, ...
                'seriesFilter', 't1_mprage_sag', ...
                'studyData', this.studyData, ...
                'preferredName', 't1_mprage_sag')
            import mlfourdfp.*;
            this.verifyTrue(FourdfpVisitor.lexist_4dfp(this.mprFqfp));
        end
        function test_sessions_to_4dfp(this)
            this.testObj.sessions_to_4dfp( ...
                'sessionFilter', 'HYGLY28*', ...
                'seriesFilter', 't1_mprage_sag', ...
                'studyData', this.studyData, ...
                'preferredName', 't1_mprage_sag');
            import mlfourdfp.*;
            this.verifyTrue(FourdfpVisitor.lexist_4dfp(this.mprFqfp));
        end
        function test_sessions_to_4dfp_CT(this)
            this.testObj.sessions_to_4dfp( ...
                'sessionFilter', 'HYGLY28_VISIT_1_CT*', ...
                'seriesFilter', {'AC_CT'}, ...
                'studyData', this.studyData, ...
                'preferredName', 'AC_CT');
            import mlfourdfp.*;
            this.verifyTrue(FourdfpVisitor.lexist_4dfp(this.ctFqfp));
        end
	end

 	methods (TestClassSetup)
		function setupDicomSorter(this)
 			import mlfourdfp.* mlraichle.*;
            this.studyData = StudyData;
            this.sessionPath = fullfile(RaichleRegistry.instance.subjectsDir, 'HYGLY28', '');
            this.sessionData = SessionData('studyData', this.studyData, 'sessionPath', this.sessionPath);
 			this.testObj_  = DicomSorter('sessionData', this.sessionData);
 		end
	end

 	methods (TestMethodSetup)
		function setupDicomSorterTest(this)
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

