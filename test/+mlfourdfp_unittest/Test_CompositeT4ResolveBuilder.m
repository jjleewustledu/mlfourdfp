classdef Test_CompositeT4ResolveBuilder < matlab.unittest.TestCase
	%% TEST_COMPOSITET4RESOLVEBUILDER 

	%  Usage:  >> results = run(mlfourdfp_unittest.Test_CompositeT4ResolveBuilder)
 	%          >> result  = run(mlfourdfp_unittest.Test_CompositeT4ResolveBuilder, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 21-Jan-2017 13:47:31
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/test/+mlfourdfp_unittest.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.  Copyright 2017 John Joowon Lee.
 	

	properties
        fvisitor
        hyglyNN = 'HYGLY28'
        ipr
        pwd0
 		studyd
        resolveTag = 'TestCT4RB'
        sessdFDG
        sessdHO
 		testObj
        theImages
        tracers = {'FDG' 'OC' 'OO' 'HO'}
        
        view = false
        quick = false
    end
    
    methods
        function c = allTracerSumtImages(this, varargin)
            % @param named typ =: ImagingContext.imagingType
            
            c = {};
            sd = this.sessdFDG;
            for it = 1:length(this.tracers)
                sd.tracer = this.tracers{it};
                c = [c sd.tracerResolvedFinalSumt(varargin{:})]; %#ok<AGROW>
            end
        end
        function h = partialTimeSumsHO(this, tinterval, varargin)
            % @param tinterval has numeric form [t0 tF]
            % @param named 'typ' =: ImagingContext.imagingType
            
            fp = sprintf('%s_t%ito%i', this.sessdHO.tracerResolvedFinalSumt('typ','fp'), tinterval(1), tinterval(2));
            if (~lexist_4dfp(fp))
                nn = mlfourd.NumericalNIfTId.load(this.sessdHO.tracerResolvedFinal);
                nn = nn.timeContracted(tinterval);
                nn.fileprefix = fp;
                nn.save;
            end
            h = this.sessdHO.fqfilenameObject(fullfile(this.pwd0, [fp '.4dfp.ifh']));
        end
    end

	methods (Test)
        function test_ctor(this)
            disp(this.testObj);
        end
        function test_finished(this)
            touchfile = fullfile(this.testObj.logPath, '.test_mlfourdfp.CompositeT4ResolveBuilder_isfinished.touch');
            if (lexist(touchfile)); delete(touchfile); end 
            this.verifyEqual(this.testObj.finished.finishedMarkerFilename, touchfile);           
            this.testObj.finished.touchFinishedMarker;
            this.verifyTrue(lexist(this.testObj.finished.finishedMarkerFilename, 'file'));
            this.verifyTrue(this.testObj.finished.isfinished);
            delete(touchfile); 
        end
        function test_resolve(this)
            tic
            this.testObj.NRevisions = 2;
            this.testObj = this.testObj.resolve('source', this.theImages);
            [e,c] = this.fvisitor.etaAndCurves( ...
                cellfun(@(x) [x '_' this.resolveTag], {this.theImages{logical([0 0 1 0 1])}}));
            this.verifyEqual(e, 0.1631, 'RelTol', 0.05);
            this.verifyEqual(c, [130 133 21 89 161 66], 'RelTol', 1);
            prod = this.testObj.product;
            entr = [3.6580 3.2309];
            for p = 1:length(prod)
                this.verifyEqual(prod{p}.niftid.entropy, entr(p), 'RelTol', 0.05);
            end            
            this.verifyTrue(this.testObj.isfinished(this.testObj.sessionData));            
            fprintf('test_resolve:  ');
            toc
            if (this.view)
                prod{p}.niftid.fslview;
            end
        end
        function test_revise(this)
            if (this.quick)
                return
            end
            this.ipr.resolved = 'testNativeCropped';
            [~,this.testObj] = this.testObj.revise(this.ipr);
            this.verifyEqual(this.testObj.sessionData.rnumber, 2);
        end
        function test_imageReg(this)
            if (this.quick)
                return
            end
            
            this.testObj.imageReg(this.ipr);
            
            in3 = 'testNativeCropped3';
            in5 = 'testNativeCropped5';
            out = sprintf('%s_on_%s', in5, in3);
            this.fvisitor.t4_mul([in3 '_to_' in5 '_t4'], [ in5 '_to_' in3 '_t4']);
            movefile([in3 '_to_' in3 '_t4'], [in3 '_to_' in3 '_t4.txt']);
            t4 = readtable([in3 '_to_' in3 '_t4'], ...
                'ReadVariableNames', 0, 'Delimiter', ' ', 'MultipleDelimsAsOne', 1, 'HeaderLines', 4, 'Format', '%f %f %f %f');
            this.verifyEqual(t4{1,1}, 1,    'AbsTol', 0.01);
            this.verifyEqual(t4{1,2}, 0.1,  'AbsTol', 0.1);
            this.verifyEqual(t4{1,3}, 0.02, 'AbsTol', 0.02);
            this.verifyEqual(t4{1,4}, 0.2,  'AbsTol', 0.2);
            
            this.fvisitor.t4img_4dfp(this.fvisitor.filenameT4(in5, in3), in5, 'out', out, 'options', ['-O' in3]);
            if (this.view)
                mlbash(sprintf('fslview %s.4dfp.img %s.4dfp.img', in3, out));
            end
        end
        function test_t4ResolveAndPaste(this)
            if (this.quick)
                return
            end
            
            this.testObj.resolveLog = loggerFilename('', 'func', 'Test_CompositeT4ResolveBuilder_test_t4ResolveAndPaste');
            [ipr_,imgFns_] = this.testObj.resolveAndPaste(this.ipr);
            this.verifyEqual(ipr_.resolved, ...
                {'testNativeCropped1_TestCT4RB' ...
                 'testNativeCropped2_TestCT4RB' ...
                 'testNativeCropped3_TestCT4RB' ...
                 'testNativeCropped4_TestCT4RB' ...
                 'testNativeCropped5_TestCT4RB'});
            this.verifyEqual(imgFns_, 'testNativeCropped3 testNativeCropped5');
            if (this.view)
                mlbash('fslview testNativeCropped3_Test_CompositeT4ResolveBuilder.4dfp.img testNativeCropped5_Test_CompositeT4ResolveBuilder.4dfp.img');
            end
        end
        
        %% UTILITY
        
        function test_lazyStageImages(this)
            this.ipr.currentIndex = nan;
            
            fqfps = this.testObj.lazyStageImages(this.ipr);
            this.verifyEqual(length(fqfps), 2);
            this.verifyEqual(mybasename(fqfps{1}), 'testNativeCropped3');
            this.verifyEqual(mybasename(fqfps{2}), 'testNativeCropped5');
            if (this.view)
                mlbash('fslview testNativeCropped.4dfp.img -b 0,2000 testNativeCropped3.4dfp.img -b 0,2000 testNativeCropped5.4dfp.img -b 0,2000');
            end
        end
        function test_lazyBlurImage(this)
            if (this.quick)
                return
            end
            
            this.ipr.currentIndex = 1;
            fqfp = this.testObj.lazyBlurImage(this.ipr, 5.5);
            this.verifyEqual(mybasename(fqfp), 'testNativeCropped1_b55');
            if (this.view)
                mlbash('fslview testNativeCropped.4dfp.img -b 0,100 testNativeCropped1_b55.4dfp.img -b 0,100');
            end
        end
        function test_lazyBlurImages(this)
            if (this.quick)
                return
            end
            
            this.ipr.currentIndex = nan;            
            fqfps = this.testObj.lazyBlurImages(this.ipr);
            this.verifyEqual(length(fqfps), 2);
            this.verifyEqual(mybasename(fqfps{1}), 'testNativeCropped3_b55');
            this.verifyEqual(mybasename(fqfps{2}), 'testNativeCropped5_b55');
            if (this.view)
                mlbash('fslview testNativeCropped.4dfp.img -b 0,2000 testNativeCropped3_b55.4dfp.img -b 0,2000 testNativeCropped5_b55.4dfp.img -b 0,2000');
            end
        end
        function test_lazyMaskForImages(this)
            if (this.quick)
                return
            end
            
            tic
            fqfp = this.testObj.lazyMaskForImages('maskForImages', 'testNativeCropped3', 'testNativeCropped5');
            this.verifyEqual(mybasename(fqfp), 'maskForImages_testNativeCropped3_testNativeCropped5');
            if (this.view)
                mlbash('fslview testNativeCropped3.4dfp.img -b 0,2000 testNativeCropped5.4dfp.img -b 0,2000 maskForImages_3_5.4dfp.img');
            end
            fprintf('test_lazyMaskForImages:  ');
            toc
        end
	end

 	methods (TestClassSetup)
		function setupCompositeT4ResolveBuilder(this)
            import mlraichle.*;
            this.studyd = SynthStudyData;
            this.sessdFDG = SynthSessionData( ...
                'studyData', this.studyd, ...
                'sessionFolder', this.hyglyNN, ...
                'vnumber', 2, ...
                'rnumber', 2, ...
                'ac', true);  
            this.sessdHO = this.sessdFDG;
            this.sessdHO.tracer = 'HO';
            
            this.pwd0 = pushd(this.sessdFDG.vLocation);        
            this.testObj_ = mlfourdfp.CompositeT4ResolveBuilder( ... 
                'sessionData', this.sessdFDG, ...  
                'NRevisions', 1, ...
                'resolveTag', this.resolveTag, ...   
                'theImages', this.theImages, ...
                'indicesLogical', this.indicesLogical, ... 
                'indexOfReference', 3);
            this.ipr = struct( ...
                'dest', [], ...
                'source', [], ...
                'maskForImages', 'none', ...
                'sourceBlur', [], ...
                'destBlur', [], ...
                'keepForensics', true, ...
                'currentIndex', nan);
            this.ipr.dest = this.theImages;
            this.ipr.source = this.theImages;
            this.ipr.sourceBlur = this.blurArray;
            this.ipr.destBlur = this.blurArray;
            this.fvisitor = mlfourdfp.FourdfpVisitor;
            this.ensureImages;
            setenv('DEBUG', ''); % cf. dbbash
 			this.addTeardown(@this.teardownClass);
 		end
	end

 	methods (TestMethodSetup)
		function setupCompositeT4ResolveBuilderTest(this)
 			this.testObj = this.testObj_;
 			this.addTeardown(@this.teardownMethod);
 		end
    end

    %% PRIVATE
    
	properties (Access = private)
 		testObj_
 	end

	methods (Access = private)
        function b = blurArray(~)
            b = zeros(1, 5);
            b(3) = 5.5;
            b(5) = 5.5;
        end
        function f = indicesLogical(~)
            f = false(1, 5);
            f(3) = true;           
            f(5) = true;
        end
        function ensureImages(this)
            if (~this.fvisitor.lexist_4dfp('testNativeCropped'))
                this.fvisitor.copy_4dfp(this.sessdHO.tracerResolvedFinal, 'testNativeCropped');
            end
            for idx = 1:5
                if (~this.fvisitor.lexist_4dfp(sprintf('testNativeCropped%i', idx)))
                    this.fvisitor.extract_frame_4dfp('testNativeCropped', idx);
                    this.fvisitor.move_4dfp(sprintf('testNativeCropped_frame%i', idx), ...
                                            sprintf('testNativeCropped%i', idx));
                end
            end
        end
		function teardownMethod(~)
 		end
		function teardownClass(this)
            if (~this.quick)
                this.testObj_.keepForensics = false;
                this.testObj_.teardownResolve(this.ipr);
                deleteExisting('*.nii');
                deleteExisting('*.nii.gz');
                if (isdir('Log')); rmdir('Log', 's'); end %#ok<*ISDIR>
                if (isdir('T4'));  rmdir('T4', 's');  end
            end
            deleteExisting('*.mat0');
            deleteExisting('*.sub');
            popd(this.pwd0);
 		end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

