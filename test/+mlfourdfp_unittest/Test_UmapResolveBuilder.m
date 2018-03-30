classdef Test_UmapResolveBuilder < matlab.unittest.TestCase
	%% TEST_UMAPRESOLVEBUILDER 

	%  Usage:  >> results = run(mlfourdfp_unittest.Test_UmapResolveBuilder)
 	%          >> result  = run(mlfourdfp_unittest.Test_UmapResolveBuilder, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 19-Jul-2016 20:08:12
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/test/+mlfourdfp_unittest.
 	%% It was developed on Matlab 8.5.0.197613 (R2015a) for MACI64.
 	

	properties
        fourdfpVisitor
 		registry
        studyd
        sessd
 		testObj
        view = true
 	end

	methods (Test)
        function test_ctor(this)
            this.verifyEqual( ...
                pwd, ...
                fullfile(mlraichle.RaichleRegistry.instance.subjectsDir, 'NP995_22', 'V1', 'FDG_V1-NAC', ''));
            this.verifyEqual(this.testObj.tracer, 'FDG');
        end
        function test_buildCTMasked(this)          
            ctm = this.testObj.buildCTMasked;
            this.verifyEqual(ctm, this.sessd.ctMasked('typ', 'fqfp'));
            this.checkImg(ctm, this.sessd.ct('typ', 'fqfp'));
        end
        function test_rescaleCT(this)
            ct = this.testObj.rescaleCT;
            this.verifyEqual(ct, this.sessd.ctRescaled('typ', 'fqfp'));
            this.checkImg(ct, this.sessd.ct('typ', 'fqfp'));
        end
        function test_alignCTToSumtResolved(this) 
            ct = this.sessd.ct('typ', 'fqfp');
            [~,ct] = this.testObj.alignCTToSumtResolved(ct);
            this.verifyEqual(ct, 'ct_on_fdgv1r1_resolved_sumt');
            this.checkImg(ct, 'fdgv1r1_resolved_sumt');
        end
        function test_buildCarneyUmap(this)
            umap = this.testObj.buildCarneyUmap(this.sessd.ctRescaled('typ', 'fqfp'));
            this.verifyTrue(lexist([umap '.4dfp.ifh'], 'file'));
            this.checkImg(umap, this.sessd.ctRescaled('typ', 'fqfp'));
        end
        function test_alignUmapToNACFrames(this)  
            umap = this.sessd.umapSynthFDG('typ', 'fqfp');
            [~,umaps] = this.testObj.alignUmapToNACFrames(umap);            
            this.verifyTrue(lexist([umap '.4dfp.ifh'], 'file'));
            this.checkImg(umaps{end}, 'ct_on_fdgv1r1_resolved_sumt.4dfp.img');
        end        
        
        function test_mpr2atl_4dfp(this)
            [~,t4] = this.testObj.mpr2atl_4dfp;
            this.verifyEqual(t4, 'mpr_to_TRIO_Y_NDC_t4');
            this.verifyTrue(lexist(t4, 'file'));
            this.checkT4( ...
                t4, ...
                this.sessd.mpr('typ', 'fqfp'), ...
                this.sessd.atlas('typ', 'fqfp'));
        end
        function test_petSumt2mpr(this)
            [~,t4] = this.testObj.petSumt2mpr('fdgv1r1_resolved_sumt');
            this.verifyEqual(t4, 'fdgv1r1_resolved_sumt_to_mpr_t4');
            this.verifyTrue(lexist(t4, 'file'));            
            this.checkT4( ...
                t4, ...
                'fdgv1r1_resolved_sumt', ...
                this.sessd.mpr('typ', 'fqfp'));
        end  
        function test_CT2mpr_4dfp(this)
            [~,t4] = this.testObj.CT2mpr_4dfp(this.sessd.ct);
            this.verifyEqual(t4, 'ct_to_mpr_t4');
            this.verifyTrue(lexist(t4, 'file'));
            this.checkT4( ...
                t4, ...
                this.sessd.ct('typ', 'fqfp'), ...
                this.sessd.mpr('typ', 'fqfp'));
        end        
        function test_mprOnPetSumt(this)
            fp = this.testObj.mprOnPetSumt(this.sessd.mpr('typ', 'fqfp'), 'fdgv1r1_resolved_sumt');
            this.verifyTrue(lexist([fp '.4dfp.ifh'], 'file'));
            this.checkImg(fp, 'fdgv1r1_resolved_sumt');
        end
        function test_ctOnPetSumt(this)
            fp = this.testObj.ctOnPetSumt(this.sessd.ctMasked('typ', 'fqfp'), 'fdgv1r1_resolved_sumt');
            this.verifyTrue(lexist([fp '.4dfp.ifh'], 'file'));
            this.checkImg(fp, 'fdgv1r1_resolved_sumt');
        end
        
        function test_buildUmaps(this)
            this.testObj.buildUmaps('frames', [0 ones(1,71)]);
        end
	end

 	methods (TestClassSetup)
		function setupUmapResolveBuilder(this)
            this.studyd = mlraichle.StudyData;
            this.sessd = mlraichle.SessionData( ...
                'studyData', this.studyd, 'sessionPath', fullfile(this.studyd.subjectsDir, 'NP995_22',  ''));            
            this.testObj_ = mlfourdfp.UmapResolveBuilder('sessionData', this.sessd, 'frames', [0 ones(1,71)]);
            this.fourdfpVisitor = mlfourdfp.FourdfpVisitor;
 		end
	end

 	methods (TestMethodSetup)
		function setupUmapResolveBuilderTest(this)
 			this.testObj = this.testObj_;
 			this.addTeardown(@this.cleanFiles);
 		end
	end

	properties (Access = private)
 		testObj_
 	end

	methods (Access = private)
		function cleanFiles(this)
%             if (lexist(this.sessd.ctMasked('typ', 'fqfn'), 'file'))
%                 delete([this.sessd.ctMasked('typ', 'fqfp') '.4dfp.*']);
%             end
        end
        function checkT4(this, t4, source, dest)
            source = myfileprefix(source);
            dest   = myfileprefix(dest);
            fp     = this.fourdfpVisitor.t4img_4dfp(t4, source, 'options', ['-O' dest]);
            ic     = mlfourd.ImagingContext([fp '.4dfp.ifh']);
            if (this.view)
                ic.view([dest '.4dfp.img']); 
            end
        end
        function checkImg(this, source, varargin)
            source = myfileprefix(source);
            ic     = mlfourd.ImagingContext([source '.4dfp.ifh']);
            if (this.view)
                if (isempty(varargin))
                    ic.view;
                else
                    ic.view([varargin{1} '.4dfp.img']);
                end
            end
        end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

