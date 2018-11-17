classdef (Abstract) AbstractT4ResolveError < mlfourdfp.AbstractSessionBuilder
	%% ABSTRACTT4RESOLVEERROR  

	%  $Revision$
 	%  was created 29-May-2018 18:15:22 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
    properties
        indicesLogical
        noiseFloorOfCounts = 200 
    end
    
    properties (Dependent)
        errMat
        resolveTag
        theImages
    end
    
    methods (Static)
        function [m,s] = meanAndStd(~, em)
            assert(isnumeric(em));
            m  = cell2mat(cellfun(@(x) mean(mean(x,[],'omitnan'),[],'omitnan'), em, 'UniformOutput', false));
            s  = cell2mat(cellfun(@(x) std( std( x,[],'omitnan'),[],'omitnan'), em, 'UniformOutput', false));            
        end
    end
    
    methods (Abstract)
        stagedImgs(this)
    end
    
	methods 
		
        %% GET
        
        function g = get.errMat(this)
            g = this.errMat_;
        end
        function g = get.theImages(this)
            g = this.theImages_;
        end
        function g = get.resolveTag(this)
            g = this.sessionData.resolveTag;
        end
        
        %%
        
        function [this,em] = estimateErr(this, simgs, varargin)
            %  @param simgs must be already screened by this.assessValidFrames.
            %  @param indicesLogical is optional.
            %  @return em, the norm-2 error of displacements & rotations.
            
            ip = inputParser;
            addOptional(ip, 'indicesLogical', this.indicesLogical, @islogical);
            addParameter(ip, 'rnumber', this.sessionData.rnumber, @isnumeric);
            parse(ip, varargin{:});
            this.indicesLogical = ip.Results.indicesLogical;
            this.sessionData_.rnumber = ip.Results.rnumber;            
            % simgs = this.ensureRNumberOfImgs(simgs); % SUSPECTED BUG
            
            this = this.updateLogging;
            this.logger.add('estimateErr: working in %s', pwd);
            len = length(simgs);
            this.errMat_ = nan(len, len);
            for m = 1:length(simgs)
                for n = 1:length(simgs)
                    if (m ~= n)
                        try
                            this.errMat_(m,n) = ...
                                this.pairedErrParser(mybasename(simgs{m}), mybasename(simgs{n}));
                        catch ME
                            dispwarning(ME);
                        end                        
                    end
                end
            end 
            em = this.errMat;
            this.logger.add('AbstractT4ResolveError.estimateErr: stagedImgs->%s\n', ...
                cell2str(simgs, 'AsRow', true));
            this.logger.add('AbstractT4ResolveError.estimateErr: this.errMat->%s\n', ...
                mat2str(this.errMat));  
            this.logger.add('saving em to %s.mat\n', this.logger.fqfileprefix);
            save([this.logger.fqfileprefix '.mat'], 'em');
            this.logger.save;
        end
        function imgs  = ensureRNumberOfImgs(this, imgs) %#ok<INUSL>
            %% imgs should end with r[0-9] to be compliant with mlfourdfp.CompositeT4ResolveBuilder,
            %  but T4ResolveBuilder has tags _framenn and no r[0-9] should follow these tags.  
            %  mlraichle.SubjectImages and other contexts.  
            %  @param imgs is char or cell.
            %  @return imgs is char or cell.
            
            %% SUSPECTED BUG            
            return
            
            % flat recursion for cell-arrays
            if (iscell(imgs))
                imgs = cellfun(@(x) this.ensureRNumberOfImgs(x), imgs, 'UniformOutput', false);
                return
            end
            
            % base case
            assert(this.sessionData.rnumber < 10);
            rstr = sprintf('r%i', this.sessionData.rnumber);
            if (~strcmp(imgs(end-1:end), rstr))
                imgs = [imgs rstr];
            end
        end
        function err   = pairedErrParser(this, obj1, obj2)
            [rmsdeg,rmsmm] = this.errParser(this.t4_resolvePair(obj1, obj2));
            err = this.errMetric(rmsdeg, rmsmm);
        end
        function [rmsdeg,rmsmm] = errParser(~, r)
            %% T4_RESOLVE_PAIR
            %  @param r is the result from a call to mlbash.
            %  @return rmsdeg is numeric error found by t4_resolve for pair.
            %  @return rmsmm  is numeric error found by t4_resolve for pair.
            % 
            % jjlee@william
            % [ /data/nil-bluearc/raichle/PPGdata/jjlee2/HYGLY14/V1 ]$ t4_resolve T1001 fdgv1r2_op_fdgv1e1to4r1_frame4_sumt -oT1001_to_op_fdgv1r2_op_fdgv1e1to4r1_frame4_sumt
            % $Id: t4_resolve.c,v 1.6 2013/09/12 00:31:12 avi Exp $
            % n=2
            % rms error averaged over   2 pairs
            %     0.2726    0.1902    0.2722    0.0986    0.2615    0.1313
            % pairs total rotation error        0.4296 (rms deg)
            % pairs total translation error     0.3088 (rms mm)
            % t4_resolve: VOI rms radius=50.0000
            % t4_resolve: begin Gauss-Newton trajectory estimation
            % estimated observational error based on   2 observations
            %     0.3722    0.2711    0.3963    0.0046    0.4277    0.2680
            % estimate total rotation error     0.6075 (rms deg)
            % estimate total translation error  0.5047 (rms mm)
            % rigid body rmserr       0.2113
            % JtJ condition number0.1638E+05
            % scale rmserr            0.0000
            % err=0.535842 err0=0.000000
            % estimated observational error based on   2 observations
            %     0.2845    0.1895    0.2715    0.0639    0.2449    0.0563
            % estimate total rotation error     0.4365 (rms deg)
            % estimate total translation error  0.2593 (rms mm)
            % rigid body rmserr       0.1330
            % JtJ condition number0.1633E+05
            % scale rmserr            0.0000
            % err=0.212365 err0=0.535842
            % estimated observational error based on   2 observations
            %     0.2846    0.1895    0.2715    0.0633    0.2449    0.0561
            % estimate total rotation error     0.4366 (rms deg)
            % estimate total translation error  0.2591 (rms mm)
            % rigid body rmserr       0.1330
            % JtJ condition number0.1633E+05
            % scale rmserr            0.0000
            % err=0.212343 err0=0.212365
            % estimated observational error based on   2 observations
            %     0.2846    0.1895    0.2715    0.0633    0.2449    0.0561
            % estimate total rotation error     0.4366 (rms deg)
            % estimate total translation error  0.2591 (rms mm)
            % rigid body rmserr       0.1330
            % JtJ condition number0.1633E+05
            % scale rmserr            0.0000
            % err=0.212343 err0=0.212343
            % estimated observational error based on   2 observations
            %     0.2846    0.1895    0.2715    0.0633    0.2449    0.0561
            % estimate total rotation error     0.4366 (rms deg)
            % estimate total translation error  0.2591 (rms mm)
            % rigid body rmserr       0.1330
            % JtJ condition number0.1633E+05
            % scale rmserr            0.0000
            % err=0.212343 err0=0.212343
            % Writing: T1001_to_op_fdgv1r2_op_fdgv1e1to4r1_frame4_sumt.sub
            % T1001	t4=T1001_to_T1001_to_op_fdgv1r2_op_fdgv1e1to4r1_frame4_sumt_t4
            % fdgv1r2_op_fdgv1e1to4r1_frame4_sumt	t4=fdgv1r2_op_fdgv1e1to4r1_frame4_sumt_to_T1001_to_op_fdgv1r2_op_fdgv1e1to4r1_frame4_sumt_t4

            rmsdeg = nan;
            rmsmm  = nan;
            
            assert(ischar(r) || isstring(r))
            r = splitlines(r);
            %r(r == "") = [];
            r = strip(r);
            foundLabel = strfind(r, 'rms error averaged over');
            ir = 0;
            while (ir <= length(foundLabel))
                ir = ir + 1;
                if (~isempty(foundLabel{ir}))
                    break
                end
            end
            if (ir >= length(foundLabel))
                return
            end
            if (contains(r(ir + 2), 'pairs total rotation error'))
                toknames = regexp(r(ir + 2), '(?<rmsdeg>\d+.\d+) \(rms deg\)', 'names');
                if (~isempty(toknames))
                    rmsdeg = str2double(toknames.rmsdeg);
                end
            end
            if (contains(r(ir + 3), 'pairs total translation error'))
                toknames = regexp(r(ir + 3), '(?<rmsmm>\d+.\d+) \(rms mm\)', 'names');
                if (~isempty(toknames))
                    rmsmm = str2double(toknames.rmsmm);
                end
            end 
        end       
        function em    = errMetric(this, rmsdeg, rmsmm)
            err = [this.rmsarc(rmsdeg) rmsmm];
            em = sqrt(norm(err(~isnan(err))));
        end
        function r     = t4_resolvePair(this, f1, f2)
            %  @param f1, f2:  filenames.
            %  @return verbose info from t4_resolve.
            
            [~,r] = this.buildVisitor.t4_resolve( ...
                this.resolveTag, [f1 ' ' f2], 'options', '-v');
        end
        function arc   = rmsarc(~, rmsdeg)
            arc = 50*pi*rmsdeg/180; % arc at 50 mm from center of mass, per Avi Snyder
        end      
        
 		function this = AbstractT4ResolveError(varargin)
 			%% ABSTRACTT4RESOLVEERROR
 			%  @param .

 			this = this@mlfourdfp.AbstractSessionBuilder(varargin{:});
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'theImages', this.sessionData.tracerRevision('typ','fqfp')); % alternative cell-array composite
            addParameter(ip, 'indicesLogical', true, @islogical);
            parse(ip, varargin{:});     
            this.theImages_ = ip.Results.theImages;       
            this.indicesLogical = ip.Results.indicesLogical;
            this.buildVisitor_ = mlfourdfp.FourdfpVisitor;
 		end
 	end 
    
    %% PROTECTED
    
    properties (Access = protected)
        errMat_
        theImages_
    end
    
    methods (Abstract, Access = protected)
        updateLogging(this)
    end
    
    methods (Access = protected)
        function this  = ensureSizeOfIndicesLogical(this, ruler)
            if (1 == length(this.indicesLogical))
                this.indicesLogical = repmat(this.indicesLogical, size(ruler));
            end
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

