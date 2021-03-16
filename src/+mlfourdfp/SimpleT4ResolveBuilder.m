classdef SimpleT4ResolveBuilder < mlfourdfp.AbstractBuilder
	%% SIMPLET4RESOLVEBUILDER avoids the complexity of AbstractT4ResolveBuilder and its implementations.

	%  $Revision$
 	%  was created 30-Dec-2020 10:44:06 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.9.0.1538559 (R2020b) Update 3 for MACI64.  Copyright 2020 John Joowon Lee.
 	
    properties
        imageRegLog
        resolveLog
    end
    
    properties (Dependent)        
        blurArg
        maskForImages
        resolveTag
        theImages
    end

	methods 
        
        %% GET, SET
        
        function g    = get.blurArg(this)
            g = this.blurArg_;
        end
        function this = set.blurArg(this, s)
            assert(isnumeric(s));
            assert(~isempty(s));
            this.blurArg_ = s;
        end
        function g    = get.maskForImages(this)
            g = this.maskForImages_;
        end
        function this = set.maskForImages(this, s)
            assert(iscell(s))
            this.maskForImages_ = s;
            this.maskForImages_ = cellfun(@(x) mlfourd.ImagingContext2(x), this.maskForImages_, 'UniformOutput', false);
            cellfun(@(x) assert(contains(x.filesuffix, '.4dfp.hdr')), this.maskForImages_, 'UniformOutput', false)
        end
        function g    = get.resolveTag(this)
            g = this.resolveTag_;
        end
        function this = set.resolveTag(this, s)
            assert(ischar(s));
            assert(~isempty(s));
            this.resolveTag_ = s;
        end
        function g    = get.theImages(this)
            g = this.theImages_;
        end
        function this = set.theImages(this, s)
            assert(iscell(s))
            this.theImages_ = s;
            this.theImages_ = cellfun(@(x) mlfourd.ImagingContext2(x), this.theImages_, 'UniformOutput', false);
            cellfun(@(x) assert(contains(x.filesuffix, '.4dfp.hdr')), this.theImages_, 'UniformOutput', false)
        end
        
        %%
		  
 		function this = SimpleT4ResolveBuilder(varargin)
 			%% SIMPLET4RESOLVEBUILDER
            %  @param blurArg default := 0.
 			%  @param maskForImages is a cell-array of objects understood by ImagingContext2; default is empty.
            %  @param resolveTag is char; default ~ 'op_theFirstImageFileprefix'.
            %  @param theImages is a cell-array of objects understood by ImagingContext2.
 			
            this = this@mlfourdfp.AbstractBuilder(varargin{:});
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'blurArg', 0, @isnumeric);
            addParameter(ip, 'maskForImages', [], @iscell);
            addParameter(ip, 'resolveTag', '', @ischar);
            addParameter(ip, 'theImages', [], @iscell);
            parse(ip, varargin{:});  
            ipr = ip.Results;
            
            this.theImages_ = ipr.theImages;
            this.theImages_ = cellfun(@(x) mlfourd.ImagingContext2(x), this.theImages_, 'UniformOutput', false);
            cellfun(@(x) assert(contains(x.filesuffix, '.4dfp.hdr')), this.theImages_, 'UniformOutput', false)  
            
            this.maskForImages_ = ipr.maskForImages;
            if isempty(this.maskForImages_)
                this.maskForImages_ = repmat({'none.4dfp.hdr'}, size(this.theImages_));
            end
            this.maskForImages_ = cellfun(@(x) mlfourd.ImagingContext2(x), this.maskForImages_, 'UniformOutput', false);
            cellfun(@(x) assert(contains(x.filesuffix, '.4dfp.hdr')), this.maskForImages_, 'UniformOutput', false)  
            
            this.blurArg_ = ipr.blurArg;
            if isscalar(this.blurArg_)
                this.blurArg_ = this.blurArg_*ones(size(this.theImages_));
            end
            
            this.resolveTag_ = ipr.resolveTag;
            if isempty(this.resolveTag_)
                this.resolveTag_ = ['op_' this.theImages_{1}.fileprefix];
            end
            
            targ = this.theImages_{1}.fileprefix;
            this.imageRegLog = [targ '_imageReg.log'];
            this.resolveLog = [targ '_resolve.log'];
        end
        
        function this = resolve(this)
            fps = this.imageReg();
            this = this.resolveAndPaste(fps); 
            this = this.finalize();
        end
        function fps = imageReg(this)
            len = length(this.theImages);
            fps = cellfun(@(x) x.fileprefix, this.theImages, 'UniformOutput', false);
            bfps = cellfun(@(x) x.fileprefix, this.blurredImages(), 'UniformOutput', false);
            mfps = cellfun(@(x) x.fileprefix, this.maskForImages, 'UniformOutput', false);
            for m = 1:len
                for n = 1:len
                    if (m ~= n)
                        t4 = this.buildVisitor.filenameT4(fps{n}, fps{m});
                        if (~this.valid_t4(t4))
                            this.buildVisitor.align_multiSpectral( ...
                                'dest',       bfps{m}, ...
                                'source',     bfps{n}, ...
                                'destMask',   mfps{m}, ...
                                'sourceMask', mfps{n}, ...
                                't4',         t4, ...
                                't4img_4dfp', false, ...
                                'log',        this.imageRegLog);
                        end
                        % t4_resolve requires an idiomatic naming convention for t4 files,
                        % based on the names of image files
                        % e. g., image1_to_image2_t4                            
                    end
                end
            end 
        end
        function this = resolveAndPaste(this, fps)
            % Must use short fileprefixes in calls to t4_resolve to avoid filenaming error by t4_resolve.  E.g.:
            % t4_resolve: /data/nil-bluearc/raichle/PPGdata/jjlee2/HYGLY28/V1/FDG_V1-NAC/E8/fdgv1e8r1_frame1_to_/data/nil-bluearc/raichle/PPGdata/jjlee2/HYGLY28/V1/FDG_V1-NAC/E8/fdgv1e8r1_frame8_t4 read error
            this.buildVisitor.t4_resolve( ...
                this.resolveTag, cell2str(fps, 'AsRow', true), ...
                'options', '-v -m -s', 'log', this.resolveLog);
            this.t4imgAll(fps); % transform ipr.dest on this.resolveTag
        end
        function this = finalize(this)
        end
    end     

    %% PROTECTED
    
    properties (Access = protected)
        blurArg_
        maskForImages_
        resolveTag_
        theImages_
    end
    
    methods (Access = protected)
        function b = blurredImages(this)
            b = cell(size(this.theImages));
            for ib = 1:length(b)
                b{ib} = this.theImages{ib}.blurred(this.blurArg(ib));
                if ~isfile(b{ib}.fqfilename)
                    b{ib}.save();
                end
            end
        end
        function this = t4imgAll(this, fps)
            tag = this.resolveTag;
            for f = 1:length(fps)
                fp = fps{f};
                if ~contains(tag, fp)
                    this.buildVisitor.t4img_4dfp( ...
                        sprintf('%s_to_%s_t4', fp, tag), ...
                        fp, ...
                        'out', sprintf('%s_%s', fp, tag), ...
                        'options', ['-O' fps{1}]);
                end
            end
        end 
        function tf = valid_t4(~, t4)
            if ~isfile(t4)
                tf = false;
                return
            end
            d = dir(t4);
            if 0 == d.bytes
                tf = false;
                return
            end
            tf = true;
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

