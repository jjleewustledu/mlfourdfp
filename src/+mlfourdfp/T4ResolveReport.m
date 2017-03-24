classdef T4ResolveReport 
	%% T4RESOLVEREPORT  

	%  $Revision$
 	%  was created 28-Feb-2016 15:06:11
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.0.0.307022 (R2016a) Prerelease for MACI64.
 	

	methods   
        function b = bar3(this, varargin)
            ip = inputParser;
            addRequired(ip, 'choice', @ischar);
            addOptional(ip, 't4rp', this.t4RParser_, @(x) isa(x, 'mlfourdfp.T4ResolveParser'));
            parse(ip, varargin{:});
            choice = ip.Results.choice;
            t4rp    = ip.Results.t4rp;
            
            switch (choice)
                case 'etas'
                    choice = '$\eta$';
                    mat = this.etas(t4rp);
                case 'curves'
                    choice = '$\partial\partial$';
                    mat = this.curves(t4rp);
                case 'z(etas)'
                    choice = sprintf('$z[\\eta]_{\\textrm{%s}}$', this.latexSafe(this.t4RParser_.imgregLogParser.filename));
                    mat = this.zEtas(t4rp);
                case 'z(curves)'
                    choice = sprintf('$z[\\partial\\partial]_{\\textrm{%s}}$', this.latexSafe(this.t4RParser_.imgregLogParser.filename));
                    mat = this.zCurves(t4rp);
                otherwise
                    error('mlfourdfp:unsupportedSwitchCase', 'T4ResolveReport.bar3.choice->%s is not supported', choice);
            end
            
            figure;
            b = bar3(mat);
            title(sprintf('%s %s', ...
                this.latexSafe(t4rp.sessionData.sessionFolder), ...
                this.latexSafe(t4rp.imgregLogParser.filename)), ...
                'Interpreter', 'latex');
            xlabel('frames');
            ylabel('frames');
            zlabel(choice, 'Interpreter', 'latex');
            this.colorbar3(b);
        end
        function p = pcolor(this, varargin)
            ip = inputParser;
            addRequired(ip, 'choice', @ischar);
            addOptional(ip, 't4rp', this.t4RParser_, @(x) isa(x, 'mlfourdfp.T4ResolveParser'));
            parse(ip, varargin{:});
            choice = ip.Results.choice;
            t4rp    = ip.Results.t4rp;
            
            switch (choice)
                case 'etas'
                    choice = '$\eta$';
                    mat = this.etas(t4rp);
                case 'curves'
                    choice = '$\partial\partial$';
                    mat = this.curves(t4rp);
                case 'z(etas)'
                    choice = sprintf('$z(\\eta)_{\\textrm{%s}}$', this.latexSafe(this.t4RParser_.imgregLogParser.filename));
                    mat = this.zEtas(t4rp);
                case 'z(curves)'
                    choice = sprintf('$z(\\partial\\partial)_{\\textrm{%s}}$', this.latexSafe(this.t4RParser_.imgregLogParser.filename));
                    mat = this.zCurves(t4rp);
                otherwise
                    error('mlfourdfp:unsupportedSwitchCase', 'T4ResolveReport.bar3.choice->%s is not supported', choice);
            end
            
            figure;
            p = pcolor(mat);
            p.EdgeColor = 'none';
            p.AmbientStrength = 0.6;
            title([sprintf('%s %s\n', ...
                this.latexSafe(t4rp.sessionData.sessionFolder), ...
                this.latexSafe(t4rp.imgregLogParser.filename)) ...
                choice], ...
                'Interpreter', 'latex');
            xlabel('frames');
            ylabel('frames');
            colorbar;
        end
        function s = surf(this, varargin)
            ip = inputParser;
            addRequired(ip, 'choice', @ischar);
            addOptional(ip, 't4rp', this.t4RParser_, @(x) isa(x, 'mlfourdfp.T4ResolveParser'));
            parse(ip, varargin{:});
            choice = ip.Results.choice;
            t4rp    = ip.Results.t4rp;
            
            switch (choice)
                case 'etas'
                    choice = '$\eta$';
                    mat = this.etas(t4rp);
                case 'curves'
                    choice = '$\partial\partial$';
                    mat = this.curves(t4rp);
                case 'z(etas)'
                    choice = sprintf('$z[\\eta]_{\\textrm{%s}}$', this.latexSafe(this.t4RParser_.imgregLogParser.filename));
                    mat = this.zEtas(t4rp);
                case 'z(curves)'
                    choice = sprintf('$z[\\partial\\partial]_{\\textrm{%s}}$', this.latexSafe(this.t4RParser_.imgregLogParser.filename));
                    mat = this.zCurves(t4rp);
                otherwise
                    error('mlfourdfp:unsupportedSwitchCase', 'T4ResolveReport.bar3.choice->%s is not supported', choice);
            end
            
            figure;
            s = surf(mat);
            %s.EdgeColor = [1 1 1];
            %s.AmbientStrength = 0.6;
            %camlight(110, 70);
            %brighten(0.6);
            s.EdgeColor = [1 1 1];
            s.AmbientStrength = 0.6;
            title(sprintf('%s %s', ...
                this.latexSafe(t4rp.sessionData.sessionFolder), ...
                this.latexSafe(t4rp.imgregLogParser.filename)), ...
                'Interpreter', 'latex');
            xlabel('frames');
            ylabel('frames');
            zlabel(choice, 'Interpreter', 'latex');
            colorbar;
        end
        function p = d(this, varargin)
            ip = inputParser;
            addRequired(ip, 'choice', @ischar);
            addOptional(ip, 't4rp',  this.t4RParser_, @(x) isa(x, 'mlfourdfp.T4ResolveParser'));
            addOptional(ip, 't4rp0', [],              @(x) isa(x, 'mlfourdfp.T4ResolveParser'));
            parse(ip, varargin{:});
            choice = ip.Results.choice;
            t4rp    = ip.Results.t4rp;
            t4rp0   = ip.Results.t4rp0;
            
            switch (choice)
                case 'etas'
                    choice = '$\mathbf{d}\eta$';
                    mat = this.etasDiff(t4rp, t4rp0);
                case 'curves'
                    choice = '$\mathbf{d}\partial\partial$';
                    mat = this.curvesDiff(t4rp, t4rp0);
                case 'z(etas)'
                    choice = sprintf('$z[\\mathbf{d}\\eta]_{\\textrm{%s}}$', this.latexSafe(this.t4RParser_.imgregLogParser.filename));
                    mat = this.zEtasDiff(t4rp, t4rp0);
                case 'z(curves)'
                    choice = sprintf('$z[\\mathbf{d}\\partial\\partial]_{\\textrm{%s}}$', this.latexSafe(this.t4RParser_.imgregLogParser.filename));
                    mat = this.zCurvesDiff(t4rp, t4rp0);
                otherwise
                    error('mlfourdfp:unsupportedSwitchCase', 'T4ResolveReport.bar3.choice->%s is not supported', choice);
            end
            
            figure;
            p = pcolor(mat);
            p.EdgeColor = 'none';
            title([sprintf('%s %s -\n%s %s:\n%s', ...
                this.latexSafe(t4rp.sessionData.sessionFolder), ...
                this.latexSafe(t4rp.imgregLogParser.filename), ...
                this.latexSafe(t4rp0.sessionData.sessionFolder), ...
                this.latexSafe(t4rp0.imgregLogParser.filename)) ...
                choice], ...
                'Interpreter', 'latex');
            xlabel('frames');
            ylabel('frames');
            colorbar;
        end
        
        function mat = etas(~, t4rp)
            assert(isa(t4rp, 'mlfourdfp.T4ResolveParser'));
            mat = nan(size(t4rp.etas));
            e   = t4rp.etas;
            for m = 1:size(mat,1)
                for n = 1:size(mat,2)
                    if (~isempty(e{m,n}))
                        mat(m,n) = e{m,n};
                    end
                end
            end
        end
        function mat = curves(~, t4rp)
            assert(isa(t4rp, 'mlfourdfp.T4ResolveParser'));
            mat = nan(size(t4rp.curves));
            c   = t4rp.curves;
            for m = 1:size(mat,1)
                for n = 1:size(mat,2)
                    if (~isempty(c{m,n}))
                        mat(m,n) = norm(c{m,n});
                    end
                end
            end
        end
        function mat = zEtas(this, t4rp)
            assert(isa(t4rp, 'mlfourdfp.T4ResolveParser'));
            mat  = nan(size(t4rp.etas));
            e    = t4rp.etas;
            mat1 = mycell2mat(this.t4RParser_.etas);
            Ee   = dipmean(mat1(~isnan(mat1)));
            Se   = dipstd( mat1(~isnan(mat1)));
            for m = 1:size(mat,1)
                for n = 1:size(mat,2)
                    if (~isempty(e{m,n}))
                        mat(m,n) = (e{m,n} - Ee)/Se;
                    end
                end
            end
        end
        function mat = zCurves(this, t4rp)
            assert(isa(t4rp, 'mlfourdfp.T4ResolveParser'));
            mat = nan(size(t4rp.curves));
            c   = t4rp.curves;
            c1  = mycell2mat(this.t4RParser_.curves);
            Ec  = dipmean(c1(~isnan(c1)));
            Sc  = dipstd( c1(~isnan(c1)));
            for m = 1:size(mat,1)
                for n = 1:size(mat,2)
                    if (~isempty(c{m,n}))
                        mat(m,n) = (norm(c{m,n}) - Ec)/Sc;
                    end
                end
            end
        end
        function mat = etasDiff(~, t4rp, t4rp0)
            assert(isa(t4rp, 'mlfourdfp.T4ResolveParser'));
            mat = nan(size(t4rp.etas));
            e   = t4rp.etas;
            e0  = t4rp0.etas;
            for m = 1:size(mat,1)
                for n = 1:size(mat,2)
                    if (~isempty(e{m,n}))
                        mat(m,n) = e{m,n} - e0{m,n};
                    end
                end
            end
        end
        function mat = curvesDiff(~, t4rp, t4rp0)
            assert(isa(t4rp, 'mlfourdfp.T4ResolveParser'));
            mat = nan(size(t4rp.curves));
            c   = t4rp.curves;
            c0  = t4rp0.curves;
            for m = 1:size(mat,1)
                for n = 1:size(mat,2)
                    if (~isempty(c{m,n}))
                        mat(m,n) = norm(c{m,n} - c0{m,n});
                    end
                end
            end
        end
        function mat = zEtasDiff(this, t4rp, t4rp0)
            assert(isa(t4rp, 'mlfourdfp.T4ResolveParser'));
            mat = nan(size(t4rp.etas));
            e    = t4rp.etas;
            e0   = t4rp0.etas;            
            mat1 = mycell2mat(this.t4RParser_.etas);
            Se   = dipstd( mat1(~isnan(mat1)));
            for m = 1:size(mat,1)
                for n = 1:size(mat,2)
                    if (~isempty(e{m,n}))
                        mat(m,n) = (e{m,n} - e0{m,n})/Se;
                    end
                end
            end
        end
        function mat = zCurvesDiff(this, t4rp, t4rp0)
            assert(isa(t4rp, 'mlfourdfp.T4ResolveParser'));
            mat = nan(size(t4rp.curves));
            c   = t4rp.curves;
            c0  = t4rp0.curves;            
            c1  = mycell2mat(this.t4RParser_.curves);
            Sc  = dipstd( c1(~isnan(c1)));
            for m = 1:size(mat,1)
                for n = 1:size(mat,2)
                    if (~isempty(c{m,n}))
                        mat(m,n) = (norm(c{m,n} - c0{m,n}))/Sc;
                    end
                end
            end
        end
        function colorbar3(~, b)
            colorbar;
            for k = 1:length(b)
                zdata = b(k).ZData;
                b(k).CData = zdata;
                b(k).FaceColor = 'interp';
            end
        end
		  
 		function this = T4ResolveReport(t4rp)
 			%% T4RESOLVEREPORT
 			%  @param t4rp is an instance of mlfourdfp.T4ResolveParser; it sets the baseline sample from which
            %  the mean and std are drawn for z-scores.

            assert(isa(t4rp, 'mlfourdfp.T4ResolveParser'));
            this.t4RParser_ = t4rp;
 		end
    end 

    %% PROTECTED
    
    properties (Access = protected)
        t4RParser_
    end
    
    methods (Access = protected)
        function s = latexSafe(~, s)
            s = strrep(s, '_', '\_');
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

