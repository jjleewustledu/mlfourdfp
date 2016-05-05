classdef T4ResolveReport 
	%% T4RESOLVEREPORT  

	%  $Revision$
 	%  was created 28-Feb-2016 15:06:11
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.0.0.307022 (R2016a) Prerelease for MACI64.
 	

	properties 		
 	end

	methods 
        
        function b = bar3(this, choice, t4r)
            switch (choice)
                case 'etas'
                    choice = '$\eta$';
                    mat = this.etas(t4r);
                case 'curves'
                    choice = '$\partial\partial$';
                    mat = this.curves(t4r);
                case 'z(etas)'
                    choice = sprintf('$z[\\eta]_{\\textrm{%s}}$', this.latexSafe(this.t4resolve_.imgregLogParser.filename));
                    mat = this.zEtas(t4r);
                case 'z(curves)'
                    choice = sprintf('$z[\\partial\\partial]_{\\textrm{%s}}$', this.latexSafe(this.t4resolve_.imgregLogParser.filename));
                    mat = this.zCurves(t4r);
                otherwise
                    error('mlfourdfp:unsupportedSwitchCase', 'T4ResolveReport.bar3.choice->%s is not supported', choice);
            end
            
            figure;
            b = bar3(mat);
            title(sprintf('%s %s', ...
                this.latexSafe(t4r.sessionData.sessionFolder), ...
                this.latexSafe(t4r.imgregLogParser.filename)), ...
                'Interpreter', 'latex');
            xlabel('frames');
            ylabel('frames');
            zlabel(choice, 'Interpreter', 'latex');
            this.colorbar3(b);
        end
        function p = pcolor(this, choice, t4r)
            switch (choice)
                case 'etas'
                    choice = '$\eta$';
                    mat = this.etas(t4r);
                case 'curves'
                    choice = '$\partial\partial$';
                    mat = this.curves(t4r);
                case 'z(etas)'
                    choice = sprintf('$z(\\eta)_{\\textrm{%s}}$', this.latexSafe(this.t4resolve_.imgregLogParser.filename));
                    mat = this.zEtas(t4r);
                case 'z(curves)'
                    choice = sprintf('$z(\\partial\\partial)_{\\textrm{%s}}$', this.latexSafe(this.t4resolve_.imgregLogParser.filename));
                    mat = this.zCurves(t4r);
                otherwise
                    error('mlfourdfp:unsupportedSwitchCase', 'T4ResolveReport.bar3.choice->%s is not supported', choice);
            end
            
            figure;
            p = pcolor(mat);
            p.EdgeColor = 'none';
            p.AmbientStrength = 0.6;
            title([sprintf('%s %s\n', ...
                this.latexSafe(t4r.sessionData.sessionFolder), ...
                this.latexSafe(t4r.imgregLogParser.filename)) ...
                choice], ...
                'Interpreter', 'latex');
            xlabel('frames');
            ylabel('frames');
            colorbar;
        end
        function s = surf(this, choice, t4r)
            switch (choice)
                case 'etas'
                    choice = '$\eta$';
                    mat = this.etas(t4r);
                case 'curves'
                    choice = '$\partial\partial$';
                    mat = this.curves(t4r);
                case 'z(etas)'
                    choice = sprintf('$z[\\eta]_{\\textrm{%s}}$', this.latexSafe(this.t4resolve_.imgregLogParser.filename));
                    mat = this.zEtas(t4r);
                case 'z(curves)'
                    choice = sprintf('$z[\\partial\\partial]_{\\textrm{%s}}$', this.latexSafe(this.t4resolve_.imgregLogParser.filename));
                    mat = this.zCurves(t4r);
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
                this.latexSafe(t4r.sessionData.sessionFolder), ...
                this.latexSafe(t4r.imgregLogParser.filename)), ...
                'Interpreter', 'latex');
            xlabel('frames');
            ylabel('frames');
            zlabel(choice, 'Interpreter', 'latex');
            colorbar;
        end
        function p = d(this, choice, t4r, t4r0)
            switch (choice)
                case 'etas'
                    choice = '$\mathbf{d}\eta$';
                    mat = this.etasDiff(t4r, t4r0);
                case 'curves'
                    choice = '$\mathbf{d}\partial\partial$';
                    mat = this.curvesDiff(t4r, t4r0);
                case 'z(etas)'
                    choice = sprintf('$z[\\mathbf{d}\\eta]_{\\textrm{%s}}$', this.latexSafe(this.t4resolve_.imgregLogParser.filename));
                    mat = this.zEtasDiff(t4r, t4r0);
                case 'z(curves)'
                    choice = sprintf('$z[\\mathbf{d}\\partial\\partial]_{\\textrm{%s}}$', this.latexSafe(this.t4resolve_.imgregLogParser.filename));
                    mat = this.zCurvesDiff(t4r, t4r0);
                otherwise
                    error('mlfourdfp:unsupportedSwitchCase', 'T4ResolveReport.bar3.choice->%s is not supported', choice);
            end
            
            figure;
            p = pcolor(mat);
            p.EdgeColor = 'none';
            title([sprintf('%s %s -\n%s %s:\n%s', ...
                this.latexSafe(t4r.sessionData.sessionFolder), ...
                this.latexSafe(t4r.imgregLogParser.filename), ...
                this.latexSafe(t4r0.sessionData.sessionFolder), ...
                this.latexSafe(t4r0.imgregLogParser.filename)) ...
                choice], ...
                'Interpreter', 'latex');
            xlabel('frames');
            ylabel('frames');
            colorbar;
        end
        
        function mat = etas(~, t4r)
            assert(isa(t4r, 'mlfourdfp.T4ResolveParser'));
            mat = nan(size(t4r.etas));
            e   = t4r.etas;
            for m = 1:size(mat,1)
                for n = 1:size(mat,2)
                    if (~isempty(e{m,n}))
                        mat(m,n) = e{m,n};
                    end
                end
            end
        end
        function mat = curves(~, t4r)
            assert(isa(t4r, 'mlfourdfp.T4ResolveParser'));
            mat = nan(size(t4r.curves));
            c   = t4r.curves;
            for m = 1:size(mat,1)
                for n = 1:size(mat,2)
                    if (~isempty(c{m,n}))
                        mat(m,n) = norm(c{m,n});
                    end
                end
            end
        end
        function mat = zEtas(this, t4r)
            assert(isa(t4r, 'mlfourdfp.T4ResolveParser'));
            mat = nan(size(t4r.etas));
            e   = t4r.etas;
            Ee  = dipmean(cell2mat(this.t4resolve_.etas));
            Se  = dipstd( cell2mat(this.t4resolve_.etas));
            for m = 1:size(mat,1)
                for n = 1:size(mat,2)
                    if (~isempty(e{m,n}))
                        mat(m,n) = (e{m,n} - Ee)/Se;
                    end
                end
            end
        end
        function mat = zCurves(this, t4r)
            assert(isa(t4r, 'mlfourdfp.T4ResolveParser'));
            mat = nan(size(t4r.curves));
            c   = t4r.curves;
            Ec  = dipmean(cell2mat(this.t4resolve_.curves));
            Sc  = dipstd( cell2mat(this.t4resolve_.curves));
            for m = 1:size(mat,1)
                for n = 1:size(mat,2)
                    if (~isempty(c{m,n}))
                        mat(m,n) = (norm(c{m,n}) - Ec)/Sc;
                    end
                end
            end
        end
        function mat = etasDiff(~, t4r, t4r0)
            assert(isa(t4r, 'mlfourdfp.T4ResolveParser'));
            mat = nan(size(t4r.etas));
            e   = t4r.etas;
            e0  = t4r0.etas;
            for m = 1:size(mat,1)
                for n = 1:size(mat,2)
                    if (~isempty(e{m,n}))
                        mat(m,n) = e{m,n} - e0{m,n};
                    end
                end
            end
        end
        function mat = curvesDiff(~, t4r, t4r0)
            assert(isa(t4r, 'mlfourdfp.T4ResolveParser'));
            mat = nan(size(t4r.curves));
            c   = t4r.curves;
            c0  = t4r0.curves;
            for m = 1:size(mat,1)
                for n = 1:size(mat,2)
                    if (~isempty(c{m,n}))
                        mat(m,n) = norm(c{m,n} - c0{m,n});
                    end
                end
            end
        end
        function mat = zEtasDiff(this, t4r, t4r0)
            assert(isa(t4r, 'mlfourdfp.T4ResolveParser'));
            mat = nan(size(t4r.etas));
            e   = t4r.etas;
            e0  = t4r0.etas;
            Se  = dipstd( cell2mat(this.t4resolve_.etas));
            for m = 1:size(mat,1)
                for n = 1:size(mat,2)
                    if (~isempty(e{m,n}))
                        mat(m,n) = (e{m,n} - e0{m,n})/Se;
                    end
                end
            end
        end
        function mat = zCurvesDiff(this, t4r, t4r0)
            assert(isa(t4r, 'mlfourdfp.T4ResolveParser'));
            mat = nan(size(t4r.curves));
            c   = t4r.curves;
            c0  = t4r0.curves;
            Sc  = dipstd( cell2mat(this.t4resolve_.curves));
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
		  
 		function this = T4ResolveReport(t4r)
 			%% T4RESOLVEREPORT
 			%  @param t4r is an instance of mlfourdfp.T4ResolveParser; it sets the baseline sample from which
            %  the mean and std are drawn for z-scores.

            assert(isa(t4r, 'mlfourdfp.T4ResolveParser'));
            this.t4resolve_ = t4r;
 		end
    end 

    %% PROTECTED
    
    properties (Access = protected)
        t4resolve_
    end
    
    methods (Access = protected)
        function s = latexSafe(~, s)
            s = strrep(s, '_', '\_');
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

