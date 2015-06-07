classdef mapfun < unbndfun

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAPFUN Class Description:
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function obj = mapfun(op, data, pref)
            
            obj.domain = [-1, 1];
            
            % Parse inputs.
            if ( nargin < 1 )
                % No input arguments; return an empty MAPFUN.
                obj.onefun = [];
                obj.mapping = [];
                return
            end

            if ( (nargin < 2) || isempty(data) )
                data = struct();
            end

            if ( (nargin < 3) || isempty(pref) )
                pref = chebfunpref();
            else
                pref = chebfunpref(pref);
            end

            data = parseDataInputs(data, pref);
            
            if ( ~isempty(data.mapping) )
                op = @(x) op(data.mapping.For(x));
                obj.mapping = data.mapping;
            end

            obj.onefun = smoothfun.constructor(op, data, pref);
            
        end
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% NON-STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false ) 
        
        function f = flipud(f)
            f.onefun = flipud(f.onefun);
        end
        
        function f = plus(f, g)
            if ( isnumeric(f) )
                f = plus(g, f);
            elseif ( isnumeric(g) )
                f.onefun = f.onefun + g;
            else
                error('not supported');
            end
        end
        
        function f = minus(f, g)
            f = plut(f, -g);
        end
        
    end
    
end

function data = parseDataInputs(data, pref)
%PARSEDATAINPUTS   Parse inputs from the DATA structure and assign defaults.

if ( isa(data, 'mapping') )
    data = struct('mapping', data);
elseif ( ~isfield(data, 'mapping') )
    data.mapping = [];
end

if ( isempty(data.mapping) && isfield(pref.mapPrefs) )
    data.mapping = pref.mapPrefs.mapping;
end

end
