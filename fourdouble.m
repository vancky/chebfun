classdef fourdouble < chebdouble
%FOURDOUBLE   Fourier double class. For example, DIFF means Fourier difference.
%
%   See the CHEBDOUBLE class for details.
%
%   This class in intended solely as a worker-class for PDE15s.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    
    methods
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%  CONSTRUCTOR  %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = fourdouble(varargin)
            
            % Call the CHEBDOUBLE constructor:
            obj = obj@chebdouble(varargin{:});
            
        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  DIFF  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function u = diff(u, k)
            %DIFF   Compute the k-th derivative of u using Fourier
            % differentiation matrices defined by diffmat.
            
            % Store the diffmat D as a persistent variable to allow speeding up
            % if we work with the same discretization at multiple time steps.
            % Note that the matrix D is independent of the domain, since it is
            % scaled separately below.
            persistent D
            
            % Assume first-order derivative
            if ( nargin == 1 )
                k = 1;
            end
            
            N = length(u.values);
            
            % Construct D if we don't match a previous discretization.
            if ( isempty(D) || numel(D) < k || size(D{k}, 1) ~= N )
                D{k} = fourtech.diffmat(N, k); % Diffmat
            end
            
            % Interval scaling
            c = 2*pi/diff(u.domain);     
            
            % Muliplying by the kth-order differentiation matrix
            u.values = c^k*(D{k}*u.values);
            
            % Update the difforder:
            u.diffOrder = u.diffOrder + k;
            
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SUM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The differential operators
        function I = sum(u, a, b)
            %SUM  Compute the integral of u.
            
            persistent W
            
            if ( isempty(W) )
                W = {};
            end
            
            % Extract the data:
            N = length(u.values);
            
            % Deal with the 3 args case. This can be integrating a sub-domain or
            % indefinite integration. (Or integrating the whole domain...)
            if ( nargin == 3 )
                x = fourpts(N, u.domain);
                if ( length(b) > 1 )
                    if ( ~all(b == x) )
                        error('CHEBFUN:pde15s:sumb', ...
                            ['Limits in sum must be scalars or the ', ...
                            'independent space variable (typically ''x'').']);
                    elseif ( a < x(1) )
                        error('CHEBFUN:pde15s:sumint', ...
                            'Limits of integration outside of domain.');
                    end
                    I = cumsum(u);
                    I = I - feval(I, a);
                    return
                elseif ( length(a) > 1 )
                    if ( ~all(a == x) )
                        error('CHEBFUN:pde15s:suma', ...
                            ['Limits in sum must be scalars or the ', ...
                            'independent space variable (typically ''x'').']);
                    elseif ( b > x(end) )
                        error('CHEBFUN:pde15s:sumint', ...
                            'Limits of integration outside of domain.');
                    end
                    I = cumsum(u);
                    I = feval(I, b) - I;
                    return
                elseif ( a ~= x(1) || b ~= x(end) )
                    if ( a < x(1) || b > x(end) )
                        error('CHEBFUN:pde15s:sumint', ...
                            'Limits of integration outside of domain.');
                    end
                    I = cumsum(u);
                    I = feval(I, b) - feval(I, a);
                    return
                else
                    % Sum(u, a, b) is the same as below!
                end
            end
            
            % Retrieve or compute weights::
            if ( N > 5 && numel(W) >= N && ~isempty(W{N}) )
                % Weights are already in storage!
            else
                c = diff(u.domain)/2; % Interval scaling.
                W{N} = c*fourtech.quadwts(N);
            end
            
            % Find the sum by muliplying by the weights vector:
            I = W{N}*u;
        end
        
        function I = integral(varargin)
            I = sum(varargin{:});
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  CUMSUM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The differential operators
        function u = cumsum(u)
            %CUMSUM   Compute the indefinite integral of the Chebyshev
            %         interpolant to u.
            
            % TODO: Add support.
            error('No support here.')
            
            persistent C

            % Extract the data:
            N = length(u.values);
            c = diff(u.domain)/2; % Interval scaling.
            
            % Compute cumsum matrix:
            if ( numel(C) ~= N )
                C = colloc2.cumsummat(N);
            end
            
            % Compute the indefinite integral:
            u.values = c*(C*u.values);
            
            % Update the difforder:
            u.diffOrder = u.diffOrder - 1;
            
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FRED  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function u = fred(K, u)
            %FRED  Fredholm operator with kernel K.
            %   FRED(K, U) computes the action of the Fredholm operator with
            %   kernel K on the Chebyshev interpolant to the points in the
            %   vector U.
            
            % TODO: Add support.
            error('No support here.')
            
            persistent X W
            if ( isempty(W) )
                X = {};
                W = {};
            end
            
            % Extract the data:
            N = length(u);
            
            % Retrieve or compute weights::
            if ( N > 5 && numel(W) >= N && ~isempty(W{N}) )
                % Weights are already in storage!
            else
                [X{N}, W{N}] = chebpts(N, u.domain);
            end
            
            % The Fredholm operator:
            [xx, yy] = ndgrid(X{N});
            u = K(xx, yy) * (W{N}.'.*u.values);
            
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  VOLT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function u = volt(K, u)
            %VOLT  Volterra operator with kernel K.
            %   VOLT(K, U) computes the action of the Volterra operator with
            %   kernel K on the Chebyshev interpolant to the points in the
            %   vector U.
            
            % TODO: Add support.
            error('No support here.')
            
            persistent X C
            
            % Extract the data:
            N = length(u.values);
            c = diff(u.domain)/2; % Interval scaling.
            
            % Compute cumsum matrix:
            if ( numel(C) ~= N )
                X = chebpts(N, u.domain);
                C = cumsummat(N);
            end
            
            % The Fredholm operator:
            [xx, yy] = ndgrid(X{N});
            u = K(xx, yy) * C * (c*u.values);
            
        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FEVAL  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function out = feval(u, y)
            %FEVAL  Evaluate polynomial interpolant of data {X_four, U} at a
            % point y using barycentric interpolation.
            
            % TODO: Add support.
            error('No support here.')
            
            [x, w, v] = fourpts(length(u.values), u.domain);
            out = fourtech.bary(y, u.values, x, v);
        end
                
    end
    
end


