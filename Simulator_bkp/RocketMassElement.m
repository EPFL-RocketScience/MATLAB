classdef RocketMassElement
    % RocketMassElement
    %
    % PROPERTIES
    % id - (str)        name of part
    % mass - (real)     mass of part, f(t)
    % Ix - (real)       cylinder moment of inertia around x axis (longitudinal)
    % Ir - (real)       cylinder moment of inertia around r axis
    % massRate - (real) rate of mass loss, f(t)
    % xPos - (real)     position of CG along x axis
    % rPos - (real)     position on CG along r axis
    % D - (real)        Diameter of part
    % L - (real)        Length of part 

    properties(SetAccess = private, GetAccess = public)
        id
        mass
        Ix
        Ir
        xPos
        rPos
        D
        L
    end
    
    methods (Access = public)
        
        % CONSTRUCTOR
        function obj = RocketMassElement(id, mass, D, L, xPos, rPos)
            obj.id = id;
            obj.mass = mass;
            obj.D = D;
            obj.L = L;
            obj.Ix = 0.5*(D/2)^2*mass;
            obj.Ir = mass/12*(3*(D/2)^2+ L^2);
            obj.xPos = xPos;
            obj.rPos = rPos;
        end
        
        % HELPERS
        function str = getStr(obj, propName)
            prop = get(obj, propName);
            if isa(prop, 'function_handle')
                str = func2str(prop);
            elseif isa(prop, 'char')
                str = prop;
            else
                str = num2str(prop);
            end
        end
        
    end
    
    methods (Access = private)
        
        function bool = isFunction(obj, h)
            % Check if value is a function handle
            % INPUT : 
            %   - h : value or function handle
            % RETURN :
            %   - bool : h is a function handle (bool = 1), else (bool = 0)
            
            bool = 0;
            if isa(h, 'function_handle')
                bool = 1;
            end
        end
        
        function h = makeFunction(obj, val)
            % If value is not a function, turn it into one. A value becomes
            % a constant function of time.
            % INPUT :
            %   - val : value to be tested and transformed if needed.
            % RETURN :
            %   - h : function handle
            
            if obj.isFunction(val)
                h = val;
            else
                h = @(t) val;
            end
        end
        
    end

end