classdef RocketAero < handle

    properties (SetAccess = private, GetAccess = public)
        x_CP        % longitudinal position of center of pressure 
    end
    
    properties (SetAccess = private, GetAccess = public)
        % private properties
        % Aerodynamic coefficients
        CP_nose; CNa_nose;
        CP_boattail; CNa_boattail;
        CP_fins; CNa_fins;
        
        % properties that are accessible through getters
        RG          % Rocket Geometry 
    end
    
    methods
        
        % CONSTRUCTOR
        function obj = RocketAero(RocketGeometry)
            if ~isa(RocketGeometry, 'RocketGeometry')
                error('RocketAero constructor requires an argument of type RocketGeometry.')
            end
            obj.RG = RocketGeometry;
        end
        
    end
    
    methods (Access = public)
        
        % Interface functions
        function update(obj)
            % update
            % update aerodynamic coefficients CP and Cna, calculate CP
            
            obj.CP_nose = obj.xCP_noseCone(obj.RG.LN);
            obj.CNa_nose = obj.cna_noseCone();
            obj.CP_boattail = obj.xCP_transition(obj.RG.LT, obj.RG.D, obj.RG.DT);
            obj.CNa_boattail = obj.cna_transition(obj.RG.D, obj.RG.D, obj.RG.DT);
            obj.CP_fins = obj.xCP_fins(obj.RG.XRT,obj.RG.CR, obj.RG.CT);
            obj.CNa_fins = obj.cna_fins(obj.RG.NF, obj.RG.D/2, obj.RG.S,...
                obj.RG.D, obj.RG.CR, obj.RG.CT, obj.RG.LF);
            
            obj.x_CP = (... 
                        obj.CP_nose*obj.CNa_nose...
                        + (obj.CP_boattail+obj.RG.XT)*obj.CNa_boattail...
                        + (obj.CP_fins+obj.RG.XF)*obj.CNa_fins...
                        ) / (obj.CNa_nose+obj.CNa_boattail+obj.CNa_fins);
        end
        
    end
    
    methods (Access = private, Static = true)
    
        % POSITIONS OF CENTER OF PRESSURE and NORMAL AERODYNAMIC
        % COEFFICIENT DERIVATIVES (NACD)
        % - body element (boattail, shoulder)
        % - nose with a conical shape
        % - nose with an ogive shape
        % - fins of trapezoidal geometry
        
        function xcp = xCP_transition(L, d1, d2)
            % xCP_transition
            % calculate position of cp on a body element (body of revolution)
            % from tip of element 
            % INPUT:
            %   - L     : length of part [m]
            %   - d1    : diameter before transition
            %   - d2    : diameter after transition
            % RETURN:
            %   - xcp   : center of pressure position along x axis [m]
            xcp = L/3*(1+(1-d1/d2)/(1-(d1/d2)^2));
        end
        
        function  xcp = xCP_noseCone(L)
            % xCP_noseCone
            % calculate cp position of a conical nose from tip of nose
            % INPUT:
            %   - L     : nose length [m]
            % RETURN:
            %   - xcp   : center of pressure position along x axis [m]
            xcp = 2/3*L;
        end
        
        function xcp = xCP_noseOgive(L)
            % xCP_noseOgive
            % calculate cp position of an ogive shaped nose from tip
            % INPUT:
            %   -L      : nose length [m]
            % RETURN:
            %   - xcp   : center of pressure position along x axis [m]
            xcp = 0.466*L;
        end
        
        function xcp = xCP_fins(Xt, Cr, Ct)
            % xCP_fins
            % calculate cp position of any type of trapezoidal fins
            % The position is the intersection of the mean aerodynamic
            % chord and the quarter cord line
            % INPUT: 
            %   - Xt    : longitudinal distance between leading edge of
            %             root chord and leading edge of tip chord. [m]
            %   - Cr    : root chord length [m]
            %   - Ct    : tip chord length [m]
            % RETURN:
            %   - xcp   : center of pressure position along x axis [m]
            xcp = Xt/3*(Cr+2*Ct)/(Cr+Ct) + 1/6*((Cr+Ct)-(Cr*Ct)/(Cr+Ct));
        end    
        
        % NORMAL AERODYNAMIC COEFFICIENT DERIVATIVES (NACD)
        
        function cna = cna_noseCone()
            % cna_noseCone
            % calculate NACD of nose with a conical or ogive shape
            % RETRUN:
            %   - cna   : ACD [1/rad]
            cna = 2;
        end
        
        function cna = cna_transition(d, d1, d2)
            % cna_transition
            % calculate NACD of a body element (shoulder, boattail, tube)
            % INPUTS:
            %   - d     : reference diameter (diameter at base of cone)
            %   - d1    : diameter before transition
            %   - d2    : diameter after transition
            % RETRUN:
            %   - cna   : ACD [1/rad]      
            cna = 2/d^2*(d2^2-d1^2);    
        end
        
        function cna = cna_fins(N, rt, S, d, Cr, Ct, l)
            % cna_fins
            % calculate NACD of N (3 or 4) trapezoidal fins.
            % INPUTS:
            %   - N     : number of fins
            %   - rt    : radius of body at fin attachement
            %   - S     : span length of fin
            %   - d     : reference diameter (diameter at base of cone)
            %   - Cr    : fin root chord
            %   - Ct    : fin tip chord
            %   - Xt    : longitudinal distance between leading edge of
            %             root chord and leading edge of tip chord. [m]
            
            if (N<3 || N>4)
                error('cna_fins can only compute the ACD for N = 3 or 4 fins.')
            end
            
            % NACD with correction factor for the fins in the presence of a
            % body.
            cna = (1+rt/(S+rt))*(N*4*(S/d)^2)/(1+sqrt(1+(2*l/(Cr+Ct))^2));
        end
        
    end 
    
end