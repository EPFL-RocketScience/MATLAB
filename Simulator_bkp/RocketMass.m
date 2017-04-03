classdef RocketMass < handle
    % RocketMass
    % Calculate CG(t), Ir(t) and Ix(t) for given mass elements
    
    properties(SetAccess = private, GetAccess = public)
        MassElements = [];
        TotalMass = 0;
        CGx = 0;
        CGr = 0;
    end
    
    methods(Access = public)
        % CONSTRUCTOR
        function obj = RocketMass()
        end
        
        % EDIT methods
        function obj = addMass(obj, RME)
            % add
            % INPUT
            %   - RME : rocket mass element object
            obj.MassElements = [obj.MassElements, RME];
        end
        
        % CHECK methods
        function printAll(obj)
            % printAll
            % Print all mass elements
            display('MASS ELEMENTS : ')
            display('******')
            for massElement = obj.MassElements
                display(['id   : ' massElement.id])
                display(['mass : ' num2str(massElement.mass) ' kg'])
                display(['xPos : ' num2str(massElement.xPos) ' m'])
                display('******')
            end
        end
        
        % CG calculator
        function update(obj)
            obj.TotalMass = 0;
            relativeMassX = 0;
            relativeMassR = 0;
            for RME = obj.MassElements
                obj.TotalMass = obj.TotalMass + RME.mass;
                relativeMassX = relativeMassX + RME.mass*RME.xPos;
                relativeMassR =relativeMassR + RME.mass*RME.rPos;
            end
            obj.CGx = relativeMassX/obj.TotalMass;
            obj.CGr = relativeMassR/obj.TotalMass;
        end
    end
    
end
