classdef RocketGeometry < handle
    % This class defines outer rocket geometry (fin sizes, body length,...)
    
    properties (SetAccess = public, GetAccess = public)
        CONE    % conical (0) or ogive (1)
        LN      % nose length
        LB      % body length
        LT      % tail length
        DT      % rear tail diameter
        XT      % fore boattail position from nose tip
        D       % body diameter
        XF      % leading edge fin position from nose tip
        CR      % fin root chord
        CT      % fin tip chord
        LF      % length of fin mid-chord line
        XRT     % longitudinal length of fin root to fin tip along leading edge
        S       % fin semi-span
        NF      % number of fins
    end
    
    methods (Access = public)
        % CONSTRUCTOR
        
        function obj = RocketGeometry()
        end
        
        % PROPERTY ASSIGNEMENT FUNCTIONS
        
        function setGeometryFromFile(obj, filepath)
            % Create rocket geometry object from parameter file 
            % INPUT : 
            %   - filepath (string) path to parameter file 
            
            CADFILE = importdata(filepath);
            
            obj.LN  = CADFILE.data(13)  ;    %length of nose  
            obj.D  = CADFILE.data(8)   ; %diameter at base of nose  
            obj.LB = CADFILE.data(9)    ; %length of Body tube 
            obj.DT  =  CADFILE.data(12)   ;  %diameter at rear of transition  
            obj.LT  = CADFILE.data(11)    ;%length of transition  
            obj.XT  =  CADFILE.data(9)+obj.LN  ;   %distance from tip of nose to front of transition (=tube lenght + nose as transition at the end) 
            obj.CR  =  CADFILE.data(1)    ; %fins root chord  
            obj.CT  =  CADFILE.data(2) ;    %fins tip chord  
            obj.S  =  (obj.CR-obj.CT)/(1/tand(CADFILE.data(4))+1/tand(CADFILE.data(3)))   ;%fins semispan  
            obj.LF  =  (obj.S^2+(obj.CR/2-obj.CT/2-obj.S*tand(CADFILE.data(3)))^2)^.5  ;   %length of fin mid-chord line  
            obj.XRT  = obj.S/tand(CADFILE.data(3))     ; %distance between fin root leading edge and fin tip leading edge parallel to body  
            obj.XF  = obj.XT-obj.CR-CADFILE.data(7)     ;%distance from nose tip to fin root chord leading edge  
            obj.NF  = CADFILE.data(5)      ;%number of fins  
        end
        
        function setGeometry(obj, LN, D, LB, DT, LT, XT, CR, CT, S, LF, XRT, XF, NF)
        % Create rocket geometry object from manually entered parameters
        % INPUT :
        %   - ...
            obj.LN  = LN;    %length of nose  
            obj.D  = D; %diameter at base of nose  
            obj.LB = LB; %length of Body tube 
            obj.DT  =  DT;  %diameter at rear of transition  
            obj.LT  = LT;%length of transition  
            obj.XT  =  XT;   %distance from tip of nose to front of transition (=tube lenght + nose as transition at the end) 
            obj.CR  =  CR; %fins root chord  
            obj.CT  =  CT;    %fins tip chord  
            obj.S  =  S;%fins semispan  
            obj.LF  =  LF;   %length of fin mid-chord line  
            obj.XRT  = XRT; %distance between fin root leading edge and fin tip leading edge parallel to body  
            obj.XF  =  XF;%distance from nose tip to fin root chord leading edge  
            obj.NF  = NF;%number of fins  
        end
        
        
    end
    
end