function TC_clean = cleanThrustCurve(TC_orig)
    % TC_clean = cleanThrustCurve(Tc_orig), cleans original thrust curve
    % TC_orig from any NaN values and time steps or thrust values that are
    % smaller than 0.
    % INPUT : 
    %   - TC_orig   :   original thrust curve
    % OUTPUT :
    %   - TC_clean  :   cleaned thrust curve
    
    % positive or nul timesteps and thrust values
    TC_clean = TC_orig(TC_orig>=0);
    
    % no NaN values
    TC_clean = TC_clean(~isnan(TC_clean));
    
    TC_clean = reshape(TC_clean, [], 2);
    
    % check if time grows monotonly
    if(~issorted(TC_clean(:,1)))
       warning('Thrust curve time steps are not monotonicaly ascending.'); 
    end
end

