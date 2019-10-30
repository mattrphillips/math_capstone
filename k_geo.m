function k = k_geo(x,y,x_steps,y_steps,k1,R)

% Constructs a room with R-value of R (converted to k2) within air of k=k1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%   OWL/OWT OWT     OWT     OWT     OWT     OWT/OWR             %
%   OWL     IWL/IWT IWT     IWT     IWT/IWR     OWR             %
%   OWL     IWL                         IWR     OWR             %
%   OWL     IWL                         IWR     OWR             %
%   OWL     IWL                         IWR     OWR             %
%   OWL     IWL/IWB IWB     IWB     IWB/IWR     OWR             %
%   OWL/OWB OWB     OWB     OWB     OWB     OWB/OWR             %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4m x 4m room w/ 1m thick walls
OWL = -2;
OWR = 2;
OWT = 2;
OWB = -2;

IWL = -1.9;
IWR = 1.9;
IWT = 1.9;
IWB = -1.9;

thickness = abs(OWL - IWL);
k2 = thickness / R;

k = zeros(y_steps,x_steps);
for i = 1:x_steps
    for j = 1:y_steps
        if x(i) >= OWL && x(i) <= OWR
            if y(j) >= OWB && y(j) <= OWT
                if x(i) >= IWL && x(i) <= IWR
                    if y(j) >= IWB && y(j) <= IWT
                        k(j,i) = k1;
                    else
                        k(j,i) = k2;
                    end
                else
                    k(j,i) = k2;
                end
            else
                k(j,i) = k1;
            end
        else
            k(j,i)=k1;
        end
    end
end

end