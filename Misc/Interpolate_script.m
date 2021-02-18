% Interpolate

% 100Hz
low=xx;
high=xx;
        % Reu The
y=interp1([X(low);X(high)],[Y(low,:);Y(high,:)],X(low:high),'cubic');
Y(low:high,:)=y;