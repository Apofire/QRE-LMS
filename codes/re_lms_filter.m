function [y,error,wts,gamma,store_mu] = re_lms_filter(x,d,w,Mu,Mu_nlms,a,ITER)

   % initialisations
   N     = length(x);     % input signal length
   W     = length(w);     % convolution window size
   x_w   = zeros(1,W);    % convolution buffer
   y     = zeros(length(Mu),N);    % filtered output
   e     = zeros(1,N);    % error signal
   error = zeros(1,ITER); % rms error after every iteration
   wts   = zeros(ITER,W); % matrix of weights at every iteration
   
   % modified lms algorithm    
   for iter = 1:ITER
       wts(iter,:) = w;
%         if iter == 1
%            mu = Mu;
%         else
%             f = 0;
%             g = 0;
% %                 for j = 1:length(e)
% %                     for i = 1:length(e)
% %                     f = f + abs(e(i)*y(j)/(a+y1(j)))/((length(e))^2);
% %                     g = g + (e(i)*y(j)/(a+y1(j)))^2/((length(e))^2);
% %                     end
% %                 end
%                   for i = 1:length(e)
%                      f = f + abs(e(i)*y(i)/(a+y1(i)))/(length(e));
%                      g = g + (e(i)*y(i)/(a+y1(i)))^2/(length(e));
%                   end
%            gamma(iter) = f/g; % store the output of the bound
%            mu = mu + mu*1.5*1/(exp(iter^2))*(f/g); % store the step size in every iteration
%            store_mu(iter) = mu;
%         end
    if iter < N

        for i = 1:N
            x_w  = [x(i),x_w(1:end-1)]; % define the input signal window for convolution
            y(i) = w*x_w';              % find y(i) by performing the convolution
            e(i) = d(i) - y(i);         % define the error signal
            y1(i) = (w.^2)*(x_w.^2)';
            %            w = w + (mu*e(i)/(a + y1(i)).* ((w*w')*x_w));
            w = w + (Mu*e(i)/(a + y1(i)).* ((w*w')*x_w));
        end
    else
        for i = 1:N
           x_w  = [x(i),x_w(1:end-1)]; % define the input signal window for convolution
           y(k,i) = w*x_w';              % find y(i) by performing the convolution
           e(i) = d(i) - y(k,i);         % define the error signal
           w = w + (Mu_nlms/(a + norm(x_w,2)^2))*e(i).*x_w;  % update the filter weights           
       end

       % define the error as the normed error      
       error(iter) = norm(e,2);   
   end
end
