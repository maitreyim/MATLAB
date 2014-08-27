function [F,J] = impliedVolatility_BAW_fsolve(x, isCall, marketPrice, S, X, r, b, T)

    try 
        
        % Set initial starting values
        S_Star = x(1);
        Sigma = x(2);
        

        if (isCall)

            % 1.) Set up variables needed to calculate f and g
            M = (2*r)/(Sigma^2);
            N = (2*b)/(Sigma^2);
            K = 1 - exp(-r*T);
            q2 = 0.5 * ( -(N-1) + sqrt( (N-1)^2 + (4*M)/K ) );

            d1_star = (log(S_Star/X) + (b+0.5*Sigma^2)*T) / (Sigma * sqrt(T));
            d2_star = d1_star - Sigma * sqrt(T);
            N_d1_Star = normcdf(d1_star,0,1);
            n_d1_Star = normpdf(d1_star,0,1);
            n_d2_Star = normpdf(d2_star,0,1);

            d1 = (log(S/X) + (b+0.5*Sigma^2)*T) / (Sigma * sqrt(T));
            d2 = d1 - Sigma * sqrt(T);
            n_d2 = normpdf(d2,0,1);

            A2 = (S_Star/q2) * ( 1 - exp((b-r)*T) * N_d1_Star );

            c_S_Star = EurpeanOptionOnCommodity(1, S_Star, X, r, b, T, Sigma);
            c = EurpeanOptionOnCommodity(1, S, X, r, b, T, Sigma);

            
            % 2.) Calculate f(x) = 0 and g(x) = 0
            f = c_S_Star + S_Star/q2 * ( 1 - exp((b-r)*T) * N_d1_Star ) - (S_Star - X);
            g = c + A2 * ((S/S_Star)^q2) - marketPrice;

            F = [f g];
            
            
            % 3.) Calculate Jacobian Matrix of f and g
                
            % Substitution variables
            y1 = 1 - exp((b-r)*T) * N_d1_Star;
            y2 = n_d1_Star / (S_Star * Sigma * sqrt(T));
            y3 = exp((b-r)*T);
            y4 = b/(Sigma^4) - ( ( (N-1)*b + 2*r/K ) / ( Sigma^4 * sqrt((N-1)^2 + 4*M/K) ) );
            y5 = X * sqrt(T) * exp(-r*T) * n_d2 / (2*Sigma);
            y6 = -0.5 * ( log(S_Star/X) + b*T ) / ( (Sigma^3 * sqrt(T)) / (4*Sigma) );
            y7 = X * sqrt(T) * exp(-r*T) * n_d2_Star / (2*Sigma);

            % Partial derivatives of 2x2Jacobi matrix
            dg_dS_Star = ((S/S_Star)^q2) * ( (1/q2)*y1 - (S_Star/q2)*y2*y3 - q2*A2*(1/S_Star) );
            df_dS_Star = y3 * N_d1_Star + (y1/q2) - (S_Star/q2)*(y2*y3) - 1;
            dg_dSigmaSqr = y5 + ((S/S_Star)^q2) * ( A2*log(S/S_Star)*y4 - (S_Star/q2)*y3*y6*n_d1_Star - y1*y4*(S_Star/(q2^2)) );
            df_dSigmaSqr = y7 - y1*y4*(S_Star/(q2^2)) - (S_Star/q2)*y3*y6*n_d1_Star;

            J = [dg_dS_Star dg_dSigmaSqr ; df_dS_Star df_dSigmaSqr];

        else

            % 1.) Set up variables needed to calculate f and g
            M = (2*r)/(Sigma^2);
            N = (2*b)/(Sigma^2);
            K = 1 - exp(-r*T);
            q1 = ( -(N-1) - sqrt( (N-1)^2 + (4*M)/K ) ) / 2;

            d1_star = (log(S_Star/X) + (b+0.5*Sigma^2)*T) / (Sigma * sqrt(T));
            d2_star = d1_star - Sigma * sqrt(T);
            N_d1_Star = normcdf(-d1_star,0,1);
            n_d1_Star = normpdf(-d1_star,0,1);
            n_d2_Star = normpdf(-d2_star,0,1);

            d1 = (log(S/X) + (b+0.5*Sigma^2)*T) / (Sigma * sqrt(T));
            d2 = d1 - Sigma * sqrt(T);
            n_d2 = normpdf(-d2,0,1);

            A1 = -(S_Star/q1) * ( 1 - exp((b-r)*T) * N_d1_Star ); 

            p_S_Star = EurpeanOptionOnCommodity(0, S_Star, X, r, b, T, Sigma);
            p = EurpeanOptionOnCommodity(0, S, X, r, b, T, Sigma);

            
            % 2.) Calculate f(x) = 0 and g(x) = 0
            f = p_S_Star - S_Star/q1 * ( 1 - exp((b-r)*T) * N_d1_Star ) - (X - S_Star);
            g = p + A1*((S/S_Star)^q1) - marketPrice;

            F = [f ; g];
            
            
            % 3.) Calculate Jacobian Matrix of f and g
            % Substitution variables
            z1 = 1 - exp((b-r)*T) * N_d1_Star;
            z2 = -n_d1_Star * exp((b-r)*T) / (S_Star * Sigma * sqrt(T));
            z3 = ( 0.5*(log(S_Star/X) + b*T)/(Sigma^3*sqrt(T)) ) - (sqrt(T)/(4*Sigma));
            z4 = ( b + ( ((N-1)*b + 2*r/K) / sqrt((N-1)^2 + 4*M/K) ) ) / (Sigma^4);
            z5 = (S_Star/q1) * ( exp((b-r)*T) * n_d1_Star * z3 ) + (S_Star/(q1^2))*z1*z4;

            % Partial derivatives of 2x2Jacobi matrix
            dg_dS_Star = ((S/S_Star)^q1) * ( (S_Star/q1)*z2 - (1/q1)*z1 - q1*A1*(1/S_Star) );
            df_dS_Star = z1 * (1 - 1/q1) + (S_Star/q1)*z2;
            dg_dSigmaSqr = ( n_d2*X*sqrt(T)*exp(-r*T) / (2*Sigma) ) +  ((S/S_Star)^q1) * ( A1*log(S/S_Star)*z4+z5 );
            df_dSigmaSqr = ( n_d2_Star*X*sqrt(T)*exp(-r*T) / (2*Sigma) ) + (S_Star/q1)*exp((b-r)*T)*n_d1_Star*z3 + (S_Star/(q1^2))*z1*z4;

            J = [dg_dS_Star dg_dSigmaSqr ; df_dS_Star df_dSigmaSqr];
            
        end
    catch
        F = [NaN NaN];
        J = [NaN NaN ; NaN NaN];
    end

end