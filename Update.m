function [x,ab] = Update(x,G,T,Pbest,Gbest,fitness,Pfitness,Gfitness,c2,c1,k,ITERmax,ab,MinT,n) 
 
 
   r=rand;
   p= rand;
   dim = G*T;
   r1 = (fitness<Pfitness); 
   r2 = (fitness<Gfitness); 
   
   if r >= 0.5
        if p > 0.18
            Time(n,:) = rand(1, dim);
            sps = x ./ Time(n,:);
            Distance_S_Travel(n,:) = sps .* Time(n,:);
            Distance_Fox_Rat(n,:) = 0.5 .* Distance_S_Travel(n,:);
            tt = sum(Time(n,:)) / dim;
            t = tt / 2;
            Jump = 0.5 * 9.81 * t^2;
            Z = Jump .* c1/ITERmax;
        elseif p <= 0.18
            Time(n,:) = rand(1, dim);
            sps = x ./ Time(n,:);
            Distance_S_Travel(n,:) = sps .* Time(n,:);
            Distance_Fox_Rat(n,:) = 0.5 .* Distance_S_Travel(n,:);
            tt = sum(Time(n,:)) / dim;
            t = tt / 2;
            Jump = 0.5 * 9.81 * t^2;
            Z = Jump .* c2/ITERmax;
        end
        if MinT > tt
            MinT = tt;
        end 
    elseif r < 0.5
        % Random walk
        Z = c2-(c2-c1)* k/ITERmax;
    end



   DZ = Z*(r1*(Pbest-x)+r2*(Gbest-x));  
   
   
    for g = 1:G 
       for t = 1:T 
             ab(:,g+G*(t-1)) = [cos(DZ(g+G*(t-1))),-sin(DZ(g+G*(t-1)));sin(DZ(g+G*(t-1))),cos(DZ(g+G*(t-1)))]*ab(:,g+G*(t-1)); 
             x(g+G*(t-1)) = (rand()<ab(2,g+G*(t-1))^2); 
       end 
    end
    