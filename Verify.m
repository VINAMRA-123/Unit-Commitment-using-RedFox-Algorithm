function [x,TOFF] = Verify(x,G,T,x0,MUT,MDT,InitialTON,InitialTOFF,PMAX,PD,SR,h) 

TON = zeros(G,T); 
TOFF = zeros(G,T); 
 
for t = 1:T 

    if t == 1 
        TON(:,t) = (InitialTON.*x(G*(t-1)+1:G*t)+x(G*(t-1)+1:G*t))'; 
        TOFF(:,t) = InitialTOFF'.*~TON(:,t)+~TON(:,t); 
    else 
        flag = (TOFF(:,t-1) == 0&TON(:,t-1) < MUT')|(TON(:,t-1) == 0&TOFF(:,t-1) < MDT'); 
        x(G*(t-1)+1:G*t) = x(G*(t-2)+1:G*(t-1)).*flag'+x(G*(t-1)+1:G*t).*~flag'; 
        TON(:,t) = TON(:,t-1).*x(G*(t-1)+1:G*t)'+x(G*(t-1)+1:G*t)'; 
        TOFF(:,t) = TOFF(:,t-1).*~TON(:,t)+~TON(:,t); 
    end 

    P = PMAX*x(G*(t-1)+1:G*t)'; 
    if P < PD(t)+SR(t) 
        for j = 1:G 
            g = h(j); 
            if x(g+G*(t-1)) == 1 
                continue; 
            else 
                x(g+G*(t-1)) = 1; 
                if TOFF(g,t) > MDT(g) 
                        if t == 1 
                            TON(g,t) = InitialTON(g)+1; 
                       else 
                            TON(g,t) = TON(g,t-1)+1; 
                        end 
                      TOFF(g,t) = 0; 
                else 
                        l = t-TOFF(g,t)+1; 
                        if l <= 0 
                            l = 1; 
                        end 
                         
                        x(g+G*(l-1:t-1)) = 1;    
                        if l == 1 
                            TON(g,l) = InitialTON(g)+1; 
                            TON(g,l+1:t) = TON(g,l)+[1:t-l]; 
                        else 
                            TON(g,l:t) = TON(g,l-1)+[1:t-l+1]; 
                        end 
                        TOFF(g,l:t) = 0; 
               end 
               P = PMAX*x(G*(t-1)+1:G*t)'; 
                if P >= (PD(t)+SR(t)) 
                   break; 
                end 
            end 
        end 
    end
    
 for i = 1:G 
     g = h(G+1-i); 
     P1 = PMAX*x(G*(t-1)+1:G*t)'; 
     if x(g+G*(t-1)) == 1 
        if P1-PMAX(g) >= PD(t)+SR(t) 
           if TON(g,t) > MUT(g)||(TON(g,t) == 1) 
               x(g+G*(t-1))=0; 
                TON(g,t)=0; 
                if t==1 
                    TOFF(g,t)=InitialTOFF(g)+1; 
                else 
                    TOFF(g,t)=TOFF(g,t-1)+1; 
                end 
           else 
               continue; 
           end 
        else  
             break; 
        end 
    end 
 end 
end 