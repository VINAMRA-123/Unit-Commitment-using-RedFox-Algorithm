 
function x = initial(G,T) 

x = ones(1,G*T);
ab = 1/sqrt(2)*ones(2,G*T); 
x = rand(1,G*T)<ab(2,:).^2;


