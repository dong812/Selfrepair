S = zeros(1,10000);
SIC = zeros(1,10000);
tca = zeros(1,10000);
tca(1,1)=1;
dt  = 0.001;



for i=2:10000
    ds = -1*S(1,i-1)/0.1 + 20*tca(1,i-1);
    ds = ds*dt;
    S(1,i)=S(1,i-1) + ds;
    dsic = (-1*SIC(1,i-1)+20*S(1,i))*dt/0.0375 ;
    SIC(1,i) = SIC(1,i-1) + dsic;
end