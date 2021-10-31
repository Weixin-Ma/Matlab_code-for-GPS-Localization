%************Caculate position for valid satellite************
function XYZ=Sate_pose(eph_para,pesudo_range)
%%%%%%%%%This function taks ephemeris data, pesudo range between satellite
%%%%%%%%%and user, and satellite id as input; and output the 

%GPS Constants:
mu=3.986005e+14;           %"mu" is Earth's universal gravitation constant
%dt_r=(3.8e-5)/86400;      %"dt_r" is Time correction due to relativistic effects
rot_rate=7.2921151467e-5; %"rot_rate" is value of earth's rotation rate 
c=299792458.0;    %"c" is the speed of light(m/s)
if (isempty(eph_para)==1)||(numel(find(isnan(eph_para)))>=1)
    error('The input ephemeris data is empty or missing parameters!');
else
    rcvr_t=eph_para(1); %"rcvr_t" is the receiver time;
    svid=eph_para(2);   %"svid" is the satellite PRN number; 
    toc=eph_para(3);    %"toc" is the reference time of clock parameters;
    toe=eph_para(4);    %"toe" is the reference time of ephemeris parameters;
    af0=eph_para(5);    %"af0" is clock bias;
    af1=eph_para(6);    %"af1" is clock drift;
    af2=eph_para(7);    %"af2" is frequency drift;
    e=eph_para(9);      %"e" is the eccentricity;
    sqrta=eph_para(10); %"sqrta" is square root of semi-major axis of the orbit plane;
    dn=eph_para(11);    %"dn" is the mean motion correction;
    m0=eph_para(12);    %"m0" is the mean anomaly at reference time;
    w=eph_para(13);     %"w" is the argument of perigee;
    omg0=eph_para(14);  %"omg0" is the right ascension (r);
    i0=eph_para(15);    %"i0" is inclination angle at reference time (r);
    odot=eph_para(16);  %"odot" is the rate of right ascension (r/s)
    idot=eph_para(17);  %"idot" is rate of inclination angle (r/s)
    cus=eph_para(18);   %"cus" is argument of latitude correction, sine;
    cuc=eph_para(19);   %"cuc" is argument of latitude correction, cosine;
    cis=eph_para(20);   %"cis" is inclination correction, sine
    cic=eph_para(21);   %"cic" is inclination correction, cosine
    crs=eph_para(22);   %"crs" is radius correction, cosine;
    crc=eph_para(23);   %"crc" is radius correction, cosine;
    
    %Caculate mean motion
    n=sqrt(mu)/(sqrta^3)+dn;
    
    %Calculate the transmit time and clock correction for each satellite
    ts=rcvr_t-pesudo_range/c;
    d_t=af0+af1*(ts-toc)+af2*(ts-toc)^2;
    %t_k=t_corrected-toe;
    t_k=ts-toe;
    
    %Mean anomaly
    M=m0+n*t_k;
    if M<0
        M=M+2*pi;
    end
    if M>2*pi
        M=M-2*pi;
    end
    
    %Caculate eccentric anomaly E
    E_old = M;
    E_new = M + e*sin(E_old);
    i = 1;
    while abs(E_new - E_old)>1e-8
        %print("i={},E={}".format(i, E_new))
        E_old = E_new;
        E_new = M + e*sin(E_old);
        i = i+1;
        if (i>10)
            break
        end
    end
    E= E_new;

    %Caculate true anomaly v
    v=atan2(sqrt(1-e^2)*sin(E),(cos(E)-e));
    
    %Caculate the correcntion for the argument of latitude
    AoL=v+w; % argument of latitude before correction
    d_AoL=cuc*cos(2*AoL)+cus*sin(2*AoL); % correction item for argument of latitude
    AoL_correct=AoL+d_AoL;  %the corrected argument of latitude
    
    %Caculate the correcntion for inclination
    d_incl=cis*sin(2*AoL)+cic*cos(2*AoL);
    incl_correct=i0+d_incl+idot*t_k;
    
    %Caculate the correction for radius
    d_r=crs*sin(2*AoL)+crc*cos(2*AoL);
    r_correct=(sqrta.^2)*(1-e*cos(E))+d_r;
    
    %Caculate the latitude for the ascension at the emission time 
    L=omg0+(odot-rot_rate)*t_k-rot_rate*toe;
    
    %Caculate the statellite coordinate in Orbital plane
    x=r_correct*cos(AoL_correct);
    y=r_correct*sin(AoL_correct);
    
    %Caculate the statellite coordinate in ECEF
    X=x*cos(L)-y*cos(incl_correct)*sin(L);
    Y=x*sin(L)+y*cos(incl_correct)*cos(L);
    Z=y*sin(incl_correct);
    XYZ=[svid;X;Y;Z;d_t;pesudo_range];
    
end
return;


