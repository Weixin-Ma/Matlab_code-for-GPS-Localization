function xyz=user_pose(pos_sate_with_dt_range,initial_pos,threshold_residual,intera_limit)
%********* The input of this function is a nx7 matri containing satellite position,
%********* relative pseudorange, and satellite clock correction, as well
%********* as initial user position, threshold for residual and interaction times. 
%********* The output of this function is the estimated user's position, the residual, 
%********* and the final iteraction times.

%GPS Constants:
c=299792458.0;    %"c" is the speed of light(m/s)


num_data=size(pos_sate_with_dt_range,1);
if num_data<4
    error("The number of availabe satellite is less than 4!")
else
    
    %Initializing the user position for the iteraction
    XYZ0=initial_pos;
    
    %Initializing the user reciver clock bias
    d_t=0;
    
    %Define the residual between the result of ith and (i+1)th iteraction
    res=1000;
    
    %Define the iteraction times
    inter_t=0;
    
    %Define the pseudorange received by the receiver
    Pr=ones(1,num_data);
    for j=1:num_data
            Pr(1,j)=pos_sate_with_dt_range(j,6)+pos_sate_with_dt_range(j,5)*c;
    end
    %Define the approx pseudorange during the iteraction
    Pr_a=ones(1,num_data);
    
    %Define accuracy thershold and epoch for iteraction 
    res_thershold=threshold_residual;
    epoch_limit=intera_limit;
    
    %Start iteraction
    while (res>res_thershold)  && (inter_t<epoch_limit)
        inter_t=inter_t+1;
        for j=1:num_data
             Pr_a(1,j)= sqrt((pos_sate_with_dt_range(j,2)-XYZ0(1))^2+(pos_sate_with_dt_range(j,3)-XYZ0(2))^2+(pos_sate_with_dt_range(j,4)-XYZ0(3))^2) + c*d_t;
             %Pr_a(1,j)= sqrt((sate_pose(j,1:3)-XYZ0)*(sate_pose(j,1:3)-XYZ0)') + c*d_t;
        end
        
        %Caculate the diference for the linearization equation
        %d_Pr=Pr-Pr_a;
        d_Pr=Pr_a-Pr;
        
        %Caculate the Matrix H
        H=ones(num_data,4);
        for n=1:num_data
           H(n,1)=(pos_sate_with_dt_range(n,2)-XYZ0(1))/Pr_a(1,n);  
           H(n,2)=(pos_sate_with_dt_range(n,3)-XYZ0(2))/Pr_a(1,n); 
           H(n,3)=(pos_sate_with_dt_range(n,4)-XYZ0(3))/Pr_a(1,n); 
        end
        
        %Update the user position XYZ0, d_t, and residual
        d_X=inv(H'*H)*H'*d_Pr';
        %d_X=((H'*H)\H')*d_Pr';
        %d_X=H\d_Pr';
        XYZ0=XYZ0+d_X(1:3,1)'; 
        d_t=d_X(4,1)/(-c);
        res=max(abs(d_X(1:3,1)));
    end
end
xyz=[XYZ0(1); XYZ0(2); XYZ0(3);res;inter_t];
return