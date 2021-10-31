file_rcvr_in='.\data\rcvr.dat';
file_eph_in='.\data\eph.dat';
data_rcvr=load(file_rcvr_in);
data_eph=load(file_eph_in);

%************Caculate position for valid satellite************
num_sate_in_eph=size(data_eph,1);
num_sate_in_rcvr=size(data_rcvr,1);
pos_sate_with_dt_range=ones(num_sate_in_eph,6)*nan;
if num_sate_in_eph~=num_sate_in_rcvr
    error('Data is incomplete, missing data in ephemeris document or receiver data!');
end

for i=1:num_sate_in_eph
    eph_para=data_eph(i,:);
    sate_id=eph_para(2);
    k=find(data_rcvr(:,2)==sate_id);
    aa=data_rcvr(k,3);
    pos_sate_with_dt_range(i,:)=Pos_Sate_new(eph_para,data_rcvr(k,3));
end

%**********Calculate the user's position**************
%Define initial user position, residual threshold, epoch_litmit
XYZ0=[-2694685.473 -4293642.366 3857878.924];
res_threshold=10e-5;
epoch_limit=1000;

%Calculation
user_pose_result=user_pose(pos_sate_with_dt_range,XYZ0,res_threshold,epoch_limit);
