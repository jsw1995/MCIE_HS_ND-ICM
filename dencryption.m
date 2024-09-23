 function [ rim,rtcip11 ] = dencryption( cip,cover,dnkey,mmr,TT,cr,emb_type,im_shape,rC )
%   DENCRYPTION 解密
%   此处显示详细说明

m=im_shape(1); n=im_shape(2);k=im_shape(3);

[m1,n1,k1]=size(cover);

a1=dnkey(1);b1=dnkey(2);x01=dnkey(3);y01=dnkey(4);
a2=dnkey(5);b2=dnkey(6);x02=dnkey(7);y02=dnkey(8);z02=dnkey(9);w02=dnkey(10);
a3=dnkey(11);b3=dnkey(12);x03=dnkey(13);y03=dnkey(14);z03=dnkey(15);w03=dnkey(16);

if emb_type == 1
    a=4;b=4;c=4;d=4;
elseif emb_type == 2
%     a=3;b=5;c=3;d=6;
    a=4;b=4;c=4;d=4;
else
    print('只有两种类型分位方式')
end

% 提取
[ x,y,z,w ] = IICM_4d( [x03,y03,z03,w03],[a3,b3],(m1+k1*n1)*0.5 );
if emb_type == 1
   rtcip = extract( cip,cover,a,b,c,d,x,y,z,w,[m/2,n/2,k],rC );
elseif emb_type == 2
   rtcip = extract2( cip,cover,a,b,c,d,x,y,z,w,[m/2,n/2,k],rC );
end

rtcip11=rtcip;

% 反向置乱
if k==3
   rtcip = [rtcip(:,:,1),rtcip(:,:,2),rtcip(:,:,3)];
end

ca = floor(rtcip/(b*c*d));
ch = floor(mod(rtcip,(b*c*d))/(c*d));
cv = floor(mod(rtcip,(c*d))/d);
cd = mod(rtcip,d);
[ x,y,z,w ] = IICM_4d( [x02,y02,z02,w02],[a2,b2], (m+k*n)*cr );
ca1 = dscram( ca,x,2 );
ch1 = dscram( ch,y,2 );
cv1 = dscram( cv,z,2 );
cd1 = dscram( cd,w,2 );
rx3 = ca1*(b*c*d) + ch1*(c*d) + cv1*d + cd1;

% 直方图移位
rx3 = mod(rx3+TT,256);

% 重构
[x,y]= IICM_2d( [x01,y01],[a1,b1], n*m*k*cr );
if k==1
    max_x3=mmr(1);min_x3=mmr(2);
    rim = refact( rx3,x,y,cr,max_x3, min_x3 );
else
    max_x3r=mmr(1);min_x3r=mmr(2);
    rimr = refact( rx3(:,1:n/2),x(1:n*m*cr),y(1:n*m*cr),cr,max_x3r, min_x3r );
    max_x3g=mmr(3);min_x3g=mmr(4);
    rimg = refact( rx3(:,n/2+1:2*n/2),x(n*m*cr+1:2*n*m*cr),y(n*m*cr+1:2*n*m*cr),cr,max_x3g, min_x3g );
    max_x3b=mmr(5);min_x3b=mmr(6);
    rimb = refact( rx3(:,2*n/2+1:3*n/2),x(2*n*m*cr+1:3*n*m*cr),y(2*n*m*cr+1:3*n*m*cr),cr,max_x3b, min_x3b );
    rim = cat(3,rimr,rimg,rimb);
end

 end

 
function [ rim ] = refact( image,x,y,cr,max_x3, min_x3 )
%   REFACT 此处显示有关此函数的摘要
%   此处显示详细说明

x = 1-2*mod(x * 10^4, 1);
y = 1-2*mod(y * 10^4, 1);

X3 = image / 255 * (max_x3 - min_x3) + min_x3;
% X3 = log((image/255 -1 )) / (-(max_x3-min_x3)) + ((max_x3+min_x3)/2);

[m,n] = size(image);
rows = int16(m/cr);
columns = int16(n/cr);

R1 = reshape(x, [int16(rows * cr), rows]);
R2 = reshape(y, [int16(columns * cr), columns]);

% 优化测量矩阵
R1(:, 1:int16(cr * rows)) = R1(:, 1:int16(cr * rows)) * 5000;
R1 = orth(R1')';
R2(:, 1:int16(cr * columns)) = R2(:, 1:int16(cr * columns)) * 5000;
R2 = orth(R2');

% 重构
[ rec ] = nsl0_2d(X3, R1, R2);

rim = idct2(rec);
rim = double(uint8(rim));

 end

function [ rim ] = nsl0_2d( y,A,B )
%   NSL0_2D 此处显示有关此函数的摘要
%   此处显示详细说明

sigma_min = 0.01;
sigma_decrease_factor = 0.05;  %
ksai = 0.01;

A_pinv = pinv(A);
B_pinv = pinv(B);

s = A_pinv * y * B_pinv;

sigma = 4 * max(max(abs(s)));
r = 0;
r0 = y - A * s * B;

while (sigma>sigma_min)

    if sum(sum((r-r0).^2)) < ksai
        
        d = -(sigma^2 * s) ./ (s.*s + sigma^2);
        s = s + d;
        s = s - A_pinv * (A * s * B - y) * B_pinv;
        r0 = y - A * s * B;
        
    end

    sigma = sigma * sigma_decrease_factor;
 
end

rim = s;

end

function [ rim ] = sl0_2d( y,A,B )
%   SL0_2D 此处显示有关此函数的摘要
%   此处显示详细说明

sigma_min = 0.001;  % 0.1
sigma_decrease_factor = 0.02;  % 0.02
L = 4;

A_pinv = pinv(A);
B_pinv = pinv(B);

s = A_pinv * y * B_pinv;

sigma = 2 * max(max(abs(s)));

while sigma>sigma_min

    for i =1:L
        delta = -(sigma^2 * s)./(sigma^2 + s.^2);
        s = s + delta;
        s = s - A_pinv * (A * s * B - y) * B_pinv;
    end
    
    sigma = sigma*sigma_decrease_factor;
    
end

rim = s;

end

function [ x,y,z,w ] = IICM_4d( init,para,L )
%   IICM 
%   此处显示详细说明

x = []; y= []; z = []; w= [];
x(1)=init(1);y(1)=init(2);z(1)=init(3);w(1)=init(4);
a=para(1);b=para(2);

for i=1:L+1000
    x(i+1)=sqrt(1-b*y(i)^2)*sin(a/(x(i)^2));
    y(i+1)=sqrt(1-b*z(i)^2)*sin(a/(y(i)^2));
    z(i+1)=sqrt(1-b*w(i)^2)*sin(a/(z(i)^2));
    w(i+1)=sqrt(1-b*x(i)^2)*sin(a/(w(i)^2));
end

x=x(1001:L+1000);
y=y(1001:L+1000);
z=z(1001:L+1000);
w=w(1001:L+1000);

end

function [ x,y ] = IICM_2d( init,para,L )
%   IICM 
%   此处显示详细说明

x = []; y= [];
x(1)=init(1);y(1)=init(2);
a=para(1);b=para(2);

for i=1:L+1000
    x(i+1)=sqrt(1-b*y(i)^2)*sin(a/(x(i)^2));
    y(i+1)=sqrt(1-b*x(i)^2)*sin(a/(y(i)^2));
end

x=x(1001:L+1000);
y=y(1001:L+1000);

end

function [ rim ] = extract( cip,cover,a,b,c,d,x,y,z,w,im_shape,rC )
%   EXTRACT 自己的嵌入操作(空域)
%   此处显示详细说明

m1=im_shape(1); n1=im_shape(2);k1=im_shape(3);
[m,n,k] = size(cip);
if k==1
    CA3 = cip(1:m/2,1:n/2);
    CH3 = cip(m/2+1:m,1:n/2);
    CV3 = cip(1:m/2,n/2+1:n);
    CD3 = cip(m/2+1:m,n/2+1:n);
    
    CA31 = cover(1:m/2,1:n/2);
    CH31 = cover(m/2+1:m,1:n/2);
    CV31 = cover(1:m/2,n/2+1:n);
    CD31 = cover(m/2+1:m,n/2+1:n);
else 
    CA3 = [cip(1:m/2,1:n/2,1),cip(1:m/2,1:n/2,2),cip(1:m/2,1:n/2,3)];
    CH3 = [cip(m/2+1:m,1:n/2,1),cip(m/2+1:m,1:n/2,2),cip(m/2+1:m,1:n/2,3)];
    CV3 = [cip(1:m/2,n/2+1:n,1),cip(1:m/2,n/2+1:n,2),cip(1:m/2,n/2+1:n,3)];
    CD3 = [cip(m/2+1:m,n/2+1:n,1),cip(m/2+1:m,n/2+1:n,2),cip(m/2+1:m,n/2+1:n,3)];
    
    CA31 = [cover(1:m/2,1:n/2,1),cover(1:m/2,1:n/2,2),cover(1:m/2,1:n/2,3)];
    CH31 = [cover(m/2+1:m,1:n/2,1),cover(m/2+1:m,1:n/2,2),cover(m/2+1:m,1:n/2,3)];
    CV31 = [cover(1:m/2,n/2+1:n,1),cover(1:m/2,n/2+1:n,2),cover(1:m/2,n/2+1:n,3)];
    CD31 = [cover(m/2+1:m,n/2+1:n,1),cover(m/2+1:m,n/2+1:n,2),cover(m/2+1:m,n/2+1:n,3)];
end

% 置乱
CA3 = scram( CA3,x,2 );
CH3 = scram( CH3,y,2 );
CV3 = scram( CV3,z,2 );
CD3 = scram( CD3,w,2 );

CA31 = scram( CA31,x,2 );
CH31 = scram( CH31,y,2 );
CV31 = scram( CV31,z,2 );
CD31 = scram( CD31,w,2 );

% 提取
CA3 = reshape(CA3,[1,m/2*n/2*k]);
CH3 = reshape(CH3,[1,m/2*n/2*k]);
CV3 = reshape(CV3,[1,m/2*n/2*k]);
CD3 = reshape(CD3,[1,m/2*n/2*k]);

CA31 = reshape(CA31,[1,m/2*n/2*k]);
CH31 = reshape(CH31,[1,m/2*n/2*k]);
CV31 = reshape(CV31,[1,m/2*n/2*k]);
CD31 = reshape(CD31,[1,m/2*n/2*k]);

CA4 = mod(mod(CA3(1:m1*n1*k1),a) - mod(CA31(1:m1*n1*k1),a), a);
CH4 = mod(mod(CH3(1:m1*n1*k1),b) - mod(CH31(1:m1*n1*k1),b), b);
CV4 = mod(mod(CV3(1:m1*n1*k1),c) - mod(CV31(1:m1*n1*k1),c), c);
CD4 = mod(mod(CD3(1:m1*n1*k1),d) - mod(CD31(1:m1*n1*k1),d), d);

% CA42 = CA4;
% CA4(CA42==2) = rC(1)-1;
% CA4(CA42==3) = rC(2)-1;
% CA4(CA42==1) = rC(3)-1;
% CA4(CA42==0) = rC(4)-1;
% 
% CH42=CH4;
% CH4(CH42==2) = rC(5)-1;
% CH4(CH42==3) = rC(6)-1;
% CH4(CH42==1) = rC(7)-1;
% CH4(CH42==0) = rC(8)-1;
% 
% CV42=CV4;
% CV4(CV42==2) = rC(9)-1;
% CV4(CV42==3) = rC(10)-1;
% CV4(CV42==1) = rC(11)-1;
% CV4(CV42==0) = rC(12)-1;
% 
% CD42=CD4;
% CD4(CD42==2) = rC(13)-1;
% CD4(CD42==3) = rC(14)-1;
% CD4(CD42==1) = rC(15)-1;
% CD4(CD42==0) = rC(16)-1;


rim = CA4*(b*c*d) + CH4*(c*d) + CV4*d + CD4;
rim = reshape(rim,[int16(m1),int16(n1),k1]);

end

function [ rim ] = extract2( cip,cover,a,b,c,d,x,y,z,w,im_shape,rC )
%   EXTRACT 自己的嵌入操作(小波变换域)
%   此处显示详细说明

m1=im_shape(1); n1=im_shape(2);k1=im_shape(3);
[m,n,k] = size(cip);

% 修改封面 (整数小波才需要修改)
max1=252;min1=3;
cover(cover<min1)=min1;
cover(cover>max1)=max1;

% 提升的小波变换
LS=liftwave('haar','Int2Int');
cip=double(cip);
cover=double(cover);
if k==1
    [CA3,CH3,CV3,CD3]=lwt2(cip,LS);
    [CA31,CH31,CV31,CD31]=lwt2(cover,LS);
else
    [CA3r,CH3r,CV3r,CD3r]=lwt2(cip(:,:,1),LS);
    [CA3g,CH3g,CV3g,CD3g]=lwt2(cip(:,:,2),LS);
    [CA3b,CH3b,CV3b,CD3b]=lwt2(cip(:,:,3),LS);
    CA3 = cat(3,CA3r,CA3g,CA3b);
    CH3 = cat(3,CH3r,CH3g,CH3b);
    CV3 = cat(3,CV3r,CV3g,CV3b);
    CD3 = cat(3,CD3r,CD3g,CD3b);
    
    [CA31r,CH31r,CV31r,CD31r]=lwt2(cover(:,:,1),LS);
    [CA31g,CH31g,CV31g,CD31g]=lwt2(cover(:,:,2),LS);
    [CA31b,CH31b,CV31b,CD31b]=lwt2(cover(:,:,3),LS);
    CA31 = cat(3,CA31r,CA31g,CA31b);
    CH31 = cat(3,CH31r,CH31g,CH31b);
    CV31 = cat(3,CV31r,CV31g,CV31b);
    CD31 = cat(3,CD31r,CD31g,CD31b);
end
% 置乱
CA3 = scram( CA3,x,2 );
CH3 = scram( CH3,y,2 );
CV3 = scram( CV3,z,2 );
CD3 = scram( CD3,w,2 );

CA31 = scram( CA31,x,2 );
CH31 = scram( CH31,y,2 );
CV31 = scram( CV31,z,2 );
CD31 = scram( CD31,w,2 );

CA4 = mod(mod(CA3(1:m1*n1*k1),a) - mod(CA31(1:m1*n1*k1),a), a);
CH4 = mod(mod(CH3(1:m1*n1*k1),b) - mod(CH31(1:m1*n1*k1),b), b);
CV4 = mod(mod(CV3(1:m1*n1*k1),c) - mod(CV31(1:m1*n1*k1),c), c);
CD4 = mod(mod(CD3(1:m1*n1*k1),d) - mod(CD31(1:m1*n1*k1),d), d);

% CA42 = CA4;
% CA4(CA42==2) = rC(1)-1;
% CA4(CA42==3) = rC(2)-1;
% CA4(CA42==1) = rC(3)-1;
% CA4(CA42==0) = rC(4)-1;
% 
% CH42=CH4;
% CH4(CH42==2) = rC(5)-1;
% CH4(CH42==3) = rC(6)-1;
% CH4(CH42==1) = rC(7)-1;
% CH4(CH42==0) = rC(8)-1;
% 
% 
% CV42=CV4;
% CV4(CV42==2) = rC(9)-1;
% CV4(CV42==3) = rC(10)-1;
% CV4(CV42==1) = rC(11)-1;
% CV4(CV42==0) = rC(12)-1;
% 
% 
% CD42=CD4;
% CD4(CD42==2) = rC(13)-1;
% CD4(CD42==3) = rC(14)-1;
% CD4(CD42==1) = rC(15)-1;
% CD4(CD42==0) = rC(16)-1;


rim = CA4*(b*c*d) + CH4*(c*d) + CV4*d + CD4;
rim = reshape(rim,[int16(m1),int16(n1),k1]);

end

