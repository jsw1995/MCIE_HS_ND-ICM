function [ cip,dnkey,mmr,TT,tcip,rC ] = encryption( im,cover,key,cr,emb_type,tt )
%   ENCRYPTION 加密
%   此处显示详细说明

[m,n,k]=size(im);
[m1,n1,k1]=size(cover);

% 分位置乱
if emb_type == 1
    a=4;b=4;c=4;d=4;
elseif emb_type == 2
%     a=3;b=5;c=3;d=6;
    a=4;b=4;c=4;d=4;
else
    print('只有两种类型分位方式')
end

% 动态密钥生成
[ dnkey ] = dkey( im,key );
a1=dnkey(1);b1=dnkey(2);x01=dnkey(3);y01=dnkey(4);
a2=dnkey(5);b2=dnkey(6);x02=dnkey(7);y02=dnkey(8);z02=dnkey(9);w02=dnkey(10);
a3=dnkey(11);b3=dnkey(12);x03=dnkey(13);y03=dnkey(14);z03=dnkey(15);w03=dnkey(16);
% 压缩
[x,y]= IICM_2d( [x01,y01],[a1,b1], n*m*k*cr );
% figure(21)
% imhist((x+1)/2)
% figure(22)
% imhist((y+1)/2)
if k==1
    [ X3,max_x3,min_x3 ] = compre( im,x,y,cr );
    mmr = [max_x3,min_x3];
else
    [ X3r,max_x3r,min_x3r ] = compre( im(:,:,1),x(1:n*m*cr),y(1:n*m*cr),cr );
    [ X3g,max_x3g,min_x3g ] = compre( im(:,:,2),x(n*m*cr+1:2*n*m*cr),y(n*m*cr+1:2*n*m*cr),cr );
    [ X3b,max_x3b,min_x3b ] = compre( im(:,:,3),x(2*n*m*cr+1:3*n*m*cr),y(2*n*m*cr+1:3*n*m*cr),cr );
    X3 = [X3r,X3g,X3b];
    mmr = [max_x3r,min_x3r,max_x3g,min_x3g,max_x3b,min_x3b];
end

% caa = floor(X3/(b*c*d));
% cha = floor(mod(X3,(b*c*d))/(c*d));
% cva = floor(mod(X3,(c*d))/d);
% cda = mod(X3,d);
% aaa = [caa,cha,cva,cda];
% aaa0 = sum(aaa(:)==0)
% aaa1 = sum(aaa(:)==1)
% aaa2 = sum(aaa(:)==2)
% aaa3 = sum(aaa(:)==3)

[hist,~] = imhist(uint8(X3));

[~,index] = max (hist);

% 直方图移位
% figure(51)
% imhist(uint8(X3))
% [~,index] = max (hist)
% mean_hist = mean(hist)
% var_hist = var(hist)
% figure(53)
% imhist(uint8(floor(X3/16)))
% figure(54)
% imhist(uint8(mod(X3,16)))


% 第二种直方图移位方式
% X31 = floor(X3/16);
% X32 = mod(X3,16);
% X3 = mod(X31-8,16)*16 + mod(X32-7,16);


TT = index-tt;
% diyizhong
X3 = mod(X3-TT,256);
% figure(52)
% imhist(uint8(X3))

% figure(55)
% imhist(uint8(floor(X3/16)))
% figure(56)
% imhist(uint8(mod(X3,16)))


ca = floor(X3/(b*c*d));
ch = floor(mod(X3,(b*c*d))/(c*d));
cv = floor(mod(X3,(c*d))/d);
cd = mod(X3,d);

% % 查看分成四份的这几个0,1是否均匀分配
% caa = ca; caa(ca==3)=1;
% sum(caa(:))
% chh = ch; chh(ch==3)=1;
% sum(chh(:))
% cvv = cv; cvv(cv==3)=1;
% sum(cvv(:))
% cdd = cd; cdd(cd==3)=1;
% sum(cdd(:))

% aaa = [ca,ch,cv,cd];
% aaa0 = sum(aaa(:)==0)
% aaa1 = sum(aaa(:)==1)
% aaa2 = sum(aaa(:)==2)
% aaa3 = sum(aaa(:)==3)

[ x1,y1,z1,w1 ] = IICM_4d( [x02,y02,z02,w02],[a2,b2], (m+k*n)*cr );
ca1 = scram( ca,x1,2 );
ch1 = scram( ch,y1,2 );
cv1 = scram( cv,z1,2 );
cd1 = scram( cd,w1,2 );
tcip = ca1*(b*c*d) + ch1*(c*d) + cv1*d + cd1;

if k==3
    tcip = cat(3,tcip(:,1:n/2),tcip(:,n/2+1:2*n/2),tcip(:,2*n/2+1:3*n/2));
end

% 嵌入
[ x2,y2,z2,w2 ] = IICM_4d( [x03,y03,z03,w03],[a3,b3],(m1+k1*n1)*0.5 );
if emb_type == 1
   [ cip, rC ] = embed( tcip,cover,a,b,c,d,x2,y2,z2,w2 );
elseif emb_type == 2
   [ cip, rC ] = embed2( tcip,cover,a,b,c,d,x2,y2,z2,w2 );
end

end



function [ dnkey ] = bit_cor( p,a )
%   DKEY 动态密钥生成
%   p明文，key动态密钥 


% tem = p(:);
% aa = []
% for i=0:a-1
%     aa(i+1)=len(tem(tem == i))
% end

[aa,~] = imhist(uint8(p),a);

x = ones(1,256);
sha_sum = 0;
sha_p = SHA(p,'SHA-256');
for i=1:32  %hex2dec只能到2^52，所以运用循环每8位来一次，也可以其他位数
    tem = ones(1,8);
    sn = dec2bin(hex2dec(sha_p((i-1)*2+1:(i-1)*2+2)),8);   
    kn = dec2bin(hex2dec(key((i-1)*2+1:(i-1)*2+2)),8);
    tem(sn==kn) = 0;
    x((i-1)*8+1:(i-1)*8+8) = tem;
    sha_sum = sha_sum + bin2dec(num2str(tem));
end

a1=0;b1=0;x01=0;y01=0;
a2=0;b2=0;x02=0;y02=0;z02=0;w02=0;
a3=0;b3=0;x03=0;y03=0;z03=0;w03=0;
for i = 1:16
    a1=a1+x(i)*2^(-i);
    b1=b1+x(i+16)*2^(-i);
    x01=x01+x(i+32)*2^(-i);
    y01=y01+x(i+48)*2^(-i);
    
    a2=a2+x(i+64)*2^(-i);
    b2=b2+x(i+80)*2^(-i);
    x02=x02+x(i+96)*2^(-i);
    y02=y02+x(i+112)*2^(-i);
    z02=z02+x(i+128)*2^(-i);
    w02=w02+x(i+144)*2^(-i);
    
    a3=a3+x(i+160)*2^(-i);
    b3=b3+x(i+176)*2^(-i);
    x03=x03+x(i+192)*2^(-i);
    y03=y03+x(i+208)*2^(-i);
    z03=z03+x(i+224)*2^(-i);
    w03=w03+x(i+240)*2^(-i);
end

a1=mod(a1*sha_sum,11)+20;b1=mod(b1*sha_sum,0.11)+0.89;x01=x01*sha_sum/8192;y01=y01*sha_sum/8192;
a2=mod(a2*sha_sum,11)+20;b2=mod(b2*sha_sum,0.11)+0.89;x02=x02*sha_sum/8192;y02=y02*sha_sum/8192;z02=z02*sha_sum/8192;w02=w02*sha_sum/8192;
a3=mod(a3*sha_sum,11)+20;b3=mod(b3*sha_sum,0.11)+0.89;x03=x03*sha_sum/8192;y03=y03*sha_sum/8192;z03=z03*sha_sum/8192;w03=w03*sha_sum/8192;
dnkey = [a1,b1,x01,y01,a2,b2,x02,y02,z02,w02,a3,b3,x03,y03,z03,w03];

end


function [ dnkey ] = dkey( p,key )
%   DKEY 动态密钥生成
%   p明文，key动态密钥 

x = ones(1,256);
sha_sum = 0;
sha_p = SHA(p,'SHA-256');
for i=1:32  %hex2dec只能到2^52，所以运用循环每8位来一次，也可以其他位数
    tem = ones(1,8);
    sn = dec2bin(hex2dec(sha_p((i-1)*2+1:(i-1)*2+2)),8);   
    kn = dec2bin(hex2dec(key((i-1)*2+1:(i-1)*2+2)),8);
    tem(sn==kn) = 0;
    x((i-1)*8+1:(i-1)*8+8) = tem;
    sha_sum = sha_sum + bin2dec(num2str(tem));
end

a1=0;b1=0;x01=0;y01=0;
a2=0;b2=0;x02=0;y02=0;z02=0;w02=0;
a3=0;b3=0;x03=0;y03=0;z03=0;w03=0;
for i = 1:16
    a1=a1+x(i)*2^(-i);
    b1=b1+x(i+16)*2^(-i);
    x01=x01+x(i+32)*2^(-i);
    y01=y01+x(i+48)*2^(-i);
    
    a2=a2+x(i+64)*2^(-i);
    b2=b2+x(i+80)*2^(-i);
    x02=x02+x(i+96)*2^(-i);
    y02=y02+x(i+112)*2^(-i);
    z02=z02+x(i+128)*2^(-i);
    w02=w02+x(i+144)*2^(-i);
    
    a3=a3+x(i+160)*2^(-i);
    b3=b3+x(i+176)*2^(-i);
    x03=x03+x(i+192)*2^(-i);
    y03=y03+x(i+208)*2^(-i);
    z03=z03+x(i+224)*2^(-i);
    w03=w03+x(i+240)*2^(-i);
end

a1=mod(a1*sha_sum,11)+20;b1=mod(b1*sha_sum,0.11)+0.89;x01=x01*sha_sum/8192;y01=y01*sha_sum/8192;
a2=mod(a2*sha_sum,11)+20;b2=mod(b2*sha_sum,0.11)+0.89;x02=x02*sha_sum/8192;y02=y02*sha_sum/8192;z02=z02*sha_sum/8192;w02=w02*sha_sum/8192;
a3=mod(a3*sha_sum,11)+20;b3=mod(b3*sha_sum,0.11)+0.89;x03=x03*sha_sum/8192;y03=y03*sha_sum/8192;z03=z03*sha_sum/8192;w03=w03*sha_sum/8192;
dnkey = [a1,b1,x01,y01,a2,b2,x02,y02,z02,w02,a3,b3,x03,y03,z03,w03];

end


function [ X3,max_x3,min_x3 ] = compre( image,x,y,cr )
%   COMPRE 此处显示有关此函数的摘要
%   此处显示详细说明

x = 1-2*mod(x * 10^4, 1);
y = 1-2*mod(y * 10^4, 1);

[rows, columns] = size(image);
X1 = dct2(image);

R1 = reshape(x, [int16(rows * cr), rows]);
R2 = reshape(y, [int16(columns * cr), columns]);

% 优化测量矩阵
R1(:, 1:uint16(cr * rows)) = R1(:, 1:uint16(cr * rows)) * 5000;
R1 = orth(R1')';
R2(:, 1:uint16(cr * columns)) = R2(:, 1:uint16(cr * columns)) * 5000;
R2 = orth(R2')';

X3 = R1 * X1 * R2';

max_x3 = max(max(X3));
min_x3 = min(min(X3));

X3 = round((X3 - min_x3) / (max_x3 - min_x3) * 255);

% ca2 = max_x3-min_x3; 
% ca3 = (max_x3+min_x3)/2; 
% X3 = round(255* (1 + exp(1)^( -ca2*( X3-ca3 )) ) );

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


function [ cip, rC ] = embed( im,cover,a,b,c,d,x,y,z,w )
%   embed 自己的嵌入操作(空间域)
%   此处显示详细说明

[m, n, k] = size(cover);
[m1,n1,k1] = size(im);
im = reshape(im,[1,m1*n1*k1]);

% % 修改封面 (整数小波才需要修改)
% max1=250;min1=5;
% cover(cover<min1)=min1;
% cover(cover>max1)=max1;

if k==1
    CA = cover(1:m/2,1:n/2);
    CH = cover(m/2+1:m,1:n/2);
    CV = cover(1:m/2,n/2+1:n);
    CD = cover(m/2+1:m,n/2+1:n);
else
    CA = [cover(1:m/2,1:n/2,1),cover(1:m/2,1:n/2,2),cover(1:m/2,1:n/2,3)];
    CH = [cover(m/2+1:m,1:n/2,1),cover(m/2+1:m,1:n/2,2),cover(m/2+1:m,1:n/2,3)];
    CV = [cover(1:m/2,n/2+1:n,1),cover(1:m/2,n/2+1:n,2),cover(1:m/2,n/2+1:n,3)];
    CD = [cover(m/2+1:m,n/2+1:n,1),cover(m/2+1:m,n/2+1:n,2),cover(m/2+1:m,n/2+1:n,3)];
end

% 置乱
CA = scram( CA,x,2 );
CH = scram( CH,y,2 );
CV = scram( CV,z,2 );
CD = scram( CD,w,2 );

% 嵌入
CA1 = reshape(CA,[1,m/2*n/2*k]);
CH1 = reshape(CH,[1,m/2*n/2*k]);
CV1 = reshape(CV,[1,m/2*n/2*k]);
CD1 = reshape(CD,[1,m/2*n/2*k]);

% 第三种方式分解
% bit1 = floor(im/2^7);
% bit2 = floor(mod(im, 2^7) / 2^6);
% bit3 = floor(mod(im, 2^6) / 2^5);
% bit4 = floor(mod(im, 2^5) / 2^4);
% bit5 = floor(mod(im, 2^4) / 2^3);
% bit6 = floor(mod(im, 2^3) / 2^2);
% bit7 = floor(mod(im, 2^2) / 2^1);
% bit8 = mod(im,2);
% sum(bit1)
% sum(bit2)
% sum(bit3)
% sum(bit4)
% sum(bit5) 
% sum(bit6)
% sum(bit7)
% sum(bit8)
% 
% 
% CA11 = bit1*2 + bit5;
% CH11 = bit2*2 + bit6;
% CV11 = bit3*2 + bit7;
% CD11 = bit4*2 + bit8;


% 第一二种方式分解
CA11 = floor(im/(b*c*d));
CH11 = floor(mod(im,(b*c*d))/(c*d));
CV11 = floor(mod(im,(c*d))/d);
CD11 = mod(im,d);

% % 第二三种方式
% CA11_0 = length(find(CA11==0));
% CA11_1 = length(find(CA11==1));
% CA11_2 = length(find(CA11==2));
% CA11_3 = length(find(CA11==3));
% [~,rCA] = sort([CA11_0,CA11_1,CA11_2,CA11_3]);
% CA12 = CA11;
% CA11(CA12==rCA(1)-1) = 2;
% CA11(CA12==rCA(2)-1) = 3;
% CA11(CA12==rCA(3)-1) = 1;
% CA11(CA12==rCA(4)-1) = 0;
% 
% CH11_0 = length(find(CH11==0));
% CH11_1 = length(find(CH11==1));
% CH11_2 = length(find(CH11==2));
% CH11_3 = length(find(CH11==3));
% [~,rCH] = sort([CH11_0,CH11_1,CH11_2,CH11_3]);
% CH12 = CH11;
% CH11(CH12==rCH(1)-1) = 2;
% CH11(CH12==rCH(2)-1) = 3;
% CH11(CH12==rCH(3)-1) = 1;
% CH11(CH12==rCH(4)-1) = 0;
% 
% CV11_0 = length(find(CV11==0));
% CV11_1 = length(find(CV11==1));
% CV11_2 = length(find(CV11==2));
% CV11_3 = length(find(CV11==3));
% [~,rCV] = sort([CV11_0,CV11_1,CV11_2,CV11_3]);
% CV12 = CV11;
% CV11(CV12==rCV(1)-1) = 2;
% CV11(CV12==rCV(2)-1) = 3;
% CV11(CV12==rCV(3)-1) = 1;
% CV11(CV12==rCV(4)-1) = 0;
% 
% CD11_0 = length(find(CD11==0));
% CD11_1 = length(find(CD11==1));
% CD11_2 = length(find(CD11==2));
% CD11_3 = length(find(CD11==3));
% [~,rCD] = sort([CD11_0,CD11_1,CD11_2,CD11_3]);
% CD12 = CD11;
% CD11(CD12==rCD(1)-1) = 2;
% CD11(CD12==rCD(2)-1) = 3;
% CD11(CD12==rCD(3)-1) = 1;
% CD11(CD12==rCD(4)-1) = 0;
% 
% rC = [rCA,rCH,rCV,rCD];

% 第一种方式
rC = [0,0,0,0];


CA1(1:m1*n1*k1) = mod(CA11 + mod(CA1(1:m1*n1*k1),a), a) + floor(CA1(1:m1*n1*k1)/a)*a;
CH1(1:m1*n1*k1) = mod(CH11 + mod(CH1(1:m1*n1*k1),b), b) + floor(CH1(1:m1*n1*k1)/b)*b;
CV1(1:m1*n1*k1) = mod(CV11 + mod(CV1(1:m1*n1*k1),c), c) + floor(CV1(1:m1*n1*k1)/c)*c;
CD1(1:m1*n1*k1) = mod(CD11 + mod(CD1(1:m1*n1*k1),d), d) + floor(CD1(1:m1*n1*k1)/d)*d;

CA1 = reshape(CA1,[m/2,k*n/2]);
CH1 = reshape(CH1,[m/2,k*n/2]);
CV1 = reshape(CV1,[m/2,k*n/2]);
CD1 = reshape(CD1,[m/2,k*n/2]);

% 修改
CA2 = correction(CA,CA1,a);
CH2 = correction(CH,CH1,b);
CV2 = correction(CV,CV1,c);
CD2 = correction(CD,CD1,d);

% 反向置乱
CA2 = dscram( CA2,x,2 );
CH2 = dscram( CH2,y,2 );
CV2 = dscram( CV2,z,2 );
CD2 = dscram( CD2,w,2 );

% 合并
if k==3
    CA2 = cat(3,CA2(:,1:n/2),CA2(:,n/2+1:2*n/2),CA2(:,2*n/2+1:3*n/2));
    CH2 = cat(3,CH2(:,1:n/2),CH2(:,n/2+1:2*n/2),CH2(:,2*n/2+1:3*n/2));
    CV2 = cat(3,CV2(:,1:n/2),CV2(:,n/2+1:2*n/2),CV2(:,2*n/2+1:3*n/2));
    CD2 = cat(3,CD2(:,1:n/2),CD2(:,n/2+1:2*n/2),CD2(:,2*n/2+1:3*n/2));
end

cip = [CA2,CV2;CH2,CD2];


end


function [ cip,rC ] = embed2( im,cover,a,b,c,d,x,y,z,w )
%   embed 自己的嵌入操作(小波变换域)
%   此处显示详细说明

[m, n, k] = size(cover);
[m1,n1,k1] = size(im);
im = reshape(im,[1,m1*n1*k1]);

% 修改封面 (整数小波才需要修改)
max1=252;min1=3;
cover(cover<min1)=min1;
cover(cover>max1)=max1;

% 提升的小波变换
cover=double(cover);
im=double(im);
LS=liftwave('haar','Int2Int');
if k==1
    [CA,CH,CV,CD]=lwt2(cover,LS);
else
    [CAr,CHr,CVr,CDr]=lwt2(cover(:,:,1),LS);
    [CAg,CHg,CVg,CDg]=lwt2(cover(:,:,2),LS);
    [CAb,CHb,CVb,CDb]=lwt2(cover(:,:,3),LS);
    CA = [CAr,CAg,CAb];
    CH = [CHr,CHg,CHb];
    CV = [CVr,CVg,CVb];
    CD = [CDr,CDg,CDb];
end

% 置乱
CA = scram( CA,x,2 );
CH = scram( CH,y,2 );
CV = scram( CV,z,2 );
CD = scram( CD,w,2 );


% 嵌入
CA1 = reshape(CA,[1,m/2*n/2*k]);
CH1 = reshape(CH,[1,m/2*n/2*k]);
CV1 = reshape(CV,[1,m/2*n/2*k]);
CD1 = reshape(CD,[1,m/2*n/2*k]);


CA11 = floor(im/(b*c*d));
CV11 = floor(mod(im,(b*c*d))/(c*d));
CH11 = floor(mod(im,(c*d))/d);
CD11 = mod(im,d);


% CA11_0 = length(find(CA11==0));
% CA11_1 = length(find(CA11==1));
% CA11_2 = length(find(CA11==2));
% CA11_3 = length(find(CA11==3));
% [~,rCA] = sort([CA11_0,CA11_1,CA11_2,CA11_3]);
% CA12 = CA11;
% CA11(CA12==rCA(1)-1) = 2;
% CA11(CA12==rCA(2)-1) = 3;
% CA11(CA12==rCA(3)-1) = 1;
% CA11(CA12==rCA(4)-1) = 0;
% 
% CH11_0 = length(find(CH11==0));
% CH11_1 = length(find(CH11==1));
% CH11_2 = length(find(CH11==2));
% CH11_3 = length(find(CH11==3));
% [~,rCH] = sort([CH11_0,CH11_1,CH11_2,CH11_3]);
% CH12 = CH11;
% CH11(CH12==rCH(1)-1) = 2;
% CH11(CH12==rCH(2)-1) = 3;
% CH11(CH12==rCH(3)-1) = 1;
% CH11(CH12==rCH(4)-1) = 0;
% 
% CV11_0 = length(find(CV11==0));
% CV11_1 = length(find(CV11==1));
% CV11_2 = length(find(CV11==2));
% CV11_3 = length(find(CV11==3));
% [~,rCV] = sort([CV11_0,CV11_1,CV11_2,CV11_3]);
% CV12 = CV11;
% CV11(CV12==rCV(1)-1) = 2;
% CV11(CV12==rCV(2)-1) = 3;
% CV11(CV12==rCV(3)-1) = 1;
% CV11(CV12==rCV(4)-1) = 0;
% 
% CD11_0 = length(find(CD11==0));
% CD11_1 = length(find(CD11==1));
% CD11_2 = length(find(CD11==2));
% CD11_3 = length(find(CD11==3));
% [~,rCD] = sort([CD11_0,CD11_1,CD11_2,CD11_3]);
% CD12 = CD11;
% CD11(CD12==rCD(1)-1) = 2;
% CD11(CD12==rCD(2)-1) = 3;
% CD11(CD12==rCD(3)-1) = 1;
% CD11(CD12==rCD(4)-1) = 0;
% 
% rC = [rCA,rCH,rCV,rCD];

rC = [0,0,0,0];

CA1(1:m1*n1*k1) = mod(CA11 + mod(CA1(1:m1*n1*k1),a), a) + floor(CA1(1:m1*n1*k1)/a)*a;
CH1(1:m1*n1*k1) = mod(CH11 + mod(CH1(1:m1*n1*k1),b), b) + floor(CH1(1:m1*n1*k1)/b)*b;
CV1(1:m1*n1*k1) = mod(CV11 + mod(CV1(1:m1*n1*k1),c), c) + floor(CV1(1:m1*n1*k1)/c)*c;
CD1(1:m1*n1*k1) = mod(CD11 + mod(CD1(1:m1*n1*k1),d), d) + floor(CD1(1:m1*n1*k1)/d)*d;


CA1 = reshape(CA1,[m/2,k*n/2]);
CH1 = reshape(CH1,[m/2,k*n/2]);
CV1 = reshape(CV1,[m/2,k*n/2]);
CD1 = reshape(CD1,[m/2,k*n/2]);

% 修改
CA2 = correction(CA,CA1,a);
CH2 = correction(CH,CH1,b);
CV2 = correction(CV,CV1,c);
CD2 = correction(CD,CD1,d);

% 反向置乱
CA2 = dscram( CA2,x,2 );
CH2 = dscram( CH2,y,2 );
CV2 = dscram( CV2,z,2 );
CD2 = dscram( CD2,w,2 );

% 合并
if k==1
    cip = ilwt2(CA2,CH2,CV2,CD2,LS);
else
    cipr = ilwt2(CA2(:,1:n/2),CH2(:,1:n/2),CV2(:,1:n/2),CD2(:,1:n/2),LS);
    cipg = ilwt2(CA2(:,n/2+1:2*n/2),CH2(:,n/2+1:2*n/2),CV2(:,n/2+1:2*n/2),CD2(:,n/2+1:2*n/2),LS);
    cipb = ilwt2(CA2(:,2*n/2+1:3*n/2),CH2(:,2*n/2+1:3*n/2),CV2(:,2*n/2+1:3*n/2),CD2(:,2*n/2+1:3*n/2),LS);
    cip = cat(3,cipr,cipg,cipb);
end

end


function [ cip2 ] = correction( cover,cip,e )
%   CORRECTION 修正
%   此处显示详细说明

tem = cip-cover;
cip2 = cip;
cip2(tem < -e/2) = cip(tem < -e/2) + e;
cip2(tem > e/2) = cip(tem > e/2) - e;
% 整数嵌入需要注意不要超范围
cip2(cip2 < 0) = cip(cip2 < 0);
cip2(cip2 > 255) = cip(cip2 > 255);

end




