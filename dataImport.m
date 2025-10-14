cd('E:\EEIOA2025\exiobase\IOT_2019_pxp')   % set directory
r=49;  % number of regions
s=200; % number of sectors
d=7;   % number of final demand categories

%% economic core

unit=readtable('unit.txt');
sectors=unique(unit.sector,'stable');
regions=unique(unit.region,'stable');

% x=importdata('x.txt');
% x=x.data;
Y=importdata('Y.txt');
Y=Y.data;
A=importdata('A.txt');
A=A.data;
% Z=importdata('Z.txt');
% Z=Z.data;

I=eye(r*s);   % identity matrix
L=inv(I-A);   % Leontief inverse

Yt=zeros(r*s,r); % aggregate final demand categories
for i=1:r
    Yt(:,i)=sum(Y(:,(i-1)*d+1:i*d),2);
end

x=L*sum(Y,2);
Z=A*diag(x);

% checks
test_Z=A*diag(x);
plot(sum(test_Z,1)./sum(Z,1))

test_x=sum(test_Z,2)+sum(Y,2);
plot(test_x./x)

test_Yi=sum(Y(:,(10-1)*d+1:10*d),2)./Yt(:,10);
plot(test_Yi)
disp(sum(sum(Yt))-sum(sum(Y)))

%% environmental extensions
cd('E:\EEIOA2025\exiobase\IOT_2019_pxp\employment')   % set directory
F_emp=importdata('F.txt');
emp=F_emp.textdata(4:end);
F_emp=F_emp.data;

Fhh_emp=importdata('F_Y.txt');
Fhh_emp=Fhh_emp.data;
Fhht_emp=zeros(12,r);
for i=1:r
    Fhht_emp(:,i)=sum(Fhh_emp(:,(i-1)*d+1:i*d),2);
end

f_emp=F_emp./x';
f_emp(isnan(f_emp))=0;
f_emp(isinf(f_emp))=0;

%% PBA & CBA

PBA_emp=sum(F_emp(1:6,:),1);
PBA_emp=sum(reshape(PBA_emp',[s,r]),1);

CBA_emp=sum(f_emp(1:6,:),1)*L*Yt + sum(Fhht_emp(1:6,:),1); % fLY+Fhh
CBA_emp_sl=f_emp(1:6,:)*L*Yt + Fhht_emp(1:6,:);
plot(CBA_emp,'LineWidth',2)
hold on
bar(CBA_emp_sl','stacked')

%% Emissions
cd('E:\EEIOA2025\exiobase\IOT_2019_pxp\air_emissions')   % set directory
F_air=importdata('F.txt');
air=F_air.textdata(4:end,:);
F_air=F_air.data;
index_co2=find(contains(air,'CO2'));
f_co2=sum(F_air(index_co2,:)./x',1);
f_co2(isnan(f_co2))=0;

f_co2(isinf(f_co2))=0;

CBA_co2=f_co2*L*Yt;  % we are ignoring Fhh here due to time limit





