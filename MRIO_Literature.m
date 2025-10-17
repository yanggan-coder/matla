clc
clear

%% Linking EF to Final Demand 

cd('D:\HuaweiMoveData\Users\yiziyaoyao\Desktop\exercises')
load('POPGDP.mat')
load('EXIOBASE_labs.mat')
t=2019;
pop=POP(:,t-1994);  
gdp=GDP(:,t-1994);
gdp_pp=gdp./pop.*1e6;    % convert million euros to euros
r=size(GDP,1);
s=size(labs_product,1);
d=7;

% Import data
cd('D:\HuaweiMoveData\Users\yiziyaoyao\Desktop\exercises\IOT_2019_pxp')
A=importdata('A.txt');
A=A.data;
Y=importdata('Y.txt');
Y=Y.data;
x=importdata('x.txt');
x=x.data;

cd('D:\HuaweiMoveData\Users\yiziyaoyao\Desktop\exercises\IOT_2019_pxp\air_emissions')
F_emission=importdata('F.txt');
F_emission=F_emission.data;
Fhh_emission=importdata('F_Y.txt');
Fhh_emission=Fhh_emission.data;
labs_emission=importdata('unit.txt');
labs_emission=labs_emission(2:end);  % remove the header

% Create essential MRIO variables
I=eye(r*s);
L=inv(I-A);
f_emission=F_emission./x';        % emission intensity
f_emission(isnan(f_emission))=0;
f_emission(isinf(f_emission))=0;

Yt=zeros(r*s,r);                               % aggregate final demand 'agents'
Fhht_emission=zeros(size(Fhh_emission,1),r);   % aggregate final demand 'agents'
for i=1:r
    Yt(:,i)=sum(Y(:,(i-1)*d+1:i*d),2);
    Fhht_emission(:,i)=sum(Fhh_emission(:,(i-1)*d+1:i*d),2);
end

% Environmental Extensions of interest
index_co2=find(contains(labs_emission,'CO2'));
labs_emission(index_co2) % display what's found by keyword
index_fc=index_co2([5 6]);
index_co2=index_co2([1 3 4 7 8 9]);
index_ch4=find(contains(labs_emission,'CH4'));
index_n2o=find(contains(labs_emission,'N2O'));

F_co2=sum(F_emission(index_co2,:),1);
F_ch4=sum(F_emission(index_ch4,:),1)*29.8;  % convert to co2 eq.
F_n2o=sum(F_emission(index_n2o,:),1)*273;   % convert to co2 eq.
F_fc=sum(F_emission(index_fc,:),1);         % already in co2 eq.

f_co2=sum(f_emission(index_co2,:),1);
f_ch4=sum(f_emission(index_ch4,:),1)*29.8;  % convert to co2 eq.
f_n2o=sum(f_emission(index_n2o,:),1)*273;   % convert to co2 eq.
f_fc=sum(f_emission(index_fc,:),1);         % already in co2 eq.

f_ghg=f_co2+f_ch4+f_n2o+f_fc;

Fhht_co2=sum(Fhht_emission(index_co2,:),1);
Fhht_ch4=sum(Fhht_emission(index_ch4,:),1)*29.8;
Fhht_n2o=sum(Fhht_emission(index_n2o,:),1)*273;
Fhht_fc=sum(Fhht_emission(index_fc,:),1);

GHG_global=[sum(F_co2)+sum(Fhht_co2) sum(F_ch4)+sum(Fhht_ch4) sum(F_n2o)+sum(Fhht_n2o) sum(F_fc)+sum(Fhht_fc)]*1e-12;   % convert from kg to Gt

% GHG footprint of nations
Fhht_ghg=Fhht_co2+Fhht_ch4+Fhht_n2o+Fhht_fc;
E_ghg=f_ghg*L*Yt + Fhht_ghg;
E_ghg_pp=E_ghg./pop'*1e-3;       % convert from kg CO2ea to tonne CO2eq per person
bar(E_ghg_pp)

% GHG footprint by country and final demand purpose (ignore Fhh here)
labs_dp=unique(labs_product(:,2),'stable');  % final demand purpose 
cc_dp=zeros(s,8);                            % concordance matrix linking 200 products to 8 final demand purposes
for j=1:size(labs_dp,1)
    cc_dp(contains(labs_product(:,2),labs_dp(j)),j)=1;
end
%%
E_ghg_s=zeros(r,s);
for i=1:r
    temp=f_ghg*L*diag(Yt(:,i));   % diagnolization brings out the product-country details in final demand
    E_ghg_s(i,:)=sum(reshape(temp',[s,r]),2);  %shape 'column' vector to product by region, then only keep product details 
end

E_ghg_dp=E_ghg_s*cc_dp;           % turn 200 products to final demand purpose (8 categories); ignore Fhh here
E_ghg_dp_pp=E_ghg_dp./pop*1e-3;   % ignore Fhh here
subplot(1,2,1)
bar(E_ghg_pp)
subplot(1,2,2)
bar(E_ghg_dp_pp,'stacked')

% Global GHG footprint by final demand agent and purpose （ignore Fhh here)
Y_global_ag=zeros(r*s,3);  % 3 columns: household, government, and investment 
Y_global_ag(:,1)=sum(Y(:,1:d:r*d),2);                         % Household final consumption
Y_global_ag(:,2)=sum(Y(:,2:d:r*d),2) + sum(Y(:,3:d:r*d),2);   % Government + NPISH final consumption
Y_global_ag(:,3)=sum(Y(:,4:d:r*d),2);                         % Inverstment or GFCF

E_ghg_global_ag=f_ghg*L*Y_global_ag;
E_ghg_global_ag_s=zeros(3,s);
for i=1:3
    temp=f_ghg*L*diag(Y_global_ag(:,i));
    E_ghg_global_ag_s(i,:)=sum(reshape(temp',[s,r]),2);
end
E_ghg_global_ag_dp=E_ghg_global_ag_s*cc_dp;
barh(E_ghg_global_ag_dp,'stacked');

% Carbon footprint as a function of the nations’ final expenditure (ignore Fhh here)
CF_s=zeros(r,s);
for i=1:r
    CF_s(i,:)=sum(reshape(f_co2*L*diag(Yt(:,i)),[s,r]),2);
end
CF_dp_pp=CF_s*cc_dp./pop.*1e-3;   % convert from kg to tonne
gdp_pp=gdp./pop;
for j=1:8
    subplot(3,3,j)
    scatter(log10(gdp_pp),log10(CF_dp_pp(:,j)));
    title(labs_dp(j))
end
    
%% Tracing Water Footprint Up Production Supply Chains

% Environmental extension: blue water consumption
cd('D:\HuaweiMoveData\Users\yiziyaoyao\Desktop\exercises\IOT_2019_pxp\water')
labs_water=importdata('unit.txt');
labs_water=labs_water(2:end);  % remove the header
F_water=importdata('F.txt');
F_water=F_water.data;
Fhh_water=importdata('F_Y.txt');
Fhh_water=Fhh_water.data;

index_bwc=find(contains(labs_water,'Water Consumption Blue'));
disp(labs_water(index_bwc,:))
F_bwc=sum(F_water(index_bwc,:),1);
Fhh_bwc=sum(Fhh_water(index_bwc,:),1);
f_bwc=F_bwc./x';  % unit: m3/euro
f_bwc(isnan(f_bwc))=0;
f_bwc(isinf(f_bwc))=0;

% Supply chain analysis (power series expansion, structural path analysis - SPA)
WF_spa=zeros(r,11);
for j=1:11
    if j==1
        y=Yt;
    end
    WF_spa(:,j)=f_bwc*y;
    y=A*y;   % create new y for the next production layer
end
WF_pld=zeros(r,12);  % production layer decomposition
WF_pld(:,1:11)=cumsum(WF_spa,2);
WF_pld(:,12)=f_bwc*L*Yt;
area(WF_pld')

%% Environemntal Footprint Dashboard

% Environmental extensions: carbon emissions, blue water consumption, land use, and materials extraction
cd('D:\HuaweiMoveData\Users\yiziyaoyao\Desktop\exercises\IOT_2019_pxp\land')
labs_land=importdata('unit.txt');
labs_land=labs_land(2:end);
F_land=importdata('F.txt');
F_land=sum(F_land.data,1);
Fhh_land=importdata('F_Y.txt');
Fhh_land=sum(Fhh_land.data,1);
f_land=F_land./x';
f_land(isnan(f_land))=0;
f_land(isinf(f_land))=0;

cd('D:\HuaweiMoveData\Users\yiziyaoyao\Desktop\exercises\IOT_2019_pxp\material')
labs_mat=importdata('unit.txt');
labs_mat=labs_mat(2:end);
F_mat=importdata('F.txt');
F_mat=sum(F_mat.data,1);
Fhh_mat=importdata('F_Y.txt');
Fhh_mat=sum(Fhh_mat.data,1);
f_mat=sum(F_mat,1)./x';
f_mat(isnan(f_mat))=0;
f_mat(isinf(f_mat))=0;

F_co2=sum(F_emission(index_co2,:),1);
Fhh_co2=sum(Fhh_emission(index_co2,:),1);

f4=[f_co2;f_bwc;f_land;f_mat];
F4=[F_co2;F_bwc;F_land;F_mat];
Fhh4=[Fhh_co2;Fhh_bwc;Fhh_land;Fhh_mat];

PBA=zeros(4,r);
Fhht4=zeros(4,r);
for i=1:r
    Fhht4(:,i)=sum(Fhh4(:,(i-1)*d+1:i*d),2);        % aggregate from 9800 to 49
    PBA(:,i)=sum(F4(:,(i-1)*s+1:i*s),2)+Fhht4(:,i); % aggregate from 343 to 49 + Fhh
end

CBA=f4*L*Yt+Fhht4;

PBA_3r=[sum(PBA(:,1:27),2) PBA(:,29) PBA(:,31) ]';  % EU, US, CN
CBA_3r=[sum(CBA(:,1:27),2) CBA(:,29) CBA(:,31) ]';  % EU, US, CN
Cov_3r=PBA_3r./CBA_3r;
bar(Cov_3r)
hold on
line([0 4], [1 1],'Color','k')
legend('Carbon','Water','Land','Materials')
ax=gca;
ax.XTickLabel={'EU','US','CN'};