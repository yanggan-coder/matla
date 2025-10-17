%%
cd("C:\Users\11797\OneDrive\Desktop\exercise");
load("POPGDP.mat");
load("EXIOBASE_labs.mat");
t=2019;
pop=POP(:,t-1994);
gdp=GDP(:,t-1994);
gdpppp=gdp./pop.*1e6;
r=size(GDP,1);
s=size(labs_product,1);
d=7;


%%
cd("C:\Users\11797\OneDrive\Desktop\IOT_2019_pxp");
A=importdata("A.txt");
A=A.data;
Y=importdata("Y.txt");
Y=Y.data;
X=importdata("X.txt");
X=X.data;


%%
cd("C:\Users\11797\OneDrive\Desktop\IOT_2019_pxp\air_emissions");
F_emission=importdata('F.txt');
F_emission=F_emission.data;
Fhh_emission=importdata("F_Y.txt");
Fhh_emission=Fhh_emission.data;
labs_emission=importdata("unit.txt");
labs_emission=labs_emission(2:end);


%%
I=eye(r*s);
L=inv(I-A);
f_emission=F_emission./X';
f_emission(isnan(f_emission)) = 0;
f_emission(isinf(f_emission)) = 0;
%%
Yt=zeros(r*s,r);
Fhht_emission=zeros(size(Fhh_emission,1),r);
for i=1:r
    Yt(:,i)=sum(Y(:,(i-1)*d+1:i*d),2);
    Fhht_emission(:,i)=sum(Fhh_emission(:,(i-1)*d+1:i*d),2);
end
%%
index_co2=find(contains(labs_emission,"CO2"));
index_emission=(index_co2);
index_fc=index_co2([5 6]);
index_co2=index_co2([1 3 4 7 8 9]);
index_ch4=find(contains(labs_emission,"CH4"));
index_n2o=find(contains(labs_emission,"N2O"));
%%
F_co2=sum(F_emission(index_co2,:),1);
F_ch4=sum(F_emission(index_ch4,:),1)*29.8;
F_n2o=sum(F_emission(index_n2o,:),1)*273;
F_fc=sum(F_emission(index_fc,:),1);
%%
f_co2=sum(f_emission(index_co2,:),1);
f_ch4=sum(f_emission(index_ch4,:),1)*29.8;
f_n2o=sum(f_emission(index_n2o,:),1)*273;
f_fc=sum(f_emission(index_fc,:),1);
%%
f_ghg=f_co2+f_ch4+f_n2o+f_fc;
%%
Fhht_co2=sum(Fhht_emission(index_co2,:),1);
Fhht_ch4=sum(Fhht_emission(index_ch4,:),1)*29.8;
Fhht_n2o=sum(Fhht_emission(index_n2o,:),1)*273;
Fhht_fc=sum(Fhht_emission(index_fc,:),1);
%%
GhG_global=[sum(F_co2)+sum(Fhht_co2) sum(Fhht_ch4) sum(Fhht_ch4) sum(F_n2o)+sum(Fhht_n2o) sum(F_fc)+sum(Fhht_fc)*1e-12];
%%
Fhht_ghg=Fhht_co2+Fhht_ch4+Fhht_n2o+Fhht_n2o+Fhht_fc;
E_ghg=f_ghg*L*Yt+Fhht_ghg;
E_ghg_pp=E_ghg./pop'*1e-3;
%%
labs_dp=unique(labs_product(:,2),'stable');
cc_dp=zeros(s,8);
for j=1:size(labs_dp,1)
    cc_dp(contains(labs_product(:,2),labs_dp(j)),j)=1;
end
%%
E_ghg_s=zeros(r,s);
for i=1:r
    temp=f_ghg*L*diag(Yt(:,1));
    E_ghg_s(1,:)=sum(reshape(temp',[s,r]),2);
end

E_ghg_dp=E_ghg_s*cc_dp;
E_ghg_dp_pp=E_ghg_dp./pop*1e-3;
subplot(1,2,2);
bar(E_ghg_dp_pp,"stacked");
%%
Y_global_ag=zeros(r*s,3);
Y_global_ag(:,1)=sum(Y(:,1:d:r*d),2);
Y_global_ag(:,2)=sum(Y(:,2:d:r*d),2)+sum(Y(:,3:d:r*d),2);
Y_global_ag(:,3)=sum(Y(:,3:d:r*d),4);

E_ghg_global_ag=f_ghg*L*Y_global_ag;
E_ghg_global_ag_s=zeros(3,s);
for i=1:3
    temp=f_ghg*L*diag(Y_global_ag(:,1));
    E_ghg_global_ag_s(i,:)=sum(reshape(temp',[s,r]),2);
end
E_ghg_global_ag_dp=E_ghg_global_ag_s*cc_dp;
barh(E_ghg_global_ag_dp,"stacked");









%%
cd("C:\Users\11797\OneDrive\Desktop\IOT_2019_pxp\water");
labs_water=importdata("unit.txt");
labs_water=labs_water(2:end);
F_water=importdata("F.txt");
F_water=F_water.data;
Fhh_water=importdata("F_Y.txt");
Fhh_water=Fhh_water.data;



index_bwc=find(contains(labs_water,"Water Consumption Blue"));
disp(labs_water(index_bwc,:));
F_bwc=sum(F_water(index_bwc,:),1);
Fhh_bwc=sum(Fhh_water(index_bwc,:),1);
f_bwc = F_bwc ./ X';
f_bwc(isnan(f_bwc))=0;
f_bwc(isinf(f_bwc))=0;




WF_spa=zeros(r,11);
for j=1:11
    if j==1
        y=Yt;
    end
    WF_spa(:,j)=f_bwc*y;
    y=A*y;
end
WF_pld=zeros(r,12);
WF_pld(:,1:11)=cumsum(WF_spa,2);
WF_pld(:,12)=f_bwc*L*Yt;
area(WF_pld')





%%
cd("C:\Users\11797\OneDrive\Desktop\IOT_2019_pxp\land");
labs_land=importdata("unit.txt");
labs_land=labs_land(2:end);
F_land=importdata("F.txt");
F_land=sum(F_land.data,1);
Fhh_land=importdata("F_Y.txt");
Fhh_land=sum(Fhh_land.data,1);
f_land = F_bwc ./ X';
f_land(isnan(f_land))=0;
f_land(isinf(f_land))=0;




%%
cd("C:\Users\11797\OneDrive\Desktop\IOT_2019_pxp\material");
labs_material=importdata("unit.txt");
labs_material=labs_material(2:end);
F_material=importdata("F.txt");
F_material=sum(F_material.data,1);
Fhh_material=importdata("F_Y.txt");
Fhh_material=sum(Fhh_material.data,1);
f_material = F_bwc ./ X';
f_material(isnan(f_material))=0;
f_material(isinf(f_material))=0;




%%
PBA=zeros(4,r);
Fhht4=zeros(4,r);
for i=1:r
    Fhht4(:,i)=sum(Fhh4(:,(i-1)*d+1,i*d),2);
    
