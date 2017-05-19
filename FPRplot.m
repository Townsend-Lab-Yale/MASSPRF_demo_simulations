clear
format long g
directory = 'C:\Users\dswle_000\Documents\massprfmss\macprf_mss\reveiwers\simulations\simulation_results\';
prefix={'BOTTLE_OUT\output_tables\'; 'CONST_OUT\output_tables\';'EXPAND_OUT\output_tables\'};
%directory = 'C:\Users\dswle_000\Documents\massprfmss\macprf_mss\reveiwers\simulations\simulation_results\SSD_out\output_tables\';
%prefix={'BOTTLE_OUT\'; 'CONSTANT_OUT\';'EXPAND_OUT\'};
suffix='*nsout_Table.txt';
prefix2={'nr\';'re\'};
badcount=zeros(length(prefix),length(prefix2));
condition={}
gammahigh=-10:.1:100;
gammalow=-10:.1:100;
ind=1;
for reornot=1:length(prefix2)
for demo=1:length(prefix)
direc= dir([directory,prefix{demo},prefix2{reornot},suffix]);
fileNames= {};
[fileNames{1:length(direc),1}]=deal(direc.name); 
foldersize = length(fileNames);
poses={};
pols={};
divs={};

for j =1 : foldersize
        FILENAME=[directory prefix{demo} prefix2{reornot} fileNames{j}];
        filedata=tdfread(FILENAME);
       genelength=length(filedata.Position);
       position=zeros(1,genelength);
       gamma=zeros(1,genelength);
       lowCI=zeros(1,genelength);
       hiCI=zeros(1,genelength);
    
        position(1,:)=filedata.Position(:,1)';
        gamma(1,:)=filedata.Gamma(:,1)';
        lowCI(1,:)=filedata.Lower_CI_Gamma(:,1)';
        hiCI(1,:)=filedata.Upper_CI_Gamma(:,1)';       
    

    poses{j,1}=position;

       while length(lowCI)<900
lowCI=[lowCI NaN];
hiCI=[hiCI NaN];
gamma=[gamma NaN];
badcount(demo,reornot)=badcount(demo,reornot)+1;
       end
    lowCIs{demo,reornot}(j,:)=lowCI;
    hiCIs{demo,reornot}(j,:)=hiCI;
    gammas{demo,reornot}(j,:)=gamma;
             
    pol=[filedata.PolymorphismMutationStatus];
    div=[filedata.DivergenceMutationStatus];
    pol(pol==' ')='';
    div(div==' ')='';
    pols{j,1}=pol;
    divs{j,1}=div;

  gammamat=[gammas{demo,reornot}];
  gammavec=reshape(gammamat,1,size(gammamat,1)*size(gammamat,2));
 %  figure(1)
 % subplot(2,3,(demo-1)+reornot.^2)
end
%   for i = 1:genelength
 %       gamma(i)=filedata.Gamma(i,:);
 %   end
%  bins=-10:5:20;
 %[h,x]= hist(gammavec,bins)
%bar(x,h./sum(h),1.0)
%xlim([min(bins) max(bins)])
%ylim([0 1])
%gca.Tick=[-10:20:40]
 lowCImat=[lowCIs{demo,reornot}];
 hiCImat=[hiCIs{demo,reornot}];
 bpnum=size(lowCImat,1)*size(lowCImat,2);
%neutral{demo,reornot}=1-length(union(find(lowCImat>=0),find(hiCImat<=0)))/bpnum;
%nearlyneutral(demo,reornot)=1- length(union(find(hiCImat<=-1),find(lowCImat>=4)))/bpnum;
%meansgamma(demo,reornot)=mean(gammavec);
%varsgamma(demo,reornot)=var(gammavec,1);
for ii=1:length(gammahigh)
FPRnegative{demo,reornot}(ii)=size(find(hiCImat<=gammalow(ii)),1)/bpnum;
FPRpositive{demo,reornot}(ii)=size(find(lowCImat>=gammahigh(ii)),1)/bpnum;
end
figure(1)
plot(gammalow,FPRnegative{demo,reornot})
title('FPR for negative selection')
xlabel('\gamma threshold')
ylabel('False positive rate')
condition{ind}=[fileNames{1}(1:end-20)];
ind=ind+1;
legend()
hold on
figure(2)
plot(gammahigh,FPRpositive{demo,reornot})
title('FPR for positive selection')
xlabel('\gamma threshold')
ylabel('False positive rate')
hold on
end
end
figure(1) 
condition={'Bottleneck without recombination','Constant without recombination','Expansion without recombination','Bottleneck with recombination','Constant with recombination','Expansion with recombination'};
legend(condition)
figure(2) 
legend(condition)