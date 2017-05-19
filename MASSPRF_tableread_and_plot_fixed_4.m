%this section reads the output table files and writes them into cell arrays  it
%is currently designed to work with the simulated massprf sequences. 
clear
outfileNames= {};
prefixout='*Table.txt';
outdirectory='/Users/zimingzhao/Dropbox/manuscripts/macprf_mss/SSD_massprf_humanGenes/Supplement69Genes/massprf_output/output_tables/'
 outfileNames={};
outdirec= dir([ outdirectory prefixout]);
[outfileNames{1:length(outdirec),1}]=deal(outdirec.name); 
foldersize = length(outfileNames);
h=zeros(1,foldersize);
p=zeros(1,foldersize);
Pr=zeros(foldersize,1);
Ps=zeros(foldersize,1);
Dr=zeros(foldersize,1);
Ds=zeros(foldersize,1);
gammas={};
DItimes={};
lowCIs={};
hiCIs={};
for i = 1:foldersize
    FILENAME=[outdirectory outfileNames{i}]; 
filedata=tdfread(FILENAME);
        poses{i,1}=filedata.Position;  
        genelength(i)=length([poses{i,1}])-1;
        genelengths(i,1)=max(filedata.Position);
if ischar(filedata.Gamma)==1
  tmpgam=filedata.Gamma;
 tmplowCI=filedata.Lower_CI_Gamma;
tmphiCI=filedata.Upper_CI_Gamma;
tmpDIT=filedata.DivergentTime;
 gammas{i,1}=zeros(length(tmpgam),1);
lowCIs{i,1}=zeros(length(tmpgam),1);
  hiCIs{i,1}=zeros(length(tmpgam),1);
  DItimes{i,1}=zeros(length(tmpgam),1);
for j = 1 :length(tmpgam)
    if tmpgam(j,1)=='N'
            gammas{i,1}(j)=NaN;
                        lowCIs{i,1}(j)=NaN;
                        hiCIs{i,1}(j)=NaN;
                        DItimes{i,1}(j)=NaN;
        else
       gammas{i,1}(j)=str2num([tmpgam(j,:)]);
                lowCIs{i,1}(j)=str2num([tmplowCI(j,:)]);
    hiCIs{i,1}(j)=str2num([tmphiCI(j,:)]);  
    %DItimes{i,1}(j)=str2num([tmpDIT(j,:)]);
        end
    end
else
                lowCIs{i,1}=filedata.Lower_CI_Gamma;
    hiCIs{i,1}=filedata.Upper_CI_Gamma;
    gammas{i,1}=filedata.Gamma;
    DItimes{i,1}=filedata.DivergentTime;
end
    Ps(i,1)=length(find(filedata.PolymorphismMutationStatus=='S'));
    Pr(i,1)=length(find(filedata.PolymorphismMutationStatus=='R'));
    Ds(i,1)=length(find(filedata.DivergenceMutationStatus=='S'));
    Dr(i,1)=length(find(filedata.DivergenceMutationStatus=='R'));
[h(i) p(i)]=fishertest([Ps(i) Pr(i); Ds(i) Dr(i)]);
gamvars(i)=var(gammas{i,1});
end
maxvars=(gamvars==max(gamvars))
minvars=(gamvars==min(gamvars))
NI=(Pr./Ps)./(Dr./Ds);
%%

%this section will output traces with the original numbers of sites. 
%This script takes output files of abridged MASSPRF (preprocessing to
%generate consensus sequences) and outputs SCALED fasta files for input
%into MASSPRF. It can also be used to perform the McDonald-Kreitman test. 

%Daniel Lee, July 2016. Email dslee.sm@gmail.com with issues or questions. 
condirectory='/Users/zimingzhao/Dropbox/manuscripts/macprf_mss/SSD_massprf_humanGenes/Supplement69Genes/preprocessed_ConsensusSequence/'
l=0;
folder_size=length(outfileNames)
plotoutdirectory='/Users/zimingzhao/Dropbox/manuscripts/macprf_mss/SSD_massprf_humanGenes/Supplement69Genes/SVG_plots/'
for i = 1:folder_size

FILE=[outfileNames{i,1}];
genename=FILE(1:min(find(FILE=='_'))-1);
conprefix=[genename '*']
condirec= dir([condirectory conprefix]);
[consensusfileNames{1:length(condirec),1}]=deal(condirec.name); 
oldFILENAME=[consensusfileNames{1}]
disp([genename ' ' num2str(i) ' out of '  num2str(folder_size)])
    text1=fileread([condirectory oldFILENAME]);
%read the file in 
text1(text1==' ') = '';
newline=sprintf('\n');
text1(text1==newline)='';
starindexes=find(text1=='*');
Polstart=min(find(text1=='*'));
Polend=max(find(text1=='D')-1);
Polseq=text1(Polstart:Polend);
genelength=length(Polseq);
Divstart=Polend+12;
Divend=Polend+12+genelength-1;
Divseq=text1(Divstart:Divend);
origlength=length(Divseq);
%Find Dr Ds Ps Pr
Pslocs=find(Polseq=='S');
Prlocs=find(Polseq=='R');
Dslocs=find(Divseq=='S');
Drlocs=find(Divseq=='R');
Ds(i)=length(Dslocs);
Dr(i)=length(Drlocs);
Pr(i)=length(Prlocs);
Ps(i)=length(Pslocs);
[h(i) p(i)]=fishertest([Ps(i) Pr(i) ; Ds(i) Dr(i)]);  %Perform fisher test 
%for significance. Should comment out if not using it; the function is VERY
%slow. 
%fig=figure(i);
 % set(fig, 'Position', [100, 100, 1920, 1080]);
    gam=gammas{i}';
hiCI=hiCIs{i}';
lowCI=lowCIs{i}';
pos=[poses{i}];
pos=pos*origlength/length(pos);
%l=l+max(pos)
fig=figure(i)
myshadedErrorBar(pos,gam,[lowCI; hiCI],[],'1');
hold on
plot(zeros(1,origlength),'-k','LineWidth',1.3) 
xlabel('Nucleotide sequence position','FontSize',16,'FontName','Arial')
ylabel('Selection intensity (\gamma)','FontSize',16,'FontName','Arial')

xloc=.75*origlength;
yloc=ceil(max(hiCI));
IPsstr =['P_{s}'];
IPrstr= ['P_{r}'];
IDsstr =['D_{s}'];
IDrstr =['D_{r}'];

Psstr =['P_{s} = ' num2str(Ps(i))];
Prstr= ['P_{r} = ' num2str(Pr(i))];
Dsstr =['D_{s} = ' num2str(Ds(i))];
Drstr =['D_{r} = ' num2str(Dr(i))];
pprestr=['P  '];
if p(i)>.01
pstr =[' =  ' num2str(p(i),1)];
else 
   expnum=ceil(log10(p(i)));
    pstr =['< 10^{' num2str(expnum) '} '];    
end
outfigfilename = ([plotoutdirectory genename '_MASSPRFtrace']);
fig.Position=[0 0 1920 1080];
xlim([1 origlength])
ymin=floor(min(lowCI));
ymax=ceil(max(hiCI));
ylim([ymin ymax ])
text(xloc,(ymax-ymin)*.9+ymin,Psstr,'FontSize',16,'FontName','Arial')
text(xloc,(ymax-ymin)*.84+ymin,Prstr,'FontSize',16,'FontName','Arial')
text(xloc,(ymax-ymin)*.78+ymin,Dsstr,'FontSize',16,'FontName','Arial')
text(xloc,(ymax-ymin)*.72+ymin,Drstr,'FontSize',16,'FontName','Arial')
text(xloc,(ymax-ymin)*.60+ymin,pprestr,'FontSize',16,'FontName','Arial','FontAngle','italic')
text(xloc+.0151*origlength,(ymax-ymin)*.60+ymin,pstr,'FontSize',16,'FontName','Arial')
title(genename,'FontSize',20,'FontName','Arial','FontWeight','bold','FontAngle','italic')
print(fig,outfigfilename,'-dpdf','-painters')
close all
end
[FDR Q Pi0]=mafdr(p); %calculate FDR

