%this parses tree data outputs with recombination. Requires indelible
%output and ms input. 
clear
numseqs=101;
directory = 'C:\Users\dswle_000\Documents\massprfmss\macprf_mss\reveiwers\simulations\INDELibleV1.03\bin\';
prefix='Expan_RE*' %can loop this for multiple inputs 
msfile=fileread(['C:\Users\dswle_000\Documents\massprfmss\macprf_mss\reveiwers\simulations\msdir\' prefix(1:end-1) '.out']);
        treedirec='C:\Users\dswle_000\Documents\massprfmss\macprf_mss\reveiwers\simulations\INDELibleV1.03\bin\';
outdirec= 'C:\Users\dswle_000\Documents\massprfmss\macprf_mss\reveiwers\simulations\INDELibleV1.03\bin\recomb_fas\';
direc= dir([directory,prefix]);
fileNames= {};
mdls=['_or'; '_ns']
for m=1:2
[fileNames{1:length(direc),1}]=deal(direc.name); 
foldersize = length(fileNames);
%FILENAME=[fileNames{i,1}];
    checklength=[];
slashsites=find(msfile=='/');
    startsites=slashsites(find(msfile(slashsites+1)=='/'));
    semisites=find(msfile==';');
   numtrees=length(startsites);
   for j = 1:numtrees
       if j<numtrees
    subtreesstops=semisites(intersect(find(semisites>startsites(j)),find(semisites<startsites(j+1))));
       else
           subtreesstops=semisites(find(semisites>startsites(j)));   
       end
       subtreestarts=[startsites(j)+9 subtreesstops(1:end-1)+8];
    numsubtrees=length(subtreesstops);
    
    subtreename={};
    subtree={};
    codonadjust(j,1)=0;
    tree =[];
allseqs=[];
    for k = 1:length(subtreesstops)    
        treefile=[prefix(1:end-1) '_r' num2str(j) '_s' num2str(k) mdls(m,:) '.fas']
        subtreefas=fileread([treedirec treefile]);
        allseqs=[allseqs subtreefas];
        if msfile(subtreestarts(k)-1)==']'
       	basenum(k)=str2num([msfile(subtreestarts(k)-4:subtreestarts(k)-2)]);
        elseif msfile(subtreestarts(k)-2)==']'
            basenum(k)=str2num([msfile(subtreestarts(k)-4:subtreestarts(k)-3)]);
        elseif msfile(subtreestarts(k)-3)==']'
            basenum(k)=str2num([msfile(subtreestarts(k)-4)]);
        end
        if mod(basenum+codonadjust(k),3)>0
        codonadjust(j,k+1)=3-mod(basenum(k)+codonadjust(j,k),3);
               else 
        codonadjust(j,k+1)=0;     
        end
        codon_num=(codonadjust(j,k)+basenum(k)+codonadjust(j,k+1))/3; 
    end     
    %excise spaces, numbers, newlines
    origseqs=allseqs;
    origstarts=find(origseqs=='>');
    allseqs(allseqs==' ')='';
    allseqs=regexprep(allseqs,'\n','');
    allseqs=regexprep(allseqs,'[0-9]','');
    allstarts=find(allseqs=='>');
    fidout=fopen([outdirec prefix(1:end-1) '_r' num2str(j) mdls(m,:) '.fas'],'w');
    for l=1:numseqs
        seqstring=[];
        for k=1:length(subtreesstops)
            str=[];   
            str=allseqs(allstarts(l+(k-1)*numseqs)+codonadjust(k)+2:allstarts(l+(k-1)*numseqs)+basenum(k)+1+codonadjust(k));            %I'm not quite sure why this works; I'm rather apprehensive of the fact that it does
      seqstring=[seqstring str];        
        end        
        fprintf(fidout,[origseqs(origstarts(l):origstarts(l)+2) '\n' seqstring '\n']);
        checklength=[length(seqstring) checklength];
    end
    fclose(fidout)
    end
    
end