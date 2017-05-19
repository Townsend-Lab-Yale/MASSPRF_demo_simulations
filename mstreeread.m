%%
%This parses for tree data output by Hudson's ms without recombination and writes it into a
%control file that can be used as an input for indelible. 
clear
directory = 'C:\Users\dswle_000\Documents\massprfmss\macprf_mss\reveiwers\simulations\msdir\';
outdirectory = 'C:\Users\dswle_000\Documents\massprfmss\macprf_mss\reveiwers\simulations\INDELibleV1.03\bin\'
suffix='*NR.out'
direc= dir([directory, suffix]);
fileNames= {};
[fileNames{1:length(direc),1}]=deal(direc.name); 
foldersize = length(fileNames);
fid=fopen([[outdirectory] 'control.txt'],'w')'
nummodels=2;
fprintf(fid,'[TYPE] CODON 1 \n')
fprintf(fid, '[MODEL] neutralmodel \n [submodel] 4.5 0.5 0 1 \n')
fprintf(fid, '[MODEL] oneratio \n [submodel] 4.5 0.5 \n')
scaleparam=100  ; %set distance scaling factor
for i = 1:foldersize
    FILENAME=[fileNames{i,1}];
    text=fileread([directory FILENAME]);
    slashsites=find(text=='/');
    startsites=find(text(slashsites+1)=='/');
    semisites=find(text==';');  
   numtrees=length(startsites);
   for j = 1:numtrees
 startloc=(slashsites(startsites(j))+4);
 endloc=(semisites(j));
 treename{j,i}=[FILENAME(1:end-4) '_r' num2str(j) ' '];
   tree{j,i}=['[TREE] ' [treename{j,i}] text(startloc:endloc) '\n'];
   colonlocs=find([tree{j,i}]==':');
   for k = 1:length(colonlocs)
        colonlocs=find([tree{j,i}]==':');
   periodlocs=find([tree{j,i}]=='.');
   endparlocs=find([tree{j,i}]==')');
   commalocs=find([tree{j,i}]==',');
       branchlengthstr=tree{j,i}(colonlocs(k)+1:colonlocs(k)+5);
        newBL=sprintf('%1.5f',(str2num(branchlengthstr)/scaleparam));
        str1=[tree{j,i}(1:colonlocs(k))];
        str2=[tree{j,i}(colonlocs(k)+6:end)];
        tree{j,i}=[str1,newBL,str2];
    end
   
   fprintf(fid, [tree{j,i}]);
   end 
   end

for i = 1:foldersize
    for j = 1:numtrees
    fprintf(fid, ['[PARTITIONS] ' [treename{j,i}(1:end-1) 'ns'] ' [' [treename{j,i}] 'neutralmodel 300] \n']);
        fprintf(fid, ['[PARTITIONS] ' [treename{j,i}(1:end-1) 'ps'] ' [' [treename{j,i}] 'oneratio 300] \n']);
    end
end
fprintf(fid,'[EVOLVE] \n')
for i = 1:foldersize
    for j = 1:numtrees
    fprintf(fid,[[treename{j,i}(1:end-1) 'ns'] ' 1 ' [treename{j,i}(1:end-1) 'ns']  '\n']);
    fprintf(fid,[[treename{j,i}(1:end-1) 'ps'] ' 1 ' [treename{j,i}(1:end-1) 'or']  '\n']);
    end
end
fclose(fid)

%%
%this parses tree data output with recombination and outputs control.txt
clear
directory = 'C:\Users\dswle_000\Documents\massprfmss\macprf_mss\reveiwers\simulations\msdir\';
suffix='*RE.out'
direc= dir([directory, suffix]);
fileNames= {};
[fileNames{1:length(direc),1}]=deal(direc.name); 
foldersize = length(fileNames);
fid=fopen('C:\Users\dswle_000\Documents\massprfmss\macprf_mss\reveiwers\simulations\INDELibleV1.03\bin\control.txt','w')'
fprintf(fid,'[TYPE] CODON 1 \n')
fprintf(fid, '[MODEL] neutralmodel \n [submodel] 4.5 0.5 0 1 \n')
fprintf(fid, '[MODEL] oneratiomodel \n [submodel] 4.5 0.5 \n')
treelist={}
scaleparam=100;
for i = 1:foldersize
    FILENAME=[fileNames{i,1}];
    text=fileread([directory FILENAME]);
    slashsites=find(text=='/');
    startsites=slashsites(find(text(slashsites+1)=='/'));
    semisites=find(text==';');
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
    lastadd=0;
    for k = 1:length(subtreesstops)
        if text(subtreestarts(k)-1)==']'
       	basenum=str2num([text(subtreestarts(k)-4:subtreestarts(k)-2)]);
        truestart=subtreestarts(k);
        elseif text(subtreestarts(k)-2)==']'
            basenum=str2num([text(subtreestarts(k)-4:subtreestarts(k)-3)]);
            truestart=subtreestarts(k)-1;
        elseif text(subtreestarts(k)-3)==']'
            basenum=str2num([text(subtreestarts(k)-4)]);
        truestart=subtreestarts(k)-2;
        end
        if mod(basenum+lastadd,3)>0
        nextadd=3-mod(basenum+lastadd,3);
               else 
            nextadd=0;                 
        end
        codon_num=(nextadd+basenum+lastadd)/3; 
subtreename{j,k}=[FILENAME(1:end-4) '_r' num2str(j) '_s' num2str(k) ' '];
    %if text(subtreestarts(k)-5)=='['
   subtree{j,k}=['[TREE] ' [subtreename{j,k}] text(truestart:subtreesstops(k)) '\n'];
   % elseif text(subtreestarts(k)-4)=='['
    %    subtree{j,k}=['[TREE] ' [subtreename{j,k}] text(subtreestarts(k)-1:subtreesstops(k)) '\n'];
    %elseif text(subtreestarts(k)-3)=='['
    %        subtree{j,k}=['[TREE] ' [subtreename{j,k}]
    %        text(subtreestarts(k)-2:subtreesstops(k)) '\n']; ignore this
    %        stuff because it's accounted for on top
       colonlocs=find([subtree{j,k}]==':');
   for m = 1:length(colonlocs)
        colonlocs=find([subtree{j,k}]==':');
   periodlocs=find([subtree{j,k}]=='.');
   endparlocs=find([subtree{j,k}]==')');
   commalocs=find([subtree{j,k}]==',');
       branchlengthstr=subtree{j,k}(colonlocs(m)+1:colonlocs(m)+5);
        newBL=sprintf('%1.5f',(str2num(branchlengthstr)/scaleparam));
        str1=[subtree{j,k}(1:colonlocs(m))];
        str2=[subtree{j,k}(colonlocs(m)+6:end)];
        subtree{j,k}=[str1,newBL,str2];
   end
   treelist=[treelist subtreename{j,k}];
   fprintf(fid, [subtree{j,k}])
    fprintf(fid, ['[PARTITIONS] ' [subtreename{j,k}(1:end-1) '_ns'] ' [' [subtreename{j,k}] 'neutralmodel ' num2str([codon_num]) '] \n']);
        fprintf(fid, ['[PARTITIONS] ' [subtreename{j,k}(1:end-1) '_or'] ' [' [subtreename{j,k}] 'oneratiomodel ' num2str([codon_num]) '] \n']);
    lastadd=nextadd;
    end
    end
end
fprintf(fid,'[EVOLVE] \n')
for i = 1:length(treelist)
    fprintf(fid,[[treelist{1,i}(1:end-1) '_ns'] ' 1 ' [treelist{1,i}(1:end-1) '_ns']  '\n']);
       fprintf(fid,[[treelist{1,i}(1:end-1) '_or'] ' 1 ' [treelist{1,i}(1:end-1) '_or']  '\n']);
end
fclose(fid)

