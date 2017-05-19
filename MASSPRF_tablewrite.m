clear
indirectory = 'C:\Users\dswle_000\Documents\MASSPRF_genes15_Aug\MASSPRF_toalign\results\';
suffix='*out.txt'
direc= dir([indirectory, suffix]);
outdirectory=[indirectory 'output_tables\']; %use this to change output directory 
%outdirectory='C:\Users\dswle_000\Documents\massprfmss\macprf_mss\reveiwers\simulations\simulation_results\Tables\'
fileNames= {};
[fileNames{1:length(direc),1}]=deal(direc.name); 
foldersize = length(fileNames);
%%  
%this section separates the table from the rest of the MASSPRF output file.
for k = 1:foldersize
    clear inputfile
    clear pared
    FILENAME=[indirectory fileNames{k,1}];
inputfile=fileread(FILENAME);
genename=fileNames{k,1}(1:end-4);
    Plocs=find(inputfile=='P');
    badfile(k)=1;
    for i=1:length(Plocs)
        if inputfile(Plocs(i)+1)=='o'  
        if inputfile(Plocs(i)+2)=='s'
        if inputfile(Plocs(i)+3)=='i'
        if inputfile(Plocs(i)+4)=='t'
        if inputfile(Plocs(i)+5)=='i'
                if inputfile(Plocs(i)+6)=='o'
        pared=inputfile(Plocs(i):end);
               badfile(k)=0;
                    
                end
        end
        end
                        end
                end
        end
    end
    if badfile(k)==0
     % plocations= find(inputfile=='p')
 %       for m = 1:length(plocations)-1
%if inputfile(plocations(m)+1)=='o'
%if inputfile(plocations(m)+2)=='l'
%if inputfile(plocations(m)+3)=='y'
%if inputfile(plocations(m)+4)=='.'
%if inputfile(plocations(m)+5)=='t'
%if inputfile(plocations(m)+6)=='x'
%if inputfile(plocations(m)+7)=='t'    
   % if inputfile(plocations(m)-2)=='_'
   %     scalefactor(k)=str2num(inputfile(plocations(m)-1));
%else
       % scalefactor(k)=str2num([inputfile(plocations(m)-2) inputfile(plocations(m)-1)]);
  %  end
%end
%end
%end
%end
%end
%end
       % end
     %   end
    Alocs=find(pared=='A');
        for i = 1:length(Alocs)
        if pared(Alocs(i)+1)=='b'
                    fintab=pared(1:Alocs(i)-1);
                    startparlocs=find(inputfile=='(');
                    endparlocs=find(inputfile==')');
time=inputfile(startparlocs(end):endparlocs(end));
colonlocs=find(time==':');
if length(colonlocs)==2
mins=str2num(time(colonlocs(1)+1:colonlocs(2)-1));
secs=str2num(time(colonlocs(2)+1:end-1));
totaltime(k,1) = 60*mins+secs;
elseif length(colonlocs)==3
hrs=str2num(time(colonlocs(1)+1:colonlocs(2)-1));
mins=str2num(time(colonlocs(2)+1:colonlocs(3)-1));
secs=str2num(time(colonlocs(3)+1:end-1));
totaltime(k,1)=3600*hrs+60*mins+secs;
        end
        end
           end
      tablefid=fopen([outdirectory genename '_Table.txt'],'w');
       fprintf(tablefid,fintab);
        fclose(tablefid);
        
     
 
    
    end
end
    %save([outdirectory 'scalefactors.dat'],'scalefactor','-ascii');
    %save([outdirectory 'totaltime.dat'],'totaltime','-ascii');
    save('badfile.dat','badfile','-ascii')