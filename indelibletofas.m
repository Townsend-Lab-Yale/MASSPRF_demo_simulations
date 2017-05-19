clear
directory = 'C:\Users\dswle_000\Documents\massprfmss\macprf_mss\reveiwers\simulations\INDELibleV1.03\bin\';
outdirectory='C:\Users\dswle_000\Documents\massprfmss\macprf_mss\reveiwers\simulations\Simulation_sequences_1\';
genelength=900;
suffix='*.fas'
direc= dir([directory, suffix]);
fileNames= {};
[fileNames{1:length(direc),1}]=deal(direc.name); 
foldersize = length(fileNames);
%This parses for tree data output by Hudson's ms without recombination and writes it into a
%control file that can be used as an input for indelible. 
for i = 1:foldersize
    FILENAME=[fileNames{i,1}];
    text=fileread([directory FILENAME]);
    text(text==' ')='';
   fiddiv=fopen([outdirectory FILENAME(1:end-4) 'div.fas'],'w');
fidpol=fopen([outdirectory FILENAME(1:end-4) 'pol.fas'],'w');
    startsites=find(text=='>');
%    divseqname=[text((slashsites(1)):(slashsites(1)+1))];
divseqname=['>' FILENAME(1:end-4) '_div'] 
divseq=[text((startsites(1)):(startsites(1)+genelength+3))];
fprintf(fiddiv,[divseq]);
fprintf(fiddiv,'\n')
fclose(fiddiv)
for j=2:length(startsites)-1
%polseqname=[text((slashsites(j)):(slashsites(j)+1))];
%polseqname=['>' FILENAME(1:end-4) '_pol' text((startsites(j)+1):(startsites(j+1)-4))];
polseq=[text((startsites(j)):(startsites(j+1)-1))];
fprintf(fidpol,[polseq]);
end
fprintf(fidpol, [text(startsites(end):end) '\n']);
fclose(fidpol)
end
%%
%this half of the script writes a .pbs file for use on the Yale cluster
clear
directory='C:\Users\dswle_000\Documents\massprfmss\macprf_mss\reveiwers\simulations\Simulation_sequences\';
divsuffix='*psdiv.fas';
divdirec= dir([directory, divsuffix]);
[divfileNames{1:length(divdirec),1}]=deal(divdirec.name); 
polsuffix='*pspol.fas';
poldirec= dir([directory, polsuffix]);
[polfileNames{1:length(poldirec),1}]=deal(poldirec.name); 
foldersize = length(divfileNames);
%This parses for tree data output by Hudson's ms without recombination and writes it into a
%control file that can be used as an input for indelible. 
fid=fopen([directory 'simtestpos.pbs'],'w');
fprintf(fid,'#!/bin/bash \n#PBS -N MASSPRF_simtest \n#PBS -q fas_high \n#PBS -l nodes=1:ppn=8 \n#PBS -l walltime=00:23:59:59 \n#PBS -k oe \n#PBS -m abe \n')
fprintf(fid, 'echo ''***\\ndate******\\n'' \nWORKDIR=~/MASS-PRF-v1.4/MASS-PRF_14June2016/bin/ \ncd $WORKDIR \n');
for i = 1:foldersize/4
string=['./massprf -p ' [polfileNames{i}] ' -o 1 -d ' [divfileNames{i}] ' -ic 0 -ci_m 1 -s 1 -r 1 -ci_r 1 -t 6 -NI 1 > ' [divfileNames{i}(1:end-7)] '_MASSPRFv1.4.txt \n']   
fprintf(fid,string)
end
fprintf(fid, 'echo ''***\\ndate******\\n''')
fclose(fid)
%fidgrace=fopen([directory 'simtest.sh'],'w');
%fprintf(fidgrace,'#BSUB -J MASSPRF_yeast\n#BSUB -n 16\n#BSUB -q week\n#BSUB -W 167:59\n#BSUB -B')

    
