#!/usr/bin/python

"""
  __    __     ______     ______     ______     __         ______     ______
 /\ "-./  \   /\  __ \   /\  ___\   /\  ___\   /\ \       /\  __ \   /\  == \
 \ \ \-./\ \  \ \ \/\ \  \ \___  \  \ \___  \  \ \ \____  \ \  __ \  \ \  __<
  \ \_\ \ \_\  \ \_____\  \/\_____\  \/\_____\  \ \_____\  \ \_\ \_\  \ \_____\
   \/_/  \/_/   \/_____/   \/_____/   \/_____/   \/_____/   \/_/\/_/   \/_____/

                                   COBRETTI
               "Does the job nobody wants" -Cobra (1986) trailer

                Contact: Snake Peterson - jakemp@iastate.edu

Automation of redundant RNA bioinformatic tasks, saving valuable time and resources

Required modules for Pronto/HPC users:
module load py-biopython
module load python

Other modules will be added during shell script builds, see 'shell_build_start' function to verify
Program locations will also need to be verified prior to the first run

Use:
python cobretti.py -stage # -email email@iastate.edu
Valid stage numbers: 1A, 1B, 1C, 2A, 2B, 2C, 2D, 3A, 3B

STAGE 1: SECONDARY STRUCTURE PREDICTION
Place fasta ('.fa') files of interest in current working directory, or provide list of accession numbers and add "-i list.txt" to command line
    1A: Finished when all ScanFold runs ('scanfold#.sh') and BLAST searches ('blast#.sh') are complete
    1B: Run once stage 1A complete
        Finished when all cm-builder runs ('cmbuilder#.sh') are complete
    1C: Run once stage 1B complete
        Finished when R-Scape run ('rscape.sh','All.Rscape.pdf') is complete and working directory is cleaned up

STAGE 2: TERTIARY STRUCTURE AND POCKET PREDICTION
Copy motifs ('.dbn') of interest into current working directory
    2A: Finished when all SimRNA runs ('simrna#.sh') are complete
    2B: Run once stage 2A complete
        Finished when QRNAS refinements ('qrnas#.sh') are complete
    2C: Run once stage 2B complete
        Finished when ARES run ('ares.sh','ares.csv') and fpocket run is complete

STAGE 3: VIRTUAL SCREENING - WORK IN PROGRESS
Copy structures ('.pdb') of interest into working directory
    3A: Finished when all DOCK 6 ru#!/usr/bin/pythonns ('dock#.sh') are complete
    3B: Run once stage 3B complete
        Finished when AnnapuRNA run ('annapurna.txt') is complete

Substages are also available for troubleshooting, and custom stages can be built as needed.

"""
import os,sys,argparse,shutil,random,subprocess



"""FUNCTIONS"""
def folder_make(directory,subdirectory): #Check if a subdirectory exists; if not, make directory
    if os.path.exists(os.path.join(directory,subdirectory)) == False:
        os.mkdir(os.path.join(directory,subdirectory))

def folder_check(directory,subdirectory): # Check if subdirectory exists; if so, return as path; if not, return directory as path
    if os.path.exists(os.path.join(directory,subdirectory)) == True:
        return(os.path.join(directory,subdirectory))
    else:
        return(directory)

def shell_build_start(filename,job,email,time=3,nodes=1,mem=0,tasks=1,notify='ALL'):
    with open(filename,'w',newline='\n') as writefile:
        writefile.writelines('#!/bin/bash -l\n')
        writefile.writelines('#SBATCH --partition=biocrunch\n') #Partition to submit to
        writefile.writelines('#SBATCH --time='+str(time)+'-00:00:00\n') #Time limit for this job
        writefile.writelines('#SBATCH --nodes='+str(nodes)+'\n') #Nodes to be used for this job during runtime
        if mem != 0:
            writefile.writelines('#SBATCH --mem='+str(mem)+'G\n') #Memory allocation may or may not be needed based on motif and database size
        writefile.writelines('#SBATCH --ntasks-per-node='+str(tasks)+'\n')
        writefile.writelines('#SBATCH --job-name='+job+'\n') #Name of this job in work queue
        writefile.writelines('#SBATCH --mail-user='+email+'\n') #Email to send notifications to
        writefile.writelines('#SBATCH --mail-type='+notify+'\n') #Email notification type (BEGIN, END, FAIL, ALL)
        writefile.writelines('\n')
        if job.startswith('scanfold'):
            writefile.writelines('module load py-biopython/1.70-py3-wos466g\n')
            writefile.writelines('module load python/3.6.5-fwk5uaj\n')
        elif job.startswith('blast'):
            writefile.writelines('module load ncbi-blast-db/latest\n')
            writefile.writelines('echo ${NCBI_BLAST_DB_PATH}: /work/LAS/BioDatabase/ncbi/blast-db/latest\n')
            writefile.writelines('module load blast-plus/2.12.0\n')
            writefile.writelines('export BLASTDB=${NCBI_BLAST_DB_PATH}\n')
        elif job.startswith('cmbuilder'):
            writefile.writelines('module load infernal\n')
        elif job.startswith('rscape'):
            writefile.writelines('module load infernal\n')
            writefile.writelines('module load gcc/7.3.0-xegsmw4\n')
            writefile.writelines('module load ghostscript\n')
            writefile.writelines('module load gnuplot\n')
            writefile.writelines('module load py-biopython\n')
            writefile.writelines('module load python\n')
        elif job.startswith('qrnas'):
            writefile.writelines('module load gcc\n')
        elif job.startswith('ares'):
            writefile.writelines('module load miniconda3\n')
        else:
            writefile.writelines('module load py-biopython\n')
            writefile.writelines('module load python\n')
        writefile.writelines('\n')


def fasta_build(filename,sequence_directory,email): #read a text file of accession numbers, download associated fasta files
    with open(filename,'r') as readfile:
        lines = readfile.readlines()
        if os.path.exists(sequence_directory) == False:
            os.mkdir(sequence_directory)
        for line in lines:
            with open(os.path.join(sequence_directory,line+'.fa'),'w',newline='\n') as writefile:
                writefile.writelines('>'+line+'\n')
                writefile.writelines(fasta_find(line,email).replace('T','U')+'\n')

def fasta_find(accession_number,email): #Download fasta sequences via Entrez
    from Bio import Entrez
    Entrez.email = email #Required for Entrez
    Entrez_search = Entrez.efetch(db='nucleotide',id=accession_number,rettype='fasta')
    Entrez_record = SeqIO.read(Entrez_search,'fasta')
    return(Entrez_record.seq)

def ScanFold_prep(sequence_directory,scanfold_directory,email): #Makes shell scripts for ScanFold runs
    all_files = os.listdir(sequence_directory)
    count = 1
    for filename in all_files:
        if filename.endswith('.fa') or filename.endswith('.fasta'):
            with open(os.path.join(sequence_directory,filename),'r') as readfile:
                lines = readfile.readlines()
                name = lines[0].rstrip().replace('.','-').strip('>').split(' ')[0]
            shell_build_start('scanfold'+str(count)+'.sh','scanfold_'+name,email,notify='END,FAIL')
            with open('scanfold'+str(count)+'.sh','a',newline='\n') as writefile:
                writefile.writelines('python %s %s --name %s --global_refold &\n' % (scanfold_directory,os.path.join(sequence_directory,filename),name))
                writefile.writelines('wait;\n')
            count += 1

def ScanFold_run(): #Run ScanFold on all fasta files
    current_directory = os.getcwd()
    all_files = os.listdir(current_directory)
    for filename in all_files:
        if filename.startswith('scanfold') and filename.endswith('.sh'):
            os.system('sbatch '+filename)

def blast_prep(sequence_directory,db_directory,email): #BLAST search using local databases, output to database directory 
    from Bio import SeqIO
    all_files = os.listdir(sequence_directory)
    count = 1
    for filename in all_files:
        if filename.endswith('.fa') or filename.endswith('.fasta'):
            for record in SeqIO.parse(os.path.join(sequence_directory,filename),'fasta'):
                name = filename.split('.')[0]
                blast_tmp = name+'_db_raw.txt'
                seq_file = os.path.join(sequence_directory,filename)
                tmp_file = os.path.join(db_directory,blast_tmp)
                shell_build_start('blast'+str(count)+'.sh','blast_'+name,email,tasks=8,notify='END,FAIL')
                with open('blast'+str(count)+'.sh','a',newline='\n') as writefile:
                    writefile.writelines('blastn -task blastn -db nt -query %s  -out %s -outfmt "6 sacc sseq" -max_target_seqs 2500 -num_threads 8 -max_hsps 1;\n' % (seq_file,tmp_file))
                    writefile.writelines('wait;')
                count += 1

def blast_run():
    current_directory = os.getcwd()
    all_files = os.listdir(current_directory)
    for filename in all_files:
        if filename.startswith('blast') and filename.endswith('.sh'):
            os.system('sbatch '+filename)

def blast_cleanup(db_directory,dbn_directory):
    current_directory = os.getcwd()
    all_files = os.listdir(current_directory)
    sh_directory = os.path.join(current_directory,'sh')
    if os.path.exists(sh_directory) == False:
        os.mkdir(sh_directory)
    out_directory = os.path.join(current_directory,'out')
    if os.path.exists(out_directory) == False:
        os.mkdir(out_directory)
    for filename in all_files:
        if (filename.startswith('scanfold') or filename.startswith('blast')) and filename.endswith('.sh'):
            try:
                shutil.move(os.path.join(current_directory,filename),os.path.join(sh_directory,filename))
            except:
                pass
        elif filename.startswith('slurm') and filename.endswith('.out'):
            try:
                shutil.move(os.path.join(current_directory,filename),os.path.join(out_directory,filename))
            except:
                pass
    db_files = os.listdir(db_directory)
    for filename in db_files:
        if filename.endswith('_db_raw.txt'):
            name = filename.split('_db_raw.txt')[0]
            db_filename = name+'_db.fa'
            with open(os.path.join(db_directory,db_filename),'w',newline='\n') as writefile:
                with open(os.path.join(db_directory,filename),'r') as readfile:
                    lines = readfile.readlines()
                    for line in lines:
                        accession = line.rstrip().split('\t')[0]
                        sequence = line.rstrip().split('\t')[1]
                        writefile.writelines('>'+accession+'\n')
                        writefile.writelines(sequence.replace('T','U')+'\n')
#            if os.path.isfile(os.path.join(db_directory,db_filename)) == True:
#                os.remove(os.path.join(db_directory,filename))
#                os.remove(os.path.join(db_directory,single_line_filename))
    current_directory = os.getcwd()
    all_files = os.listdir(current_directory)
    for root, dirs, files in os.walk(current_directory):
        for name in files:
            if root.endswith('motifs'):
                pass
            else:
                try:
                    if name.endswith('.dbn') and name.split('_')[-2] == 'motif':
                        shutil.copy(os.path.join(root, name),os.path.join(dbn_directory,name))
                except:
                    pass

def extend_motifs(sequence_directory,dbn_directory,dbn_writefile='extended.dbn'):
    open(dbn_writefile,'w',newline='\n')
    all_files = os.listdir(dbn_directory)
    for filename in all_files:
        try:
            if filename.endswith('.dbn') and filename.split('_')[-2] == 'motif':
                dbn_extend(dbn_directory,filename,sequence_directory,dbn_writefile)
        except:
            continue

def dbn_extend(dbn_directory,dbn_readfile,sequence_directory,dbn_writefile): #Extends ScanFold DBN motifs by 30 nts on each side
    from Bio import SeqIO
    with open(os.path.join(dbn_directory,dbn_readfile),'r') as readfile:
        dbn_lines = readfile.readlines()
        tag = dbn_readfile.split('_')[0]
        dbn_header = str(dbn_lines[0]).rstrip()
        dbn_sequence  = str(dbn_lines[1]).rstrip()
        dbn_structure = str(dbn_lines[2]).rstrip()
    all_files = os.listdir(sequence_directory)
    for filename in all_files:
        if filename.startswith(tag) and (filename.endswith('.fa') or filename.endswith('.fasta')):
            for record in SeqIO.parse(os.path.join(sequence_directory,filename),'fasta'):
                fasta_sequence = str(record.seq).replace('T','U')
                seq_start = fasta_sequence.find(dbn_sequence)
                seq_length = len(dbn_sequence)
                seq_end = seq_start + seq_length
                if seq_start < 30:
                    dbn_start = 0
                    five_prime_filler = seq_start
                else:
                    dbn_start = seq_start - 30
                    five_prime_filler = 30
                if seq_end + 30 > len(fasta_sequence):
                    dbn_end = len(fasta_sequence)
                    three_prime_filler = dbn_end - seq_end
                else:
                    dbn_end = seq_end + 30
                    three_prime_filler = 30
    with open(dbn_writefile,'a+',newline='\n') as writefile:
        writefile.writelines(dbn_header+'\n')
        writefile.writelines(fasta_sequence[dbn_start:dbn_end]+'\n')
        writefile.writelines('.'*five_prime_filler+dbn_structure+'.'*three_prime_filler+'\n')

def pk_fold(knotty_program,hfold_program,dbn_readfile='extended.dbn',pk_writefile='tmppk.txt'): #Fold dbn motifs and output to single textfile
    with open(pk_writefile,'w',newline='\n') as writefile, open(dbn_readfile,'r') as readfile:
        i = 0
        lines = readfile.readlines()
        for line_index, line in enumerate(lines):
            if (line_index + 3) >= len(lines):
                break
            #print("first character in line: " + str(line[0]))
            if line[0] =='>':
                dbn_header = str(lines[i]).rstrip()
                dbn_sequence  = str(lines[i+1]).rstrip()
                dbn_structure = str(lines[i+2]).rstrip()
                left_pos = 0
                right_pos = len(dbn_sequence)
                nucleotides = ('A','C','G','U')
                left_brackets = ('(','[','{','<')
                right_brackets = (')',']','}','>')
                nucleotide_codes = {
                    "N": ('A','C','G','U'),
                    "B": ('C','G','U'),
                    "D": ('A','G','U'),
                    "H": ('A','C','U'),
                    "V": ('A','C','G'),
                    "K": ('G','U'),
                    "M": ('A','C'),
                    "R": ('A','G'),
                    "S": ('C','G'),
                    "W": ('A','U'),
                    "Y": ('C','U')
                                   }
                #print("dbn_sequence:  " + dbn_sequence)
                #print("dbn_structure: " + dbn_structure)
                if dbn_sequence[0] == '>' or dbn_structure[0] == '>':
                    continue
                if not dbn_sequence[0].isalnum():
                    continue
                if not (dbn_structure[0] == '.' or dbn_structure[0] in left_brackets or dbn_structure[0] in right_brackets):
                    continue
                for k,j in enumerate(dbn_sequence):
                    #print("index: " + str(k))
                    #print("item:  " + j)
                    if j in nucleotides or dbn_structure[k] in left_brackets:
                        left_pos = k
                        print('k left = '+str(left_pos))
                        break

                for k,j in enumerate(dbn_sequence[::-1]):
                    if j in nucleotides or dbn_structure[k] in right_brackets:
                        right_pos -= k
                        print('k right = '+str(right_pos))
                        break
                dbn_sequence = dbn_sequence[left_pos:right_pos]
                dbn_structure = dbn_structure[left_pos:right_pos]
                print("final left position: "+str(left_pos))
                print("final right position: "+str(right_pos))
                print("dbn_sequence:  " + dbn_sequence)
                print("dbn_structure: " + dbn_structure)
                new_sequence = ""
                for k in dbn_sequence:
                    if k in nucleotide_codes:
                        new_sequence += random.choice(nucleotide_codes[k])

                    else:
                        new_sequence += k
                dbn_sequence = new_sequence

                writefile.writelines(dbn_header+' ScanFold\n')
                writefile.writelines(dbn_sequence+'\n')
                writefile.writelines(dbn_structure+'\n')
                writefile.writelines(dbn_header+' Knotty\n')
                #writefile.writelines(dbn_sequence+'\n') #Knotty writes own sequence line
                knotty_input = str('%s %s' % (knotty_program,dbn_sequence))
                #print(knotty_input)
                writefile.flush()
                knotty = subprocess.run(knotty_input,stdout=writefile,shell=True)
                writefile.writelines(dbn_header+' HFold\n')
                #writefile.writelines(dbn_sequence+'\n') #HFold writes own sequence line
                blank_structure = str('_'*len(dbn_sequence)) #HFold requires '_' instead of '.' to fold properly
                hfold_blank_input = str('%s --s "%s" --r "%s"' % (hfold_program,dbn_sequence,blank_structure))
                #print(hfold_blank_input)
                writefile.flush()
                hfold_only = subprocess.run(hfold_blank_input,stdout=writefile,shell=True)                       
                writefile.writelines(dbn_header+' HFold ScanFold\n')
                # writefile.writelines(dbn_sequence+'\n') #HFold writes own sequence line
                hfold_input = str('%s --s "%s" --r "%s"' % (hfold_program,dbn_sequence,dbn_structure.replace('.','_'))) #HFold requires '_' instead of '.' to fold properly
                #print(hfold_input)
                writefile.flush()
                hfold = subprocess.run(hfold_input,stdout=writefile,shell=True)
                i += 1
            else:
                i += 1
                
                
    
def replace_nucleotide_codes_with_random_nucleotide(sequence, codes):
    new_sequence = ""
    for k in sequence:
        if k in codes:
            new_sequence += random.choice(codes[k])

        else:
            new_sequence += k
    return new_sequence
                

def pk_cleanup(pk_readfile='tmppk.txt',pk_writefile='pk_clean.txt'): #Clean up pseudoknot file for easier visual comparison and breakdown
    i = 0
    with open(pk_readfile,'r') as readfile:
        with open (pk_writefile,'w',newline='\n') as writefile:
            lines = readfile.readlines()
            for line in lines:
                if line.rstrip().endswith(' HFold ScanFold'):
                    writefile.writelines(lines[i].rstrip()+'\n')
                    try:
                        sequence = lines[i+5].rstrip().split(' ')
                        writefile.writelines(sequence[1]+'\n')
                    except:
                        writefile.writelines('HFold ScanFold Sequence ERROR\n')
                    try:
                        structure = lines[i+6].rstrip().split(' ')
                        writefile.writelines(structure[1]+'\n')
                    except:
                        writefile.writelines('HFold ScanFold Structure ERROR\n')
                    writefile.writelines('\n')
                    i += 1
                elif line.rstrip().endswith(' Knotty'):
                    writefile.writelines(lines[i].rstrip()+'\n')
                    try:
                        sequence = lines[i+1].rstrip().split(' ')
                        writefile.writelines(sequence[1]+'\n')
                    except:
                        writefile.writelines('Knotty Sequence ERROR\n')
                    try:
                        structure = lines[i+2].rstrip().split(' ')
                        writefile.writelines(structure[1]+'\n')
                    except:
                        writefile.writelines('Knotty Structure ERROR\n')
                    i += 1
                elif line.rstrip().endswith(' HFold'):
                    writefile.writelines(lines[i].rstrip()+'\n')
                    try:
                        sequence = lines[i+5].rstrip().split(' ')
                        writefile.writelines(sequence[1]+'\n')
                    except:
                        writefile.writelines('HFold Sequence ERROR\n')
                    try:
                        structure = lines[i+6].rstrip().split(' ')
                        writefile.writelines(structure[1]+'\n')
                    except:
                        writefile.writelines('HFold Structure ERROR\n')
                    i += 1
                elif line.rstrip().endswith(' ScanFold'):
                    writefile.writelines(lines[i].rstrip()+'\n')
                    try:
                        writefile.writelines(lines[i+1].rstrip()+'\n')
                    except:
                        writefile.writelines('ScanFold Sequence ERROR\n')
                    try:
                        writefile.writelines(lines[i+2].rstrip()+'\n')
                    except:
                        writefile.writelines('ScanFold Structure ERROR\n')
                    i += 1
                else:
                    i += 1

def pk_breakdown(pk_directory,pk_readfile='pk_clean.txt'): #Takes clean dbn file (pk_cleanup) and converts into separate dbn files - breaks pseudoknots into individual non-pseudoknotted motifs
    i = 0
    with open(pk_readfile,'r') as readfile:
        lines = readfile.readlines()
        for line in lines:
            if line[0] == '>':
                dbn_header = lines[i].rstrip()
                if lines[i+1].endswith('ERROR') or lines[i+2].endswith('ERROR'):
                    pk_filename = dbn_header.strip('>').split('_coordinates')
                    pk_header = dbn_header
                    pk_filename = pk_filename[0]+'_error.dbn'
                elif dbn_header.endswith(' HFold ScanFold'):
                    pk_filename = dbn_header.strip('>').split('_coordinates')
                    pk_header = '>'+pk_filename[0]+' HFold ScanFold'
                    pk_filename = pk_filename[0]+'_hfold_scanfold.dbn'
                elif dbn_header.endswith(' Knotty'):
                    pk_filename = dbn_header.strip('>').split('_coordinates')
                    pk_header = '>'+pk_filename[0]+' Knotty'
                    pk_filename = pk_filename[0]+'_knotty.dbn'
                elif dbn_header.endswith(' HFold'):
                    pk_filename = dbn_header.strip('>').split('_coordinates')
                    pk_header = '>'+pk_filename[0]+' HFold'
                    pk_filename = pk_filename[0]+'_hfold.dbn'
                elif dbn_header.endswith(' ScanFold'):
                    pk_filename = dbn_header.strip('>').split('_coordinates')
                    pk_header = '>'+pk_filename[0]+' ScanFold'
                    pk_filename = pk_filename[0]+'_scanfold.dbn'
                pk_sequence = lines[i+1].rstrip()
                pk_structure = lines[i+2].rstrip()
                if pk_sequence.endswith('ERROR') or pk_sequence.endswith('source') or pk_structure.endswith('ERROR') or pk_structure.endswith('HFold'):
                    with open(os.path.join(pk_directory,pk_filename),'w',newline='\n') as writefile:
                        writefile.writelines(pk_header+'\n')
                        writefile.writelines(pk_sequence+'\n')
                        writefile.writelines(pk_structure+'\n') 
                    i += 1
                elif pk_istrue(pk_structure) == False:
                    pk_structure = pk_structure.replace('[','(').replace(']',')').replace('{','(').replace('}',')').replace('<','(').replace('>',')') #Pseudoknot software may use multiple brackets when no pseudoknot is present, this defaults all non-pseudoknot structures to open/close parentheses
                    with open(os.path.join(pk_directory,pk_filename),'w',newline='\n') as writefile:
                        left_pos = 0
                        right_pos = len(pk_sequence)
                        left_brackets = ('(','[','{','<')
                        right_brackets = (')',']','}','>')
                        while True:
                            if left_pos == right_pos:
                                left_pos = 0
                                break
                            elif pk_structure[left_pos] in left_brackets:
                                break
                            else:
                                left_pos += 1
                        while True:
                            if left_pos == right_pos:
                                right_pos = len(pk_sequence)
                                break
                            elif pk_structure[right_pos - 1] in right_brackets:
                                break
                            else:
                                right_pos -= 1
                        writefile.writelines(pk_header+'\n')
                        writefile.writelines(pk_sequence[left_pos:right_pos]+'\n')
                        writefile.writelines(pk_structure[left_pos:right_pos]+'\n') 
                    i += 1
                else:
                    pk_splitter(pk_directory,pk_filename,pk_header,pk_sequence,pk_structure)
                    i += 1
            else:
                i += 1

def pk_istrue(structure): #Will check for presence of pseudoknots, return either True or False - able to search 4 levels deep (),[],{},<> which is assumed to be enough depth for <200 nt structures
    bracket1 = structure.count('(')
    bracket2 = structure.count('[')
    bracket3 = structure.count('{')
    bracket4 = structure.count('<')
    if bracket1 >= 0 and bracket2 == 0 and bracket3 == 0 and bracket4 == 0: #Any sequence with only one bracket type can't be a pseudoknot structure
        return(False)
    elif bracket1 == 0 and bracket2 > 0 and bracket3 == 0 and bracket4 == 0:
        return(False)
    elif bracket1 == 0 and bracket2 == 0 and bracket3 > 0 and bracket4 == 0:
        return(False)
    elif bracket1 == 0 and bracket2 == 0 and bracket3 == 0 and bracket4 > 0:
        return(False)
    else:
        ijcoords = {}
        icoords = []
        jcoords = []
        klcoords = {}
        kcoords = []
        lcoords = []
        mncoords = {}
        mcoords = []
        ncoords = []
        opcoords = {}
        ocoords = []
        pcoords = []
        column = 0
        for character in structure: #Build lists of coordinates (i, j, etc.)
            if character == '(':
                icoords.append(column)
                column += 1
            elif character == ')':
                jcoords.append(column)
                column += 1
            elif character == '[':
                kcoords.append(column)
                column += 1
            elif character == ']':
                lcoords.append(column)
                column += 1
            elif character == '{':
                mcoords.append(column)
                column += 1
            elif character == '}':
                ncoords.append(column)
                column += 1
            elif character == '<':
                ocoords.append(column)
                column += 1
            elif character == '>':
                pcoords.append(column)
                column += 1
            else:
                column += 1
        for jcoord in jcoords: #Finds matching i,j coordinates by taking the left-most (lowest) j value and pairing it with the right-most (highest) i value that is less than j, then remove that value from the list
            for icoord in icoords[::-1]:
                if int(icoord) < int(jcoord):
                    ijcoords.update({icoord:jcoord})
                    icoords.remove(icoord)
                    break
        for lcoord in lcoords:
            for kcoord in kcoords[::-1]:
                if int(kcoord) < int(lcoord):
                    klcoords.update({kcoord:lcoord})
                    kcoords.remove(kcoord)
                    break
        for ncoord in ncoords:
            for mcoord in mcoords[::-1]:
                if int(mcoord) < int(ncoord):
                    mncoords.update({mcoord:ncoord})
                    mcoords.remove(mcoord)
                    break
        for pcoord in pcoords:
            for ocoord in ocoords[::-1]:
                if int(ocoord) < int(pcoord):
                    opcoords.update({ocoord:pcoord})
                    ocoords.remove(ocoord)
                    break
        for i,j in ijcoords.items(): #Checks for presence of pseudoknots in all possible base pair bracket configurations
            for k,l in klcoords.items():
                if i<k<j<l:
                    return(True)
                elif k<i<l<j:
                    return(True)
                else:
                    continue
            for m,n in mncoords.items():
                if i<m<j<n:
                    return(True)
                elif m<i<n<j:
                    return(True)
                else:
                    continue
            for o,p in opcoords.items():
                if i<o<j<p:
                    return(True)
                elif o<i<p<j:
                    return(True)
                else:
                    continue
        for k,l in klcoords.items():
            for m,n in mncoords.items():
                if k<m<l<n:
                    return(True)
                elif m<k<n<l:
                    return(True)
                else:
                    continue
            for o,p in opcoords.items():
                if k<o<l<p:
                    return(True)
                elif o<k<p<l:
                    return(True)
                else:
                    continue
        for m,n in mncoords.items():
            for o,p in opcoords.items():
                if m<o<n<p:
                    return(True)
                elif o<m<p<n:
                    return(True)
                else:
                    continue
        return(False)

def pk_splitter(dbn_directory,dbn_filename,dbn_header,dbn_sequence,dbn_structure): #Takes dbn data and breaks it down into multiple non-pseudoknot structures based on bracket notations
    bracket1 = dbn_structure.count('(')
    bracket2 = dbn_structure.count('[')
    bracket3 = dbn_structure.count('{')
    bracket4 = dbn_structure.count('<')
    count = 0
    dbn_filename = dbn_filename.rsplit('.dbn')
    dbn_filename = dbn_filename[0]
    if bracket1 > 0:
        count += 1
        dbn_filename1 = str(dbn_filename)+'_'+str(count)+'.dbn'
        dbn_structure1 = dbn_structure.replace('[','.').replace(']','.').replace('{','.').replace('}','.').replace('<','.').replace('>','.') #Defaults all extracted structures to open/close parentheses
        with open(os.path.join(dbn_directory,dbn_filename1),'w',newline='\n') as writefile:
            writefile.writelines(dbn_header+' '+str(count)+'\n')
            left_pos = 0
            right_pos = len(dbn_sequence)
            left_brackets = ('(')
            right_brackets = (')')
            while True:
                if left_pos == right_pos:
                    left_pos = 0
                    break
                elif dbn_structure1[left_pos] in left_brackets:
                    break
                else:
                    left_pos += 1
            while True:
                if left_pos == right_pos:
                    right_pos = len(dbn_sequence)
                    break
                elif dbn_structure1[right_pos - 1] in right_brackets:
                    break
                else:
                    right_pos -= 1
            dbn_sequence1 = dbn_sequence[left_pos:right_pos]
            dbn_structure1 = dbn_structure1[left_pos:right_pos]
            writefile.writelines(dbn_sequence1+'\n')
            writefile.writelines(dbn_structure1)
    if bracket2 > 0:
        count += 1
        dbn_filename2 = str(dbn_filename)+'_'+str(count)+'.dbn'
        dbn_structure2 = dbn_structure.replace('(','.').replace(')','.').replace('{','.').replace('}','.').replace('<','.').replace('>','.')
        dbn_structure2 = dbn_structure2.replace('[','(').replace(']',')')
        with open(os.path.join(dbn_directory,dbn_filename2),'w',newline='\n') as writefile:
            writefile.writelines(dbn_header+' '+str(count)+'\n')
            left_pos = 0
            right_pos = len(dbn_sequence)
            left_brackets = ('(')
            right_brackets = (')')
            while True:
                if left_pos == right_pos:
                    left_pos = 0
                    break
                elif dbn_structure2[left_pos] in left_brackets:
                    break
                else:
                    left_pos += 1
            while True:
                if left_pos == right_pos:
                    right_pos = len(dbn_sequence)
                    break
                elif dbn_structure2[right_pos - 1] in right_brackets:
                    break
                else:
                    right_pos -= 1
            dbn_sequence2 = dbn_sequence[left_pos:right_pos]
            dbn_structure2 = dbn_structure2[left_pos:right_pos]
            writefile.writelines(dbn_sequence2+'\n')
            writefile.writelines(dbn_structure2)
    if bracket3 > 0:
        count += 1
        dbn_filename3 = str(dbn_filename)+'_'+str(count)+'.dbn'
        dbn_structure3 = dbn_structure.replace('[','.').replace(']','.').replace('(','.').replace(')','.').replace('<','.').replace('>','.')
        dbn_structure3 = dbn_structure3.replace('{','(').replace('}',')')
        with open(os.path.join(dbn_directory,dbn_filename3),'w',newline='\n') as writefile:
            writefile.writelines(dbn_header+' '+str(count)+'\n')
            left_pos = 0
            right_pos = len(dbn_sequence)
            left_brackets = ('(')
            right_brackets = (')')
            while True:
                if left_pos == right_pos:
                    left_pos = 0
                    break
                elif dbn_structure3[left_pos] in left_brackets:
                    break
                else:
                    left_pos += 1
            while True:
                if left_pos == right_pos:
                    right_pos = len(dbn_sequence)
                    break
                elif dbn_structure3[right_pos - 1] in right_brackets:
                    break
                else:
                    right_pos -= 1
            dbn_sequence3 = dbn_sequence[left_pos:right_pos]
            dbn_structure3 = dbn_structure3[left_pos:right_pos]
            writefile.writelines(dbn_sequence3+'\n')
            writefile.writelines(dbn_structure3)
    if bracket4 > 0:
        count += 1
        dbn_filename4 = str(dbn_filename4)+'_'+str(count)+'.dbn'
        dbn_structure4 = dbn_structure.replace('[','.').replace(']','.').replace('{','.').replace('}','.').replace('(','.').replace(')','.')
        dbn_structure4 = dbn_structure4.replace('<','(').replace('>',')')
        with open(os.path.join(dbn_directory,dbn_filename4),'w',newline='\n') as writefile:
            writefile.writelines(dbn_header+' '+str(count)+'\n')
            left_pos = 0
            right_pos = len(dbn_sequence)
            left_brackets = ('(')
            right_brackets = (')')
            while True:
                if left_pos == right_pos:
                    left_pos = 0
                    break
                elif dbn_structure4[left_pos] in left_brackets:
                    break
                else:
                    left_pos += 1
            while True:
                if left_pos == right_pos:
                    right_pos = len(dbn_sequence)
                    break
                elif dbn_structure4[right_pos - 1] in right_brackets:
                    break
                else:
                    right_pos -= 1
            dbn_sequence4 = dbn_sequence[left_pos:right_pos]
            dbn_structure4 = dbn_structure4[left_pos:right_pos]
            writefile.writelines(dbn_sequence4+'\n')
            writefile.writelines(dbn_structure4)
            
def cmbuilder_prep(seq_directory,db_directory,dbn_directory,cmbuilder_program,perl_program,rnaframework_directory,email,rscape_program,cobretti_program,cm_writefile='cmbuilder'): #Reads all dbn files and outputs a shell script for cm-builder, based on Vans code
    count = 1
    current_size = 0
    max_size = 10
    max_db_size = 8000000000 #Files over 8GB will be run separately to aid with HPC OOM errors
    current_directory = os.getcwd()
    all_seq_files = os.listdir(seq_directory)
    all_db_files = os.listdir(db_directory)
    all_dbn_files = os.listdir(dbn_directory)
    for dbn_filename in all_dbn_files:
        if dbn_filename.endswith('_error.dbn'):
            print('Error with pseudoknot motif file(s), check DBN files in /pk_motifs directory and then run stage 1BC')
            sys.exit()
    for dbn_filename in all_dbn_files:
        if dbn_filename.endswith('.dbn'):
            gene = dbn_filename.split('_')[0]
            dbn_file = str(os.path.join(dbn_directory,dbn_filename))
            if current_size == 0:
                shell_build_start(cm_writefile+str(count)+'.sh',cm_writefile+str(count),email,mem=100,tasks=20,notify='END,FAIL')
                with open(cm_writefile+str(count)+'.sh','a',newline='\n') as writefile:
                    writefile.writelines('export PERL5LIB='+perl_program+'\n')
                    writefile.writelines('export PERL5LIB='+rnaframework_directory+'\n')
            with open(cm_writefile+str(count)+'.sh','a',newline='\n') as writefile:
                for seq_filename in all_seq_files:
                    if seq_filename.startswith(gene) and seq_filename.endswith('.fa'):
                        seq_file = str(os.path.join(seq_directory,seq_filename))
                        for db_filename in all_db_files:
                            if db_filename.startswith(gene) and db_filename.endswith('_db.fa'):
                                db_file = str(os.path.join(db_directory,db_filename))
                                db_size = int(os.path.getsize(db_file))
                                if db_size > max_db_size:
                                    oversize_db_count = 1
                                    while os.path.isfile(cm_writefile+'_'+gene+str(oversize_db_count)+'.sh'):
                                        oversize_db_count += 1
                                    shell_build_start(cm_writefile+'_'+gene+str(oversize_db_count)+'.sh',cm_writefile+'_'+gene+str(oversize_db_count),email,mem=100,tasks=20,notify='END,FAIL')
                                    with open(cm_writefile+'_'+gene+str(oversize_db_count)+'.sh','a',newline='\n') as oversize_writefile:
                                        oversize_writefile.writelines('perl %s -s %s -m %s -d %s -c 4 -T ./tmp%s -k -t 1 &\n' % (cmbuilder_program,seq_file,dbn_file,db_file,count))
                                        oversize_writefile.writelines('wait;\n')
                                else:
                                    writefile.writelines('perl %s -s %s -m %s -d %s -c 4 -T ./tmp%s -k -t 1 &\n' % (cmbuilder_program,seq_file,dbn_file,db_file,count))
                                    current_size += 1
            if current_size >= max_size:
                with open(cm_writefile+str(count)+'.sh','a',newline='\n') as writefile:
                    writefile.writelines('wait;\n')
                current_size = 0
                count += 1
    with open(cm_writefile+str(count)+'.sh','a',newline='\n') as writefile:
        writefile.writelines('wait;\n')
    shell_build_start('rscape.sh','rscape',email,time=1,notify='END,FAIL')
    with open('rscape.sh','a',newline='\n') as writefile:
        writefile.writelines('export PERL5LIB='+perl_program+'\n')
        writefile.writelines('export PERL5LIB='+rnaframework_directory+'\n')
        rscape_runs=10
        writefile.writelines('for f in *.stockholm; do '+rscape_program+' -s --ntree '+str(rscape_runs)+' $f; done\n')
        writefile.writelines('gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=All.Rscape.pdf -dBATCH *.R2R.sto.pdf\n')
        writefile.writelines('wait;\n')
        writefile.writelines('python '+cobretti_program+' -stage 1CA -email '+email+' &\n')
        writefile.writelines('wait;\n')

def cmbuilder_run(): #Console command to sbatch all the cmbuilder shell scripts in current directory
    current_directory = os.getcwd()
    all_files = os.listdir(current_directory)
    for filename in all_files:
        if filename.startswith('cmbuilder') and filename.endswith('.sh'):
            os.system('sbatch '+filename)

def cmbuilder_cleanup(): #Clean up working directory, run R-Scape
    current_directory = os.getcwd()
    folders = ['cm','dbn','out','sh']
    for folder in folders:
        folder_make(current_directory,folder)
    all_files = os.listdir(current_directory)
    for filename in all_files:
        for folder in folders:
            if filename != 'rscape.sh':
                try:
                    if os.path.isdir(os.path.join(current_directory,filename)):
                        if filename.startswith('tmp'):
                            shutil.rmtree(filename)
                    else:
                        if filename.endswith(folder):
                            shutil.move(os.path.join(current_directory,filename),os.path.join(current_directory,folder,filename))
                except:
                    pass
    os.system('sbatch rscape.sh') #Automatically starts cmbuilder_cleanup2 upon completion of shell script

def cmbuilder_cleanup2(cm_writefile='covariance.txt'): #Move remaining files into folders, write covariance power and output to file
    current_directory = os.getcwd()
    covariance_read(cm_writefile)
    current_directory = os.getcwd()
    folders = ['cm','cov','dbn','out','pdf','power','ps','sh','sto','stockholm','surv','svg']
    for folder in folders:
        folder_make(current_directory,folder)
    all_files = os.listdir(current_directory)
    for filename in all_files:
        for folder in folders:
            if not filename.endswith('Rscape.pdf'):
                if not os.path.isdir(os.path.join(current_directory,filename)):
                    if filename.endswith(folder):
                        try:
                            shutil.move(os.path.join(current_directory,filename),os.path.join(current_directory,folder,filename))
                        except:
                            pass

def covariance_read(cm_writefile): #Read power files and output all power results to single file - Adapted from Collin's script
    with open(cm_writefile,'w',newline='\n') as writefile:
        writefile.writelines(f'Name\tBPs\tavg_substitutions\texpected BPs\tSTD_DEV\tObserved BPs\t# of BP power 0-0.1\t# of BP power 0.1-0.25\t# of BP power >=0.25\t0-0.1 BP info\t0.1-0.25 BP info\t>=0.25 BP info\n')
        current_directory = os.getcwd()
        all_files = os.listdir(current_directory)
        for filename in all_files:
            if filename.endswith('.power'):
                with open(filename,'r') as readfile:
                    name = filename.replace("_1.power","")
                    low_cv = []
                    mid_cv = []
                    high_cv = []
                    lines = readfile.readlines()
                    for line in lines:
                        if len(line) <= 1:
                            pass
                        elif line.startswith('# BPAIRS observed to covary'):
                            obv_bp = line.rstrip().split(' ')[5]
                        elif line.startswith('# BPAIRS expected to covary'):
                            ex_bp = line.rstrip().split(' ')[5]
                            std_dev = line.rstrip().split(' ')[7]
                        elif line.startswith('# avg substitutions per BP'):
                            subs = line.rstrip().split(' ')[6]
                        elif line.startswith('# BPAIRS'):
                            bps =  line.rstrip().split(' ')[2]
                        elif line.startswith('     *'):
                            power = line.rsplit()
                            if float(power[-1]) >= 0.25:
                                high_cv.append((int(power[1]),int(power[2]),float(power[-1])))
                            elif float(power[-1]) >= 0.1:
                                mid_cv.append((int(power[1]),int(power[2]),float(power[-1])))
                            else:
                                low_cv.append((int(power[1]),int(power[2]),float(power[-1])))
                        else:
                            pass
                    low_num = len(low_cv)
                    mid_num = len(mid_cv)
                    high_num = len(high_cv)
                    writefile.writelines(f'{name}\t{bps}\t{subs}\t{ex_bp}\t{std_dev}\t{obv_bp}\t{low_num}\t{mid_num}\t{high_num}\t{low_cv}\t{mid_cv}\t{high_cv}\n')

def SimRNA_prep(SimRNA_directory,email,config_file='None',instances=8,replicas=10,simrna_writefile='simrna'): #Writes shell scripts for optimized SimRNA runs
    SimRNA_program = os.path.join(SimRNA_directory,'SimRNA')
    SimRNA_data = os.path.join(SimRNA_directory,'data')
    SimRNA_clustering = os.path.join(SimRNA_directory,'clustering')
    SimRNA_trafl2pdb = os.path.join(SimRNA_directory,'SimRNA_trafl2pdbs')
    if config_file == 'None':
        SimRNA_config = os.path.join(SimRNA_directory,'config.dat')
    else:
        SimRNA_config = config_file
    current_directory = os.getcwd()
    data_link = os.path.join(current_directory,'data')
    if os.path.exists(data_link) != True:
          os.symlink(SimRNA_data, data_link)
    file_count = 1
    all_files = os.listdir(current_directory)
    for filename in all_files:
        if filename.endswith('.dbn'):
            dbn_headers = os.path.splitext(filename)[0]
            sequence_tempfile = dbn_headers+'.sequence' #Split sequence and structure into separate inputs for SimRNA
            structure_tempfile = dbn_headers+'.structure'
            rmsd = [5.0,7.0,10.0]
            with open(filename,'r') as readfile:
                with open(sequence_tempfile,'w',newline='\n') as writefile:
                    line = readfile.readlines()
                    writefile.write(line[1])
                    sequence_length = len(line[1])
                    rmsd.append(float(sequence_length/10.0))
                    #print(rmsd)
                with open(structure_tempfile,'w',newline='\n') as writefile:
                    pk_count = line[2].count('[')+line[2].count('{')+line[2].count('<')
                    if pk_count == 0:
                        writefile.write(line[2])
                    else:
                        newline = line[2].replace('[','.').replace(']','.').replace('{','.').replace('}','.').replace('<','.').replace('>','.').rstrip()
                        writefile.write(newline+'\n')
                        if line[2].count('[') != 0:
                            newline = line[2].replace('(','.').replace(')','.').replace('{','.').replace('}','.').replace('<','.').replace('>','.').replace('[','(').replace(']',')').rstrip()
                            writefile.write(newline+'\n')
                        if line[2].count('{') != 0:
                            newline = line[2].replace('(','.').replace(')','.').replace('[','.').replace(']','.').replace('<','.').replace('>','.').replace('{','(').replace('}',')').rstrip()
                            writefile.write(newline+'\n')
                        if line[2].count('<') != 0:
                            newline = line[2].replace('(','.').replace(')','.').replace('{','.').replace('}','.').replace('[','.').replace(']','.').replace('<','(').replace('>',')').rstrip()
                            writefile.write(newline+'\n')
            shell_build_start(simrna_writefile+str(file_count)+'.sh',simrna_writefile+str(file_count),email,tasks=50,notify='END,FAIL')
            with open(simrna_writefile+str(file_count)+'.sh','a',newline='\n') as writefile:
                if instances == 1:
                    writefile.writelines('%s -c %s -E %s -s %s -S %s -o %s -R %s >& %s.log &\n'  % (SimRNA_program,SimRNA_config,replicas,sequence_tempfile,structure_tempfile,dbn_headers,random.randint(0,9999999999),dbn_headers))
                elif instances < 1:
                        print('error in number of instances: less than one')
                elif instances > 100:
                        print('error in number of instances: too big')
                else:
                    for instance_count in range(instances):
                        writefile.writelines('%s -c %s -E %s -s %s -S %s -o %s'  % (SimRNA_program,SimRNA_config,replicas,sequence_tempfile,structure_tempfile,dbn_headers))
                        if instance_count < 10:
                            writefile.writelines('_0'+str(instance_count)+' -R '+str(random.randint(0,9999999999))+' >& '+dbn_headers+'_0'+str(instance_count)+'.log &\n')
                        else:
                            writefile.writelines('_'+str(instance_count)+' -R '+str(random.randint(0,9999999999))+' >& '+dbn_headers+'_'+str(instance_count)+'.log &\n')
                writefile.writelines('wait;\n')
                if instances != 1:
                    writefile.writelines('cat '+dbn_headers+'_??')
                    if replicas != 1:
                            writefile.writelines('_??')
                    else:
                            writefile.writelines('_119')
                    writefile.writelines('.trafl > '+dbn_headers+'_all.trafl &\n')
                    writefile.writelines('wait;\n')
                    writefile.writelines(SimRNA_clustering+' '+dbn_headers+'_all.trafl 0.01 ')
                    for value in rmsd:
                        writefile.writelines(str(value)+' ')
                    writefile.writelines(' >& '+dbn_headers+'_clust.log &\n')
                    writefile.writelines('wait;\n')
                    for value in rmsd:
                        for rank in range(1,6):
                            writefile.writelines('%s %s_01_01-000001.pdb %s_all_thrs%s0A_clust0%s.trafl 1 AA &\n' % (SimRNA_trafl2pdb,dbn_headers,dbn_headers,value,rank))
                    writefile.writelines('wait;\n')
                else:
                    if replicas != 1:
                        writefile.writelines('cat '+dbn_headers+'_??.trafl > '+dbn_headers+'_all.trafl &\n')
                        writefile.writelines('wait;\n')
                    writefile.writelines(SimRNA_clustering+' '+dbn_headers)
                    if replicas != 1:
                        writefile.writelines('_all')
                    else:
                        writefile.writelines('_119')
                    writefile.writelines('.trafl 0.01 ')
                    for value in rmsd:
                        writefile.writelines(value+' ')
                    writefile.writelines(' >& '+dbn_headers+'_clust.log &\n')
                    writefile.writelines('wait;\n')
                    for value in rmsd:
                        for rank in range(1,6):
                            writefile.writelines('%s %s_119-000001.pdb %s_119_thrs%s0A_clust0%s.trafl 1 AA &\n' % (SimRNA_trafl2pdb,dbn_headers,dbn_headers,value,rank))
                    writefile.writelines('wait;\n')
            file_count += 1

def SimRNA_run():
    current_directory = os.getcwd()
    all_files = os.listdir(current_directory)
    for filename in all_files:
        if filename.startswith('simrna') and filename.endswith('.sh'):
            os.system('sbatch '+filename)

def SimRNA_cleanup():
    current_directory = os.getcwd()
    folders = ['out','sh','log','trafl','motifs','bonds','ss_detected','models','cluster_models']
    for folder in folders:
        folder_make(current_directory,folder)
    all_files = os.listdir(current_directory)
    for filename in all_files:
        if not os.path.isdir(os.path.join(current_directory,filename)):
            if filename.endswith('_clust.log') or filename.endswith('_AA.pdb') or (filename.endswith('.ss_detected') and filename.__contains__('clust')):
                try:
                    shutil.move(os.path.join(current_directory,filename),os.path.join(current_directory,'cluster_models',filename))
                except:
                    pass
            elif (filename.endswith('.structure') or filename.endswith('.sequence') or filename.endswith('.dbn')):
                try:
                    shutil.move(os.path.join(current_directory,filename),os.path.join(current_directory,'motifs',filename))
                except:
                    pass
            elif filename.endswith('.pdb'):
                try:
                    shutil.move(os.path.join(current_directory,filename),os.path.join(current_directory,'models',filename))
                except:
                    pass
            for folder in folders:
                if filename.endswith(folder):
                    try:
                        shutil.move(os.path.join(current_directory,filename),os.path.join(current_directory,folder,filename))
                    except:
                        pass

def QRNAS_prep(qrnas_program,qrnas_ff_directory,email):
    count = 1
    current_size = 0
    max_size = 10
    current_directory = os.getcwd()
    cluster_directory = folder_check(current_directory,'cluster_models')
    all_files = os.listdir(cluster_directory)
    for filename in all_files:
        if filename.endswith('_AA.pdb'):
            output_file = str(filename.strip('.pdb')+'_QRNAS.pdb')
            log_file = str(filename.strip('.pdb')+'_QRNAS.log')
            if current_size == 0:
                shell_build_start('qrnas'+str(count)+'.sh','qrnas'+str(count),email,tasks=10,notify='END,FAIL')
                with open('qrnas'+str(count)+'.sh','a',newline='\n') as writefile:
                    writefile.writelines('export QRNAS_FF_DIR='+qrnas_ff_directory+'\n')
            with open('qrnas'+str(count)+'.sh','a',newline='\n') as writefile:
                writefile.writelines('%s -i %s -o %s 2>&1> %s &\n' % (qrnas_program,os.path.join(cluster_directory,filename),output_file,log_file))
                current_size += 1
            if current_size >= max_size:
                with open('qrnas'+str(count)+'.sh','a',newline='\n') as writefile:
                    writefile.writelines('wait;\n')
                    count += 1
                    current_size = 0
    with open('qrnas'+str(count)+'.sh','a',newline='\n') as writefile:
        writefile.writelines('wait;\n')

def QRNAS_run():
    current_directory = os.getcwd()
    all_files = os.listdir(current_directory)
    for filename in all_files:
        if filename.startswith('qrnas') and filename.endswith('.sh'):
            os.system('sbatch '+filename)

def QRNAS_cleanup():
    current_directory = os.getcwd()
    folders = ['out','sh','log','qrnas_models']
    for folder in folders:
        folder_make(current_directory,folder)
    all_files = os.listdir(current_directory)
    for filename in all_files:
        if not os.path.isdir(os.path.join(current_directory,filename)):
            if filename.endswith('_AA_QRNAS.pdb'):
                try:
                    shutil.move(os.path.join(current_directory,filename),os.path.join(current_directory,'qrnas_models',filename))
                except:
                    pass
            for folder in folders:
                if filename.endswith(folder):
                    try:
                        shutil.move(os.path.join(current_directory,filename),os.path.join(current_directory,folder,filename))
                    except:
                        pass
def ARES_prep(ares_directory,ares_environment,email):
    current_directory = os.getcwd()
    qrnas_directory = folder_check(current_directory,'qrnas_models') # ARES runs on ALL PDBs within the directory provided - including subdirectories. If running ARES alone, isolate PDB files in a subdirectory to avoid running unoptimized models
    shell_build_start('ares.sh','ares',email,nodes=8,notify='END,FAIL')
    with open('ares.sh','a',newline='\n') as writefile:
        writefile.writelines('conda activate '+ares_environment+'\n')
        writefile.writelines('cd '+ares_directory+'\n')
        writefile.writelines('python -m ares.predict %s data/epoch=0-step=874.ckpt %s/ares.csv -f pdb --nolabels\n' % (qrnas_directory,current_directory))

def ARES_run():
    os.system('sbatch ares.sh')
    
def fpocket_run(fpocket_program):
    current_directory = os.getcwd()
    qrnas_directory = folder_check(current_directory,'qrnas_models')
    all_files = os.listdir(qrnas_directory)
    for filename in all_files:
        if filename.endswith('_QRNAS.pdb'):
            os.system(fpocket_program+' -f '+str(os.path.join(qrnas_directory,filename)))

def fpocket_cleanup():
    current_directory = os.getcwd()
    folder_make(current_directory,'pockets')
    pocket_directory = os.path.join(current_directory,'pockets')
    qrnas_directory = folder_check(current_directory,'qrnas_models')
    all_files = os.listdir(qrnas_directory)
    for filename in all_files:
        if os.path.isdir(os.path.join(qrnas_directory,filename)) and filename != 'pockets':
            try:
                shutil.move(os.path.join(qrnas_directory,filename),os.path.join(pocket_directory,filename))
            except:
                pass

def fpocket_read():
    current_directory = os.getcwd()
    with open('pockets.txt','w',newline='\n') as writefile:
        writefile.writelines('Model Fold RMSD Cluster Pocket Score Druggability\n')
        pocket_directory = folder_check(current_directory,'pockets')
        all_models = os.listdir(pocket_directory)
        for model_directory in all_models:
            all_files = os.listdir(os.path.join(pocket_directory,model_directory))
            for filename in all_files:
                if filename.endswith('_info.txt'):
                    model = ''
                    fold = ''
                    rmsd = ''
                    cluster = ''
                    header = filename.split('_')
                    model = header[0]
                    for i in range(1,20):
                        try:
                            if header[i] == 'all':
                                fold = ''
                                for  j in range(1,i):
                                    fold += header[j]
                            if header[i].startswith('thrs'):
                                rmsd = header[i].split('thrs')[1]
                                cluster = header[i+1].split('-')[0]
                        except:
                            pass
                    with open(os.path.join(pocket_directory,model_directory,filename),'r') as readfile:
                        lines = readfile.readlines()
                        for line in lines:
                            if line.startswith('Pocket'):
                                pocket_number = line.split()[1]
                            if line.startswith('	Score'):
                                pocket_score = line.split()[2]
                            if line.startswith('	Druggability'):
                                drug_score = line.split()[3]
                                writefile.writelines('%s %s %s %s %s %s %s\n' % (model,fold,rmsd,cluster,pocket_number,pocket_score,drug_score))
                                pocket_number = ''
                                pocket_score = ''
                                drug_score = ''

def dock6_prep():
    #read SimRNA pdb and pocket pdb (might need dummy molecule?), convert to mol2 format, add hydrogens and charges, make surface file (Google "Chimera Programmer's Guide")
    #find and center on pocket
    #location of wmoss-lab ZINC15 database(s), copy database of choice to working directory? Or just reference
    #copy /in files
    #write shell script(s)
    #write INSPH file
    print('WIP')

def dock6_run():
    current_directory = os.getcwd()
    all_files = os.listdir(current_directory)
    for filename in all_files:
        if filename.startswith('dock') and filename.endswith('.sh'):
            os.system('sbatch '+filename)

def AnnapuRNA_run():
    #AnnapuRNA directory, read current directory, find dock results
    #run AnnapuRNA on results
    #module load rnaframe/2.7.2
    #annapurna_env (might not need for shell)
    #run_annapurna_env annapurna.py -r 1AJU.pdb -l ARG.sdf -m kNN_modern -o output --groupby
    print('WIP')


def main():

    
    current_directory = os.getcwd()
    
    
    
    """PROGRAM LOCATIONS"""
    cobretti_program = '/work/LAS/wmoss-lab/scripts/cobretti.py'
    ScanFold_program = '/work/LAS/wmoss-lab/scripts/ScanFold.py'
    cmbuilder_program = '/work/LAS/wmoss-lab/scripts/labtools/cm-builder'
    RScape_program = '/work/LAS/wmoss-lab/programs/rscape_v2.0.0.k/bin/R-scape'
    Perl_directory = '/work/LAS/wmoss-lab/programs/lib/perl5/'
    RNAFramework_directory = '/work/LAS/wmoss-lab/programs/RNAFramework/lib/'
    Knotty_program = '/work/LAS/wmoss-lab/programs/knotty/knotty'
    HFold_program = '/work/LAS/wmoss-lab/programs/hfold/HFold_iterative'
    SimRNA_directory = '/work/LAS/wmoss-lab/programs/SimRNA'
    QRNAS_program = '/work/LAS/wmoss-lab/programs/qrnas/QRNA'
    QRNAS_ff = '/work/LAS/wmoss-lab/programs/qrnas/forcefield'
    ARES_program = '/work/LAS/wmoss-lab/programs/ares'
    ARES_environment = '/work/LAS/wmoss-lab/programs/envs/ares'
    fpocket_program = '/work/LAS/wmoss-lab/programs/fpocket2/bin/fpocket'
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-stage',type=str,default='0',help='input stage number to run')
    parser.add_argument('-email',type=str,default='',help='input email address')
    
    parser.add_argument('-i',type=str,default='',help='input list of fasta accession numbers (stage 1A)')
    parser.add_argument('-seq',type=str,default=current_directory,help='input location of fasta sequences')
    parser.add_argument('-db',type=str,default=current_directory,help='input location of BLAST sequence files')
    parser.add_argument('-dbn',type=str,default=current_directory,help='input location of fasta sequences')
    parser.add_argument('-pk',type=str,default=current_directory,help='input location of pseudoknot motifs')
    parser.add_argument('-set_dbn_extend',type=bool,default=True,help='extend pseudoknot motifs')
    
    parser.add_argument('-cobretti',type=str,default=cobretti_program,help='input location of cobretti.py')
    parser.add_argument('-sf',type=str,default=ScanFold_program,help='input location of ScanFold.py')
    parser.add_argument('-cmb',type=str,default=cmbuilder_program,help='input location of cm-builder')
    parser.add_argument('-rs',type=str,default=RScape_program,help='input location of R-Scape')
    parser.add_argument('-perl',type=str,default=Perl_directory,help='input location of Perl directory (perl5) (for R-Scape)')
    parser.add_argument('-rf',type=str,default=RNAFramework_directory,help='input location of RNAFramework directory (RNAFramework/lib/) (for R-Scape)')
    parser.add_argument('-ky',type=str,default=Knotty_program,help='input location of Knotty (knotty)')
    parser.add_argument('-hf',type=str,default=HFold_program,help='input location of Iterative HFold (HFold_iterative)')
    parser.add_argument('-sim',type=str,default=SimRNA_directory,help='input location of SimRNA directory')
    parser.add_argument('-qrnas',type=str,default=QRNAS_program,help='input location of QRNAS')
    parser.add_argument('-qrnasff',type=str,default=QRNAS_ff,help='input location of QRNAS force field')
    parser.add_argument('-ares',type=str,default=ARES_program,help='input location of ARES')
    parser.add_argument('-aresenv',type=str,default=ARES_environment,help='input location of ARES conda environment')
    parser.add_argument('-fpocket',type=str,default=fpocket_program,help='input location of fpocket')
    
    args = parser.parse_args()
    fasta_list = args.i
    stage = args.stage
    email = args.email
    seq_dir = args.seq
    db_dir = args.db
    dbn_dir = args.dbn
    pk_dir = args.pk
    set_dbn_extend = args.set_dbn_extend
    cobretti_prog = args.cobretti
    scanfold_prog = args.sf
    cmbuilder_prog = args.cmb
    rscape_prog = args.rs
    perl_prog = args.perl
    rnaframework_dir = args.rf
    knotty_prog = args.ky
    hfold_prog = args.hf
    simrna_dir = args.sim
    qrnas_prog = args.qrnas
    qrnas_ff = args.qrnasff
    ares_prog = args.ares
    ares_env = args.aresenv
    fpocket_prog = args.fpocket

    """STAGES"""
    if email == '':
        print('No email provided, use "-email" to define string. Exiting...')
        sys.exit()
    if stage == '0':
        print('No stage specified, use "-stage" to define string (e.g., "1A"). Exiting...')
        sys.exit()
    
    if stage == '1A':
        if seq_dir == current_directory:
            folder_make(seq_dir,'sequences')
            seq_dir = os.path.join(seq_dir,'sequences')
        if fasta_list != '':
            fasta_build(fasta_list,seq_dir,email)
        else:
            all_files = os.listdir(current_directory)
            for filename in all_files:
                if filename.endswith('.fa') or filename.endswith('.fasta'):
                    try:
                        shutil.move(os.path.join(current_directory,filename),os.path.join(seq_dir,filename))
                    except:
                        pass
        ScanFold_prep(seq_dir,scanfold_prog,email)
        ScanFold_run()
        if db_dir == current_directory:
            folder_make(db_dir,'databases')
            db_dir = os.path.join(db_dir,'databases')
        blast_prep(seq_dir,db_dir,email)
        blast_run()
    elif stage == '1AA':
        print('1AA: '+stage)
        if seq_dir == current_directory:
            folder_make(seq_dir,'sequences')
            seq_dir = os.path.join(seq_dir,'sequences')
        if fasta_list != '':
            fasta_build(fasta_list,seq_dir,email)
        else:
            all_files = os.listdir(current_directory)
            for filename in all_files:
                if filename.endswith('.fa') or filename.endswith('.fasta'):
                    try:
                        shutil.move(os.path.join(current_directory,filename),os.path.join(seq_dir,filename))
                    except:
                        pass
        ScanFold_prep(seq_dir,scanfold_prog,email)
        ScanFold_run()
    elif stage == '1AB':
        if seq_dir == current_directory:
            seq_dir = folder_check(seq_dir,'sequences')
        if db_dir == current_directory:
            folder_make(db_dir,'databases')
            db_dir = os.path.join(db_dir,'databases')
        blast_prep(seq_dir,db_dir,email)
        blast_run()
        
    elif stage == '1B':
        if seq_dir == current_directory:
            seq_dir = folder_check(seq_dir,'sequences')
        if db_dir == current_directory:
            db_dir = folder_check(db_dir,'databases')
        if dbn_dir == current_directory:
            folder_make(dbn_dir,'motifs')
            dbn_dir = os.path.join(dbn_dir,'motifs')
        blast_cleanup(db_dir,dbn_dir)
        if pk_dir == current_directory:
            folder_make(pk_dir,'pk_motifs')
            pk_dir = os.path.join(pk_dir,'pk_motifs')
        if set_dbn_extend == True:
            extend_motifs(seq_dir,dbn_dir)
        pk_fold(knotty_prog,hfold_prog)
        pk_cleanup()
        pk_breakdown(pk_dir)
        cmbuilder_prep(seq_dir,db_dir,pk_dir,cmbuilder_prog,perl_prog,rnaframework_dir,email,rscape_prog,cobretti_prog)
        cmbuilder_run()
    elif stage == '1BA':
        if seq_dir == current_directory:
            seq_dir = folder_check(seq_dir,'sequences')
        if db_dir == current_directory:
            db_dir = folder_check(db_dir,'databases')
        if dbn_dir == current_directory:
            folder_make(dbn_dir,'motifs')
            dbn_dir = os.path.join(dbn_dir,'motifs')
        blast_cleanup(db_dir,dbn_dir)
    elif stage == '1BB':
        if seq_dir == current_directory:
            seq_dir = folder_check(seq_dir,'sequences')
        if dbn_dir == current_directory:
            dbn_dir = folder_check(dbn_dir,'motifs')
        if pk_dir == current_directory:
            folder_make(pk_dir,'pk_motifs')
            pk_dir = os.path.join(pk_dir,'pk_motifs')
        if set_dbn_extend == True: #set_dbn_extend to False if 'extended.dbn' was made correctly but pks failed
            extend_motifs(seq_dir,dbn_dir)
        pk_fold(knotty_prog,hfold_prog)
        pk_cleanup()
        pk_breakdown(pk_dir)
    elif stage == '1BC': #Substages in case of run failure after PK files are made (usually due to mismatch in db/dbn/fasta file locations/naming convention)
        if seq_dir == current_directory:
            seq_dir = folder_check(seq_dir,'sequences')
        if db_dir == current_directory:
            db_dir = folder_check(db_dir,'databases')
        if pk_dir == current_directory:
            pk_dir = folder_check(pk_dir,'pk_motifs')
        cmbuilder_prep(seq_dir,db_dir,pk_dir,cmbuilder_prog,perl_prog,rnaframework_dir,email,rscape_prog,cobretti_prog)
        cmbuilder_run()
    
    elif stage == '1C':
        cmbuilder_cleanup()
    elif stage == '1CA': #Substage to finish cleanup upon conclusion of R-Scape run
        cmbuilder_cleanup2()
    
    elif stage == '2A':
        SimRNA_prep(simrna_dir,email)
        SimRNA_run()
    elif stage == '2AA':
        SimRNA_prep(simrna_dir,email)
    elif stage == '2AB':
        SimRNA_run()
    
    elif stage == '2B':
        SimRNA_cleanup()
        QRNAS_prep(qrnas_prog,qrnas_ff,email)
        QRNAS_run()
    elif stage == '2BA':
        SimRNA_cleanup()
    elif stage == '2BB':
        QRNAS_prep(qrnas_prog,qrnas_ff,email)
    elif stage == '2BC':
        QRNAS_run()
    
    elif stage == '2C':
        QRNAS_cleanup()
        ARES_prep(ares_prog,ares_env,email)
        ARES_run()
    elif stage == '2CA':
        QRNAS_cleanup()
    elif stage == '2CB':
        ARES_prep(ares_prog,ares_env)
        ARES_run()
    
    elif stage == '2D': #Separated out to avoid ARES reading fpocket .pdb files in subdirectories and to allow for final stage cleanup
        fpocket_run(fpocket_prog)
        fpocket_cleanup()
        fpocket_read()
    
    elif stage == '3A':
        dock6_prep()
        dock6_run()
        
    elif stage == '3B':
        AnnapuRNA_run()
    
    else:
        print('Incorrect stage specified, use "-stage" to define string (e.g., "1A"). Exiting...')
        sys.exit()

            
if __name__ == "__main__":
    main()
