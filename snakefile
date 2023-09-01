localrules:make_list

'''
resources set for 80k samples
'''

#FILES
OUT_DIR='/cluster/work/pausch/naveen/SNPDATA2/HOL'
vcf='/cluster/work/pausch/naveen/SNPDATA2/CLEAN/HOL/hd_imputed/CHR{chr}/combined.vcf.gz'


#PROGRAMS
GCTA='/cluster/home/nkadri/PROGRAMS/GCTA/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 '

#WILDCARDS
chromosomes=range (1,30)

rule all:
    input:
        expand(OUT_DIR + '/GENOME/genome.{ext}',ext=['eigenvec', 'eigenvalues'])

rule make_bed:
    input:
        vcf=vcf
    output:
        binaries=expand(OUT_DIR + '/CHR{{chr}}/geno.{ext}', ext=['fam', 'bim','bed'])
    params:
        prefix=lambda wc, output:output.binaries[0][:-4]
    resources:
        mem_mb=16000,
        walltime='04:00:00'
    envmodules:
        'gcc/8.2.0',
        'plink/1.9-beta6.18'
    shell:
        '''
        plink --cow \
        --vcf {input.vcf} \
        --make-bed \
        --out {params.prefix} \
        --keep-allele-order
        '''
        
rule make_grm:
    input:
        binaries=rules.make_bed.output.binaries
    output:
        grm=expand(OUT_DIR + '/CHR{{chr}}/geno.{ext}', ext=['grm.id','grm.N.bin','grm.bin'])
    params:
        prefix=lambda wc, input:input.binaries[0][:-4],
        maf=0.05
    resources:
        mem_mb=12000,
        walltime='04:00:00'
    threads:
        4    
    shell:
        '''
        {GCTA} --autosome-num 29 \ 
        --bfile {params.prefix} \
        --maf {params.maf} \
        --make-grm \
        --out {params.prefix} \
        --threads {threads}
        '''

rule make_list:
    input:
        grm=expand(OUT_DIR + '/CHR{chr}/geno.{ext}', ext=['grm.id','grm.N.bin','grm.bin'], chr=chromosomes)
    output:
        outfile=OUT_DIR + '/grm_list.txt'
    run:
        with open (output.outfile, 'w') as out:
            for myfile in input.grm:
                if 'grm.id' in myfile:
                    out.write(f'{myfile[:-7]}\n')
    

rule merge_grms:
    input:
        grm=expand(OUT_DIR + '/CHR{chr}/geno.{ext}', ext=['grm.id','grm.N.bin','grm.bin'], chr=chromosomes),
        grm_list=rules.make_list.output.outfile
    output:
        grm=expand(OUT_DIR + '/GENOME/genome.{ext}',ext=['grm.id','grm.N.bin','grm.bin'])
    params:
        prefix=lambda wc,output:output.grm[0][:-7]
    resources:
        mem_mb=180000,
        walltime='24:00:00'
    threads:
        1    
    shell:
        '''
        {GCTA}  --autosome-num 29 \ 
        --mgrm {input.grm_list} \
        --make-grm \
        --out {params.prefix} \
        --threads {threads}
        '''
        
rule pca:
    input:
        grm=rules.merge_grms.output.grm
    params:
        prefix=lambda wc,input:input.grm[0][:-7],
        npc=20
    output:
        pca=expand(OUT_DIR + '/GENOME/genome.{ext}',ext=['eigenvec', 'eigenvalues'])
    resources:
        mem_mb=8000,
        walltime='24:00:00'
    threads:
        20    
    shell:
        '''
        {GCTA}  --autosome-num 29 \ 
        --grm {params.prefix} \
        --pca {params.npc} \
        --out {params.prefix} \
        --threads {threads} \
        '''

