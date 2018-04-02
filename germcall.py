"""
script to call germline variants
Dec 26 2017 Chris Yoon (cjyoon@kaist.ac.kr)

Tools to use
* GATK best practice HaplotypeCaller with gvcf mode
* FreeBayes
* Strelka2

Specify option to run it on WGS or WXS

"""

import subprocess
import sys
import os
import re
import argparse
import yaml
import shlex
import pysam

def argument_parser():
    """parses argument passed on from command line"""
    parser = argparse.ArgumentParser(description='Run germline variation calling using various pipelines')
    parser.add_argument('-b', '--bams', required=False, nargs='+',help='Bam files to run germline calls')
    parser.add_argument('-c', '--config', required=False, default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.yml'))
    parser.add_argument('-e', '--exomeBed', required=False, default=0, help='Exome capture BED file, supply only if Bam files specified in -b arguments are whole exome sequencing')
    parser.add_argument('-l', '--list_of_tools', required=False, nargs='+', default=['strelka', 'freebayes', 'gatk'],
                        help='List of tools to run, by default runs all available tools. Curently supported tools: strelka, freebayes, gatk. If running multiple tools, simply append additional tool name separated with a space. ex) -l strelka gatk')
    parser.add_argument('-o', '--output_dir', required=False, default=os.getcwd(), help='output directory to write results')
    parser.add_argument('-d', '--dryrun', required=False, type=int, default=1,
                        help='If set to 1, will print out shell commands that ar run. If set to 0 will actually run the commands')


    args = vars(parser.parse_args())
    bams = args['bams']
    config = os.path.abspath(args['config'])
    if args['exomeBed'] != 0:
        exomeBed = os.path.abspath(args['exomeBed'])
    else:
        exomeBed = 0
    list_of_tools = args['list_of_tools']
    output_dir = os.path.abspath(args['output_dir'])
    dryrun = args['dryrun']

    # get absolute paths to all the bams
    absBamPaths = []
    for bam in bams:
        absBamPaths.append(os.path.abspath(bam))

    return absBamPaths, config, exomeBed, list_of_tools, output_dir, dryrun

def execute(commandString, dryrun):
    """executes the given commandString. 
    If dryrun == 1, then will print out the commands instead of executing them. """
    if dryrun != 1:
        executeSubprocess = subprocess.Popen(shlex.split(commandString))
        executeSubprocess.wait()
    else:
        print(commandString)

    return 0

def run_bgzip(bgzip_path, afile, dryrun):
    """create a bgzip file for a given vcf file"""
    bgzipCMD = f'{bgzip_path} -f {afile}'
    execute(bgzipCMD, dryrun)

    return os.path.abspath(afile + '.gz')


def run_tabix(tabix_path, gzfile, filetype, dryrun):
    """create a tabix file from a .gz file"""
    tabixCMD = f'{tabix_path} -p {filetype} {gzfile}'
    execute(tabixCMD, dryrun)

    return gzfile


def bgzip_tabix(bgzip_path, tabix_path, afile, filetype, dryrun):
    """one step creation of bgzip and tabix file from a vcf file"""
    gzfile  = run_bgzip(bgzip_path, afile, dryrun)
    tabix = run_tabix(tabix_path, gzfile, filetype, dryrun)

    return gzfile


def sampleNameBam(bamFile):
    """get @RG SM: information as sample name from BAM header"""
    bam = pysam.AlignmentFile(bamFile)
    name = bam.header['RG'][0]['SM']
    return name



def create_outputFolders(output_dir, list_of_tools, dryrun):
    """creates a folder for each tool under output_dir"""

    outputdir_tools = dict()
    for tool in list_of_tools:
        folder = os.path.join(output_dir, tool)
        if not os.path.isdir(folder):
            execute('mkdir ' + folder, dryrun)

        outputdir_tools.update({tool: folder})
    return outputdir_tools


def configure(configPath):
    """configures the path to softwares and input files from a yaml file"""
    with open(configPath, 'r') as f:
        pathMaps = yaml.safe_load(f)

    return pathMaps

 
def run_strelka(python2_path,  strelka_germline_path, tabix_path, bgzip_path, output_dir, reference_path, bamList, exome_bed, dryrun):
    """runs strelka germline with python2, creates a runnable script to execute
    Need to run manta first to use its output"""
    # create strelka configuration file
    bamString = ''
    for bam in bamList:
        bamString += ' --bam ' + bam


    if exome_bed != 0:
        # unlike freebayes, strelka requires BED file to be gzipped with tabix. 
        # thus need to create a new copy of BED file for gzipping and tabixing. 
        strelka_bed = f'{exome_bed}.strelka.bed'
        strelka_bed_gz = f'{exome_bed}.strelka.bed.gz' 
        if os.path.isfile(strelka_bed_gz): 
            # if strelka bed.gz file is alreay present, then skip 
            pass
        else:
            execute(f'cp {exome_bed} {strelka_bed}', dryrun)
            # once copy of BED file is created then, create .bed.gz and .bed.gz.tbi files
            bgzip_tabix(bgzip_path, tabix_path, strelka_bed, 'bed', dryrun)

        exomeString = f' --callRegions {strelka_bed_gz}'
    else:
        exomeString = ''


    strelka_configCMD = f'{python2_path} {strelka_germline_path} {bamString} --referenceFasta {reference_path} --runDir {output_dir} {exomeString}'
    execute(strelka_configCMD, dryrun)

    # now run the runWorkflow.py script generated from'runWorkflow.py' above
    # command
    strelka_script = os.path.join(output_dir, 'runWorkflow.py')
    strelka_runCMD = f'{python2_path} {strelka_script} -m local -j 8 --quiet'
    execute(strelka_runCMD, dryrun)

    return 0


def run_gatk_hayplotypecaller(java_path, gatk_path, output_dir, reference_path, bamPath, dbsnp, dryrun):
    """run gatk haplotypeCaller according to best practices
    dbsnp is only used for annotation, not used in calculation
    """
    output_gvcf = os.path.join(output_dir, sampleNameBam(bamPath)+ '.g.vcf')
    
    haplotypeCaller_runCMD = f"{java_path} -jar {gatk_path} -I HaplotypeCaller -R {reference_fasta} --emitRefConfidence GVCF --dbsnp {dbsnp} -o {output_gvcf}"
    execute(haplotypeCaller_runCMD, dryrun)

    return output_gvcf

def run_gatk_genotypegvcfs(java_path, output_dir, reference_path, gvcf_list, dbsnp, dryrun):
    """once all g.vcf files are created, can do joint genotyping"""
    output_joint_genotypevcf = os.path.join(output_dir, 'joint_genotyped.vcf')
    # create a gvcfstring to be used for the command 
    gvcfString = '' 
    for gvcf in gvcf_list:
        gvcfString = gvcfString +  ' --variant ' + gvcf 

    # command line to run
    genotypegvcf_runCMD = f"{java_path} -jar {gatk_path} -I GenotypeGVCFs -R {refernce_fasta} --dbsnp {dbsnp} -o {output_joint_genotypegvcf} {gvcfString} "
    execute(genotypegvcf_runCMD, dryrun)

    return 0


def run_freebayes(freebayes_path, tabix_path, bgzip_path, output_dir, reference_path, bamPath, exome_bed, dryrun):
    """runs freebayes
    note freebayes works with BED target file, instead of the .bed.gz file like strelka"""
    freebayes_output_dir = os.path.join(output_dir, sampleNameBam(bamPath))
    if os.path.isdir(freebayes_output_dir):
        pass
    else:
        execute(f'mkdir {freebayes_output_dir}', dryrun)
    
    # specify output file 
    output_vcf = os.path.join(freebayes_output_dir, sampleNameBam(bamPath) + '.freebayes.vcf')

    if exome_bed != 0:
        exomeString = f' --target {exome_bed}'
    else:
        exomeString = ''

    freebayes_runCMD = f'{freebayes_path} -f {reference_path} -v {output_vcf} --standard-filters {bamPath} {exomeString}'
    execute(freebayes_runCMD, dryrun)

    # create bgzipped and tabixed file
    bgzip_tabix(bgzip_path, tabix_path, output_vcf, 'vcf', dryrun)

    return 0


def run_freebayes_naive(freebayes_path, tabix_path, bgzip_path, output_dir, reference_path, bamPath, exome_bed, dryrun):
    """runs freebayes
    note freebayes works with BED target file, instead of the .bed.gz file like strelka"""
    freebayes_output_dir = os.path.join(output_dir, sampleNameBam(bamPath))
    if os.path.isdir(freebayes_output_dir):
        pass
    else:
        execute(f'mkdir {freebayes_output_dir}', dryrun)

    # specify output file
    output_vcf = os.path.join(freebayes_output_dir, sampleNameBam(bamPath) + '.freebayes.vcf')

    if exome_bed != 0:
        exomeString = f' --target {exome_bed}'
    else:
        exomeString = ''

    freebayes_runCMD = f'{freebayes_path} -0 --haplotype-length 0 --min-alternate-count 2 --min-alternate-fraction 0.005 -f {reference_path} -v {output_vcf} --standard-filters {bamPath} {exomeString}' # -0 argument is stringent input base and mapping quality filter -m 30 -q 20 -R 0 -S 0
    execute(freebayes_runCMD, dryrun)

    # create bgzipped and tabixed file
    bgzip_tabix(bgzip_path, tabix_path, output_vcf, 'vcf', dryrun)

    return 0

def create_outputFolders(output_dir, list_of_tools, dryrun):
    """creates a folder for each tool under output_dir"""

    outputdir_tools = dict()
    for tool in list_of_tools:
        folder = os.path.join(output_dir, tool)
        if not os.path.isdir(folder):
            execute('mkdir ' + folder, dryrun)

        outputdir_tools.update({tool: folder})
    return outputdir_tools


def main():
    bams, config, exomeBed, list_of_tools, output_dir, dryrun = argument_parser()
    print(exomeBed)
    path = configure(config)
    output_dir_tools = create_outputFolders(output_dir, list_of_tools, dryrun)

    if 'strelka' in list_of_tools:
        run_strelka(path['python2'], path['strelka_germline'], path['tabix'], path['bgzip'], output_dir_tools['strelka'], path['reference'], bams, exomeBed, dryrun)
    
    if 'freebayes' in list_of_tools:
        for bam in bams:
            run_freebayes(path['freebayes'], path['tabix'], path['bgzip'], output_dir_tools['freebayes'], path['reference'], bam, exomeBed, dryrun)
    
    if 'freebayes_naive' in list_of_tools:
        for bam in bams:
            run_freebayes_naive(path['freebayes'], path['tabix'], path['bgzip'], output_dir_tools['freebayes_naive'], path['reference'], bam, exomeBed, dryrun)
 
    if 'gatk' in list_of_tools:
        pass

    return 0

if __name__ == '__main__':
    main()
