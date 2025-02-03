import os
import pdb
import re
import sys
import gzip
import time
import tskit
import subprocess
import pickle
import copy
import msprime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix
import umap.umap_ as umap


# ======================================================================================================================
# generic utils
def exeCmd(cmd, logFile=None, loggerObj=None):  # todo: update the chiang lab version
    if loggerObj:
        loggerObj.info(f"Begining execution: {cmd}")
    else:
        print(f"[{time.strftime('%R:%S-%D')}]Begining execution: {cmd}")
    if logFile:
        with open(logFile, 'w') as fw:
            exe = subprocess.Popen(cmd, stderr=subprocess.STDOUT, close_fds=True,
                                   stdout=fw, universal_newlines=True, shell=True, bufsize=1)
            exe.communicate()  # Wait for process to terminate and set the returncode attribute
            status = exe.returncode
            if status:
                print(f"[{time.strftime('%R:%S-%D')} - ERROR] The detail is recorded in {logFile}")
    else:
        exe = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = exe.communicate()
        status = exe.returncode
        if status:
            if loggerObj:
                loggerObj.error(f"{stderr}")
            else:
                print(f"[{time.strftime('%R:%S-%D')} - ERROR]{stderr}")
    if loggerObj:
        loggerObj.info(f"Finished execution: {cmd}")
    else:
        print(f"[{time.strftime('%R:%S-%D')}]Finished execution: {cmd}")


def eigendecomposition_on_grm(X, num_dims=100):
    """
    PCA on genetic relationship matrix(grm)
    eGRM is a kind of grm
    Copied from egrm-master/manuscript/simulate

    Args:
        X:
        num_dims:

    Returns:

    """
    w, v = np.linalg.eig(X)  # w is eigenvalues, v is eigenvectors
    idx = w.argsort()[::-1]
    w = w[idx]
    v = v[:, idx]
    ndim = min(num_dims, (w >= 0).sum())
    PCs = np.dot(v[:, :ndim], np.diag(np.power(w[:ndim], 0.5)))
    PCs = PCs.real
    return PCs


# ======================================================================================================================
# VCF utils
def createOutVcfFile(outVcfFile_tmp, outVcfFile):  # todo: update the chiang lab version
    if not outVcfFile_tmp.endswith('.gz'):
        exeCmd(f'bgzip {outVcfFile_tmp}')
        outVcfFile_tmp = outVcfFile_tmp + '.gz'
    exeCmd(f'tabix {outVcfFile_tmp}')
    outVcfFile_tmp_index = outVcfFile_tmp + '.tbi'
    outVcfFile_index = outVcfFile + '.tbi'
    if not os.path.exists(outVcfFile_tmp_index):
        exeCmd(f'tabix -C {outVcfFile_tmp}')
        outVcfFile_tmp_index = outVcfFile_tmp + '.csi'
        outVcfFile_index = outVcfFile + '.csi'
    if os.path.exists(outVcfFile):
        os.remove(outVcfFile)
    if os.path.exists(outVcfFile_index):
        os.remove(outVcfFile_index)
    os.rename(outVcfFile_tmp, outVcfFile)
    os.rename(outVcfFile_tmp_index, outVcfFile_index)


def rm_mono_sites(inVcfFile, outVcfFile=None):
    """
    Remove monomorphic sites from a VCF file
    # Args
        inVcfFile: str, input VCF file ending with .gz
        outVcfFile: str, output VCF file ending with .gz, if None, the result will be saved to inVcfFile
        overWrite: bool, if True and the outVcfFile exists, the outVcfFile will be re-generated
    """
    print(f'Removing monomorphic sites from {inVcfFile} ...')
    # replace = False
    if outVcfFile is None:
        outVcfFile = inVcfFile
        # replace = True
    if not outVcfFile.endswith('.gz'):
        outVcfFile += '.gz'
    outVcfFile_tmp = re.sub('.gz', '.tmp', outVcfFile)
    fr = gzip.open(inVcfFile, 'rt')
    with open(outVcfFile_tmp, 'w') as fw:
        row = fr.readline()
        posStart = False
        while row:
            if posStart:
                cols = row.split('\t')
                genotypes = ' '.join(cols[9:])
                genotypes = genotypes.rstrip('\n')
                genotypes = re.sub('\|', ' ', genotypes)  # phased
                genotypes = re.sub('/', ' ', genotypes)  # unphased
                genotypes = genotypes.split(' ')
                if len(list(set(genotypes))) > 1:
                    fw.write(row)
            else:
                if row[:6] == '#CHROM':
                    posStart = True
                fw.write(row)
            row = fr.readline()
    fr.close()
    createOutVcfFile(outVcfFile_tmp, outVcfFile)


def subset_by_samples(inVcfFile, outVcfFile, sampleIDs, force_samples=False, rm_mono=False, suppress_format_error=False,
                      generateLogFile=True):
    """
    subset vcf file by sample IDs
    # Args
        inVcfFile: input vcf file, ending with ".vcf.gz"
        outVcfFile: output vcf file, ending with ".vcf.gz"
        sampleIDs: list, sample IDs to be subsetted
        force_samples: False or True, whether to suppress the error when some of the sample IDs are not found in the input vcf file. If True, how many sample IDs are not found will be recorded in the .subset_by_sampleIDs.log file
        rm_mono: False or True, whether to remove the monomorphic sites produced by the subsetting
        suppress_format_error: when GT field is not consistent with other fields, for example, REF or ALT is "." but GT has no ".", "bcftools view" will raise error. This option will suppress the error
        suppress_format_error: False or True, whether to suppress the error due to the inconsistency between the GT field and other fields, for example, REF or ALT is "." but GT has no "."
    """
    assert outVcfFile.endswith('.vcf.gz')
    os.makedirs(os.path.dirname(outVcfFile), exist_ok=True)
    cmd = 'bcftools view '
    sampleIDsFile = None
    if len(sampleIDs) < 100:
        cmd += f" -s {','.join(sampleIDs)} "
    else:
        sampleIDsFile = outVcfFile.rstrip('.vcf.gz') + '.samples.tmp'
        with open(sampleIDsFile, 'w') as fw:
            fw.write('\n'.join(sampleIDs))
        cmd += f' -S {sampleIDsFile} '
    cmd += f' {inVcfFile} -Oz -o {outVcfFile}'
    if force_samples:
        cmd += ' --force-samples'
    if suppress_format_error:
        cmd += ' -I'
    exeCmd(cmd)
    exeCmd(f'tabix {outVcfFile}')
    if not os.path.exists(outVcfFile + '.tbi'):
        exeCmd(f'tabix -C {outVcfFile}')

    # remove monomorphic sites
    if rm_mono:
        rm_mono_sites(outVcfFile)

    if force_samples and generateLogFile:
        samplesNumFile = outVcfFile.rstrip('.vcf.gz') + '.subset_by_samples.log'
        cmd = f'bcftools query -l {outVcfFile}'
        samples = [i.rstrip('\n') for i in os.popen(cmd).readlines()]
        with open(samplesNumFile, 'w') as fw:
            fw.write(f"# Number of samples\ninput\toutput\n{len(sampleIDs)}\t{len(samples)}")
    if sampleIDsFile is not None:
        os.remove(sampleIDsFile)


def sort_positions(inVcfFile, outVcfFile=None):
    """
    Sort positions in a VCF file
    """
    if outVcfFile is None:
        outVcfFile = inVcfFile
    if not outVcfFile.endswith('.gz'):
        outVcfFile += '.gz'
    outVcfFile_tmp = re.sub('.gz', '.tmp.gz', outVcfFile)
    cmd = f"bcftools sort -Oz -o {outVcfFile_tmp} {inVcfFile}"
    exeCmd(cmd)
    createOutVcfFile(outVcfFile_tmp, outVcfFile)


def rm_dup_positions(inVcfFile, outVcfFile=None):
    """
    Remove duplicated positions from a VCF file
    """
    if outVcfFile is None:
        outVcfFile = inVcfFile
    if not outVcfFile.endswith('.gz'):
        outVcfFile += '.gz'
    outVcfFile_tmp = re.sub('.gz', '.tmp.gz', outVcfFile)
    cmd = f"bcftools norm -d any {inVcfFile} -Oz -o {outVcfFile_tmp}"
    exeCmd(cmd)
    createOutVcfFile(outVcfFile_tmp, outVcfFile)


def read_samples_positions(vcfFile, read_samples, read_positions, multi_chr=False):
    rt = []
    if read_samples:
        cmd = f'bcftools query -l {vcfFile}'
        samples = [i.rstrip('\n') for i in os.popen(cmd).readlines()]
        rt.append(samples)
    if read_positions:
        cmd = f'bcftools query -f "%CHROM:%POS\\n" {vcfFile}'
        positions = [i.rstrip('\n') for i in os.popen(cmd).readlines()]
        if not multi_chr:
            positions = [int(i.split(':')[1]) for i in positions]
        rt.append(positions)
    if len(rt) == 1:
        rt = rt[0]
    return rt


# ======================================================================================================================
# utils for local ancestry inference
def calc_genetic_map_for_sim_data(positions, recombRate):
    posPre, mapPre = positions[0], 0
    pos_map_pairs = [[posPre, mapPre]]
    for p in positions[1:]:
        map = mapPre + (p - posPre) * recombRate
        pos_map_pairs.append([p, map])
        posPre, mapPre = p, map
    return pos_map_pairs


def run_rfmix(queryFile, refFile, sampleMapFile, recombMapFile, outFile, chrN, logFile=None):
    # To avoid the error "/lib64/libstdc++.so.6: version `CXXABI_1.3.9' not found(required by /project/jitang_1167/software/rfmix-master/rfmix)"
    # exeCmd('export PATH=$PATH:/project/jitang_1167/.conda/envs/hawaii/lib')
    # exeCmd('export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/project/jitang_1167/.conda/envs/hawaii/lib')

    os.makedirs(os.path.dirname(outFile), exist_ok=True)
    os.makedirs(os.path.dirname(logFile), exist_ok=True)

    # cmd = '/project/jitang_1167/software/rfmix-master/rfmix '
    cmd = PATH_TO_RFMix  # '/project/chia657_28/programs/rfmix/rfmix '
    cmd += ' -f ' + queryFile
    cmd += ' -r ' + refFile
    cmd += ' -m ' + sampleMapFile
    cmd += ' -g ' + recombMapFile
    cmd += ' -o ' + outFile
    # cmd += ' --chromosome=chr' + chrN
    cmd += ' --chromosome=' + str(chrN)

    # print(cmd)
    exeCmd(cmd, logFile=logFile)


def run_rfmix_wrapper(rfmix_input_output_path,
                      query_pops,
                      ref_pops,
                      sampleIDbyPopNameDict,
                      trees=None,
                      hapMat=None,
                      geneticMapFile=None,
                      recombRate=1e-8,
                      positions=None,
                      positionIDs=None,
                      chrN=None,
                      label=None,
                      genotypes_query_file=None,
                      genotypes_ref_file=None):
    assert hapMat is None or trees is None, 'Only one of hapMat and trees can be provided.'
    outFile = f'{rfmix_input_output_path}/rfmixOut'
    if label is not None:
        outFile += f'_{label}'
    if os.path.exists(outFile + '.msp.tsv'):
        return
    os.makedirs(rfmix_input_output_path, exist_ok=True)

    if chrN is None: chrN = '1'

    # generate sample map file of reference populations
    refSamples = []
    sampleMapFile = f'{rfmix_input_output_path}/ref_samples_map'
    if label is not None:
        sampleMapFile += f'_{label}'
    sampleMapFile += '.tsv'
    refSamples_haps = []
    with open(sampleMapFile, 'w') as fw:
        for pop in ref_pops:
            sampleStart, sampleEnd = sampleIDbyPopNameDict[pop]
            for sample in range(sampleStart, sampleEnd):
                fw.write(f'tsk_{sample}\t{pop}\n')
                refSamples.append(f'tsk_{sample}')
                refSamples_haps += [sample * 2, sample * 2 + 1]

    # get query and ref files from vcf or hapMat
    if genotypes_query_file is None or genotypes_ref_file is None:
        querySamples = []
        querySamples_haps = []
        for pop in query_pops:
            sampleStart, sampleEnd = sampleIDbyPopNameDict[pop]
            for sample in range(sampleStart, sampleEnd):
                querySamples.append(f'tsk_{sample}')
                querySamples_haps += [sample * 2, sample * 2 + 1]
        genotypes_query_file = f'{rfmix_input_output_path}/genotypes_query'
        genotypes_ref_file = f'{rfmix_input_output_path}/genotypes_ref'
        genotypes_all_file = f'{rfmix_input_output_path}/genotypes_all'
        if label is not None:
            genotypes_query_file += f'_{label}'
            genotypes_ref_file += f'_{label}'
            genotypes_all_file += f'_{label}'
        genotypes_query_file += '.vcf'
        genotypes_ref_file += '.vcf'
        genotypes_all_file += '.vcf'
        # if hapMat is not None:
        #     hap2vcf(hapMat[querySamples_haps, :], chrN, positions, positionIDs, genotypes_query_file, querySamples)
        #     hap2vcf(hapMat[refSamples_haps, :], chrN, positions, positionIDs, genotypes_ref_file, refSamples)
        if trees is not None:
            with open(genotypes_all_file, 'w') as fw:
                trees.write_vcf(fw)
            chrN = '1'

            # sort by positions for vcf files
            cmd = f"bcftools sort -o {genotypes_all_file.rstrip('.vcf')}_sorted.vcf {genotypes_all_file}"
            exeCmd(cmd)

            # remove duplicated positions
            # cmd = f"bcftools norm -m+any -Ov {genotypes_all_file.rstrip('.vcf')}_sorted.vcf | awk '!uniq[$2]++' > {genotypes_all_file.rstrip('.vcf')}_rmdup.vcf"
            cmd = f"bcftools norm -d any {genotypes_all_file.rstrip('.vcf')}_sorted.vcf -o {genotypes_all_file.rstrip('.vcf')}_rmdup.vcf"
            exeCmd(cmd)

            os.remove(genotypes_all_file.rstrip('.vcf') + '_sorted.vcf')
            os.remove(genotypes_all_file)
            os.rename(genotypes_all_file.rstrip('.vcf') + '_rmdup.vcf', genotypes_all_file)

            # extract sub-set in terms of query and reference samples
            cmd = f"bcftools view -s {','.join(querySamples)} " + genotypes_all_file + " -o " + genotypes_query_file
            exeCmd(cmd)
            cmd = f"bcftools view -s {','.join(refSamples)} " + genotypes_all_file + " -o " + genotypes_ref_file
            exeCmd(cmd)

    # generate genetic map file
    if geneticMapFile is None:
        if positions is None:
            cmd = "bcftools query -f '%POS,' " + genotypes_query_file
            positions = os.popen(cmd).readlines()[0].rstrip(',').split(',')
            positions = [int(i) for i in positions]
        geneticMapFile = f'{rfmix_input_output_path}/genetic_map'
        if label is not None:
            geneticMapFile += f'_{label}'
        geneticMapFile += '.txt'
        pos_map_pairs = calc_genetic_map_for_sim_data(positions, recombRate)
        with open(geneticMapFile, 'w') as fw:
            for pos, map in pos_map_pairs:
                # fw.write(f'chr1\t{pos}\t{map}\n')
                fw.write(f'{chrN}\t{pos}\t{map}\n')

    # run rfmix
    logFile = f'{outFile}.log'
    run_rfmix(genotypes_query_file, genotypes_ref_file, sampleMapFile, geneticMapFile, outFile, chrN=chrN,
              logFile=logFile)


# ======================================================================================================================
# utils for ARG reconstruction
def vcf2relateHap(outHapFile, outSampleFile, inVCFfile):
    """
    PATH_TO_RELATE / bin / RelateFileFormats \
    - -mode
    ConvertFromVcf \
    - -haps
    example.haps \
    - -sample
    example.sample \
    - i
    example
    Returns:

    """

    cmd = PATH_TO_RELATE + '/bin/RelateFileFormats'
    cmd += ' --mode ConvertFromVcf'
    cmd += ' --haps ' + outHapFile
    cmd += ' --sample ' + outSampleFile
    cmd += ' --input ' + inVCFfile
    exeCmd(cmd)


def FlipHapsUsingAncestor(inHapFile, inSampleFile, inAncestalAllele, outFile):
    cmd = PATH_TO_RELATE + '/bin/RelateFileFormats '
    cmd += ' --mode FlipHapsUsingAncestor'
    cmd += ' --haps ' + inHapFile
    cmd += ' --sample ' + inSampleFile
    cmd += ' --ancestor ' + inAncestalAllele
    cmd += ' -o ' + outFile
    exeCmd(cmd)


def run_relate(hapFile, sampleFile, recombMapFile, Ne, mutRate, outFileName, RelateOutPath,
               memory=None,  # only number(unit: gigabytes)
               logFile=None, idx=None,
               totalN=None, step=None, overWrite=False):
    if not overWrite:
        if os.path.exists(RelateOutPath + '/' + outFileName + '.anc') and \
                os.path.exists(RelateOutPath + '/' + outFileName + '.mut'):
            return
    if idx is not None and totalN is not None:
        if step is None:
            step = 10000
        if idx % step == 0:
            print('Progress of hap2RelateARG: %i/%i' % (idx, totalN))
    cmd = PATH_TO_RELATE + '/bin/Relate --mode All'
    if memory:
        cmd += ' --memory ' + memory
    cmd += ' --haps ' + hapFile
    cmd += ' --sample ' + sampleFile
    cmd += ' --map ' + recombMapFile
    cmd += ' -N ' + str(Ne)
    cmd += ' -m ' + str(mutRate)
    cmd += ' -o ' + outFileName
    cmd += ' --seed 1 '

    current_dir = os.getcwd()
    os.chdir(RelateOutPath)
    exeCmd(cmd, logFile=logFile)
    os.chdir(current_dir)


def run_relate_wrapper(relate_input_output_path, Ne, mutRate, recombRate=None, recombMapFile=None,
                       vcfFile=None, chrN=None, fileName=None, treesFile=None,
                       positions=None, memory='16', logFile=None):
    # if vcfFile is not None:
    inputType = 'vcf'
    if fileName is None:
        fileName = vcfFile.split('/')[-1].rstrip('.vcf').rstrip('.vcf.gz')
    outFileName = f'{fileName}.relate'

    os.makedirs(relate_input_output_path, exist_ok=True)
    if treesFile is None:
        treesFile = f"{relate_input_output_path}/{outFileName}.trees"
    if not os.path.exists(treesFile):
        print(f"generating relate trees for {outFileName} ...", flush=True)
        hapFile = f"{relate_input_output_path}/{outFileName}.haps"
        sampleFile = f"{relate_input_output_path}/{outFileName}.sample"
        if inputType == 'vcf':
            vcf2relateHap(hapFile, sampleFile, vcfFile.rstrip('.vcf').rstrip('.vcf.gz'))
            # confirm chromosome is integer, for example 1 instead of chr1
            hapFile_tmp = hapFile + '.tmp'
            with open(hapFile, 'r') as f:
                with open(hapFile_tmp, 'w') as f2:
                    for line in f:
                        if line.startswith(f'chr{chrN}'):
                            line = line.lstrip('chr')
                        f2.write(line)
            exeCmd(f'rm -f {hapFile}')
            exeCmd(f'mv {hapFile_tmp} {hapFile}')
        if positions is None:
            positions = []
            with open(hapFile, 'r') as fr:
                for i in fr.readlines():
                    positions.append(int(i.split()[2]))
        exeCmd(f'gzip -f {hapFile} > {hapFile}.gz && rm -f {hapFile}')
        hapFile = hapFile + '.gz'
        if recombMapFile is None:
            recombMapFile = f"{relate_input_output_path}/{outFileName}.map"
            pos_map_pairs = calc_genetic_map_for_sim_data(positions, recombRate)
            with open(recombMapFile, "w") as fw:
                fw.write("pos COMBINED_rate Genetic_Map\n")
                for pos, map in pos_map_pairs:
                    string = str(pos) + " " + str(recombRate * 1000000) + " "
                    string = string + str(map) + "\n"
                    fw.write(string)

        run_relate(hapFile, sampleFile, recombMapFile, Ne, mutRate, outFileName=f'{outFileName}',
                   RelateOutPath=relate_input_output_path, logFile=logFile, memory=memory)
        os.system(
            PATH_TO_RELATE + "/bin/RelateFileFormats --mode ConvertToTreeSequence "
            + "-i "
            + f"{relate_input_output_path}/{outFileName}"
            + " -o "
            + f"{re.sub('.trees', '', treesFile)}"
        )
    else:
        print('File exists: ' + treesFile)

    trees_relate = tskit.load(treesFile)
    return trees_relate


# ======================================================================================================================
# utils for demographic model
def step_mig_mat_pairs(pops, n_rows, n_cols, symmetric):
    assert len(pops) == n_rows * n_cols
    pmat = np.array(pops).reshape(n_rows, n_cols)
    mig_mat_pairs = []
    for i in range(n_rows):
        for j in range(n_cols):
            if i > 0:
                mig_mat_pairs.append((pmat[i - 1, j], pmat[i, j]))  # with top
            if i < n_rows - 1:
                mig_mat_pairs.append((pmat[i + 1, j], pmat[i, j]))  # with bottom
            if j > 0:
                mig_mat_pairs.append((pmat[i, j - 1], pmat[i, j]))  # with left
            if j < n_cols - 1:
                mig_mat_pairs.append((pmat[i, j + 1], pmat[i, j]))  # with right
    if symmetric:
        mig_mat_pairs_symm = []
        for p1, p2 in mig_mat_pairs:
            if (p1, p2) not in mig_mat_pairs_symm and (p2, p1) not in mig_mat_pairs_symm:
                mig_mat_pairs_symm.append((p1, p2))
        return mig_mat_pairs_symm
    else:
        return mig_mat_pairs


def defineDemog(args, sampleAnc=False):
    demography = msprime.Demography()
    popIDbyNameDict = {}
    sampleIDbyPopNameDict = {}
    if 'twoAncTwoDer' in args['demogName']:
        # seems the root has to be one population, can not be two populations or more
        a2_Ne = args['Ne']
        if 'a2-Ne' in args:
            a2_Ne = args['a2-Ne']
        demography.add_population(name="d1", initial_size=args['Ne'], growth_rate=0)
        demography.add_population(name="d2", initial_size=args['Ne'], growth_rate=0)
        demography.add_population(name="a1", initial_size=args['Ne'], growth_rate=0)
        demography.add_population(name="a2", initial_size=a2_Ne, growth_rate=0)
        popIDbyNameDict.update({'d1': 0, 'd2': 1, 'a1': 2, 'a2': 3})
        sampleSize = args['sample_size']
        if sampleAnc:
            sampleIDbyPopNameDict = {'d1': (0, sampleSize), 'd2': (sampleSize, sampleSize * 2),
                                     'a1': (sampleSize * 2, sampleSize * 3), 'a2': (sampleSize * 3, sampleSize * 4)}
        else:
            sampleIDbyPopNameDict = {'d1': (0, sampleSize), 'd2': (sampleSize, sampleSize * 2)}
        if args['demogName'] in ['twoAncTwoDer6', 'twoAncTwoDer7', 'twoAncTwoDer8', 'twoAncTwoDer9', 'twoAncTwoDer10',
                                 'twoAncTwoDer11']:
            demography.add_population(name="a1_1", initial_size=args['Ne'], growth_rate=0)
            demography.add_population(name="a1_2", initial_size=args['Ne'], growth_rate=0)
            popIDbyNameDict.update({'a1_1': 4, 'a1_2': 5})
            demography.add_admixture(time=args['t_admix'],
                                     derived="d1",
                                     # ancestral=["d1", "a2"], # the ancestral and the derived must be different
                                     ancestral=["a1_1", "a2"],
                                     proportions=[1 - args['a2_to_d1_prop'], args['a2_to_d1_prop']])
            demography.add_admixture(time=args['t_admix'],
                                     derived="d2",
                                     ancestral=["a1_2", "a2"],
                                     proportions=[1 - args['a2_to_d2_prop'], args['a2_to_d2_prop']])
            if args['migration_rate'] > 0:
                demography.set_symmetric_migration_rate(populations=[popIDbyNameDict["d1"], popIDbyNameDict["d2"]],
                                                        rate=args['migration_rate'])

                demography.set_symmetric_migration_rate(populations=[popIDbyNameDict["a1_1"], popIDbyNameDict["a1_2"]],
                                                        rate=0)  #
                if args['migration_stop_t'] > args['t_admix']:
                    demography.add_symmetric_migration_rate_change(time=args['t_admix'],
                                                                   populations=[popIDbyNameDict["a1_1"],
                                                                                popIDbyNameDict["a1_2"]],
                                                                   rate=args['migration_rate'])
                    demography.add_symmetric_migration_rate_change(time=args['migration_stop_t'],
                                                                   populations=[popIDbyNameDict["a1_1"],
                                                                                popIDbyNameDict["a1_2"]],
                                                                   rate=0)
                else:
                    demography.add_symmetric_migration_rate_change(time=args['migration_stop_t'],
                                                                   populations=[popIDbyNameDict["d1"],
                                                                                popIDbyNameDict["d2"]],
                                                                   rate=0)
            demography.add_population_split(time=args['t_split_2'], derived=["a1_1", "a1_2"], ancestral="a1")
            demography.add_population_split(time=args['t_split_1'], derived=["a2"], ancestral="a1")
        else:
            demography.add_admixture(time=args['t_admix'],
                                     derived="d1",
                                     ancestral=["a1", "a2"],
                                     proportions=[args['a1_to_d1_prop'], 1 - args['a1_to_d1_prop']])
            demography.add_admixture(time=args['t_admix'], derived="d2", ancestral=["a1", "a2"],
                                     proportions=[1 - args['a2_to_d2_prop'], args['a2_to_d2_prop']])
            demography.add_population_split(time=args['t_split'], derived=["a2"], ancestral="a1")
    if 'GridDer' in args['demogName']:
        popID = -1
        sampleSize = args['sample_size']
        mig_mat_pops = []
        for i in range(args['n_rows']):
            for j in range(args['n_cols']):
                popID += 1
                mig_mat_pops.append(popID)
                demography.add_population(name=f"d{i}{j}", initial_size=args['Ne'], growth_rate=0)
                popIDbyNameDict.update({f"d{i}{j}": popID})
                sampleIDbyPopNameDict.update({f"d{i}{j}": (sampleSize * popID, sampleSize * (popID + 1))})

        demography.add_population(name="a1", initial_size=args['Ne'], growth_rate=0)
        popID += 1
        popIDbyNameDict.update({'a1': popID})
        if sampleAnc:
            sampleIDbyPopNameDict.update({'a1': (sampleSize * popID, sampleSize * (popID + 1))})
        demography.add_population(name="a2", initial_size=args['Ne'], growth_rate=0)
        popID += 1
        popIDbyNameDict.update({'a2': popID})
        if sampleAnc:
            sampleIDbyPopNameDict.update({'a2': (sampleSize * popID, sampleSize * (popID + 1))})
        a1_d_list = []
        mig_mat_anc_pops = []
        for i in range(args['n_rows']):
            for j in range(args['n_cols']):
                popID += 1
                mig_mat_anc_pops.append(popID)
                demography.add_population(name=f"a1_d{i}{j}", initial_size=args['Ne'], growth_rate=0)
                popIDbyNameDict.update({f"a1_d{i}{j}": popID})
                a1_d_list.append(f"a1_d{i}{j}")
                prop_from_a2 = args['prop_from_a2_mat'][i][j]
                demography.add_admixture(time=args['t_admix'],
                                         derived=f"d{i}{j}",
                                         ancestral=[f"a1_d{i}{j}", "a2"],
                                         proportions=[1 - prop_from_a2, prop_from_a2])

        # set up migration matrix, should be equal to "egrm-master->manuscript->simulation->step_mig_mat()"
        if args['migration_rate'] > 0:
            mig_pairs = step_mig_mat_pairs(mig_mat_pops, args['n_rows'], args['n_cols'], symmetric=True)
            for p1, p2 in mig_pairs:
                demography.set_symmetric_migration_rate(populations=[p1, p2], rate=args['migration_rate'])
            mig_pairs_anc = step_mig_mat_pairs(mig_mat_anc_pops, args['n_rows'], args['n_cols'], symmetric=True)
            for p1, p2 in mig_pairs_anc:
                demography.set_symmetric_migration_rate(populations=[p1, p2], rate=0)

            if args['migration_stop_t'] > args['t_admix']:
                for p1, p2 in mig_pairs_anc:
                    demography.add_symmetric_migration_rate_change(time=args['t_admix'],
                                                                   populations=[p1, p2],
                                                                   rate=args['migration_rate'])
                for p1, p2 in mig_pairs_anc:
                    demography.add_symmetric_migration_rate_change(time=args['migration_stop_t'],
                                                                   populations=[p1, p2],
                                                                   rate=0)
            else:
                for p1, p2 in mig_pairs:
                    demography.add_symmetric_migration_rate_change(time=args['migration_stop_t'],
                                                                   populations=[p1, p2],
                                                                   rate=0)
        demography.add_population_split(time=args['t_split_2'], derived=a1_d_list, ancestral="a1")
        demography.add_population_split(time=args['t_split_1'], derived=["a2"], ancestral="a1")
    if 'gLikeFig6A' in args['demogName']:
        sampleSize = args['sample_size']
        params_latino = [25, 353, 1018, 2094, 0.107, 0.442, 0.0, 83157, 9973, 26682, 10000, 146339, 6185, 5895, 5693,
                         # 0.078
                         0.132]
        t1, t2, t3, t4, r1, r2, r3, N_admix, N_afr, N_eur, N_asia, N_pol, N_aa, N_ooa, N_anc, gr = params_latino
        # migration_rate = 0.01
        # demography.add_population(name="admix", initial_size=N_admix, growth_rate=gr)
        d1_d2_Ne = args['d1_d2_Ne']
        pol_d1_d2_Ne = args['pol_d1_d2_Ne']
        latino_growth_rate = 0
        if args['latino_growth']:
            latino_growth_rate = gr
        # define d1
        demography.add_population(name="d1", initial_size=d1_d2_Ne, growth_rate=latino_growth_rate)
        popID = 0
        popIDbyNameDict.update({f"d1": popID})
        sampleIDbyPopNameDict.update({f"d1": (sampleSize * popID, sampleSize * (popID + 1))})

        # define d2
        demography.add_population(name="d2", initial_size=d1_d2_Ne, growth_rate=latino_growth_rate)
        popID += 1
        popIDbyNameDict.update({f"d2": popID})
        sampleIDbyPopNameDict.update({f"d2": (sampleSize * popID, sampleSize * (popID + 1))})

        # stop growth if needed
        if args['latino_growth']:
            if args['latino_growth_stop_t'] > 0:
                demography.add_population_parameters_change(time=args['latino_growth_stop_t'],
                                                            initial_size=args['Ne_after_growth_stop'],
                                                            growth_rate=0,
                                                            population="d1")
                demography.add_population_parameters_change(time=args['latino_growth_stop_t'],
                                                            initial_size=args['Ne_after_growth_stop'],
                                                            growth_rate=0,
                                                            population="d2")

        demography.add_population(name="afr", initial_size=N_afr)
        popID += 1
        popIDbyNameDict.update({f"afr": popID})
        sampleIDbyPopNameDict.update({f"afr": (sampleSize * popID, sampleSize * (popID + 1))})

        demography.add_population(name="eur", initial_size=N_eur)
        popID += 1
        popIDbyNameDict.update({f"eur": popID})
        sampleIDbyPopNameDict.update({f"eur": (sampleSize * popID, sampleSize * (popID + 1))})

        demography.add_population(name="pol", initial_size=N_pol)
        popID += 1
        popIDbyNameDict.update({f"pol": popID})
        sampleIDbyPopNameDict.update({f"pol": (sampleSize * popID, sampleSize * (popID + 1))})

        demography.add_population(name="pol_d1", initial_size=pol_d1_d2_Ne, growth_rate=0)
        popID += 1
        popIDbyNameDict.update({f"pol_d1": popID})
        demography.add_population(name="pol_d2", initial_size=pol_d1_d2_Ne, growth_rate=0)
        popID += 1
        popIDbyNameDict.update({f"pol_d2": popID})
        demography.add_population(name="asia", initial_size=N_asia)

        # demography.add_admixture(time=t1, derived="admix", ancestral=["afr", "eur", "asia", "pol"],
        #                          proportions=[r1, r2, r3, 1 - r1 - r2 - r3])
        afr_d1_prop, eur_d1_prop = args['d1_anc_prop']
        demography.add_admixture(time=args['t_admix'],
                                 derived="d1",
                                 ancestral=["afr", "eur", "pol_d1"],
                                 # proportions=[r1, r2, 1 - r1 - r2]
                                 proportions=[afr_d1_prop, eur_d1_prop, 1 - afr_d1_prop - eur_d1_prop]
                                 )
        afr_d2_prop, eur_d2_prop = args['d2_anc_prop']
        demography.add_admixture(time=args['t_admix'],
                                 derived="d2",
                                 ancestral=["afr", "eur", "pol_d2"],
                                 proportions=[afr_d2_prop, eur_d2_prop, 1 - afr_d2_prop - eur_d2_prop])

        if args['latino_migration_rate'] > 0:
            demography.set_symmetric_migration_rate(populations=[popIDbyNameDict["d1"], popIDbyNameDict["d2"]],
                                                    rate=args['latino_migration_rate'])

            demography.set_symmetric_migration_rate(populations=[popIDbyNameDict["pol_d1"], popIDbyNameDict["pol_d2"]],
                                                    rate=0)
            if args['latino_migration_stop_t'] > 0:
                if args['latino_migration_stop_t'] > args['t_admix']:
                    demography.add_symmetric_migration_rate_change(time=args['t_admix'],
                                                                   populations=[popIDbyNameDict["pol_d1"],
                                                                                popIDbyNameDict["pol_d2"]],
                                                                   rate=args['latino_migration_rate'])
                    demography.add_symmetric_migration_rate_change(time=args['latino_migration_stop_t'],
                                                                   populations=[popIDbyNameDict["pol_d1"],
                                                                                popIDbyNameDict["pol_d2"]],
                                                                   rate=0)
                else:
                    demography.add_symmetric_migration_rate_change(time=args['latino_migration_stop_t'],
                                                                   populations=[popIDbyNameDict["d1"],
                                                                                popIDbyNameDict["d2"]],
                                                                   rate=0)

        demography.add_population_split(time=args['t_split'], derived=["pol_d1", "pol_d2"], ancestral="pol")
        demography.add_population_split(time=t2, derived=["pol"], ancestral="asia")
        demography.add_population_parameters_change(time=t2, initial_size=N_aa, growth_rate=0, population="asia")

        demography.add_population_split(time=t3, derived=["asia"], ancestral="eur")
        demography.add_population_parameters_change(time=t3, initial_size=N_ooa, growth_rate=0, population="eur")

        demography.add_population_split(time=t4, derived=["eur"], ancestral="afr")
        demography.add_population_parameters_change(time=t4, initial_size=N_anc, growth_rate=0, population="afr")

    demography.sort_events()
    demography.debug()
    return demography, popIDbyNameDict, sampleIDbyPopNameDict


def schemesForSim(simScheme, sampleAnc, sample_size=None):
    def get_args(simScheme, sample_size):
        argsList = []
        if 'GridDer' in simScheme:
            n_cols, n_rows = 0, 0
            if '5x5GridDer' in simScheme:
                n_cols, n_rows = 5, 5
            if '3x3GridDer' in simScheme:
                n_cols, n_rows = 3, 3
            if sample_size is None:
                sample_size = 100
            # t_admix_list = [90, 50, 20]
            # t_admix_list = [60, 20]
            t_admix_list = [25]
            t_split_1_list = [2000]
            recombination_rate = 1e-8
            mutation_rate = 1e-8
            Ne = 500  # used when growth rate is 0
            t_split_list = []
            prop_from_a2_mats = []
            migration_rate = 0
            migration_stop_t = ''
            if '3x3GridDer3' in simScheme:
                t_split_list = [50, 100, 300]
                prop_from_a2_mats = [
                    ([[0.2, 0.3, 0.4]] * n_rows, '0.2-0.3-0.4'),
                    ([[0.3, 0.4, 0.5]] * n_rows, '0.3-0.4-0.5'),
                    ([[0.4, 0.5, 0.6]] * n_rows, '0.4-0.5-0.6'),
                ]
                migration_rate = 0.01
                migration_stop_t = 10
            for t_split_1 in t_split_1_list:
                for prop_from_a2_mat, prop_from_a2_mat_label in prop_from_a2_mats:
                    for t in t_split_list:
                        for t_admix in t_admix_list:
                            argsList.append({'n_cols': n_cols,
                                             'n_rows': n_rows,
                                             'prop_from_a2_mat': prop_from_a2_mat,
                                             'prop_from_a2_mat_label': prop_from_a2_mat_label,
                                             't_admix': t_admix,
                                             't_split_1': t_split_1,
                                             't_split_2': t,
                                             # 'migration_rate': m,
                                             'recombination_rate': recombination_rate,
                                             'sample_size': sample_size,
                                             'mutation_rate': mutation_rate,
                                             'migration_stop_t': migration_stop_t,
                                             'Ne': Ne,
                                             'chrLength': 10000000,
                                             # 'chrLength': 100000,
                                             'demogName': simScheme,
                                             'migration_rate': migration_rate
                                             })
        if 'gLikeFig6A' in simScheme:
            t_split_list = [50, 100, 300]
            t_admix_list = [25]
            # t_admix_list = [20]
            # Ne should use 10000 instead of 146339, because it is about the sub-populations
            d1_d2_Ne = 10000
            pol_d1_d2_Ne = 10000
            latino_growth = False
            latino_migration = False
            Ne = 10000
            recombination_rate = 1e-8
            mutation_rate = 1e-8
            if sample_size is None:
                sample_size = 100
            latino_migration_rates = [0]
            afr_prop, eur_prop = 0.107, 0.442
            d2_anc_props = [(afr_prop, eur_prop)]
            d1_anc_prop = (0.107, 0.442)

            # generate argsList
            for t in t_split_list:
                for t_admix in t_admix_list:
                    for mr in latino_migration_rates:
                        for d2_anc_prop in d2_anc_props:
                            argsList.append({'t_admix': t_admix,
                                             't_split': t,
                                             'recombination_rate': recombination_rate,
                                             'sample_size': sample_size,
                                             'mutation_rate': mutation_rate,
                                             'd1_d2_Ne': d1_d2_Ne,
                                             'pol_d1_d2_Ne': pol_d1_d2_Ne,
                                             'Ne': Ne,
                                             'chrLength': 10000000,
                                             'demogName': simScheme,
                                             'latino_growth': latino_growth,
                                             # 'latino_migration': latino_migration,
                                             'latino_migration_rate': mr,
                                             'latino_migration_stop_t': t,
                                             'd2_anc_prop': d2_anc_prop,
                                             'd1_anc_prop': d1_anc_prop,
                                             })
        return argsList

    if simScheme == 'twoAncTwoDer8':
        argsList = get_args(simScheme='twoAncTwoDer6', sample_size=sample_size)
        argsList_1 = []
        for i in argsList:
            i['recombination_rate'] = 1e-8
            i['migration_rate'] = 0.01
            i['migration_stop_t'] = 10
            i['demogName'] = 'twoAncTwoDer8'
            argsList_1.append(i)
        argsList = argsList_1
    elif simScheme == 'gLikeFig6A-6':
        argsList = get_args(simScheme='gLikeFig6A-4', sample_size=sample_size)
        argsList_1 = []
        for i in argsList:
            i['latino_migration_rate'] = 0.01
            i['latino_migration_stop_t'] = 10
            i['latino_growth_stop_t'] = 10
            i['Ne_after_growth_stop'] = 10000
            i['demogName'] = 'gLikeFig6A-6'
            argsList_1.append(i)
        argsList = argsList_1
    else:
        argsList = get_args(simScheme=simScheme, sample_size=sample_size)

    if 'twoAncTwoDer' in simScheme:
        for i in range(len(argsList)):
            sample_size = argsList[i]['sample_size']
            if sampleAnc:
                sampleSizeByPopDict = {'d1': sample_size, 'd2': sample_size, 'a1': sample_size, 'a2': sample_size}
                totalSampleSize = sample_size * 4
                argsList[i]['rfmix_query_pops'] = ['d1', 'd2']
                argsList[i]['rfmix_ref_pops'] = ['a1', 'a2']
            else:
                sampleSizeByPopDict = {'d1': sample_size, 'd2': sample_size}
                totalSampleSize = sample_size * 2,
            argsList[i]['sampleSizeByPopDict'] = sampleSizeByPopDict
            argsList[i]['totalSampleSize'] = totalSampleSize
            totalQuerySampleSize = sample_size * 2
            totalQuerySampleSize_dip = int(totalQuerySampleSize / 2)
            argsList[i]['totalQuerySampleSize'] = totalQuerySampleSize
            argsList[i]['totalQuerySampleSize_dip'] = totalQuerySampleSize_dip
    if 'gLikeFig6A' in simScheme:
        for i in range(len(argsList)):
            sample_size = argsList[i]['sample_size']
            sampleSizeByPopDict = {'d1': sample_size, 'd2': sample_size, 'afr': sample_size, 'eur': sample_size,
                                   'pol': sample_size}
            totalSampleSize = sample_size * 5
            totalQuerySampleSize = sample_size * 2
            totalQuerySampleSize_dip = int(totalQuerySampleSize / 2)
            argsList[i]['rfmix_query_pops'] = ['d1', 'd2']
            argsList[i]['rfmix_ref_pops'] = ['afr', 'eur', 'pol']
            argsList[i]['sampleSizeByPopDict'] = sampleSizeByPopDict
            argsList[i]['totalSampleSize'] = totalSampleSize
            argsList[i]['totalQuerySampleSize'] = totalQuerySampleSize
            argsList[i]['totalQuerySampleSize_dip'] = totalQuerySampleSize_dip
    if 'GridDer' in simScheme:
        for i_args in range(len(argsList)):
            sample_size = argsList[i_args]['sample_size']
            sampleSizeByPopDict = {}
            totalSampleSize = 0
            argsList[i_args]['rfmix_query_pops'] = []
            for i in range(argsList[i_args]['n_rows']):
                for j in range(argsList[i_args]['n_cols']):
                    sampleSizeByPopDict.update({f"d{i}{j}": sample_size})
                    argsList[i_args]['rfmix_query_pops'].append(f"d{i}{j}")
                    totalSampleSize += sample_size
            if sampleAnc:
                argsList[i_args]['rfmix_ref_pops'] = ['a1', 'a2']
                for i in argsList[i_args]['rfmix_ref_pops']:
                    sampleSizeByPopDict.update({i: sample_size})
                    totalSampleSize += sample_size
            argsList[i_args]['sampleSizeByPopDict'] = sampleSizeByPopDict
            argsList[i_args]['totalSampleSize'] = totalSampleSize
            if '3x3GridDer' in argsList[i_args]['demogName']:
                totalQuerySampleSize = sample_size * 9
                totalQuerySampleSize_dip = int(totalQuerySampleSize / 2)
                argsList[i_args]['totalQuerySampleSize'] = totalQuerySampleSize
                argsList[i_args]['totalQuerySampleSize_dip'] = totalQuerySampleSize_dip

    # add sample size of diploid mode
    for i in range(len(argsList)):
        sampleSizeByPopDict_ = copy.deepcopy(argsList[i]['sampleSizeByPopDict'])
        totalSampleSize_dip = 0
        for pop in sampleSizeByPopDict_:
            sampleSizeByPopDict_[pop] = int(sampleSizeByPopDict_[pop] / 2)
            totalSampleSize_dip += sampleSizeByPopDict_[pop]
        argsList[i]['sampleSizeByPopDict_dip'] = sampleSizeByPopDict_
        argsList[i]['totalSampleSize_dip'] = totalSampleSize_dip
    return argsList


def get_samples_by_pop(samplesFile=None, pops=None, sampleIDbyPopNameDict=None, diploid=False):
    if not os.path.exists(samplesFile):
        samples = []
        for pop in pops:
            sampleStart, sampleEnd = sampleIDbyPopNameDict[pop]
            for sample in range(sampleStart, sampleEnd):
                samples.append(sample)
        if diploid:
            samples_dip = []
            for i in range(0, len(samples), 2):
                samples_dip.append(int(samples[i] / 2))
            samples = samples_dip
        with open(samplesFile, 'w') as fw:
            fw.write(str(samples))
    else:
        with open(samplesFile, 'r') as fr:
            samples = eval(fr.read())
    return samples


# ======================================================================================================================
# plotting utils
def CalcAncProp(simScheme, chrLen, samples, local_anc_calls, anc_id_by_name_dict=None):
    print('CalcAncProp ...')
    ancPropBySampleDict = {}
    if type(local_anc_calls) != dict:
        if 'twoAnc' in simScheme:
            lens_from_a1 = np.zeros(len(samples))
            a1_id = anc_id_by_name_dict['a1']
            for start, end, anc in local_anc_calls:
                for i in np.where(np.array(anc) == a1_id)[0]:
                    lens_from_a1[i] += end - start
            for i in samples:
                a1_prop = lens_from_a1[i] / chrLen
                ancPropBySampleDict[i] = {'a1': a1_prop, 'a2': 1 - a1_prop}
        if 'gLikeFig6A' in simScheme:
            for i in samples:
                ancPropBySampleDict[i] = {}
            for p in ['afr', 'eur', 'pol']:
                lens_from_p = np.zeros(len(samples))
                p_id = anc_id_by_name_dict[p]
                for start, end, anc in local_anc_calls:
                    for i in np.where(np.array(anc) == p_id)[0]:
                        lens_from_p[i] += end - start
                for i in samples:
                    ancPropBySampleDict[i].update({p: lens_from_p[i] / chrLen})
    print('CalcAncProp done!')
    return ancPropBySampleDict


def average_anc_props(ancPropBySampleDict_list, samples, ancPops):
    ancPropBySampleDict = {}
    for i in samples:
        ancPropBySampleDict[i] = {}
    for p in ancPops:
        for i in samples:
            ancPropBySampleDict[i].update({p: 0})
        for ancPropBySampleDict_i in ancPropBySampleDict_list:
            for i in samples:
                ancPropBySampleDict[i][p] += ancPropBySampleDict_i[i][p]
        for i in samples:
            ancPropBySampleDict[i][p] /= len(ancPropBySampleDict_list)
    return ancPropBySampleDict


def separation_index(true_labels, x, dmat=None):
    if dmat is None:
        dmat = distance_matrix(x, x)
    omat = dmat.argsort(axis=1)

    scores = []
    for i in np.arange(true_labels.shape[0]):
        fellows = np.where(true_labels == true_labels[i])[0]
        neighbors = omat[i, : fellows.shape[0]]
        intersects = np.intersect1d(neighbors, fellows)
        score = intersects.shape[0] / fellows.shape[0]
        scores.append(score)

    scores = np.array(scores)
    return scores


def colorSampleByPop(pops, sampleSizeByPopDict):
    colors = [[228, 26, 28],
              [55, 126, 184],
              [77, 175, 74],
              [152, 78, 163],
              [255, 127, 0],
              [255, 255, 51],
              [166, 86, 40],
              [247, 129, 191],
              [153, 153, 153]]
    assert len(pops) <= len(colors)
    sampleRGBs = []
    for pop, color in zip(pops, colors):
        sampleRGBs += [color] * sampleSizeByPopDict[pop]
    sampleRGBs = np.array(sampleRGBs) / 255
    return sampleRGBs


def colorSampleByAncProp(ancPropBySampleDict, samples, colorLabel):
    sampleRGBs = []
    for sample in samples:
        if colorLabel == 'gLikeFig6A':
            sampleRGBs.append([ancPropBySampleDict[sample]['pol'], ancPropBySampleDict[sample]['afr'],
                               ancPropBySampleDict[sample]['eur']])
    sampleRGBs = np.array(sampleRGBs)
    return sampleRGBs


def defineSampleShape(pops, sampleSizeByPopDict):
    # letters = list(string.ascii_lowercase)# can not be colored
    letters = ['o', 'v', '^', '<', '>',
               '1', '2', '3', '4', '8',
               's', 'p', 'P', '*', 'h',
               'H', '+', 'x', 'X', 'D',
               'd', '|', '_', '$...$', 'o']
    assert len(pops) <= len(letters)
    sampleShapes = []
    for i, pop in enumerate(pops):
        # sampleShapes += [f'${i}$'] * sampleSizeByPopDict[pop] # hard to tell 1 and 2 and 12 # can not be colored
        # sampleShapes += [f'${letters[i]}$'] * sampleSizeByPopDict[pop] # can not be colored
        sampleShapes += [f'{letters[i]}'] * sampleSizeByPopDict[pop]
    return sampleShapes


def plot_pca_umap(PCs, true_labels, pcaFig, umapFigPath, pcaFigTitle, umapFigTitle, sampleRGBs, sampleShapes,
                  pca_on_dip=False, plot4pub=False):
    def Plot(outfig, components, figTitle=None, axisName=None, textBoxAnnot=None, sampleRGBs=None,
             sampleShapes=None, edgeColors=None, markerSizes=30, no_axis_label=False):
        component1 = components[:, 0]
        component2 = components[:, 1]
        sampleSize = len(component1)

        xLabel, yLabel = None, None
        if axisName is not None:
            xLabel, yLabel = axisName + '1', axisName + '2'

        if sampleRGBs is None:
            sampleRGBs = np.array([[30, 144, 255]] * sampleSize) / 255
        if sampleShapes is None:
            sampleShapes = ['o'] * sampleSize
        if markerSizes is None:
            markerSizes = [3] * sampleSize
        if type(markerSizes) == int or type(markerSizes) == float:
            markerSizes = [markerSizes] * sampleSize

        # plot
        plt.close()
        # figsize = [6.4, 4.8]
        figsize = [6, 4]
        fig, ax = plt.subplots(dpi=500, figsize=figsize)
        for i in range(sampleSize):
            plt.scatter(x=component1[i], y=component2[i],
                        c=sampleRGBs[i].reshape(1, 3),
                        marker=sampleShapes[i], edgecolors=edgeColors, s=markerSizes[i])

        if textBoxAnnot is not None:
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            ax.text(0.02, 0.98, textBoxAnnot, transform=ax.transAxes,
                    fontsize=12, weight="bold",
                    verticalalignment='top', bbox=props)
        if no_axis_label:
            ax.set_xticks([])
            ax.set_yticks([])
        plt.legend('', frameon=False)
        if xLabel:
            plt.xlabel(xLabel)
        if yLabel:
            plt.ylabel(yLabel)
        if figTitle:
            plt.title(figTitle)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        print('Plotting to ' + outfig)
        plt.savefig(outfig)
        plt.close()

    pcSI = round(separation_index(true_labels, PCs[:, :2]).mean(), 2)
    textBoxAnnot = f'SI:{pcSI}'
    if pca_on_dip:
        pcaFig = re.sub('.png', '_dip.png', pcaFig)
    axisName = 'PC'
    if plot4pub:
        axisName = None
        pcaFigTitle = None
    edgeColors = None
    markerSizes = 30
    Plot(pcaFig, axisName=axisName, components=PCs, sampleRGBs=sampleRGBs, sampleShapes=sampleShapes,
         # explainedVars=exp_var_pca,
         figTitle=pcaFigTitle, edgeColors=edgeColors, textBoxAnnot=textBoxAnnot, markerSizes=markerSizes)
    os.makedirs(umapFigPath, exist_ok=True)
    umapSIs = []
    umap_metrics = ['euclidean']
    umap_minDists = [0.1]
    umap_neighbors = [15]
    umap_pcNums = [20, 50]
    for umap_PCnum in umap_pcNums:
        # PCs = PCs[:, :umap_PCnum]
        PCs_umap = PCs[:, :umap_PCnum]
        for umap_metric in umap_metrics:
            for umap_minDist in umap_minDists:
                for umap_neighbor in umap_neighbors:
                    figNameLabel = f'metric-{umap_metric}_neighbors{umap_neighbor}_pc{umap_PCnum}_minDist{umap_minDist}'
                    umapFig = f'{umapFigPath}/{figNameLabel}.png'
                    if pca_on_dip:
                        umapFig = re.sub('.png', '_dip.png', umapFig)
                    umapArray = umap.UMAP(  # random_state=1,
                        n_components=3, min_dist=umap_minDist, metric=umap_metric,
                        n_neighbors=umap_neighbor).fit_transform(PCs_umap)
                    umapSI = round(separation_index(true_labels, umapArray[:, :2]).mean(), 2)
                    umapSIs.append(umapSI)
                    textBoxAnnot = f'SI:{umapSI}'
                    axisName = 'UMAP'
                    if plot4pub:
                        umapFigTitle = None
                        axisName = None
                    Plot(umapFig, axisName=axisName, components=umapArray, sampleRGBs=sampleRGBs,
                         sampleShapes=sampleShapes,
                         figTitle=umapFigTitle, edgeColors=edgeColors, textBoxAnnot=textBoxAnnot,
                         markerSizes=markerSizes, no_axis_label=True)


### simulation
def asegrm_on_sim(simScheme, task, taskIdx):
    assert task in ['sim', 'calc', 'plot']
    overwrite = False
    asegrm_pca = 'pca'  # standard PCA

    LACname = 'rfmixLAC'
    ARGname = 'relateARG'

    sampleAnc = True
    plot4pub = True

    methodName = 'as-eGRM'

    chrNums = list(range(1, 11))
    chrNums_plot = list(range(1, 11))

    specifiedSampleSize = 500
    ploidy_sim = 2
    pca_on_dip = True

    tasks = []
    chrN = None
    # if not returnArgs:
    if task == 'calc':
        # for methodName in methodNames:
        # for LACname, ARGname in LAC_ARG_list:
        # for simScheme in simSchemes:
        argsList = schemesForSim(simScheme, sampleAnc=sampleAnc, sample_size=specifiedSampleSize)
        for args in argsList:
            for chrN in chrNums:
                tasks.append([methodName, simScheme, args, chrN, LACname, ARGname])
        method, simScheme, args, chrN, LACname, ARGname = tasks[taskIdx]
    else:
        # for methodName in methodNames:
        # for LACname, ARGname in LAC_ARG_list:
        # for simScheme in simSchemes:
        argsList = schemesForSim(simScheme, sampleAnc=sampleAnc, sample_size=specifiedSampleSize)
        for args in argsList:
            tasks.append([methodName, simScheme, args, LACname, ARGname])
        method, simScheme, args, LACname, ARGname = tasks[taskIdx]
    # pdb.set_trace()

    args.update({'simScheme': simScheme})
    schemePath = f'{os.getcwd()}/simulation/{simScheme}'
    os.makedirs(schemePath, exist_ok=True)

    args['schemePath'] = schemePath

    # record_migrations = False if LACname is True else True
    record_migrations = True
    colorLabel = '9bins'
    ancPops, tgtAnc = [], ''
    if 'gLikeFig6A' in simScheme:
        LACname, ARGname = 'rfmixLAC', 'relateARG'
        tgtAnc = 'pol'
        ancPops = ['pol', 'afr', 'eur']
        colorLabel = 'gLikeFig6A'
    if 'twoAnc' in simScheme:
        tgtAnc = 'a1'
        ancPops = ['a1', 'a2']

    if 'GridDer' in simScheme:
        colorBy = 'pop'
    else:
        colorBy = 'ancProp'

    true_labels = []
    if 'twoAncTwoDer' in simScheme or 'gLikeFig6A' in simScheme:
        true_labels = np.repeat([0, 1], [args['sample_size'], args['sample_size']])
    if 'GridDer' in simScheme:
        ns = [args['sample_size']] * args['n_rows'] * args['n_cols']
        ns_array = np.array(ns)
        true_labels = np.repeat(np.arange(len(ns)), ns_array)

    label = ''
    if 'gLikeFig6A' in simScheme:
        label = f"tSplit-{args['t_split']}_tAdmix{args['t_admix']}_lm{args['latino_migration_rate']}-{args['latino_migration_stop_t']}_" \
                f"d1AncProp{args['d1_anc_prop'][0]}-{args['d1_anc_prop'][1]}_d2AncProp{args['d2_anc_prop'][0]}-{args['d2_anc_prop'][1]}"
    else:
        if 'twoAncTwoDer' in simScheme:
            label = f"tSplit1-{args['t_split_1']}_tSplit2-{args['t_split_2']}_tAdmix{args['t_admix']}_m1-{args['a2_to_d1_prop']}-m2-{args['a2_to_d2_prop']}"
        if 'GridDer' in simScheme:
            label = f"tSplit1-{args['t_split_1']}_tSplit2-{args['t_split_2']}_tAdmix{args['t_admix']}_propMat-{args['prop_from_a2_mat_label']}"
        if args['migration_rate'] > 0:
            label += f"_mig{args['migration_rate']}-{args['migration_stop_t']}"

    chrLength_all = args['chrLength'] * len(chrNums)
    label += f"_{int(chrLength_all / 1000000)}Mb"

    genotypesPath = f'{schemePath}/genotypes'
    os.makedirs(genotypesPath, exist_ok=True)
    rfmix_input_output_path = f"{args['schemePath']}/rfmix_input_output"
    anc_prop_path = f"{schemePath}/anc_prop_{LACname}"
    os.makedirs(anc_prop_path, exist_ok=True)
    anc_prop_file_forPlot = f"{anc_prop_path}/{label}_chr{chrNums_plot[0]}-chr{chrNums_plot[-1]}.p"

    samplesFile_query = f'{schemePath}/samples_query.txt'
    samplesFile_ref = f'{schemePath}/samples_ref.txt'
    samplesFile_query_dip = f'{schemePath}/samples_query_diploid.txt'
    samplesFile_ref_dip = f'{schemePath}/samples_ref_diploid.txt'

    treeSeqPath = f'{schemePath}/tree_sequence'
    os.makedirs(treeSeqPath, exist_ok=True)
    treeSeqFile_all = f'{treeSeqPath}/{label}'
    treeSeqFile_all += '_record_migrations'

    LAClabel = LACname
    ARGlabel = ARGname
    methodPath = f"{schemePath}/{method}_{ARGlabel}_{LAClabel}"
    os.makedirs(methodPath, exist_ok=True)
    method_output_path = f'{methodPath}/{method}_output'
    os.makedirs(method_output_path, exist_ok=True)
    demography, popIDbyNameDict, sampleIDbyPopNameDict = defineDemog(args, sampleAnc=sampleAnc)
    args.update({'popIDbyNameDict': popIDbyNameDict, 'sampleIDbyPopNameDict': sampleIDbyPopNameDict})
    sampleSizeByPopDict_ = args['sampleSizeByPopDict']
    if ploidy_sim == 2:
        sampleSizeByPopDict_ = args['sampleSizeByPopDict_dip']

    method_output_path_label = f"{method_output_path}/{label}"
    os.makedirs(method_output_path_label, exist_ok=True)

    # generate samples
    querySamples = get_samples_by_pop(samplesFile_query, args['rfmix_query_pops'], sampleIDbyPopNameDict,
                                      diploid=False)
    refSamples = get_samples_by_pop(samplesFile_ref, args['rfmix_ref_pops'], sampleIDbyPopNameDict,
                                    diploid=False)
    querySamples_dip = get_samples_by_pop(samplesFile_query_dip, args['rfmix_query_pops'],
                                          sampleIDbyPopNameDict, diploid=True)
    refSamples_dip = get_samples_by_pop(samplesFile_ref_dip, args['rfmix_ref_pops'], sampleIDbyPopNameDict,
                                        diploid=True)

    popNameBySampleDict = {}
    for i in querySamples:
        for popName, (start, end) in sampleIDbyPopNameDict.items():
            if start <= i < end:
                popNameBySampleDict[i] = popName
                break
    if task == 'sim':
        if not os.path.exists(treeSeqFile_all):
            print('simulation starts ...')
            trees = msprime.sim_ancestry(samples=sampleSizeByPopDict_,
                                         sequence_length=chrLength_all,
                                         recombination_rate=args['recombination_rate'],
                                         demography=demography,
                                         ploidy=ploidy_sim,
                                         record_migrations=record_migrations
                                         )
            trees = msprime.sim_mutations(trees, rate=args['mutation_rate'], discrete_genome=False)
            trees.dump(treeSeqFile_all)

    elif task in ['calc']:
        label_chr = f'{label}_chr{chrN}'
        treeSeqFile_chr = treeSeqFile_all + '_chr' + str(chrN)
        if not os.path.exists(treeSeqFile_chr):
            trees = tskit.load(treeSeqFile_all)
            chunks = []
            steps = [i for i in range(0, chrLength_all + 1, args['chrLength'])]
            for chunk_start, chunk_end in zip(steps[:-1], steps[1:]):
                chunks.append([chunk_start, chunk_end])
            chunk = chunks[chrN - 1]
            simplify = False if record_migrations else True
            trees = trees.keep_intervals(np.array([chunk]), simplify=simplify)
            # todo: the interval of the last tree in the subset tree sequence is not correct,
            #  for example, after applying chunk=[0,10000000] to the tree sequence with a interval=[0,50000000], the interval of the last tree
            #  in the subset tree sequence is [10000000, 50000000], but it should be not larger than 10000000.
            #  Therefore, it is better to drop the last tree in the subset tree sequence.
            trees.dump(treeSeqFile_chr)
        else:
            print('loading existing tree sequence ...')
            trees = tskit.load(treeSeqFile_chr)

        # generate genotypes of vcf format
        genotypes_query_file = f'{genotypesPath}/genotypes_query_{label_chr}.vcf.gz'
        genotypes_ref_file = f'{genotypesPath}/genotypes_ref_{label_chr}.vcf.gz'
        genotypes_all_file = f'{genotypesPath}/genotypes_all_{label_chr}.vcf.gz'
        if LACname != 'trueLAC' or ARGname != 'trueARG':
            if not os.path.exists(genotypes_all_file):
                with open(genotypes_all_file.rstrip('.gz'), 'w') as fw:
                    trees.write_vcf(fw)
                exeCmd(f"bgzip {genotypes_all_file.rstrip('.gz')}")
                # sorting positions and removing duplicated postions are also necessary for the following sub-setting
                sort_positions(genotypes_all_file)
                rm_dup_positions(genotypes_all_file)
            if not os.path.exists(genotypes_query_file):
                querySamples_ = querySamples_dip if ploidy_sim == 2 else querySamples
                querySamples_ = [f'tsk_{i}' for i in querySamples_]
                subset_by_samples(genotypes_all_file, genotypes_query_file, querySamples_, rm_mono=True)
            if not os.path.exists(genotypes_ref_file):
                refSamples_ = refSamples_dip if ploidy_sim == 2 else refSamples
                refSamples_ = [f'tsk_{i}' for i in refSamples_]
                subset_by_samples(genotypes_all_file, genotypes_ref_file, refSamples_, rm_mono=True)

        # infer local ancestry calls
        rfmix_msp_file = f"{rfmix_input_output_path}/rfmixOut_{label_chr}.msp.tsv"
        if not os.path.exists(rfmix_msp_file):
            sampleIDbyPopNameDict_ = copy.deepcopy(args['sampleIDbyPopNameDict'])
            for pop in sampleIDbyPopNameDict_:
                sampleIDbyPopNameDict_[pop] = (
                    int(sampleIDbyPopNameDict_[pop][0] / 2), int(sampleIDbyPopNameDict_[pop][1] / 2))
            chrN_ = '1' if trees is not None else chrN
            run_rfmix_wrapper(rfmix_input_output_path=rfmix_input_output_path,
                              query_pops=args['rfmix_query_pops'],
                              ref_pops=args['rfmix_ref_pops'],
                              sampleIDbyPopNameDict=sampleIDbyPopNameDict_,
                              trees=trees,
                              recombRate=args['recombination_rate'],
                              chrN=chrN_,
                              label=label_chr,
                              genotypes_query_file=genotypes_query_file,
                              genotypes_ref_file=genotypes_ref_file
                              )

        # compute ancestry proportions
        ancPropsFile = f"{anc_prop_path}/{label_chr}.p"
        if not os.path.exists(ancPropsFile):
            local_anc_calls, anc_id_by_name_dict = [], {}
            with open(rfmix_msp_file, 'r') as fr:
                rows = fr.readlines()
                for i in re.sub('#Subpopulation order/codes: ', '', rows[0].rstrip('\n')).split('\t'):
                    anc_name, anc_id = i.split('=')
                    anc_id_by_name_dict[anc_name] = int(anc_id)
                for row in rows[2:]:
                    cols = row.rstrip('\n').split('\t')
                    segStart, segEnd = cols[1:3]
                    segStart, segEnd = int(segStart), int(segEnd)
                    localAncBySampleOrderArray = [int(i) for i in cols[6:]]
                    local_anc_calls.append([segStart, segEnd, localAncBySampleOrderArray])
            local_anc_calls.sort(key=lambda x: x[0], reverse=False)
            ancPropBySampleDict = CalcAncProp(simScheme, args['chrLength'], samples=querySamples,
                                              local_anc_calls=local_anc_calls,
                                              anc_id_by_name_dict=anc_id_by_name_dict
                                              )
            with open(ancPropsFile, 'wb') as fw:
                pickle.dump(ancPropBySampleDict, fw)

        # run ARG method
        # if ARGname != 'trueARG':
        # specify the directory for running ARG method and the file to save trees
        ARGmethod = ARGname.rstrip('ARG')
        ARGmethod_path = f"{schemePath}/{ARGmethod}_input_output"
        hapMat, positions, positionIDs = None, None, None
        genotypes_query_file_ = genotypes_query_file
        treesFile = f"{ARGmethod_path}/{label_chr}.{ARGmethod}.trees"
        # treesFile = f"{ARGmethod_path}/relateOut_{label_chr}.trees" # tmp
        if not os.path.exists(treesFile):
            os.makedirs(ARGmethod_path, exist_ok=True)
            if ARGname == 'relateARG':
                run_relate_wrapper(relate_input_output_path=ARGmethod_path,
                                   recombRate=args['recombination_rate'],
                                   Ne=args['Ne'],
                                   mutRate=args['mutation_rate'],
                                   fileName=label_chr,
                                   treesFile=treesFile,
                                   vcfFile=genotypes_query_file_,
                                   positions=positions,
                                   )

        # run asegrm
        asegrm_file = f'{method_output_path_label}/{label_chr}.trees.{tgtAnc}.asegrm.npy'
        asegrm_mu_file = f'{method_output_path_label}/{label_chr}.trees.{tgtAnc}.asegrm.mu.npy'
        if not os.path.exists(asegrm_file) and not os.path.exists(asegrm_mu_file):
            # prepare leaf ids and genetic map
            samples, positions = read_samples_positions(genotypes_query_file_, read_samples=True, read_positions=True)
            leaf_ids_file = f"{method_output_path_label}/leaf_ids_chr{chrN}.txt"
            genetic_map_file = f"{method_output_path_label}/genetic_map_chr{chrN}.txt"
            recombRate = args['recombination_rate']
            pos_map_pairs = calc_genetic_map_for_sim_data(positions, recombRate)
            with open(genetic_map_file, "w") as fw:
                fw.write("pos COMBINED_rate Genetic_Map\n")
                for pos, map in pos_map_pairs:
                    string = str(pos) + " " + str(recombRate) + " "
                    string = string + str(map) + "\n"
                    fw.write(string)
            with open(leaf_ids_file, 'w') as fw:
                content = []
                for i in samples:
                    content.append(f'{i}.0')
                    content.append(f'{i}.1')
                fw.write('\n'.join(content))

            cmd = 'asegrm compute '
            cmd += f' --trees {treesFile}'
            cmd += f' --leaf_ids {leaf_ids_file}'
            cmd += f' --local_ancestry {rfmix_msp_file}'
            cmd += f' --target_ancestry {tgtAnc}'
            cmd += f' --genetic_map {genetic_map_file}'
            cmd += f' --output_path {method_output_path_label}'
            exeCmd(cmd)

    elif task == 'plot':
        plotPath = f"{methodPath}/plot_chr{chrNums_plot[0]}-chr{chrNums_plot[-1]}"
        os.makedirs(plotPath, exist_ok=True)
        with open(samplesFile_query, 'r') as f:
            samples = eval(f.read())

        # generate ancestry proportions on multiple chromosomes
        ancPropBySampleDict_chrs = {}
        if colorBy == 'ancProp':
            if overwrite or (not os.path.exists(anc_prop_file_forPlot)):
                ancPropBySampleDict_chrs = []
                for chrN in chrNums_plot:
                    ancPropsFile = f"{anc_prop_path}/{label}_chr{chrN}.p"
                    with open(ancPropsFile, "rb") as f:
                        ancPropBySampleDict_chrs.append(pickle.load(f))
                ancPropBySampleDict_chrs = average_anc_props(ancPropBySampleDict_chrs, samples, ancPops)
                with open(anc_prop_file_forPlot, "wb") as f:
                    pickle.dump(ancPropBySampleDict_chrs, f)
            else:
                with open(anc_prop_file_forPlot, "rb") as f:
                    ancPropBySampleDict_chrs = pickle.load(f)

        # generate method result on multiple chromosomes
        if pca_on_dip:
            methodResultFile = f"{method_output_path_label}/merged.{tgtAnc}.asegrm.diploid.npy"
        else:
            methodResultFile = f"{method_output_path_label}/merged.{tgtAnc}.asegrm.haploid.npy"
        if overwrite or (not os.path.exists(methodResultFile)):
            cmd = f'asegrm merge {method_output_path_label}'
            exeCmd(cmd)
        K = np.load(methodResultFile)

        print('computing PCA-UMAP and plotting ...')
        # specify samplesShapes and sampleColors
        sampleShapes, sampleRGBs = [], []
        if colorBy == 'ancProp':
            sampleRGBs = colorSampleByAncProp(ancPropBySampleDict_chrs, samples, colorLabel)
            sampleShapes = defineSampleShape(args['rfmix_query_pops'], args['sampleSizeByPopDict'])
        if colorBy == 'pop':
            for i, pop in enumerate(args['rfmix_query_pops']):
                sampleShapes += ['o'] * args['sampleSizeByPopDict'][pop]
            sampleRGBs = colorSampleByPop(args['rfmix_query_pops'], args['sampleSizeByPopDict'])

        if pca_on_dip:
            maternals = np.array(range(0, len(sampleRGBs), 2))
            sampleRGBs = sampleRGBs[maternals]
            sampleShapes = [sampleShapes[i] for i in maternals]
            true_labels = true_labels[maternals]
        PCs = eigendecomposition_on_grm(K, num_dims=100)
        pcaFig = f'{plotPath}/pcaPlot'
        pcaFig += f'_{label}.png'
        methodLabel = f'{method}_tgtAnc-{tgtAnc}' if 'as-' in method else method
        pcaFigTitle = f'{simScheme}+{methodLabel}+PCA'
        umapFigPath = f'{plotPath}/{asegrm_pca}-umapPlot_{label}'
        plot_pca_umap(PCs, true_labels, pcaFig, umapFigPath, pcaFigTitle=pcaFigTitle,
                      umapFigTitle=f'{pcaFigTitle}+UMAP', sampleRGBs=sampleRGBs, sampleShapes=sampleShapes,
                      pca_on_dip=pca_on_dip, plot4pub=plot4pub)


if __name__ == '__main__':
    if sys.argv[2] == 'calc':
        PATH_TO_RELATE = sys.argv[4]  # '/project/jitang_1167/software/relate_v1.2.0_x86_64_static'
        PATH_TO_RFMix = sys.argv[5]  # '/project/chia657_28/programs/rfmix/rfmix '
    asegrm_on_sim(simScheme=sys.argv[1], task=sys.argv[2], taskIdx=int(sys.argv[3]))
