import os


output_path = './output'
cmd = 'asegrm compute '
cmd += ' --trees chr22.part-02.relate.trees'
cmd += ' --leaf_ids leaf_ids.txt'
cmd += ' --local_ancestry chr22.msp.tsv'
cmd += ' --target_ancestry amr'
cmd += ' --genetic_map chr22.map'
cmd += f' --output_path {output_path}'
os.system(cmd)

cmd = f'asegrm merge {output_path}'
os.system(cmd)