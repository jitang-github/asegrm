import os


output_path = './output'
for chrN in ['1', '2']:
    cmd = 'asegrm compute '
    cmd += f' --trees chr{chrN}.trees'
    cmd += ' --leaf_ids leaf_ids.txt'
    cmd += f' --local_ancestry chr{chrN}.msp.tsv'
    cmd += ' --target_ancestry a1'
    cmd += f' --genetic_map chr{chrN}.map'
    cmd += f' --output_path {output_path}'
    print(cmd)
    os.system(cmd)

cmd = f'asegrm merge {output_path}'
os.system(cmd)