from pathlib import Path
from glob import glob
import os
import yaml
import shutil

class PrepareMapping:
    
    def prepare_ids(self, mode):
        
        output_directory = self.config_dict["general"]["output_directory"]
        fastq_info = self.config_dict["general"]["fastq_info"]

        # Results directory
        results_directory = os.path.join(output_directory, "results")
        Path(results_directory).mkdir(parents=True, exist_ok=True)

        # Snakemake directory
        snakemake_path = os.path.join(Path(__file__).parent.resolve(), "snakemake", mode)
        snakemake_directory = os.path.join(output_directory, "snakemake_mapping")
        shutil.copytree(snakemake_path, snakemake_directory, dirs_exist_ok=True)
        
        with open(f"{results_directory}/mapping_scripts.txt", 'w') as scripts, \
            open(fastq_info) as f:
    
            for line in f:
                line = line.strip()
                if len(line) == 0:
                    continue
                if line[0] == "#":
                    continue
                line = line.split()
                id, r1_fastq, r2_fastq = line[0], line[1], line[2]
                
                id_run_directory = os.path.join(results_directory, id)
                id_run_directory_path = Path(id_run_directory)
                id_run_directory_path.mkdir(parents=True, exist_ok=True)

                run_config = f"{id_run_directory}/run_config.csv"
                snakemake_cmd = f"{id_run_directory}/mapping_cmd.txt"
                scripts.write(snakemake_cmd + "\n")
                        
                with open(snakemake_cmd, 'w') as f:
                    
                    cmd = f"snakemake -d {id_run_directory} --profile none --workflow-profile none --sdm conda "
                    cmd += f"--snakefile {snakemake_directory}/mapping.smk "
                    cmd += f"--configfile {self.config} {self.snakemake_params}"
                    f.write(cmd + '\n')

                
                with open(run_config, 'w') as f:
                    
                    f.write(",".join(["r1", "r2"]) + '\n')
                    f.write(",".join([id, r1_fastq, r2_fastq]) + '\n')

    def __init__(self,
                 config,
                 snakemake_params
                ):
        
        with open(config) as f:
            self.config_dict = yaml.safe_load(f)

        self.config = config
        
        self.mode = self.config_dict["general"]["mode"]
        
        self.snakemake_params = snakemake_params
        
        if not self.config_dict["general"]["fastq_info"]:

            raise Exception("FASTQ info must be defined")

        if self.mode in ["bsdna", "dna", "snm3Cseq"]:
            interpret_mode = "m3C"
        elif self.mode in ["snmCTseq"]:
            interpret_mode = "mCT"
        self.prepare_ids(interpret_mode)
