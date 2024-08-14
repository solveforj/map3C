from pathlib import Path
from glob import glob
import os
import yaml
import shutil

class PrepareMapping:

    def prepare_plates(self):

        barcodes = self.config_dict["general"]["snm3Cseq_plate_input"]["barcodes"]
        plate_info = self.config_dict["general"]["snm3Cseq_plate_input"]["plate_info"]
        output_directory = self.config_dict["general"]["output_directory"]

        # Results directory
        results_directory = os.path.join(output_directory, "results")
        Path(results_directory).mkdir(parents=True, exist_ok=True)

        # Snakemake directory
        snakemake_path = os.path.join(Path(__file__).parent.resolve(), "snakemake")
        snakemake_directory = os.path.join(output_directory, "snakemake_mapping")
        shutil.copytree(snakemake_path, snakemake_directory, dirs_exist_ok=True)

        cell_ids = []
        with open(barcodes) as f:
            for line in f:
                line = line.strip()
                if line[0] == ">":
                    cell_ids.append(line[1:])
        
        with open(plate_info) as f:
            for line in f:
                line = line.strip()
                if len(line) == 0:
                    continue
                if line[0] == "#":
                    continue
                line = line.split()
                plate, fastq_directory = line[0], line[1]
                
                plate_run_directory = os.path.join(results_directory, plate)
                Path(plate_run_directory).mkdir(parents=True, exist_ok=True)
    
                with open(f"{plate_run_directory}/mapping_scripts.txt", 'w') as scripts:
                    
                    for cell in cell_ids:
                        cell_local_directory = f"{plate}_{cell}"
                        cell_run_directory = os.path.join(plate_run_directory, 
                                                          cell_local_directory)
                        Path(cell_run_directory).mkdir(parents=True, exist_ok=True)

                        r1_fastq = os.path.join(cell_run_directory, f"{cell_local_directory}_indexed_R1.fastq.gz")
                        r2_fastq = os.path.join(cell_run_directory, f"{cell_local_directory}_indexed_R2.fastq.gz")
                        
                        run_config = f"{cell_run_directory}/run_config.csv"
                        snakemake_cmd = f"{cell_run_directory}/mapping_cmd.txt"
                        scripts.write(snakemake_cmd + "\n")
                        
                        with open(snakemake_cmd, 'w') as f:
                            
                            cmd = f"snakemake -d {cell_run_directory} --profile none --workflow-profile none --sdm conda "
                            cmd += f"--snakefile {snakemake_directory}/mapping.smk "
                            cmd += f"--configfile {self.config} {self.snakemake_params}"
                            f.write(cmd + '\n')

                        with open(run_config, 'w') as f:
                            
                            f.write(",".join(["r1", "r2"]) + '\n')
                            f.write(",".join([cell_local_directory, r1_fastq, r2_fastq]) + '\n')
                    
    
    def prepare_ids(self):
        
        output_directory = self.config_dict["general"]["output_directory"]
        fastq_info = self.config_dict["general"]["general_input"]["fastq_info"]

        # Results directory
        results_directory = os.path.join(output_directory, "results")
        Path(results_directory).mkdir(parents=True, exist_ok=True)

        # Snakemake directory
        snakemake_path = os.path.join(Path(__file__).parent.resolve(), "snakemake")
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

        if ((self.config_dict["general"]["snm3Cseq_plate_input"]["plate_info"] or 
            self.config_dict["general"]["snm3Cseq_plate_input"]["barcodes"]) and
            self.config_dict["general"]["general_input"]["fastq_info"]):

            raise Exception("Cannot have both snm3C-seq plate input and general input defined")
            
        if (not self.config_dict["general"]["general_input"]["fastq_info"] and not
                (self.config_dict["general"]["snm3Cseq_plate_input"]["plate_info"] or
                     self.config_dict["general"]["snm3Cseq_plate_input"]["barcodes"])):

            raise Exception("Must have either snm3C-seq plate input or general input defined")

        # Plate-level data (demultiplexed by this package)
        if (self.config_dict["general"]["snm3Cseq_plate_input"]["plate_info"] and 
            self.config_dict["general"]["snm3Cseq_plate_input"]["barcodes"]):
            self.prepare_plates()

        # Id-level data  (already demultiplexed)
        elif self.config_dict["general"]["general_input"]["fastq_info"]:
            self.prepare_ids()
