from pathlib import Path
import os
import yaml
import shutil

class PrepareDemultiplex:

    def prepare_plates(self):

        plate_info = self.config_dict["general"]["snm3Cseq_plate_input"]["plate_info"]
        output_directory = self.config_dict["general"]["output_directory"]

        # Results directory
        results_directory = os.path.join(output_directory, "results")
        Path(results_directory).mkdir(parents=True, exist_ok=True)

        # Snakemake directory
        snakemake_path = os.path.join(Path(__file__).parent.resolve(), "snakemake")
        snakemake_directory = os.path.join(output_directory, "snakemake_demultiplex")
        shutil.copytree(snakemake_path, snakemake_directory, dirs_exist_ok=True)
                    
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

                run_config = f"{plate_run_directory}/run_config.csv"
                
                with open(run_config, "w") as f:
                    f.write(",".join(["fastq_dir"]) + "\n")
                    f.write(",".join([plate, fastq_directory]) + "\n")

                with open(f"{plate_run_directory}/demultiplex_cmd.txt", 'w') as outfile:
                    cmd = f"snakemake -d {plate_run_directory} --profile none --workflow-profile none --sdm conda "
                    cmd += f"--snakefile {snakemake_directory}/demultiplex.smk "
                    cmd += f"--configfile {self.config} {self.snakemake_params} "
                    outfile.write(cmd + '\n')

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
            
        if not (self.config_dict["general"]["snm3Cseq_plate_input"]["plate_info"] and 
            self.config_dict["general"]["snm3Cseq_plate_input"]["barcodes"]):

            raise Exception("Must have snm3C-seq plate info and barcodes defined")
            
        # Plate-level data (demultiplexed by this package)
        if (self.config_dict["general"]["snm3Cseq_plate_input"]["plate_info"] and 
            self.config_dict["general"]["snm3Cseq_plate_input"]["barcodes"]):
            
            self.prepare_plates()
