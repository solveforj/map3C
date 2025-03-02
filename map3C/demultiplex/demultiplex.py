from pathlib import Path
import os
import yaml
import shutil

class PrepareDemultiplex:

    def prepare_plates_snm3Cseq(self):

        plate_info = self.config_dict["general"]["fastq_info"]
        output_directory = self.config_dict["general"]["output_directory"]

        # Results directory
        results_directory = os.path.join(output_directory, "results")
        Path(results_directory).mkdir(parents=True, exist_ok=True)
        all_cmds = os.path.join(results_directory, "demultiplex_scripts.txt")
        
        # Snakemake directory
        snakemake_path = os.path.join(Path(__file__).parent.resolve(), "snakemake", "snm3Cseq")
        snakemake_directory = os.path.join(output_directory, "snakemake_demultiplex")
        shutil.copytree(snakemake_path, snakemake_directory, dirs_exist_ok=True)
                    
        with open(plate_info) as f, open(all_cmds, "w") as acmd:
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

                acmd.write(f"{plate_run_directory}/demultiplex_cmd.txt" + "\n")


    def prepare_plates_GAGEseq(self):

        plate_info = self.config_dict["general"]["fastq_info"]
        output_directory = self.config_dict["general"]["output_directory"]

        bc2 = self.config_dict["demultiplex_protocols"]["GAGE-seq"]["BC2"]["barcodes"]
        bc2_ids = []
        
        with open(bc2) as f:
            for line in f:
                line = line.strip()
                if line[0] == ">":
                    bc2_ids.append(line[1:])


        # Results directory
        results_directory = os.path.join(output_directory, "results")
        Path(results_directory).mkdir(parents=True, exist_ok=True)
        all_bc1_cmds = os.path.join(results_directory, "demultiplex_BC1_scripts.txt")
        all_bc2_cmds = os.path.join(results_directory, "demultiplex_BC2_scripts.txt")
        
        # Snakemake directory
        snakemake_path = os.path.join(Path(__file__).parent.resolve(), "snakemake", "GAGE-seq")
        snakemake_directory = os.path.join(output_directory, "snakemake_demultiplex")
        shutil.copytree(snakemake_path, snakemake_directory, dirs_exist_ok=True)
                    
        with open(plate_info) as f, open(all_bc1_cmds, "w") as abc1, open(all_bc2_cmds, "w") as abc2:
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

                with open(f"{plate_run_directory}/demultiplex_cmd_BC2.txt", 'w') as outfile:
                    cmd = f"snakemake -d {plate_run_directory} --profile none --workflow-profile none --sdm conda "
                    cmd += f"--snakefile {snakemake_directory}/demultiplex_BC2.smk "
                    cmd += f"--configfile {self.config} {self.snakemake_params} "
                    outfile.write(cmd + '\n')

                abc2.write(f"{plate_run_directory}/demultiplex_cmd_BC2.txt" + "\n")

                for id in bc2_ids:
                    plate_bc2_directory = os.path.join(results_directory, plate, f"{plate}-BC2-{id}")
                    Path(plate_bc2_directory).mkdir(parents=True, exist_ok=True)

                    run_config_id = f"{plate_bc2_directory}/run_config.csv"
                    
                    with open(run_config_id, "w") as f:
                        f.write(",".join(["R1", "R2"]) + "\n")
                        R1 = os.path.join(plate_bc2_directory, f"{plate}-BC2-{id}_R1.fastq.gz")
                        R2 = os.path.join(plate_bc2_directory, f"{plate}-BC2-{id}_R2.fastq.gz")
                        f.write(",".join([f"{plate}-BC2-{id}", R1, R2]) + "\n")

                    with open(f"{plate_bc2_directory}/demultiplex_cmd_BC1.txt", 'w') as outfile:
                        cmd = f"snakemake -d {plate_bc2_directory} --profile none --workflow-profile none --sdm conda "
                        cmd += f"--snakefile {snakemake_directory}/demultiplex_BC1.smk "
                        cmd += f"--configfile {self.config} {self.snakemake_params} "
                        outfile.write(cmd + '\n')

                    abc1.write(f"{plate_bc2_directory}/demultiplex_cmd_BC1.txt\n")
                
                    
    def __init__(self,
                 config,
                 snakemake_params
                ):
        
        with open(config) as f:
            self.config_dict = yaml.safe_load(f)

        self.config = config
        
        self.mode = self.config_dict["general"]["demultiplex_protocol"]
        
        self.snakemake_params = snakemake_params

        if not self.config_dict["general"]["fastq_info"]:

            raise Exception("FASTQ info must be defined")

        if self.mode in ["snm3Cseq", "snmCTseq"]:
            
            self.prepare_plates_snm3Cseq()

        elif self.mode == "GAGE-seq":

            self.prepare_plates_GAGEseq()