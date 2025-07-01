#!/bin/bash
# Run each script in the background
sbatch leafcutter_step_6_DPseq_sACC.sh
sbatch leafcutter_step_6_BPseq_sACC.sh
sbatch leafcutter_step_6_DPseq_Amygdala.sh
sbatch leafcutter_step_6_BPseq_Amygdala.sh
