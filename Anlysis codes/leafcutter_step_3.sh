#!/bin/bash
# Run each script in the background
sbatch leafcutter_step_3_Amygdala_BPseq.sh
sbatch leafcutter_step_3_sACC_BPseq.sh
sbatch leafcutter_step_3_Amygdala_DPseq.sh
sbatch leafcutter_step_3_sACC_DPseq.sh

