#!/usr/bin/env bash

#SBATCH --time=24:00:00
#SBATCH --error=./log/err
#SBATCH --output=./log/out

source /scratch/groups/hheyn/software/anaconda3/bin/activate scrublet
python 01-doublet_detection_scrublet.py BCLLATLAS_33 ifZOgenn_TpMNTvBa;
python 01-doublet_detection_scrublet.py BCLLATLAS_33 mLuLpVxi_v0fLyotc;
python 01-doublet_detection_scrublet.py BCLLATLAS_38 B20O1bh7_VmM99YZJ;
python 01-doublet_detection_scrublet.py BCLLATLAS_38 rdFRFhrU_ZdYeOZlf;
python 01-doublet_detection_scrublet.py BCLLATLAS_38 WToIzInl_LudU7hVX;
python 01-doublet_detection_scrublet.py BCLLATLAS_40 CfdzDgHe_IMOTbrIP;
python 01-doublet_detection_scrublet.py BCLLATLAS_40 LxBTpkPO_8TcNpBg4;
python 01-doublet_detection_scrublet.py BCLLATLAS_40 pseMjZsU_qgZNOhOQ;
python 01-doublet_detection_scrublet.py BCLLATLAS_40 ujxNn2kq_lG2VLlYd;
python 01-doublet_detection_scrublet.py BCLLATLAS_40 uqJAc4r9_BYScOzxA;
python 01-doublet_detection_scrublet.py BCLLATLAS_46 BZEECBEG_GXkc6Q1y;
python 01-doublet_detection_scrublet.py BCLLATLAS_46 HjqdPU0E_aGDmEY5F;
python 01-doublet_detection_scrublet.py BCLLATLAS_46 KETfaLdx_Ub1mtE13;
python 01-doublet_detection_scrublet.py BCLLATLAS_46 SOJZt9kY_qpnv20QN;
python 01-doublet_detection_scrublet.py BCLLATLAS_46 XV1SLOR2_HRF5D9A3;
