# Config file

cellranger_path = "/scratch/groups/hheyn/software/cellranger/4.0.0/cellranger"
reference_path = "/scratch/groups/hheyn/data/reference/refdata-gex-GRCh38-2020-A"
project_path = "/scratch/devel/rmassoni/tonsil_atlas/current/scRNA-seq/1-cellranger_mapping/"
REPEAT_MASK_GTF = "/scratch/groups/hheyn/data/reference/GRCh38_repeat_mask/GRCh38_repeat_mask.gtf"
REFERENCE_GTF = "/scratch/devel/rmassoni/reference/human/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf"
SUBPROJECTS = ["BCLLATLAS_05", "BCLLATLAS_06", "BCLLATLAS_07", "BCLLATLAS_16", "BCLLATLAS_19", "BCLLATLAS_21", "BCLLATLAS_22", "BCLLATLAS_25", "BCLLATLAS_34", "BCLLATLAS_41"]
ALL_TARGET_FILES = ['projects/BCLLATLAS_05/jobs/izi9unx1_8qdzhivu/izi9unx1_8qdzhivu.cmd',
                    'projects/BCLLATLAS_05/jobs/jb6vuao4_g4vi9ur0/jb6vuao4_g4vi9ur0.cmd',
                    'projects/BCLLATLAS_06/jobs/d8kwy76j_7blj9otf/d8kwy76j_7blj9otf.cmd',
                    'projects/BCLLATLAS_06/jobs/giz3qso4_783dbpu6/giz3qso4_783dbpu6.cmd',
                    'projects/BCLLATLAS_06/jobs/wkup7hvl_reo7jg84/wkup7hvl_reo7jg84.cmd',
                    'projects/BCLLATLAS_07/jobs/i5udk3x0_57gv6ncx/i5udk3x0_57gv6ncx.cmd',
                    'projects/BCLLATLAS_07/jobs/indj51br_urakqcmi/indj51br_urakqcmi.cmd',
                    'projects/BCLLATLAS_07/jobs/umt51kfr_p8ei65ms/umt51kfr_p8ei65ms.cmd',
                    'projects/BCLLATLAS_16/jobs/bw94nf57_vm85woki/bw94nf57_vm85woki.cmd',
                    'projects/BCLLATLAS_16/jobs/ggq3nifm_jkilwp1x/ggq3nifm_jkilwp1x.cmd',
                    'projects/BCLLATLAS_16/jobs/odyctre0_4qtf9h1z/odyctre0_4qtf9h1z.cmd',
                    'projects/BCLLATLAS_19/jobs/dvcbn9p8_ix0j3k8b/dvcbn9p8_ix0j3k8b.cmd',
                    'projects/BCLLATLAS_19/jobs/ff8s19u3_7e96iusr/ff8s19u3_7e96iusr.cmd',
                    'projects/BCLLATLAS_19/jobs/lf02okcv_j85nmx31/lf02okcv_j85nmx31.cmd',
                    'projects/BCLLATLAS_21/jobs/dvdzq8et_eot75su8/dvdzq8et_eot75su8.cmd',
                    'projects/BCLLATLAS_21/jobs/md651vbh_eymr91s7/md651vbh_eymr91s7.cmd',
                    'projects/BCLLATLAS_21/jobs/n1b3su0a_l7shyi35/n1b3su0a_l7shyi35.cmd',
                    'projects/BCLLATLAS_21/jobs/x739d5z1_dsamhgey/x739d5z1_dsamhgey.cmd',
                    'projects/BCLLATLAS_22/jobs/et5kifg3_0ufc6p4w/et5kifg3_0ufc6p4w.cmd',
                    'projects/BCLLATLAS_25/jobs/bz5rpwtv_kg7w108r/bz5rpwtv_kg7w108r.cmd',
                    'projects/BCLLATLAS_25/jobs/hg9au5oq_9ih1c0vd/hg9au5oq_9ih1c0vd.cmd',
                    'projects/BCLLATLAS_25/jobs/wf4su8ny_h4yj8bv7/wf4su8ny_h4yj8bv7.cmd',
                    'projects/BCLLATLAS_34/jobs/kjzv2rwx_sfomyxok/kjzv2rwx_sfomyxok.cmd',
                    'projects/BCLLATLAS_34/jobs/v8g80gtx_ps9bamz7/v8g80gtx_ps9bamz7.cmd',
                    'projects/BCLLATLAS_34/jobs/y7qn780g_p6jkgk63/y7qn780g_p6jkgk63.cmd',
                    'projects/BCLLATLAS_41/jobs/ejto2bae_y5mydeam/ejto2bae_y5mydeam.cmd',
                    'projects/BCLLATLAS_41/jobs/z3of7uaq_mzbhy4tt/z3of7uaq_mzbhy4tt.cmd'] 


ALL_VELOCYTO_FILES = ['projects/BCLLATLAS_05/jobs/izi9unx1_8qdzhivu/velocyto_izi9unx1_8qdzhivu.cmd',
                     'projects/BCLLATLAS_05/jobs/jb6vuao4_g4vi9ur0/velocyto_jb6vuao4_g4vi9ur0.cmd',
                     'projects/BCLLATLAS_06/jobs/d8kwy76j_7blj9otf/velocyto_d8kwy76j_7blj9otf.cmd',
                     'projects/BCLLATLAS_06/jobs/giz3qso4_783dbpu6/velocyto_giz3qso4_783dbpu6.cmd',
                     'projects/BCLLATLAS_06/jobs/wkup7hvl_reo7jg84/velocyto_wkup7hvl_reo7jg84.cmd',
                     'projects/BCLLATLAS_07/jobs/i5udk3x0_57gv6ncx/velocyto_i5udk3x0_57gv6ncx.cmd',
                     'projects/BCLLATLAS_07/jobs/indj51br_urakqcmi/velocyto_indj51br_urakqcmi.cmd',
                     'projects/BCLLATLAS_07/jobs/umt51kfr_p8ei65ms/velocyto_umt51kfr_p8ei65ms.cmd',
                     'projects/BCLLATLAS_16/jobs/bw94nf57_vm85woki/velocyto_bw94nf57_vm85woki.cmd',
                     'projects/BCLLATLAS_16/jobs/ggq3nifm_jkilwp1x/velocyto_ggq3nifm_jkilwp1x.cmd',
                     'projects/BCLLATLAS_16/jobs/odyctre0_4qtf9h1z/velocyto_odyctre0_4qtf9h1z.cmd',
                     'projects/BCLLATLAS_19/jobs/dvcbn9p8_ix0j3k8b/velocyto_dvcbn9p8_ix0j3k8b.cmd',
                     'projects/BCLLATLAS_19/jobs/ff8s19u3_7e96iusr/velocyto_ff8s19u3_7e96iusr.cmd',
                     'projects/BCLLATLAS_19/jobs/lf02okcv_j85nmx31/velocyto_lf02okcv_j85nmx31.cmd',
                     'projects/BCLLATLAS_21/jobs/dvdzq8et_eot75su8/velocyto_dvdzq8et_eot75su8.cmd',
                     'projects/BCLLATLAS_21/jobs/md651vbh_eymr91s7/velocyto_md651vbh_eymr91s7.cmd',
                     'projects/BCLLATLAS_21/jobs/n1b3su0a_l7shyi35/velocyto_n1b3su0a_l7shyi35.cmd',
                     'projects/BCLLATLAS_21/jobs/x739d5z1_dsamhgey/velocyto_x739d5z1_dsamhgey.cmd',
                     'projects/BCLLATLAS_22/jobs/et5kifg3_0ufc6p4w/velocyto_et5kifg3_0ufc6p4w.cmd',
                     'projects/BCLLATLAS_25/jobs/bz5rpwtv_kg7w108r/velocyto_bz5rpwtv_kg7w108r.cmd',
                     'projects/BCLLATLAS_25/jobs/hg9au5oq_9ih1c0vd/velocyto_hg9au5oq_9ih1c0vd.cmd',
                     'projects/BCLLATLAS_25/jobs/wf4su8ny_h4yj8bv7/velocyto_wf4su8ny_h4yj8bv7.cmd',
                     'projects/BCLLATLAS_34/jobs/kjzv2rwx_sfomyxok/velocyto_kjzv2rwx_sfomyxok.cmd',
                     'projects/BCLLATLAS_34/jobs/v8g80gtx_ps9bamz7/velocyto_v8g80gtx_ps9bamz7.cmd',
                     'projects/BCLLATLAS_34/jobs/y7qn780g_p6jkgk63/velocyto_y7qn780g_p6jkgk63.cmd',
                     'projects/BCLLATLAS_41/jobs/ejto2bae_y5mydeam/velocyto_ejto2bae_y5mydeam.cmd',
                     'projects/BCLLATLAS_41/jobs/z3of7uaq_mzbhy4tt/velocyto_z3of7uaq_mzbhy4tt.cmd'] 
