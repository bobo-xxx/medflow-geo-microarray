data_dir <- "/work/run/projects/bio-30/projects/8-GEO/data/"

# 查看一级文件和目录
list.files(data_dir)

# 递归查看所有文件
list.files(data_dir, recursive = TRUE)
data_dir <- "/work/run/projects/bio-30/projects/8-GEO/data/"
> 
  > # 查看一级文件和目录
  > list.files(data_dir)
[1] "GSE100155" "GSE102952" "GSE10300"  "GSE106099" "GSE106986" "GSE108000" "GSE109430" "GSE109857" "GSE111761" "GSE113212" "GSE113725" "GSE115797" "GSE120103" "GSE121618" "GSE12195" 
[16] "GSE125989" "GSE127952" "GSE129183" "GSE130588" "GSE13213"  "GSE133057" "GSE13355"  "GSE136701" "GSE139994" "GSE141864" "GSE14333"  "GSE146615" "GSE150082" "GSE153007" "GSE162102"
[31] "GSE16449"  "GSE16461"  "GSE16561"  "GSE166467" "GSE17260"  "GSE173608" "GSE17755"  "GSE180394" "GSE18606"  "GSE18850"  "GSE188944" "GSE19274"  "GSE19422"  "GSE19617"  "GSE20194" 
[46] "GSE20864"  "GSE2191"   "GSE21933"  "GSE22317"  "GSE23558"  "GSE24080"  "GSE24287"  "GSE26566"  "GSE26886"  "GSE26966"  "GSE27034"  "GSE27984"  "GSE28894"  "GSE29111"  "GSE29221" 
[61] "GSE30029"  "GSE30122"  "GSE30186"  "GSE30528"  "GSE30529"  "GSE30759"  "GSE31312"  "GSE31348"  "GSE31370"  "GSE32894"  "GSE34198"  "GSE34248"  "GSE34822"  "GSE35452"  "GSE36002" 
[76] "GSE36238"  "GSE37816"  "GSE37837"  "GSE38396"  "GSE38860"  "GSE39281"  "GSE39340"  "GSE39958"  "GSE40355"  "GSE40360"  "GSE41657"  "GSE41662"  "GSE41664"  "GSE42834"  "GSE42861" 
[91] "GSE43256"  "GSE43378"  "GSE43696"  "GSE43974"  "GSE44132"  "GSE44295"  "GSE44314"  "GSE44711"  "GSE45001"  "GSE45603"  "GSE46394"  "GSE47472"  "GSE47915"  "GSE49515"  "GSE51588" 
[106] "GSE52068"  "GSE52093"  "GSE52793"  "GSE53552"  "GSE53849"  "GSE54236"  "GSE54388"  "GSE54536"  "GSE54618"  "GSE54837"  "GSE5500"   "GSE55747"  "GSE56420"  "GSE56885"  "GSE57691" 
[121] "GSE58121"  "GSE58294"  "GSE58435"  "GSE59444"  "GSE60436"  "GSE61616"  "GSE62336"  "GSE63409"  "GSE63695"  "GSE63881"  "GSE64634"  "GSE66271"  "GSE67530"  "GSE68004"  "GSE68020" 
[136] "GSE6891"   "GSE69223"  "GSE70453"  "GSE71647"  "GSE73461"  "GSE73463"  "GSE73894"  "GSE75819"  "GSE76427"  "GSE76826"  "GSE76895"  "GSE7904"   "GSE81211"  "GSE83456"  "GSE84796" 
[151] "GSE85195"  "GSE85446"  "GSE87053"  "GSE87211"  "GSE90074"  "GSE92324"  "GSE95233"  "GSE97466"  "GSE98224"  "GSE98278"  "GSE98770"  "GSE98895"  "platforms"
> 
  > # 递归查看所有文件
  > list.files(data_dir, recursive = TRUE)
[1] "GSE100155/suppl/filelist.txt"                                                                       
[2] "GSE100155/suppl/GSE100155_Liver_Transplant_non-normalized.txt.gz"                                   
[3] "GSE100155/suppl/GSE100155_Liver_Transplant_normalized.txt.gz"                                       
[4] "GSE100155/suppl/GSE100155_RAW.tar"                                                                  
[5] "GSE100155/suppl/GSE100155_series_matrix.txt.gz"                                                     
[6] "GSE102952/suppl/filelist.txt"                                                                       
[7] "GSE102952/suppl/GSE102952_RAW.tar"                                                                  
[8] "GSE102952/suppl/GSE102952_unmethylated_and_methylated_data.txt.gz"                                  
[9] "GSE10300/suppl/filelist.txt"                                                                        
[10] "GSE10300/suppl/GSE10300_RAW.tar"                                                                    
[11] "GSE106099/suppl/filelist.txt"                                                                       
[12] "GSE106099/suppl/GSE106099_methylation_dAEC_vs_AEC.xlsx"                                             
[13] "GSE106099/suppl/GSE106099_methylation_dVEC_vs_VEC.xlsx"                                             
[14] "GSE106099/suppl/GSE106099_RAW.tar"                                                                  
[15] "GSE106986/GSE106986_series_matrix.txt.gz"                                                           
[16] "GSE106986/suppl/filelist.txt"                                                                       
[17] "GSE106986/suppl/GSE106986_log2_norm_annotated_matrix.txt.gz"                                        
[18] "GSE106986/suppl/GSE106986_RAW.tar"                                                                  
[19] "GSE108000/GSE108000_series_matrix.txt.gz"                                                           
[20] "GSE108000/suppl/filelist.txt"                                                                       
[21] "GSE108000/suppl/GSE108000_RAW.tar"                                                                  
[22] "GSE109430/suppl/filelist.txt"                                                                       
[23] "GSE109430/suppl/GSE109430_RAW.tar"                                                                  
[24] "GSE109857/suppl/filelist.txt"                                                                       
[25] "GSE109857/suppl/GSE109857_RAW.tar"                                                                  
[26] "GSE111761/suppl/filelist.txt"                                                                       
[27] "GSE111761/suppl/GSE111761_RAW.tar"                                                                  
[28] "GSE113212/suppl/filelist.txt"                                                                       
[29] "GSE113212/suppl/GSE113212_RAW.tar"                                                                  
[30] "GSE113725/suppl/filelist.txt"                                                                       
[31] "GSE113725/suppl/GSE113725_detectionP.csv.gz"                                                        
[32] "GSE113725/suppl/GSE113725_methylatedIntensities.csv.gz"                                             
[33] "GSE113725/suppl/GSE113725_RAW.tar"                                                                  
[34] "GSE113725/suppl/GSE113725_rawBetas.csv.gz"                                                          
[35] "GSE113725/suppl/GSE113725_unmethylatedIntensities.csv.gz"                                           
[36] "GSE115797/suppl/filelist.txt"                                                                       
[37] "GSE115797/suppl/GSE115797_RAW.tar"                                                                  
[38] "GSE120103/GSE120103_series_matrix.txt.gz"                                                           
[39] "GSE120103/suppl/filelist.txt"                                                                       
[40] "GSE120103/suppl/GSE120103_RAW.tar"                                                                  
[41] "GSE121618/suppl/filelist.txt"                                                                       
[42] "GSE121618/suppl/GSE121618_RAW.tar"                                                                  
[43] "GSE12195/suppl/filelist.txt"                                                                        
[44] "GSE12195/suppl/GSE12195_RAW.tar"                                                                    
[45] "GSE125989/suppl/filelist.txt"                                                                       
[46] "GSE125989/suppl/GSE125989_RAW.tar"                                                                  
[47] "GSE127952/suppl/filelist.txt"                                                                       
[48] "GSE127952/suppl/GSE127952_RAW.tar"                                                                  
[49] "GSE129183/suppl/filelist.txt"                                                                       
[50] "GSE129183/suppl/GSE129183_RAW.tar"                                                                  
[51] "GSE129183/suppl/GSE129183_signals.txt.gz"                                                           
[52] "GSE130588/suppl/filelist.txt"                                                                       
[53] "GSE130588/suppl/GSE130588_RAW.tar"                                                                  
[54] "GSE13213/suppl/filelist.txt"                                                                        
[55] "GSE13213/suppl/GSE13213_AD117_patient_info.txt"                                                     
[56] "GSE13213/suppl/GSE13213_RAW.tar"                                                                    
[57] "GSE133057/suppl/filelist.txt"                                                                       
[58] "GSE133057/suppl/GSE133057_RAW.tar"                                                                  
[59] "GSE13355/suppl/filelist.txt"                                                                        
[60] "GSE13355/suppl/GSE13355_RAW.tar"                                                                    
[61] "GSE136701/suppl/filelist.txt"                                                                       
[62] "GSE136701/suppl/GSE136701_RAW.tar"                                                                  
[63] "GSE139994/suppl/filelist.txt"                                                                       
[64] "GSE139994/suppl/GSE139994_Normalized_Data.txt.gz"                                                   
[65] "GSE139994/suppl/GSE139994_RAW.tar"                                                                  
[66] "GSE141864/suppl/filelist.txt"                                                                       
[67] "GSE141864/suppl/GSE141864_Processed_datafile_Heart_Lungs_Kidneys_Liver_Spleen_FFPE_filt_5_log2.xlsx"
[68] "GSE141864/suppl/GSE141864_RAW.tar"                                                                  
[69] "GSE14333/suppl/filelist.txt"                                                                        
[70] "GSE14333/suppl/GSE14333_RAW.tar"                                                                    
[71] "GSE146615/suppl/filelist.txt"                                                                       
[72] "GSE146615/suppl/GSE146615_non-normalized_data.txt.gz"                                               
[73] "GSE146615/suppl/GSE146615_RAW.tar"                                                                  
[74] "GSE150082/suppl/filelist.txt"                                                                       
[75] "GSE150082/suppl/GSE150082_RAW.tar"                                                                  
[76] "GSE153007/suppl/filelist.txt"                                                                       
[77] "GSE153007/suppl/GSE153007_RAW.tar"                                                                  
[78] "GSE162102/suppl/filelist.txt"                                                                       
[79] "GSE162102/suppl/GSE162102_raw_and_normalized_signal.txt.gz"                                         
[80] "GSE162102/suppl/GSE162102_RAW.tar"                                                                  
[81] "GSE16449/suppl/filelist.txt"                                                                        
[82] "GSE16449/suppl/GSE16449_RAW.tar"                                                                    
[83] "GSE16461/suppl/filelist.txt"                                                                        
[84] "GSE16461/suppl/GSE16461_RAW.tar"                                                                    
[85] "GSE16561/suppl/filelist.txt"                                                                        
[86] "GSE16561/suppl/GSE16561_RAW.tar"                                                                    
[87] "GSE16561/suppl/GSE16561_RAW.txt.gz"                                                                 
[88] "GSE166467/suppl/filelist.txt"                                                                       
[89] "GSE166467/suppl/GSE166467_non-normalized.txt.gz"                                                    
[90] "GSE166467/suppl/GSE166467_RAW.tar"                                                                  
[91] "GSE17260/suppl/filelist.txt"                                                                        
[92] "GSE17260/suppl/GSE17260_additional_clinical_information.xls.gz"                                     
[93] "GSE17260/suppl/GSE17260_clinical_information.txt.gz"                                                
[94] "GSE17260/suppl/GSE17260_RAW.tar"                                                                    
[95] "GSE173608/suppl/filelist.txt"                                                                       
[96] "GSE173608/suppl/GSE173608_RAW.tar"                                                                  
[97] "GSE17755/suppl/filelist.txt"                                                                        
[98] "GSE17755/suppl/GSE17755_RAW.tar"                                                                    
[99] "GSE180394/suppl/filelist.txt"                                                                       
[100] "GSE180394/suppl/GSE180394_RAW.tar"                                                                  
[101] "GSE18606/suppl/filelist.txt"                                                                        
[102] "GSE18606/suppl/GSE18606_RAW.tar"                                                                    
[103] "GSE18850/suppl/filelist.txt"                                                                        
[104] "GSE18850/suppl/GSE18850_RAW.tar"                                                                    
[105] "GSE188944/suppl/filelist.txt"                                                                       
[106] "GSE188944/suppl/GSE188944_RAW.tar"                                                                  
[107] "GSE19274/suppl/filelist.txt"                                                                        
[108] "GSE19274/suppl/GSE19274_clinical_annotation_codes.txt.gz"                                           
[109] "GSE19274/suppl/GSE19274_non-normalized.txt.gz"                                                      
[110] "GSE19274/suppl/GSE19274_RAW.tar"                                                                    
[111] "GSE19422/suppl/filelist.txt"                                                                        
[112] "GSE19422/suppl/GSE19422_RAW.tar"                                                                    
[113] "GSE19617/suppl/filelist.txt"                                                                        
[114] "GSE19617/suppl/GSE19617_RAW.tar"                                                                    
[115] "GSE20194/suppl/filelist.txt"                                                                        
[116] "GSE20194/suppl/GSE20194_MDACC_Sample_Info.xls.gz"                                                   
[117] "GSE20194/suppl/GSE20194_RAW.tar"                                                                    
[118] "GSE20864/suppl/filelist.txt"                                                                        
[119] "GSE20864/suppl/GSE20864_RAW.tar"                                                                    
[120] "GSE2191/GSE2191_series_matrix.txt.gz"                                                               
[121] "GSE22317/suppl/filelist.txt"                                                                        
[122] "GSE22317/suppl/GSE22317_RAW.tar"                                                                    
[123] "GSE23558/suppl/filelist.txt"                                                                        
[124] "GSE23558/suppl/GSE23558_RAW.tar"                                                                    
[125] "GSE24080/suppl/filelist.txt"                                                                        
[126] "GSE24080/suppl/GSE24080_MM_UAMS565_ClinInfo_27Jun2008_LS_clean.xls.gz"                              
[127] "GSE24080/suppl/GSE24080_RAW.tar"                                                                    
[128] "GSE24287/suppl/filelist.txt"                                                                        
[129] "GSE24287/suppl/GSE24287_RAW.tar"                                                                    
[130] "GSE26566/suppl/filelist.txt"                                                                        
[131] "GSE26566/suppl/GPL6104_Illumina_HumanRef-8_V2_0_R1_11223162_A.bgx"                                  
[132] "GSE26566/suppl/GPL6104_Illumina_HumanRef-8_V2_0_R1_11223162_A.bgx.gz"                               
[133] "GSE26566/suppl/GSE26566_non-normalized.txt.gz"                                                      
[134] "GSE26566/suppl/GSE26566_RAW.tar"                                                                    
[135] "GSE26886/suppl/filelist.txt"                                                                        
[136] "GSE26886/suppl/GSE26886_RAW.tar"                                                                    
[137] "GSE26966/suppl/filelist.txt"                                                                        
[138] "GSE26966/suppl/GSE26966_RAW.tar"                                                                    
[139] "GSE27034/suppl/filelist.txt"                                                                        
[140] "GSE27034/suppl/GSE27034_RAW.tar"                                                                    
[141] "GSE27984/suppl/filelist.txt"                                                                        
[142] "GSE27984/suppl/GSE27984_RAW.tar"                                                                    
[143] "GSE28894/suppl/filelist.txt"                                                                        
[144] "GSE28894/suppl/GSE28894_non-normalized.txt.gz"                                                      
[145] "GSE28894/suppl/GSE28894_RAW.tar"                                                                    
[146] "GSE29111/suppl/filelist.txt"                                                                        
[147] "GSE29111/suppl/GSE29111_file1_DE Day30.txt.gz"                                                      
[148] "GSE29111/suppl/GSE29111_file1_esetFull_miR_NO_MI_combv5_exprs.txt.gz"                               
[149] "GSE29111/suppl/GSE29111_file2_DEvisit4.txt.gz"                                                      
[150] "GSE29111/suppl/GSE29111_file2_esetFull_miR_NO_MI_combv4_exprs.txt.gz"                               
[151] "GSE29111/suppl/GSE29111_RAW.tar"                                                                    
[152] "GSE29221/suppl/filelist.txt"                                                                        
[153] "GSE29221/suppl/GSE29221_non-normalized.txt.gz"                                                      
[154] "GSE29221/suppl/GSE29221_RAW.tar"                                                                    
[155] "GSE30029/suppl/filelist.txt"                                                                        
[156] "GSE30029/suppl/GSE30029_non-normalized.txt.gz"                                                      
[157] "GSE30029/suppl/GSE30029_normalized.txt.gz"                                                          
[158] "GSE30029/suppl/GSE30029_RAW.tar"                                                                    
[159] "GSE30122/suppl/filelist.txt"                                                                        
[160] "GSE30122/suppl/GSE30122_RAW.tar"                                                                    
[161] "GSE30186/suppl/filelist.txt"                                                                        
[162] "GSE30186/suppl/GSE30186_non_normalized.txt.gz"                                                      
[163] "GSE30186/suppl/GSE30186_RAW.tar"                                                                    
[164] "GSE30528/suppl/filelist.txt"                                                                        
[165] "GSE30528/suppl/GSE30528_RAW.tar"                                                                    
[166] "GSE30529/suppl/filelist.txt"                                                                        
[167] "GSE30529/suppl/GSE30529_RAW.tar"                                                                    
[168] "GSE30759/suppl/filelist.txt"                                                                        
[169] "GSE30759/suppl/GSE30759_RAW.tar"                                                                    
[170] "GSE31312/suppl/filelist.txt"                                                                        
[171] "GSE31312/suppl/GSE31312_Microarray_and_clinical_data_DLBCL_475_cases_PMID_22437443.pdf.gz"          
[172] "GSE31312/suppl/GSE31312_RAW.tar"                                                                    
[173] "GSE31348/suppl/filelist.txt"                                                                        
[174] "GSE31348/suppl/GSE31348_RAW.tar"                                                                    
[175] "GSE31370/suppl/filelist.txt"                                                                        
[176] "GSE31370/suppl/GSE31370_non_normalized.txt.gz"                                                      
[177] "GSE31370/suppl/GSE31370_RAW.tar"                                                                    
[178] "GSE32894/suppl/filelist.txt"                                                                        
[179] "GSE32894/suppl/GSE32894_non-normalized_308UCsamples.txt.gz"                                         
[180] "GSE32894/suppl/GSE32894_RAW.tar"                                                                    
[181] "GSE32894/suppl/GSE32894_reps_normals_preprocess_matrix-non-normalized.txt.gz"                       
[182] "GSE32894/suppl/GSE32894_reps_normals_preprocess_metadata.txt.gz"                                    
[183] "GSE34198/suppl/filelist.txt"                                                                        
[184] "GSE34198/suppl/GSE34198_controlsTable.txt.gz"                                                       
[185] "GSE34198/suppl/GSE34198_non-normalized.txt.gz"                                                      
[186] "GSE34198/suppl/GSE34198_RAW.tar"                                                                    
[187] "GSE34248/suppl/filelist.txt"                                                                        
[188] "GSE34248/suppl/GSE34248_RAW.tar"                                                                    
[189] "GSE34822/suppl/filelist.txt"                                                                        
[190] "GSE34822/suppl/GSE34822_RAW.tar"                                                                    
[191] "GSE35452/suppl/filelist.txt"                                                                        
[192] "GSE35452/suppl/GSE35452_RAW.tar"                                                                    
[193] "GSE36002/suppl/filelist.txt"                                                                        
[194] "GSE36002/suppl/GSE36002_non_normalized.txt.gz"                                                      
[195] "GSE36002/suppl/GSE36002_RAW.tar"                                                                    
[196] "GSE36238/suppl/filelist.txt"                                                                        
[197] "GSE36238/suppl/GSE36238_RAW.tar"                                                                    
[198] "GSE37816/suppl/filelist.txt"                                                                        
[199] "GSE37816/suppl/GSE37816_RAW.tar"                                                                    
[200] "GSE37816/suppl/GSE37816_signal_intensities.txt.gz"                                                  
[201] "GSE37837/suppl/filelist.txt"                                                                        
[202] "GSE37837/suppl/GSE37837_RAW.tar"                                                                    
[203] "GSE38396/suppl/filelist.txt"                                                                        
[204] "GSE38396/suppl/GSE38396_RAW.tar"                                                                    
[205] "GSE38860/suppl/filelist.txt"                                                                        
[206] "GSE38860/suppl/GSE38860_methylated_unmethylated_for_tumor_28.txt.gz"                                
[207] "GSE38860/suppl/GSE38860_methylated_unmethylated_signals.txt.gz"                                     
[208] "GSE38860/suppl/GSE38860_RAW.tar"                                                                    
[209] "GSE39281/suppl/filelist.txt"                                                                        
[210] "GSE39281/suppl/GSE39281_all_data_by_genes.txt.gz"                                                   
[211] "GSE39281/suppl/GSE39281_all_lesions.conf_95.txt.gz"                                                 
[212] "GSE39281/suppl/GSE39281_all_thresholded.by_genes.txt.gz"                                            
[213] "GSE39281/suppl/GSE39281_amp_genes.conf_95.txt.gz"                                                   
[214] "GSE39281/suppl/GSE39281_amp_qplot.pdf.gz"                                                           
[215] "GSE39281/suppl/GSE39281_arraylistfile.txt.gz"                                                       
[216] "GSE39281/suppl/GSE39281_broad_data_by_genes.txt.gz"                                                 
[217] "GSE39281/suppl/GSE39281_broad_significance_results.txt.gz"                                          
[218] "GSE39281/suppl/GSE39281_broad_values_by_arm.txt.gz"                                                 
[219] "GSE39281/suppl/GSE39281_del_genes.conf_95.txt.gz"                                                   
[220] "GSE39281/suppl/GSE39281_del_qplot.pdf.gz"                                                           
[221] "GSE39281/suppl/GSE39281_focal_data_by_genes.txt.gz"                                                 
[222] "GSE39281/suppl/GSE39281_freqarms_vs_ngenes.pdf.gz"                                                  
[223] "GSE39281/suppl/GSE39281_GISTIC_README.txt"                                                          
[224] "GSE39281/suppl/GSE39281_raw_copy_number.pdf.gz"                                                     
[225] "GSE39281/suppl/GSE39281_RAW.tar"                                                                    
[226] "GSE39281/suppl/GSE39281_regions_track.conf_95.bed.gz"                                               
[227] "GSE39281/suppl/GSE39281_sample_cutoffs.txt.gz"                                                      
[228] "GSE39281/suppl/GSE39281_spanishaCGH_clinical.txt.gz"                                                
[229] "GSE39281/suppl/GSE39281_spanishaCGH.glad.txt.gz"                                                    
[230] "GSE39340/suppl/filelist.txt"                                                                        
[231] "GSE39340/suppl/GSE39340_non_normalized.txt.gz"                                                      
[232] "GSE39340/suppl/GSE39340_RAW.tar"                                                                    
[233] "GSE39958/suppl/filelist.txt"                                                                        
[234] "GSE39958/suppl/GSE39958_RAW.tar"                                                                    
[235] "GSE39958/suppl/GSE39958_signal_intensity.txt.gz"                                                    
[236] "GSE40355/suppl/filelist.txt"                                                                        
[237] "GSE40355/suppl/GSE40355_RAW.tar"                                                                    
[238] "GSE40360/suppl/filelist.txt"                                                                        
[239] "GSE40360/suppl/GSE40360_GenomeStudio_norm_bg.txt.gz"                                                
[240] "GSE40360/suppl/GSE40360_RAW_nonorm_nobg.txt.gz"                                                     
[241] "GSE40360/suppl/GSE40360_RAW.tar"                                                                    
[242] "GSE41657/suppl/filelist.txt"                                                                        
[243] "GSE41657/suppl/GSE41657_RAW.tar"                                                                    
[244] "GSE41662/suppl/filelist.txt"                                                                        
[245] "GSE41662/suppl/GSE41662_RAW.tar"                                                                    
[246] "GSE41664/suppl/filelist.txt"                                                                        
[247] "GSE41664/suppl/GSE41664_RAW.tar"                                                                    
[248] "GSE42834/suppl/filelist.txt"                                                                        
[249] "GSE42834/suppl/GSE42834_RAW.tar"                                                                    
[250] "GSE42861/suppl/filelist.txt"                                                                        
[251] "GSE42861/suppl/GSE42861_methylation_signal_matrix_SUBSETS.tar.gz"                                   
[252] "GSE42861/suppl/GSE42861_methylation_signal_matrix.txt.gz"                                           
[253] "GSE42861/suppl/GSE42861_non-methylated_signal_matrix_SUBSETS.tar.gz"                                
[254] "GSE42861/suppl/GSE42861_non-methylated_signal_matrix.txt.gz"                                        
[255] "GSE42861/suppl/GSE42861_processed_methylation_matrix_SUBSETS.tar.gz"                                
[256] "GSE42861/suppl/GSE42861_processed_methylation_matrix.txt.gz"                                        
[257] "GSE42861/suppl/GSE42861_RAW.tar"                                                                    
[258] "GSE42861/suppl/GSE42861_Readme.txt"                                                                 
[259] "GSE43256/suppl/filelist.txt"                                                                        
[260] "GSE43256/suppl/GSE43256_methylated_unmethylated_signals.txt.gz"                                     
[261] "GSE43256/suppl/GSE43256_RAW.tar"                                                                    
[262] "GSE43378/suppl/filelist.txt"                                                                        
[263] "GSE43378/suppl/GSE43378_RAW.tar"                                                                    
[264] "GSE43696/suppl/filelist.txt"                                                                        
[265] "GSE43696/suppl/GSE43696_RAW.tar"                                                                    
[266] "GSE43974/suppl/filelist.txt"                                                                        
[267] "GSE43974/suppl/GSE43974_non-normalized_data.txt.gz"                                                 
[268] "GSE43974/suppl/GSE43974_paired_analyses_DBD_DCD_Living.xls.gz"                                      
[269] "GSE43974/suppl/GSE43974_RAW.tar"                                                                    
[270] "GSE44132/suppl/filelist.txt"                                                                        
[271] "GSE44132/suppl/GSE44132_FinalReport_zkaminsky_PPD_GEO.txt.gz"                                       
[272] "GSE44132/suppl/GSE44132_PPDarray_detection_Pvals2.txt.gz"                                           
[273] "GSE44132/suppl/GSE44132_RAW.tar"                                                                    
[274] "GSE44295/suppl/filelist.txt"                                                                        
[275] "GSE44295/suppl/GPL6883_HumanRef-8_V3_0_R0_11282963_A.bgx"                                           
[276] "GSE44295/suppl/GPL6883_HumanRef-8_V3_0_R0_11282963_A.bgx.gz"                                        
[277] "GSE44295/suppl/GSE44295_Expression_Matrix_non-normalized.csv"                                       
[278] "GSE44295/suppl/GSE44295_Expression_Matrix_non-normalized.csv.gz"                                    
[279] "GSE44295/suppl/GSE44295_Expression_Matrix_normalized.csv.gz"                                        
[280] "GSE44295/suppl/GSE44295_RAW.tar"                                                                    
[281] "GSE44314/suppl/filelist.txt"                                                                        
[282] "GSE44314/suppl/GSE44314_RAW.tar"                                                                    
[283] "GSE44711/suppl/filelist.txt"                                                                        
[284] "GSE44711/suppl/GSE44711_non-normalized_data.txt.gz"                                                 
[285] "GSE44711/suppl/GSE44711_RAW.tar"                                                                    
[286] "GSE45001/suppl/filelist.txt"                                                                        
[287] "GSE45001/suppl/GSE45001_RAW.tar"                                                                    
[288] "GSE45603/suppl/filelist.txt"                                                                        
[289] "GSE45603/suppl/GSE45603_Discovery_rawdata.txt.gz"                                                   
[290] "GSE45603/suppl/GSE45603_RAW.tar"                                                                    
[291] "GSE46394/suppl/filelist.txt"                                                                        
[292] "GSE46394/suppl/GSE46394_methylated_unmethylated_data.txt.gz"                                        
[293] "GSE46394/suppl/GSE46394_RAW.tar"                                                                    
[294] "GSE47472/suppl/filelist.txt"                                                                        
[295] "GSE47472/suppl/GSE47472_RAW.tar"                                                                    
[296] "GSE47472/suppl/GSE47472_Rawdata_GEO_AAA_Neck.txt.gz"                                                
[297] "GSE47915/suppl/filelist.txt"                                                                        
[298] "GSE47915/suppl/GSE47915_RAW.tar"                                                                    
[299] "GSE47915/suppl/GSE47915_signal_intensities.txt.gz"                                                  
[300] "GSE49515/suppl/filelist.txt"                                                                        
[301] "GSE49515/suppl/GSE49515_RAW.tar"                                                                    
[302] "GSE51588/suppl/filelist.txt"                                                                        
[303] "GSE51588/suppl/GSE51588_RAW.tar"                                                                    
[304] "GSE52068/suppl/filelist.txt"                                                                        
[305] "GSE52068/suppl/GSE52068_RAW.tar"                                                                    
[306] "GSE52068/suppl/GSE52068_signal_intensities.txt.gz"                                                  
[307] "GSE52093/suppl/filelist.txt"                                                                        
[308] "GSE52093/suppl/GSE52093_non-normalized.txt.gz"                                                      
[309] "GSE52093/suppl/GSE52093_RAW.tar"                                                                    
[310] "GSE52793/suppl/filelist.txt"                                                                        
[311] "GSE52793/suppl/GSE52793_non-normalized.txt.gz"                                                      
[312] "GSE52793/suppl/GSE52793_RAW.tar"                                                                    
[313] "GSE53552/suppl/filelist.txt"                                                                        
[314] "GSE53552/suppl/GSE53552_RAW.tar"                                                                    
[315] "GSE53849/suppl/filelist.txt"                                                                        
[316] "GSE53849/suppl/GSE53849_methylated_unmethylated_signal_intensities.txt.gz"                          
[317] "GSE53849/suppl/GSE53849_RAW.tar"                                                                    
[318] "GSE54236/suppl/filelist.txt"                                                                        
[319] "GSE54236/suppl/GSE54236_RAW.tar"                                                                    
[320] "GSE54388/suppl/filelist.txt"                                                                        
[321] "GSE54388/suppl/GSE54388_RAW.tar"                                                                    
[322] "GSE54536/suppl/filelist.txt"                                                                        
[323] "GSE54536/suppl/GSE54536_Non-normalized_data.txt.gz"                                                 
[324] "GSE54536/suppl/GSE54536_RAW.tar"                                                                    
[325] "GSE54618/suppl/filelist.txt"                                                                        
[326] "GSE54618/suppl/GSE54618_non-normalized.txt.gz"                                                      
[327] "GSE54618/suppl/GSE54618_RAW.tar"                                                                    
[328] "GSE54837/suppl/filelist.txt"                                                                        
[329] "GSE54837/suppl/GSE54837_RAW.tar"                                                                    
[330] "GSE5500/suppl/filelist.txt"                                                                         
[331] "GSE5500/suppl/GSE5500_RAW.tar"                                                                      
[332] "GSE55747/suppl/filelist.txt"                                                                        
[333] "GSE55747/suppl/GSE55747_Non-normalized_data.txt.gz"                                                 
[334] "GSE55747/suppl/GSE55747_RAW.tar"                                                                    
[335] "GSE56420/suppl/filelist.txt"                                                                        
[336] "GSE56420/suppl/GSE56420_RAW.tar"                                                                    
[337] "GSE56420/suppl/GSE56420_unmethylated_methylated.txt.gz"                                             
[338] "GSE56885/suppl/filelist.txt"                                                                        
[339] "GSE56885/suppl/GSE56885_RAW.tar"                                                                    
[340] "GSE57691/suppl/filelist.txt"                                                                        
[341] "GSE57691/suppl/GSE57691_non-normalized_data.txt.gz"                                                 
[342] "GSE57691/suppl/GSE57691_RAW.tar"                                                                    
[343] "GSE58121/suppl/filelist.txt"                                                                        
[344] "GSE58121/suppl/GSE58121_RAW.tar"                                                                    
[345] "GSE58294/suppl/filelist.txt"                                                                        
[346] "GSE58294/suppl/GSE58294_RAW.tar"                                                                    
[347] "GSE58435/suppl/filelist.txt"                                                                        
[348] "GSE58435/suppl/GSE58435_RAW.tar"                                                                    
[349] "GSE59444/suppl/filelist.txt"                                                                        
[350] "GSE59444/suppl/GSE59444_RAW.tar"                                                                    
[351] "GSE60436/suppl/filelist.txt"                                                                        
[352] "GSE60436/suppl/GSE60436_non-normalized_data.txt.gz"                                                 
[353] "GSE60436/suppl/GSE60436_RAW.tar"                                                                    
[354] "GSE61616/suppl/filelist.txt"                                                                        
[355] "GSE61616/suppl/GSE61616_RAW.tar"                                                                    
[356] "GSE62336/suppl/filelist.txt"                                                                        
[357] "GSE62336/suppl/GSE62336_RAW.tar"                                                                    
[358] "GSE62336/suppl/GSE62336_signals.txt.gz"                                                             
[359] "GSE63409/suppl/filelist.txt"                                                                        
[360] "GSE63409/suppl/GSE63409_RAW.tar"                                                                    
[361] "GSE63695/suppl/filelist.txt"                                                                        
[362] "GSE63695/suppl/GSE63695_RAW.tar"                                                                    
[363] "GSE63695/suppl/GSE63695_unmethylated_methylated.txt.gz"                                             
[364] "GSE63881/suppl/filelist.txt"                                                                        
[365] "GSE63881/suppl/GSE63881_non-normalized.txt.gz"                                                      
[366] "GSE63881/suppl/GSE63881_RAW.tar"                                                                    
[367] "GSE64634/suppl/filelist.txt"                                                                        
[368] "GSE64634/suppl/GSE64634_RAW.tar"                                                                    
[369] "GSE66271/suppl/filelist.txt"                                                                        
[370] "GSE66271/suppl/GSE66271_RAW.tar"                                                                    
[371] "GSE67530/suppl/filelist.txt"                                                                        
[372] "GSE67530/suppl/GSE67530_RAW.tar"                                                                    
[373] "GSE67530/suppl/GSE67530_signals.txt.gz"                                                             
[374] "GSE68004/suppl/filelist.txt"                                                                        
[375] "GSE68004/suppl/GSE68004_non-normalized_ProbeRowSignalDetectionDataset_5383_143523_1.txt.gz"         
[376] "GSE68004/suppl/GSE68004_non-normalized_ProbeRowSignalDetectionDataset_5383_143529_1.txt.gz"         
[377] "GSE68004/suppl/GSE68004_RAW.tar"                                                                    
[378] "GSE68020/suppl/filelist.txt"                                                                        
[379] "GSE68020/suppl/GSE68020_non-normalized.txt.gz"                                                      
[380] "GSE68020/suppl/GSE68020_RAW.tar"                                                                    
[381] "GSE6891/suppl/filelist.txt"                                                                         
[382] "GSE6891/suppl/GSE6891_RAW.tar"                                                                      
[383] "GSE69223/suppl/filelist.txt"                                                                        
[384] "GSE69223/suppl/GSE69223_RAW.tar"                                                                    
[385] "GSE70453/suppl/filelist.txt"                                                                        
[386] "GSE70453/suppl/GSE70453_Matrix_signal_intensities.csv.gz"                                           
[387] "GSE70453/suppl/GSE70453_RAW.tar"                                                                    
[388] "GSE71647/suppl/filelist.txt"                                                                        
[389] "GSE71647/suppl/GSE71647_RAW.tar"                                                                    
[390] "GSE73461/suppl/filelist.txt"                                                                        
[391] "GSE73461/suppl/GSE73461_GEOupload_Discovery_Dataset_Normalised_Sept_15_n_459.txt.gz"                
[392] "GSE73461/suppl/GSE73461_GEOupload_Discovery_Dataset_Raw_Sept_15_n_459.txt.gz"                       
[393] "GSE73461/suppl/GSE73461_RAW.tar"                                                                    
[394] "GSE73463/suppl/filelist.txt"                                                                        
[395] "GSE73463/suppl/GSE73463_GEOupload_Validation_HT12V4_Dataset_Normalised_Sept_15_n_233.txt.gz"        
[396] "GSE73463/suppl/GSE73463_GEOupload_Validation_HT12V4_Dataset_Raw_Sept_15_n_233.txt.gz"               
[397] "GSE73463/suppl/GSE73463_RAW.tar"                                                                    
[398] "GSE73894/suppl/filelist.txt"                                                                        
[399] "GSE73894/suppl/GSE73894_Matrix_signal_intensities.txt.gz"                                           
[400] "GSE73894/suppl/GSE73894_RAW.tar"                                                                    
[401] "GSE75819/suppl/filelist.txt"                                                                        
[402] "GSE75819/suppl/GSE75819_non-normalized.txt.gz"                                                      
[403] "GSE75819/suppl/GSE75819_RAW.tar"                                                                    
[404] "GSE76427/suppl/filelist.txt"                                                                        
[405] "GSE76427/suppl/GSE76427_non-normalized.txt.gz"                                                      
[406] "GSE76427/suppl/GSE76427_RAW.tar"                                                                    
[407] "GSE76826/suppl/filelist.txt"                                                                        
[408] "GSE76826/suppl/GSE76826_RAW.tar"                                                                    
[409] "GSE76895/suppl/filelist.txt"                                                                        
[410] "GSE76895/suppl/GSE76895_DifferentialExpression-T2DvsND.txt.gz"                                      
[411] "GSE76895/suppl/GSE76895_RAW.tar"                                                                    
[412] "GSE7904/suppl/filelist.txt"                                                                         
[413] "GSE7904/suppl/GSE7904_RAW.tar"                                                                      
[414] "GSE81211/suppl/filelist.txt"                                                                        
[415] "GSE81211/suppl/GSE81211_non_normalized.txt.gz"                                                      
[416] "GSE81211/suppl/GSE81211_RAW.tar"                                                                    
[417] "GSE83456/suppl/filelist.txt"                                                                        
[418] "GSE83456/suppl/GSE83456_matrix_EPTB_unnormalized.txt.gz"                                            
[419] "GSE83456/suppl/GSE83456_RAW.tar"                                                                    
[420] "GSE84796/suppl/filelist.txt"                                                                        
[421] "GSE84796/suppl/GSE84796_RAW.tar"                                                                    
[422] "GSE85195/suppl/filelist.txt"                                                                        
[423] "GSE85195/suppl/GSE85195_RAW.tar"                                                                    
[424] "GSE85446/suppl/filelist.txt"                                                                        
[425] "GSE85446/suppl/GSE85446_processed_data.txt.gz"                                                      
[426] "GSE85446/suppl/GSE85446_RAW.tar"                                                                    
[427] "GSE87053/suppl/filelist.txt"                                                                        
[428] "GSE87053/suppl/GSE87053_RAW.tar"                                                                    
[429] "GSE87211/suppl/filelist.txt"                                                                        
[430] "GSE87211/suppl/GSE87211_RAW.tar"                                                                    
[431] "GSE90074/suppl/filelist.txt"                                                                        
[432] "GSE90074/suppl/GSE90074_RAW.tar"                                                                    
[433] "GSE92324/suppl/filelist.txt"                                                                        
[434] "GSE92324/suppl/GSE92324_Non-normalized_data.txt.gz"                                                 
[435] "GSE92324/suppl/GSE92324_RAW.tar"                                                                    
[436] "GSE92324/suppl/GSE92324_sample_probe_profile.txt.gz"                                                
[437] "GSE95233/suppl/filelist.txt"                                                                        
[438] "GSE95233/suppl/GSE95233_RAW.tar"                                                                    
[439] "GSE97466/suppl/filelist.txt"                                                                        
[440] "GSE97466/suppl/GSE97466_RAW.tar"                                                                    
[441] "GSE98224/suppl/filelist.txt"                                                                        
[442] "GSE98224/suppl/GSE98224_average_beta.csv.gz"                                                        
[443] "GSE98224/suppl/GSE98224_expr.txt.gz"                                                                
[444] "GSE98224/suppl/GSE98224_RAW.tar"                                                                    
[445] "GSE98224/suppl/GSE98224_signal_intensities.csv.gz"                                                  
[446] "GSE98278/suppl/filelist.txt"                                                                        
[447] "GSE98278/suppl/GSE98278_Non-normalized_data.txt.gz"                                                 
[448] "GSE98278/suppl/GSE98278_RAW.tar"                                                                    
[449] "GSE98770/suppl/filelist.txt"                                                                        
[450] "GSE98770/suppl/GSE98770_RAW.tar"                                                                    
[451] "GSE98895/suppl/filelist.txt"                                                                        
[452] "GSE98895/suppl/GSE98895_non-normalized.txt.gz"                                                      
[453] "GSE98895/suppl/GSE98895_RAW.tar"                                                                    
[454] "platforms/GPL10150/GPL10150.soft.gz"                                                                
[455] "platforms/GPL10558/GPL10558.soft.gz"                                                                
[456] "platforms/GPL1261/GPL1261.soft.gz"                                                                  
[457] "platforms/GPL1291/GPL1291.soft.gz"                                                                  
[458] "platforms/GPL13497/GPL13497.soft.gz"                                                                
[459] "platforms/GPL13534/GPL13534.soft.gz"                                                                
[460] "platforms/GPL1355/GPL1355.soft.gz"                                                                  
[461] "platforms/GPL14550/GPL14550.soft.gz"                                                                
[462] "platforms/GPL14951/GPL14951.soft.gz"                                                                
[463] "platforms/GPL15207/GPL15207.soft.gz"                                                                
[464] "platforms/GPL17077/GPL17077.soft.gz"                                                                
[465] "platforms/GPL17586/GPL17586.soft.gz"                                                                
[466] "platforms/GPL19983/GPL19983.soft.gz"                                                                
[467] "platforms/GPL20995/GPL20995.soft.gz"                                                                
[468] "platforms/GPL21185/GPL21185.soft.gz"                                                                
[469] "platforms/GPL570/GPL570.soft.gz"                                                                    
[470] "platforms/GPL571/GPL571.soft.gz"                                                                    
[471] "platforms/GPL6102/GPL6102.soft.gz"                                                                  
[472] "platforms/GPL6104/GPL6104.soft.gz"                                                                  
[473] "platforms/GPL6480/GPL6480.soft.gz"                                                                  
[474] "platforms/GPL6848/GPL6848.soft.gz"                                                                  
[475] "platforms/GPL6883/GPL6883.soft.gz"                                                                  
[476] "platforms/GPL6884/GPL6884.soft.gz"                                                                  
[477] "platforms/GPL6885/GPL6885.soft.gz"                                                                  
[478] "platforms/GPL6947/GPL6947.soft.gz"                                                                  
[479] "platforms/GPL7202/GPL7202.soft.gz"                                                                  
[480] "platforms/GPL8300/GPL8300.soft.gz"                                                                  
[481] "platforms/GPL8490/GPL8490.soft.gz"                                                                  
[482] "platforms/GPL96/GPL96.soft.gz" 


#Step 0：准备路径
data_dir <- "/work/run/projects/bio-30/projects/8-GEO/data"

#Step 1：找到所有 filelist.txt
filelist_paths <- list.files(
  data_dir,
  pattern = "filelist.txt",
  recursive = TRUE,
  full.names = TRUE
)
#Step 2：读取 + 整理每个 GSE 的文件信息
library(dplyr)
library(stringr)

parse_one_filelist <- function(filelist_path) {
  
  gse_id <- basename(dirname(dirname(filelist_path)))  # GSEXXXX
  
  df <- read.table(
    filelist_path,
    header = FALSE,
    sep = "",
    stringsAsFactors = FALSE,
    fill = TRUE,
    quote = ""
  )
  
  colnames(df)[1:6] <- c("entry_type", "filename", "date", "time", "size", "suffix")
  
  df$GSE <- gse_id
  
  df
}

all_files_df <- bind_rows(
  lapply(filelist_paths, parse_one_filelist)
)


#tep 3：按 GSE 级别判断数据类型（核心）
gse_summary <- all_files_df %>%
  group_by(GSE) %>%
  summarise(
    has_RAW_tar = any(grepl("_RAW\\.tar", filename, ignore.case = TRUE)),
    has_CEL     = any(grepl("\\.CEL(\\.gz)?$", filename, ignore.case = TRUE)),
    has_GPR     = any(grepl("\\.GPR(\\.gz)?$", filename, ignore.case = TRUE)),
    has_IDAT    = any(grepl("\\.idat(\\.gz)?$", filename, ignore.case = TRUE)),
    has_BPM     = any(grepl("\\.bpm", filename, ignore.case = TRUE)),
    has_matrix  = any(grepl("series_matrix", filename, ignore.case = TRUE)),
    file_types  = paste(sort(unique(suffix)), collapse = ";"),
    .groups = "drop"
  )

#Step 4：推断“数据类型 + 是否可处理”
gse_summary <- gse_summary %>%
  mutate(
    inferred_type = case_when(
      has_CEL ~ "Affymetrix_microarray",
      has_GPR ~ "Agilent_microarray",
      has_IDAT & has_BPM ~ "Methylation_array",
      has_matrix ~ "Matrix_only",
      TRUE ~ "Unknown"
    ),
    processable = inferred_type %in% c(
      "Affymetrix_microarray",
      "Agilent_microarray"
    )
  )

#step 5：导出成“总览文档”（非常重要）
write.csv(
  gse_summary,
  file = "GSE_data_type_summary.csv",
  row.names = FALSE
)


# 创建一个列表存储所有 filelist 的内容
all_filelists <- lapply(filelist_paths, readLines)

# 如果想合并成一个向量
all_files <- unlist(all_filelists)

# 读取所有 filelist.txt 内容，并合并成一个数据框
all_files <- lapply(filelist_paths, function(f) {
  df <- read.table(f, header = FALSE, stringsAsFactors = FALSE)
  df$source_file <- f  # 保留来源文件信息
  return(df)
})

all_files_df <- do.call(rbind, all_files)
head(all_files_df)

library(dplyr)
library(stringr)

all_files_df <- all_files_df %>%
  rename(
    entry_type = V1,
    filename   = V2,
    date       = V3,
    time       = V4,
    size       = V5,
    suffix     = V6
  )
head(all_files_df)


all_files_df <- all_files_df %>%
  mutate(
    GSE = str_extract(source_file, "GSE\\d+")
  )

unique(all_files_df$GSE)[1:10]

gse_summary <- all_files_df %>%
  group_by(GSE) %>%
  summarise(
    has_RAW_tar = any(suffix == "TAR"),
    has_CEL     = any(suffix == "CEL"),
    has_GPR     = any(suffix == "GPR"),
    has_IDAT    = any(suffix == "IDAT"),
    has_BPM     = any(suffix == "BPM"),
    has_TXT     = any(suffix == "TXT"),
    has_CSV     = any(suffix == "CSV"),
    has_XLSX    = any(suffix == "XLSX"),
    has_matrix  = any(grepl("series_matrix", filename, ignore.case = TRUE)),
    suffix_set  = paste(sort(unique(suffix)), collapse = ";"),
    .groups = "drop"
  )

gse_summary <- gse_summary %>%
  mutate(
    inferred_type = case_when(
      has_CEL ~ "Affymetrix_microarray",
      has_GPR ~ "Agilent_microarray",
      has_IDAT & has_BPM ~ "Methylation_array",
      has_matrix ~ "Matrix_only",
      TRUE ~ "Unknown"
    ),
    processable = inferred_type %in% c(
      "Affymetrix_microarray",
      "Agilent_microarray"
    )
  )



write.csv(
  gse_summary,
  file = "GSE_microarray_data_inventory.csv",
  row.names = FALSE
)

# 
# 读取临川信息
# # 安装包
# # if(!require("GEOquery"))
# #   BiocManager::install("GEOquery",update = F,ask = F)
# 
# #加载包
# library(GEOquery)
# library(dplyr)
# 
# 
# geoID <- "GSE2191"
# GSE <- getGEO(GEO = geoID, destdir = "raw", getGPL = F)
# #先下载好再读入
# GSE <- getGEO(filename = paste0("data/GSE2191/",geoID,"_series_matrix.txt.gz"),getGPL = F)  
# 
# #提取表型数据或临床信息数据
# clinical <- pData(GSE)    #用pData()函数获取分组信息，这里包含了所有临床数据
# 
# 
