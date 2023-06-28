set -eux

CTL_PREFIX=/home/aiokada/sandbox/nanomonsv-tutorial/control_panel/hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control/hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control

tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr10:58717444-58717494   > control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr10:58717622-58717686  >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr10:7017529-7017579    >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr10:7090881-7090934    >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr1:224458881-224458956 >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr12:72272793-72272813  >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr12:72273092-72273151  >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr15:23440440-23440485  >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr15:23462431-23462463  >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr15:23464836-23464857  >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr15:84141953-84142009  >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr20:15019944-15019998  >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr20:15019966-15019986  >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr3:26390406-26390481   >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr3:26390406-26390532   >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr3:25359548-25359595   >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr12:72273249-72273314  >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr19:17285957-17286022  >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr19:17286811-17286867  >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr1:224612373-224612438 >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr12:72273110-72273130  >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr3:25359085-25359133   >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr20:38645803-38645850  >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr15:23462432-23462464  >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr15:23464884-23464904  >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr7:151049504-151049591 >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr20:15033176-15033262  >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr20:15033185-15033205  >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr6:26193783-26193858   >> control_panel.bp_info.bed
tabix ${CTL_PREFIX}.bp_info.sorted.bed.gz chr6:26193783-26193860   >> control_panel.bp_info.bed

tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr10:58717444-58717494   > control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr10:58717622-58717686  >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr10:7017529-7017579    >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr10:7090881-7090934    >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr1:224458881-224458956 >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr12:72272793-72272813  >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr12:72273092-72273151  >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr15:23440440-23440485  >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr15:23462431-23462463  >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr15:23464836-23464857  >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr15:84141953-84142009  >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr20:15019944-15019998  >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr20:15019966-15019986  >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr3:26390406-26390481   >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr3:26390406-26390532   >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr3:25359548-25359595   >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr12:72273249-72273314  >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr19:17285957-17286022  >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr19:17286811-17286867  >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr1:224612373-224612438 >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr12:72273110-72273130  >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr3:25359085-25359133   >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr20:38645803-38645850  >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr15:23462432-23462464  >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr15:23464884-23464904  >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr7:151049504-151049591 >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr20:15033176-15033262  >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr20:15033185-15033205  >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr6:26193783-26193858   >> control_panel.deletion.bed
tabix ${CTL_PREFIX}.deletion.sorted.bed.gz chr6:26193783-26193860   >> control_panel.deletion.bed

tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr10:58717444-58717494   > control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr10:58717622-58717686  >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr10:7017529-7017579    >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr10:7090881-7090934    >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr1:224458881-224458956 >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr12:72272793-72272813  >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr12:72273092-72273151  >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr15:23440440-23440485  >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr15:23462431-23462463  >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr15:23464836-23464857  >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr15:84141953-84142009  >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr20:15019944-15019998  >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr20:15019966-15019986  >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr3:26390406-26390481   >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr3:26390406-26390532   >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr3:25359548-25359595   >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr12:72273249-72273314  >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr19:17285957-17286022  >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr19:17286811-17286867  >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr1:224612373-224612438 >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr12:72273110-72273130  >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr3:25359085-25359133   >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr20:38645803-38645850  >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr15:23462432-23462464  >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr15:23464884-23464904  >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr7:151049504-151049591 >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr20:15033176-15033262  >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr20:15033185-15033205  >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr6:26193783-26193858   >> control_panel.insertion.bed
tabix ${CTL_PREFIX}.insertion.sorted.bed.gz chr6:26193783-26193860   >> control_panel.insertion.bed

tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr10:58717444-58717494   > control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr10:58717622-58717686  >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr10:7017529-7017579    >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr10:7090881-7090934    >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr1:224458881-224458956 >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr12:72272793-72272813  >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr12:72273092-72273151  >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr15:23440440-23440485  >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr15:23462431-23462463  >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr15:23464836-23464857  >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr15:84141953-84142009  >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr20:15019944-15019998  >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr20:15019966-15019986  >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr3:26390406-26390481   >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr3:26390406-26390532   >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr3:25359548-25359595   >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr12:72273249-72273314  >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr19:17285957-17286022  >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr19:17286811-17286867  >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr1:224612373-224612438 >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr12:72273110-72273130  >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr3:25359085-25359133   >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr20:38645803-38645850  >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr15:23462432-23462464  >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr15:23464884-23464904  >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr7:151049504-151049591 >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr20:15033176-15033262  >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr20:15033185-15033205  >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr6:26193783-26193858   >> control_panel.rearrangement.bedpe
tabix ${CTL_PREFIX}.rearrangement.sorted.bedpe.gz chr6:26193783-26193860   >> control_panel.rearrangement.bedpe

sort -k1,1 -k2,2n -k3,3n control_panel.bp_info.bed > control_panel.bp_info.sorted.bed
sort -k1,1 -k2,2n -k3,3n control_panel.deletion.bed > control_panel.deletion.sorted.bed
sort -k1,1 -k2,2n -k3,3n control_panel.insertion.bed > control_panel.insertion.sorted.bed
sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n control_panel.rearrangement.bedpe > control_panel.rearrangement.sorted.bedpe

bgzip control_panel.bp_info.sorted.bed
bgzip control_panel.deletion.sorted.bed
bgzip control_panel.insertion.sorted.bed
bgzip control_panel.rearrangement.sorted.bedpe

tabix -f -p bed control_panel.bp_info.sorted.bed.gz
tabix -f -p bed control_panel.deletion.sorted.bed.gz
tabix -f -p bed control_panel.insertion.sorted.bed.gz
tabix -f -p bed control_panel.rearrangement.sorted.bedpe.gz

rm control_panel.bp_info.bed control_panel.deletion.bed control_panel.insertion.bed control_panel.rearrangement.bedpe


#chr10:58717444-58717494  
#chr10:58717622-58717686  
#chr10:7017529-7017579    
#chr10:7090881-7090934    
#chr1:224458881-224458956 
#chr12:72272793-72272813  
#chr12:72273092-72273151  
#chr15:23440440-23440485  
#chr15:23462431-23462463  
#chr15:23464836-23464857  
#chr15:84141953-84142009  
#chr20:15019944-15019998  
#chr20:15019966-15019986  
#chr3:26390406-26390481   
#chr3:26390406-26390532   
#chr3:25359548-25359595   
#chr12:72273249-72273314  
#chr19:17285957-17286022  
#chr19:17286811-17286867  
#chr1:224612373-224612438 
#chr12:72273110-72273130  
#chr3:25359085-25359133   
#chr20:38645803-38645850  
#chr15:23462432-23462464  
#chr15:23464884-23464904  
#chr7:151049504-151049591 
#chr20:15033176-15033262  
#chr20:15033185-15033205  
#chr6:26193783-26193858   
#chr6:26193783-26193860   

