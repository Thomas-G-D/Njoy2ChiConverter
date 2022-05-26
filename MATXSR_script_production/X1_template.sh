#
# Simple cross-section for O16 using mostly defaults and 
# FIXME: wrong title and missing files in repo for custom wt 
#
CWD=$PWD



#================================= Set properties here
export ENDF_ROOT=/home/tdaggie/ENDF-B-VIII.0
neutron_file="Neutron_filename"

# sab_file="tsl-graphite.endf"

output_directory="/home/tdaggie/MGXS/outputs/ENDF-B-VIII-0/172gxs/"
output_file_prefix="output_prefix"

#================================= Run NJOY
cd /home/tdaggie/MGXS/Njoy2ChiConverter/njoy_runner || exit

python3 X0_NjoyInput_gen.py \
--njoy_exec_name=njoy21 \
--path_to_neutron_endf=$ENDF_ROOT/neutrons/$neutron_file \
--temperature=293.6 \
--neutron_group_structure=22 \
--output_directory=$output_directory \
--output_filename=$output_file_prefix.njoy \
--matxsr_huse=userID 


# --path_to_neutron_endf=$ENDF_ROOT/neutrons/$neutron_file \
# --path_to_sab=$ENDF_ROOT/thermal_scatt/$sab_file \
# --inelastic_thermal_number=229 \
# --inelastic_thermal_num_atoms=1 \
# --elastic_thermal_number=230 \
# --path_to_gamma_endf= \
# --temperature=296.0 \
# --neutron_group_structure=22 \
# --neutron_weight_function=8 \
# --output_directory=$output_directory \
# --output_filename=$output_file_prefix.njoy \
# --gamma_group_structure=0 \
# --gamma_weight_function=2 \
# --custom_neutron_gs_file="" \
# --custom_gamma_gs_file="" \
# --custom_neutron_wt_file="" \
# --custom_gamma_wt_file="" \

cd "$CWD" || exit

#================================= Run converter
cd /home/tdaggie/MGXS/Njoy2ChiConverter/njoy_processor || exit

python3 njoy_processor.py \
--output_path=$output_directory \
--njoy_output_filename=$output_file_prefix.njoy \
--chixs_filename=$output_file_prefix.cxs \
## --plot

cd "$CWD" || exit

