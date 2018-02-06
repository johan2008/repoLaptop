import os
import sys
import shutil
import ConfigParser
import time

current_path = os.path.dirname(__file__)
# config file
config = ConfigParser.RawConfigParser()
config.read('%s/config.cfg'%current_path)
# bin paths
internal_bin_path = os.path.join(current_path, config.get('Bin-Directories', 'internal_bin_rel_path'))
external_bin_path = os.path.join(current_path, config.get('Bin-Directories', 'external_bin_rel_path'))
# settings
blended_settings_path = os.path.join(current_path, config.get('Bin-Directories', 'blended_settings_rel_path'))
blended_vis_settings_path = os.path.join(current_path, config.get('Bin-Directories', 'blended_vis_settings_rel_path'))
bench_settings_path = os.path.join(current_path, config.get('Bin-Directories', 'bench_settings_rel_path'))
bench_vis_settings_path = os.path.join(current_path, config.get('Bin-Directories', 'bench_vis_settings_rel_path'))
agd_settings_path = os.path.join(current_path, config.get('Bin-Directories', 'agd_settings_rel_path'))
agd_vis_settings_path = os.path.join(current_path, config.get('Bin-Directories', 'agd_vis_settings_rel_path'))
sym_settings_path = os.path.join(current_path, config.get('Bin-Directories', 'sym_settings_rel_path'))
sym_vis_settings_path = os.path.join(current_path, config.get('Bin-Directories', 'sym_vis_settings_rel_path'))

# data paths
detected_axis_path = os.path.join(current_path, config.get('Data-Directories', 'detected_axis_rel_path'))
marked_axis_path = os.path.join(current_path, config.get('Data-Directories', 'marked_axis_rel_path'))
marked_axisalign_path = os.path.join(current_path, config.get('Data-Directories', 'marked_axisalign_rel_path'))
marked_corrs_path = os.path.join(current_path, config.get('Data-Directories', 'marked_corrs_rel_path'))
