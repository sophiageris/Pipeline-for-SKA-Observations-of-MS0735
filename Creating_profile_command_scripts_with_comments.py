#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 12:02:43 2023

@author: sophiageris
"""

template_file = 'template.txt' # sample file that contains information to be copied to the new files
sz_files = 'sz_files.txt'

integration_time = '10'
for HA_replace in [0.5, 1, 2]:
    for days_replace in [1, 2]:

        # Function to get the lines from a file
        def get_lines(filename):
            lines = []
            with open(filename, 'r') as f:
                lines = f.readlines()
            return lines

        # Function to determine the type of sz_file (thermal or nonthermal)
        def sz_type(filename):
            if 'non' in filename:
                return 'nonthermal'
            else:
                return 'thermal'
            
        def types(filename):
            if 'upper' in filename:
                return 'upper'
            else:
                return 'lower'    
    
        # Function to generate duplicate filenames for the new files
        def dup_filename(sz_files, HA_replace, days_replace):
            sz_names = get_lines(sz_files)
            type_name = get_lines(sz_files)
            
            type_list = []
            SZ_type_list = []
            for name in sz_names:
                SZ_type = sz_type(name)
                SZ_type_list.append(SZ_type)
                Types = types(name)
                type_list.append(Types)
        
            # Create duplicate filenames based on HA_replace, days_replace, and SZ_type
            dup_filename_list = ['sim_' + 'ms0735_' + i + '_' + 'MUSTANG' + '_' + j + '_' + 'Jy_' + 'pix_' +  str(HA_replace) + '_' + integration_time + '_' + str(days_replace) + '.sh' for i, j in zip(SZ_type_list, type_list)]

            return dup_filename_list

        # Function to generate the path list for sz_files
        def sz_path_lists(sz_files):
            sz_names = get_lines(sz_files)
            path_list = ['./' + name for name in sz_names]
            return path_list

        # Function to modify the sz_path in the file
        def sz_path(path_name, file):
            with open(file, 'r') as f:
                contents = f.read()
            contents = contents.replace('YMAP.fits', path_name.strip())
            with open(file, 'w') as f: 
                f.write(contents)  
        
        # Function to replace variables in the file
        def replace_vars_(HA_replace, days_replace, file):
            name = file.replace('.sh', '')
            with open(file, 'r') as f:
                contents = f.read()
            contents = contents.replace('HA', str(HA_replace))
            contents = contents.replace('days', str(days_replace))
            contents = contents.replace('SIM_OBS', name, 1)
            contents = contents.replace('name', name, 1)
            contents = contents.replace('integration_time', integration_time)
            with open(file, 'w') as f: 
                f.write(contents)
        
        # Function to create a copy of the template file with the given filenames
        def _make_copy_(template_file, dup_filenames):
            with open(template_file, 'r') as rf: 
                with open(dup_filenames, 'w') as wf:
                    for line in rf:
                        wf.write(line)
            return dup_filenames

        # Generate the new filenames for the combination of HA_replace and days_replace
        new_files_names = dup_filename(sz_files, HA_replace, days_replace)
        # Generate the path list for sz_files
        sz_path_list = sz_path_lists(sz_files)

        # Loop over the new filenames and path list to create and modify the files
        for i, j in zip(new_files_names, sz_path_list):
            # Create a copy of the template file with the new filename
            file = _make_copy_(template_file, i)
            # Replace variables in the file with the corresponding values
            replace_vars_(HA_replace, days_replace, file)
            # Modify the sz_path in the file
            sz_path(j, file)

