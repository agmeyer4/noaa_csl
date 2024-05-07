'''
Module: noaa_csl_funcs.py
Author: Aaron G. Meyer (agmeyer4@gmail.com)
Last Updated: May 2, 2024
Description:
Useful functions for dealing with noaa csl data

More information can be found here: https://csl.noaa.gov/groups/csl7/measurements/2020covid-aqs/emissions/
'''

##################################################################################################################################
#Import Packages
import os
import subprocess
import shutil
import glob

##################################################################################################################################
# Define Functions
def replace_all_strs(text,dic):
    '''Replaces all strings contained in text matching the keys of dic with the values of dic
    
    Args:
    text (str) : the string we want to replace stuff in
    dic (dict) : a dictionary with keys as what should be replaced in text and values with what the key should be replaced with
    
    Returns :
    text (str): the new text with keys replaced
    '''

    for key,value in dic.items():
        text = text.replace(key,value)
    return text

def month_int_to_str(month_int):
        '''Converts a month integer to a character string like "MonthXX"
        
        Args: 
        month_int (int) : the integer of the month (1,2,12,etc)
        
        Returns (str) : a string of form "MonthXX" where XX is the two digit month (1 --> Month01)'''

        return f'Month{month_int:02d}'

def check_space(path,excep_thresh='8Tb'):
    '''Checks the amount of space on a filesystem given a path, and raise an error if the amount of space is below the threshold
    
    Args:
    path (str) : path to check how much room there would be based on the filesystem
    excep_thresh (str) : threshold below which to raise an exception, defaults '8T'
    
    Returns:
    usage_details (tuple) : output of shutil.disk_usage(path)
    
    Raises:
    MemoryError : if the amount of space on the filesystem containing the path is below the excep_thresh
    '''

    disk_usage = shutil.disk_usage(path)
    total_bytes,used_bytes,free_bytes = int(disk_usage.total), int(disk_usage.used), int(disk_usage.free)
    thresh_bytes = human_to_bytes(excep_thresh)
    if free_bytes < thresh_bytes:
        print(f'Path: {path}')
        print(f'Total space: {bytes_to_human(total_bytes)}')
        print(f'Used space: {bytes_to_human(used_bytes)}')
        print(f'Free space: {bytes_to_human(free_bytes)}')
        print(f'Required space (a user defined property): {excep_thresh}\n')
        raise MemoryError(f'There is less than {excep_thresh} in the destination filesystem. Ensure there is enough room for your data download, or change the threshold')
    else:
        print(f'Sufficient space detected in {path}')
        print(f'Found {bytes_to_human(free_bytes)} free space')
        return disk_usage

def human_to_bytes(human_readable, units=['b', 'Kb', 'Mb', 'Gb', 'Tb', 'Pb']):
    '''Makes a human readable storage size into a bytes integer
    #TODO Different definitions for this and the other direction (1000 bytes per kb vs 1024)
    
    Args:
    human_readable (str) : a human readable string representing space (something like 2.1Gb, or 5Tb)
    units (list) : list of strings defining the units 
    
    Returns:
    bytes (int) : integer representing how many bytes the human readable string is
    '''

    # Split the input string into value and unit
    value_str = ''.join([ch for ch in human_readable if ch.isdigit() or ch == '.'])
    unit_str = ''.join([ch for ch in human_readable if not ch.isdigit() and ch != '.']).strip()
    value = float(value_str)# Convert value string to a float
    
    # Find the unit in the provided units list
    try:
        unit_index = units.index(unit_str)
    except ValueError:
        raise ValueError(f"Invalid unit '{unit_str}' in input string.")
    
    multiplier = 1024 ** unit_index # Calculate the multiplier based on the unit index
    bytes = int(value * multiplier)   # Convert the human-readable value to bytes    
    return bytes

def bytes_to_human(bytes,units=['b','Kb','Mb','Gb','Tb','Pb']):
    '''Converts an integer bytes into a human readable string notation
    #TODO Different definitions for this and the other direction (1000 bytes per kb vs 1024)
    
    Args: 
    bytes (int) : integer bytes
    units (list) : unit convention for factors of 10 
    
    Returns (str) : the human readable string
    '''

    return str(bytes) + units[0] if bytes < 1024 else bytes_to_human(bytes>>10, units[1:])

