# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 13:06:59 2026

@author: arun2
"""


def excel_read_funct(column_name, file_name, file_address = 'False'):
    import pandas as pd
    if file_address == 'False':
        pass
    else:
        file_name = file_address + file_name
    excel_data = pd.read_excel(file_name)
    data_array = excel_data[column_name]
    data_array_output = []
    for data_value in data_array:
        data_array_output.append(data_value)
    return data_array_output

def fill_time_funct(ambient_temperature, tank_volume, fill_rate):
    import CoolProp as CP
    HEOS = CP.AbstractState("HEOS&BICUBIC",'NitrousOxide')
    HEOS.update(CP.QT_INPUTS, 0, ambient_temperature)
    density = HEOS.rhomass()
    time = tank_volume/(density*fill_rate)
    return time

#output = excel_read_funct('Age','Data.xlsx')
fill_time = fill_time_funct(298, 4.6/1000, 5/1000)
print(fill_time, 'seconds')