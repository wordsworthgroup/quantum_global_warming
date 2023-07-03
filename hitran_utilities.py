# Utilities for loading and manipulating HITRAN data

import pandas as pd
import numpy as np
import fortranformat as ff
import re

# see e.g., Rothmann 2004, Table 2
hitran_field_names = ['molecule number','isotopologue number','vacuum wavenumber','intensity',
                      'Einstein A-coefficient','air-broadened half-width','self-broadened half-width',
                      'lower-state energy','temperature-dependence exponent for gamma_air','air pressure-induced line shift',
                      'upper-state global quanta','lower-state global quanta','upper-state local quanta',
                      'lower-state local quanta','uncertainty indices','reference indices','line mixing flag',
                      'statistical weight of the upper state','statistical weight of the lower state']

def parse_hitran_line(line):
    
    # check line length
    if (len(line)!=161):
        print('error: unexpected line length',len(line))
        return('')
    
    # remove trailing \n
    line = line[:-1]
    
    # separate out the raw fields by character length
    # see e.g., Rothmann 2004, Table 1
    field_lengths = [2,1,12,10,10,5,5,10,4,8,15,15,15,15,6,12,1,7,7]
    
    raw_fields = []
    
    field_start_indices = np.insert(np.cumsum(field_lengths),0,0)
    
    for i in range(len(field_start_indices)-1):

        raw_fields.append(line[field_start_indices[i]:field_start_indices[i+1]])
        
        
    # format the raw fields according to fortran format identifiers
    fields = []
    hitran_line_formats = ['I2', 'I1', 'F12.6', 'E10.3', 'E10.3', 'F5.4', 'F5.4', 
                           'F10.4', 'F4.2', 'F8.6', 'A15', 'A15', 'A15', 'A15', 
                           '6I1', '6I2', 'A1', 'F7.1', 'F7.1']
    
    for i in range(len(raw_fields)):
        field_format = '('+hitran_line_formats[i]+')'
        reader = ff.FortranRecordReader(field_format)
        field = reader.read(raw_fields[i])
        
        if (len(field)==1):
            field = field[0]
        
        fields.append(field)
    
    return(fields)

def get_branch(local_quanta_field):
    split = re.findall('(\d+|[A-Za-z]+)', local_quanta_field)
    return(split[0])

def get_J(local_quanta_field):
    split = re.findall('(\d+|[A-Za-z]+)', local_quanta_field)
    J = int(split[1])
    return(J)

def get_iso_data(gas_file,isotope_number,nlines_max,status_freq=1000):
    
    with open(gas_file, 'r') as f:
        nlines = sum(1 for line in f)
        
    f = open(gas_file, 'r')
    
    # fields we want are wavenumber, line intensity, upper global quanta, lower global quanta, and branch
    df = pd.DataFrame({'wavenumber'         : pd.Series(dtype='float'),
                       'intensity'          : pd.Series(dtype='float'),
                       'Einstein_A'         : pd.Series(dtype='float'),
                       'Gamma_air'          : pd.Series(dtype='float'),
                       'Gamma_self'         : pd.Series(dtype='float'),
                       'E_gnd'              : pd.Series(dtype='float'),
                       'n_broad'            : pd.Series(dtype='float'),
                       'pres_shift'         : pd.Series(dtype='float'),
                       'upper_global'       : pd.Series(dtype='str'),
                       'lower_global'       : pd.Series(dtype='str'),
                       'branch'             : pd.Series(dtype='str'),
                       'J'                  : pd.Series(dtype='int64'),
                       'band_ID'            : pd.Series(dtype='int64')})
    count = 0
    
    # rewind to beginning of file      
    f.seek(0,0)
    
    band_dict = {}
    nbands    = 0
    
    for line in f:

        if (count < nlines_max):
            
            if (count%status_freq==0):
                print('line ',count)
            
            fields = parse_hitran_line(line)

            if (fields==''):
                print('reached empty line, assume end of file')
                break

            # only get selected isotope
            if (fields[1]==isotope_number or isotope_number==-1):

                line = {}
                line['wavenumber']   = fields[2]
                line['intensity']    = fields[3]
                line['Einstein_A']   = fields[4]
                line['Gamma_air']    = fields[5]
                line['Gamma_self']   = fields[6]
                line['E_gnd']        = fields[7]
                line['n_broad']      = fields[8]
                line['pres_shift']   = fields[9]
                line['upper_global'] = "".join(fields[10].split())                
                line['lower_global'] = "".join(fields[11].split())
                line['branch']       = get_branch(fields[13])
                line['J']            = get_J(fields[13])

                # assign line to band
                transition = line['lower_global']+(line['upper_global'])
                if (transition in band_dict.values()):
                    band_ID = list({i for i in band_dict if band_dict[i]==transition}).pop()
                else:
                    nbands = nbands + 1
                    band_ID = nbands
                    band_dict[band_ID] = transition
                
                line['band_ID'] = band_ID
                
                df = df.append(line,ignore_index=True)

                count = count + 1
                
            else:
                break
    print('')
    print('read ',count,'lines')
    print('found ',nbands,'bands')

    return(df)

# given a list of transition wavenumbers (nu_list), this function
# writes a hitran file (outfile) containing just those transitions
# selected from the infile
def write_hitran_subset(infile,outfile,nu_list,status_freq=1000):

    # check that selected wavenumbers are strictly increasing (code assumes this is the case)
    if (len(nu_list)>1):
        if (not np.all(nu_list[1:]-nu_list[:-1]>=0.)):
            print('Error: nu_select is not monotonically increasing')

    in_f  = open(infile,  'r')
    out_f = open(outfile, 'w')

    nlines_out = len(nu_list)
    print('Looking for',nlines_out,'lines')
    out_count = 0

    for line in in_f:

        fields = parse_hitran_line(line)

        nu = fields[2]

        if (nu>nu_list[-1]):
            break

        if (nu==nu_list[out_count]):
            #print('found ',nu)

            out_f.write(line)

            out_count = out_count + 1

            if (out_count%status_freq==0):
                print('Found',out_count,'lines')

        if (out_count == nlines_out):
            break

    if (out_count != nlines_out):
        print('Error: found unexpected number of lines')
    else:
        print('Found',nlines_out,'lines')

    out_f.close()
    in_f.close()
    
