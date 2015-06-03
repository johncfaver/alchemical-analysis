################################################################################
# BOSS parser for alchemical analysis 
#
# John C. Faver 2015
#
#
#This library is free software; you can redistribute it and/or
#modify it under the terms of the GNU Lesser General Public
#License as published by the Free Software Foundation; either
#version 2.1 of the License, or (at your option) any later version.
#
#This library is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#Lesser General Public License for more details.
#
#You should have received a copy of the GNU Lesser General Public
#License along with this library; if not, see <http://www.gnu.org/licenses/>.
################################################################################

from glob import iglob
import numpy as np
import os
import sys

def readDataBoss(P):
    """Read log files from boss run.
    Input:
        Set of log files which contain solute energies at each reference lambda
        e.g. for sampling at lambda=0, first column is reference, columns
        2 and 3 are the 1st and 2nd perturbed solute energies in kcal/mol
        L = 0.00          L = 0.00          L = 0.05     Etotal @ L = 0.00
       -3.26159476       -3.26159476       -3.56442970    -4972.75616341
       -4.45238655       -4.45238655       -4.36447665    -4986.76173010
       -3.92963897       -3.92963897       -4.47076849    -4983.67419519 
    
    Return values:
        lambda_values:  nparray of nparray for each lambda
                        (value of each reference lambda)
        nsnapshots:     nparray of integers for each lambda
                        (number of samples for each reference lambda)
        dhdlt[k,n,t]:   None
                        (<dH/dl>i Used for TI calculations.)
        u_klt[k,m,t]:   3D nparray, beta*energy for configuration t,
                        of state k, evaluated at lambda=m
    """
    data_dict = {}          #Accumulate all data in a dict

    beta = P.beta * 4.184   #Convert beta to inverse kcal/mol 

    logfile_names = '{}.*.{}'.format(P.prefix,P.suffix)
    logfile_path = os.path.join(P.datafile_directory,logfile_names)
    file_count = 0
    #Iterate over log files
    for logfile in iglob(logfile_path):
        file_count += 1 
        print 'Reading file: {} '.format(logfile), 
        with open(logfile) as f:
            first_line = f.readline().split()
            try:
                k0 = first_line[2] #Reference solute
                k1 = first_line[5] #1st perturbed solute
                k2 = first_line[8] #2nd perturbed solute
            except Exception:
                print 'Could not read lambda values!'
                continue

            print 'reference lambda= {} '.format(k0),
            data_dict[k0] = {   'k1'         : k1, 
                                'k2'         : k2,
                                 k0          : [], 
                                 k1          : [], 
                                 k2          : [],    
                                'nsnapshots' : 0 
                            }
            #Iterate over snapshots
            for line in f.readlines():
                energy_line = line.split()
                #Energy of reference
                u_0 = beta * float(energy_line[0])   
                data_dict[k0][k0].append(u_0)
                #Energy of 1st perturbed Solute
                if k1 != k0:
                    u_1 = beta * float(energy_line[1])
                    data_dict[k0][k1].append(u_1)
                #Energy of 2nd perturbed Solute
                if k2 != k0:
                    u_2 = beta * float(energy_line[2])
                    data_dict[k0][k2].append(u_2)
                #Increment snapshot count
                data_dict[k0]['nsnapshots'] += 1
            print 'nsnapshots={} '.format(data_dict[k0]['nsnapshots'])
    print 'Read {} files.'.format(file_count)
    
    #Done reading files. Create variables for alchemical_analysis.py
    
    #number of states
    n_lambdas = len(data_dict)
   
    #list of states - needs to be sorted by lambda
    lambda_values = sorted([ l for l in data_dict ])
    print 'sorted lambda_values: {}'.format(lambda_values)
    
    #nparray of number of samples for each reference state
    nsnapshots = np.array([ data_dict[l]['nsnapshots'] for l in lambda_values ], np.int)
    
    #nparray of each reduced potential energy (b*E) t from simulation 
    #of lambda = k calculated with lambda = l
    u_klt = np.zeros((n_lambdas, n_lambdas, max(nsnapshots)))
    for k0 in data_dict:
        k0_idx = lambda_values.index(k0)
        n = nsnapshots[k0_idx]

        u_klt[k0_idx][k0_idx] = data_dict[k0][k0]
        
        l = data_dict[k0]['k1']
        if l in lambda_values:
            l_idx = lambda_values.index(l)
            u_klt[k0_idx][l_idx][:n] = data_dict[k0][l]

        l = data_dict[k0]['k2']
        if l in lambda_values:
            l_idx = lambda_values.index(l)
            u_klt[k0_idx][l_idx][:n] = data_dict[k0][l]
    
    #convert to nparray of floats
    lambda_values = np.array([ [float(l)] for l in lambda_values])

    #Needed by alchemical analysis for some reason
    P.lv_names=''

    #Log everyting
    output_file = os.path.join(P.output_directory,'bossparser.log')
    with open(output_file,'w') as f:
        f.write('beta = {}\n'.format(beta))
        f.write('lambda_values: {} \n'.format(lambda_values))
        f.write('nsnapshots: {} \n'.format(nsnapshots))
        for ix,i in enumerate(u_klt):
            for jx,j in enumerate(i):
                f.write('u_klt[{}][{}][:10] = '.format(ix,jx))
                for k in u_klt[ix][jx][:10]:
                    f.write('{} '.format(k))
                f.write(' \n')

    return nsnapshots,lambda_values,None,u_klt
        

