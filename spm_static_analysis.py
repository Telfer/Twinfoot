import spm1d
import numpy as np
import matplotlib.pyplot as plt


## MZ data
export_array = np.zeros([101, 323])

for i in range(1, 323):
    # file paths
    fname1 = 'Z:/Projects/049 TwinFoot/DATA/Emed/sensor_dfs/mz1sensor_' + str(i) + '.csv'
    fname2 = 'Z:/Projects/049 TwinFoot/DATA/Emed/sensor_dfs/mz2sensor_' + str(i) + '.csv'
        
    # load data
    data1 = np.genfromtxt(fname1, delimiter = ',')
    data2 = np.genfromtxt(fname2, delimiter = ',')

    # run test
    t  = spm1d.stats.ttest_paired(data1, data2)
    ti = t.inference(alpha= (0.05 / 323), two_tailed=False)
    
    # plot data
    #ti.plot()
    #plt.show()
    
    # store z values in array
    export_array[:, (i - 1)] = ti.z

# save data
np.savetxt('Z:/Projects/049 TwinFoot/DATA/Emed/mz_static_results.csv', export_array, delimiter = ',')


## DY data
export_array = np.zeros([101, 323])

for i in range(1, 323):
    # file paths
    fname1 = 'Z:/Projects/049 TwinFoot/DATA/Emed/sensor_dfs/dy1sensor_' + str(i) + '.csv'
    fname2 = 'Z:/Projects/049 TwinFoot/DATA/Emed/sensor_dfs/dy2sensor_' + str(i) + '.csv'
    
    # load data
    data1 = np.genfromtxt(fname1, delimiter = ',')
    data2 = np.genfromtxt(fname2, delimiter = ',')

    # run test
    t  = spm1d.stats.ttest_paired(data1, data2)
    ti = t.inference(alpha = (0.05 / 323), two_tailed = False)
    
    # plot data
    #ti.plot()
    #plt.show()
    
    # store z values in array
    export_array[:, (i - 1)] = ti.z

# save data
np.savetxt('Z:/Projects/049 TwinFoot/DATA/Emed/dy_static_results.csv', export_array, delimiter = ',')