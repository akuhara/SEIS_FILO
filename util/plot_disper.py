import argparse
from sub import filoplot 
        
#-----------------------------------------------------------------------

if __name__ == "__main__":

    # Get argument
    parser = argparse.ArgumentParser(
        description='Make plot for surface wave inversion'
        )
    parser.add_argument('parameter_file',  
                        help='main parameter file of inversion'
                        )
    parser.add_argument('data_ID',
                        type=int,
                        help='ID of data file to be plotted')
    parser.add_argument('--vs_true',
                        default=None,
                        help='True Vs model (for synthetic test)'
                        )
    args = parser.parse_args()

    rslt = filoplot.InvResult(args.parameter_file,
                              args.vs_true)

    rslt.draw_surface_wave_set(args.data_ID)
