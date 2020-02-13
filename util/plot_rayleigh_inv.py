import sys
from sub import filoplot 



        
#-----------------------------------------------------------------------

if __name__ == "__main__":
    # Get argument
    args = sys.argv
    if len(args) != 2:
        print("USAGE: plot.py [parameter file]")
        sys.exit()
        
    param_file = args[1]
    rslt = filoplot.InvResult(param_file)
    rslt.draw_surface_wave_set()
