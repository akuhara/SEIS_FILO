import sys
from sub import filoplot 
        
#-----------------------------------------------------------------------

if __name__ == "__main__":
    # Get argument
    args = sys.argv
    if len(args) != 3:
        print("USAGE: plot.py [parameter file] [dispersion curve ID]")
        sys.exit()
        
    param_file = args[1]
    curve_id = int(args[2])
    rslt = filoplot.InvResult(param_file)
    rslt.draw_surface_wave_set(curve_id)
