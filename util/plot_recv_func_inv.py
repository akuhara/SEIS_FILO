import sys
from sub import filoplot


if __name__ == "__main__":
    # Get argument
    args = sys.argv
    if len(args) != 3:
        print("USAGE: plot.py [parameter file] [trace ID]")
        sys.exit()

    param_file = args[1]
    trace_id = int(args[2])
    rslt = filoplot.InvResult(param_file)
    rslt.draw_recv_func_set(trace_id)
    
