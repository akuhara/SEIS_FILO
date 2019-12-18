import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

eps = 1.0e-8

class InvRslt:
    def __init__(self, param_file):

        # Check if file exists or not
        if os.path.exists(param_file) == False:
            print("ERROR:" + param_file + "is not found", \
                  file=sys.stderr)
            sys.exit(1)
        self._param_file = param_file
        
        # Get parameters
        self._param = {}
        self._read_param_file()
        self._read_obs_header()
        print(self._param)

    #-------------------------------------------------------------------
    
    def _read_param_file(self):
        print("Reading parameter file " + self._param_file)
        with open(self._param_file, 'r') as f:
            for line in f:
                tmp_line = line.rstrip()
                # Dealing with comment out
                i = tmp_line.find('#')
                if i >= 0:
                    tmp_line = line[0:i]
                
                # Skip if line is null
                if len(tmp_line) == 0:
                    continue

                # Remove space
                tmp_line = tmp_line.replace(' ', '')
                
                # Get parameter
                item = tmp_line.split("=")
                if len(item) == 2:
                    (name, val) = (item[0], item[1])
                    self._param[name] = val
                else:
                    print("ERROR: erroneous line in " + \
                          self._param_file + ": " + line, \
                          file=sys.stderr)
                    sys.exit(1)
        
    #---------------------------------------------------------------
    
    def _plot_n_layers(self, fig, ax):
        file = "n_layers.ppd"
        xlabel = "# of layers"
        ylabel = "Posterior probability"
        df = pd.read_csv(file, delim_whitespace=True, header=None, \
                         names=(xlabel, ylabel))
        df.plot(x=xlabel, y=ylabel, ax=ax, kind="area", legend=None)
        ax.set_ylabel(ylabel)

    #---------------------------------------------------------------

    def _plot_vz(self, fig, ax, mode):
        param = self._param
        if mode == "vs":
            vlabel = "S wave velocity (km/s)"
            file = "vs_z.ppd"
            v_min = float(param["vs_min"])
            v_max = float(param["vs_max"])
            nbin_v = int(param["nbin_vs"])
        elif mode == "vp":
            vlabel = "P wave velocity (km/s)"
            file = "vp_z.ppd"
            v_min = float(param["vp_min"])
            v_max = float(param["vp_max"])
            nbin_v = int(param["nbin_vp"])
            
        zlabel = "Depth (km)"
        plabel = "Posterior probability"
        del_v = (v_max - v_min) / nbin_v
        z_min = 0.0
        z_max = float(param["z_max"])
        nbin_z = int(param["nbin_z"])
        del_z = (z_max - z_min) / nbin_z
        
        df = pd.read_csv(file, delim_whitespace=True, header=None, \
                         names=(vlabel, zlabel, plabel))
        
        z, v = np.mgrid[slice(z_min, z_max + eps, del_z), \
                        slice(v_min, v_max + eps, del_v)]
        
        data = df.pivot(zlabel, vlabel, plabel)
        mappable = ax.pcolormesh(v, z, data, cmap='hot_r')
        cbar = fig.colorbar(mappable, ax=ax)
        cbar.ax.set_ylabel(plabel)
        
        ax.set_xlabel(vlabel)
        ax.set_ylabel(zlabel)
        ax.set_ylim([z_max, 0])
        

    #---------------------------------------------------------------

    def _plot_disp(self, fig, ax, mode):
        param = self._param
        clabel = "Phase velocity (km/s)"
        ulabel = "Group velocity (km/s)"
        if mode == "c":
            file = "f_c.ppd"
            vlabel = clabel
        elif mode == "u":
            file = "f_u.ppd"
            vlabel = ulabel
        flabel = "Frequency (Hz)"
        plabel = "Posterior probability"
        df = pd.read_csv(file, delim_whitespace=True, header=None, \
                         names=(flabel, vlabel, plabel))
        v_min = float(param["cmin"])
        v_max = float(param["cmax"])
        del_v = float(param["dc"])
        f_min = float(param["fmin"])
        del_f = float(param["df"])
        nf    = int(param["nf"])
        f_max = f_min + nf * del_f
        
        v, f = np.mgrid[slice(v_min - 0.5 * del_v, \
                              v_max + 0.5 *  del_v, \
                              del_v), \
                        slice(f_min - 0.5 * del_f, \
                              f_max + 0.5 * del_f, \
                              del_f)]
        data = df.pivot(vlabel, flabel, plabel)
        
        mappable = ax.pcolor(f, v, data, cmap='hot_r')
        ax.set_xlabel(flabel)
        ax.set_ylabel(vlabel)
        cbar = fig.colorbar(mappable, ax=ax)
        cbar.ax.set_ylabel(plabel)
        
        # Observation
        df = pd.read_csv(param["obs_in"], delim_whitespace=True, \
                         header=0, \
                         names=(clabel, "c_err", ulabel, "u_err"))
        df[flabel] = np.arange(f_min, f_max, del_f)
        if mode == "c":
            df.plot.scatter(flabel, clabel, ax=ax)
        elif mode == "u":
            df.plot.scatter(flabel, ulabel, ax=ax)
            
    #---------------------------------------------------------------       

    def _read_obs_header(self):
        param = self._param
        file = param["obs_in"]
        with open(file, 'r') as f:
            line = f.readline()
            item = line.split()
            self._param["nf"] = item[0]
            self._param["fmin"] = item[1]
            self._param["df"] = item[2]
        
                         
    #---------------------------------------------------------------
    
    def _plot_likelihood_history(self, fig, ax):
        param = self._param
        file = "likelihood.history"
        df = pd.read_csv(file, delim_whitespace=True, header=None)
        df.plot(ax=ax, legend=None, linewidth=0.6)
        ax.set_ylim([-3000,0])
        ax.set_xlabel("Iteration #")
        ax.set_ylabel("Log-likelihood")

    #---------------------------------------------------------------

    def _plot_temp_history(self, fig, ax):
        param = self._param
        file = "temp.history"
        df = pd.read_csv(file, delim_whitespace=True, header=None)
        df.plot(ax=ax, legend=None, linewidth=0.6)
        
        ax.set_xlabel("Iteration #")
        ax.set_ylabel("Temperature")

    #---------------------------------------------------------------
    
    def draw_figure(self):
        grid_geom = (4, 4)
        fig_size = (22, 13)
        param = self._param
        sns.set()
        sns.set_style('ticks')
        fig = plt.figure(figsize=fig_size)
        
        fig.subplots_adjust(wspace=0.8, hspace=0.4)
        
        # Number of layers
        ax = plt.subplot2grid(grid_geom, (0, 0), colspan=2, fig=fig)
        self._plot_n_layers(fig, ax)
        
        # Phase velocity
        ax = plt.subplot2grid(grid_geom, (1, 0), fig=fig)
        self._plot_disp(fig, ax, "c")

        # Group velocity
        ax = plt.subplot2grid(grid_geom, (1, 1), fig=fig)
        self._plot_disp(fig, ax, "u")

        # Vs-z
        ax = plt.subplot2grid(grid_geom, (2, 0), rowspan=2, fig=fig)
        self._plot_vz(fig, ax, "vs")

        # Vp-z
        if param["solve_vp"].lower() == ".true.":
            ax = plt.subplot2grid(grid_geom, (2, 1), rowspan=2, fig=fig)
            self._plot_vz(fig, ax, "vp")
            
        # Likelihood history
        ax = plt.subplot2grid(grid_geom, (0, 2), colspan=2, rowspan=1, \
                              fig=fig)
        self._plot_likelihood_history(fig, ax)

        # Temperature hisotry
        ax = plt.subplot2grid(grid_geom, (1, 2), colspan=2, rowspan=1, \
                              fig=fig)
        self._plot_temp_history(fig, ax)

        plt.show()
        
#-----------------------------------------------------------------------

if __name__ == "__main__":
    # Get argument
    args = sys.argv
    if len(args) != 2:
        print("USAGE: plot.py [parameter file]")
        sys.exit()
        
    param_file = args[1]
    rslt = InvRslt(param_file)
    rslt.draw_figure()
    


