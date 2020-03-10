import os 
import sys
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import struct

class InvResult:
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
        lines = [] # for legend
        labels = []
        eps = 1.0e-8
        if mode == "vs":
            vlabel = "S wave velocity (km/s)"
            file = "vs_z.ppd"
            mean_file = "vs_z.mean"
            v_min = float(param["vs_min"])
            v_max = float(param["vs_max"])
            nbin_v = int(param["n_bin_vs"])
        elif mode == "vp":
            vlabel = "P wave velocity (km/s)"
            file = "vp_z.ppd"
            mean_file = "vp_z.mean"
            v_min = float(param["vp_min"])
            v_max = float(param["vp_max"])
            nbin_v = int(param["n_bin_vp"])
            
        zlabel = "Depth (km)"
        plabel = "Posterior probability"
        del_v = (v_max - v_min) / nbin_v
        z_min = 0.0
        z_max = float(param["z_max"])
        nbin_z = int(param["n_bin_z"])
        del_z = (z_max - z_min) / nbin_z
        
        df = pd.read_csv(file, delim_whitespace=True, header=None, \
                         names=(vlabel, zlabel, plabel))
        z, v = np.mgrid[slice(z_min, z_max + eps, del_z), \
                        slice(v_min, v_max + eps, del_v)]
        df_below_sea = df[df[zlabel] > float(param["z_min"])]
        p_max = df_below_sea[plabel].max()
        #print(p_max)
        data = df.pivot(zlabel, vlabel, plabel)
        mappable = ax.pcolormesh(v, z, data, cmap='hot_r', vmax=p_max)
        cbar = fig.colorbar(mappable, ax=ax)
        cbar.ax.set_ylabel(plabel)
        
        df = pd.read_csv(mean_file, delim_whitespace=True, \
                         header=None, names=(zlabel, vlabel))
        line, = ax.plot(df[vlabel], df[zlabel], color="blue")
        lines.append(line)
        labels.append("Mean model")
        ax.legend(lines, labels, loc="lower left")
        
        ax.set_xlabel(vlabel)
        ax.set_ylabel(zlabel)
        ax.set_ylim([z_max, 0])
        ax.set_xlim([v_min, v_max])
        
    #---------------------------------------------------------------
    
    def _plot_recv_func(self, fig, ax, trace_id):
        param = self._param
        ppd_file = "syn_rf" + str(trace_id).zfill(3) + ".ppd"
        obs_file = param["recv_func_in"]

        self._read_recv_func_obs(obs_file, trace_id)
        param = self._param
        
        tlabel = "Time (s)"
        alabel = "RF amp."
        plabel = "Posterior probability"
        
        df = pd.read_csv(ppd_file, delim_whitespace=True, \
                         header=None, \
                         names=(alabel, tlabel, plabel))
        amp_min = -1.0
        amp_max = 1.0
        n_bin_amp = 100
        del_amp = (amp_max - amp_min) / n_bin_amp
        
        t_min = float(param["t_start"])
        t_max = 2 * float(param["t_end"])
        del_t = float(param["delta"])
        
        a, t = np.mgrid[slice(amp_min, \
                              amp_max + del_amp, \
                              del_amp), \
                        slice(t_min - 0.5 * del_t, \
                              t_max + 0.5 * del_t, \
                              del_t)]
        data = df.pivot(tlabel, alabel, plabel)
        mappable = ax.pcolormesh(t, a, data, cmap='hot_r')
        cbar = fig.colorbar(mappable, ax=ax)
        cbar.ax.set_ylabel(plabel) 
        ax.set_xlim([t_min, t_max / 2])

        # plot observation
        b, delta, t, data = self._read_sac(param["syn_rf_file"])
        ax.plot(t, data, color="black")
        
        
    #---------------------------------------------------------------
    
    def _read_recv_func_obs(self, obs_file, trace_id):
        with open(obs_file, 'r') as f:
            count = 0
            for line in f:
                tmp_line = line.rstrip()
                
                # Dealing with comment out
                i = tmp_line.find('#')
                if i >= 0:
                    tmp_line = line[0:i]
                
                # Skip if line is null
                if len(tmp_line) == 0:
                    continue
                
                count += 1
                
                # Skip a line of N_trc
                if count == 1:
                    continue
                
                if (count - 2) // 6 + 1 == trace_id:
                    iloc = (count - 2) % 6
                    if iloc == 0:
                        self._param["syn_rf_file"] = tmp_line
                    elif iloc == 1:
                        item = tmp_line.split(" ")
                        self._param["delta"] = item[2]
                    elif iloc == 3:
                        item = tmp_line.split(" ")
                        self._param["t_start"] = item[0]
                        self._param["t_end"] = item[1]

    #---------------------------------------------------------------

    def _read_sac(self, sac_file):

        endian = "<"

        # read entire part of SAC file
        f = open(sac_file, 'rb')
        sac_data = f.read()

        # check endian
        header = struct.unpack_from(endian + 'i', sac_data, offset=4*76)
        if (header[0] != 6):
            endian = ">"
            
        # get headers
        header = struct.unpack_from(endian + 'f', sac_data, offset=0)
        delta = header[0]
        header = struct.unpack_from(endian + 'f', sac_data, offset=4*5)
        b = header[0]
        header = struct.unpack_from(endian + 'i', sac_data, offset=4*79)
        npts = header[0]
        
        # get data
        data = struct.unpack_from(endian + 'f' * npts, sac_data, offset=4*158)
        
        # make time index
        t = np.arange(b, b + npts * delta, delta)
        
        return b, delta, t, data
    
    #---------------------------------------------------------------

    def _plot_dispersion(self, fig, ax, mode, curve_id):
        param = self._param
        obs_file = param["disper_in"]
        self._read_dispersion_obs(obs_file, curve_id)

        clabel = "Phase velocity (km/s)"
        ulabel = "Group velocity (km/s)"
        if mode == "c":
            ppd_file = "syn_phase" + str(curve_id).zfill(3) + ".ppd"
            vlabel = clabel
        elif mode == "u":
            ppd_file = "syn_group" + str(curve_id).zfill(3) + ".ppd"
            vlabel = ulabel
        flabel = "Frequency (Hz)"
        plabel = "Posterior probability"
        print(ppd_file)
        df = pd.read_csv(ppd_file, delim_whitespace=True, header=None, \
                         names=(flabel, vlabel, plabel))
        v_min = float(param["cmin"])
        v_max = float(param["cmax"])
        del_v = float(param["dc"])
        f_min = float(param["fmin"])
        del_f = float(param["df"])
        nf    = int(param["nf"])
        f_max = f_min + (nf -1) * del_f
        print(v_min, v_max, del_v, f_min, del_f, nf, f_max)
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
        df = pd.read_csv(param["obs_disper_file"],\
                         delim_whitespace=True, \
                         header=None, \
                         names=(clabel, "c_err", ulabel, "u_err"), \
                         comment='#')
        print(f_min, f_max, del_f)
        print(np.arange(f_min, f_max, del_f))
        df[flabel] = np.arange(f_min, f_max + 0.5 * del_f, del_f)

        df.plot.scatter(flabel, vlabel, ax=ax, s=5, marker=".", \
                        c="blue")

    #---------------------------------------------------------------       
    def _read_dispersion_obs(self, obs_file, curve_id):
        with open(obs_file, 'r') as f:
            count = 0
            for line in f:
                tmp_line = line.rstrip()
                
                # Dealing with comment out
                i = tmp_line.find('#')
                if i >= 0:
                    tmp_line = line[0:i]
                
                # Skip if line is null
                if len(tmp_line) == 0:
                    continue
                
                count += 1
                
                # Skip a line of N_disp
                if count == 1:
                    continue

                if (count - 2) // 4 + 1 == curve_id:
                    iloc = (count - 2) % 4
                    if iloc == 0:
                        self._param["obs_disper_file"] = tmp_line
                    elif iloc == 2:
                        item = tmp_line.split(" ")
                        self._param["nf"]   = item[0]
                        self._param["fmin"] = item[1]
                        self._param["df"]   = item[2]
                    elif iloc == 3:
                        item = tmp_line.split(" ")
                        self._param["cmin"] = item[0]
                        self._param["cmax"] = item[1]
                        self._param["dc"] = item[2]


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
    
    def _plot_proposal_count(self, fig, ax):
        param = self._param
        file = "proposal.count"
        df = pd.read_csv(file, delim_whitespace=True, header=None, \
                         index_col=0)
    
        if param["solve_vp"].lower() == ".true.":
            df.index = ['Birth', 'Death', 'Depth', 'Vs', 'Vp']
        elif param["solve_vp"].lower() == "t":
            df.index = ['Birth', 'Death', 'Depth', 'Vs', 'Vp']
        else:
            df.index = ['Birth', 'Death', 'Depth', 'Vs']
        
        df.columns = ['Proposed', 'Accepted']
        df.plot(kind='bar', ax=ax)
        

    #---------------------------------------------------------------

    
    def draw_surface_wave_set(self, curve_id):
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
        self._plot_dispersion(fig, ax, "c", curve_id)

        # Group velocity
        ax = plt.subplot2grid(grid_geom, (1, 1), fig=fig)
        self._plot_dispersion(fig, ax, "u", curve_id)

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

        # Proposal count
        ax = plt.subplot2grid(grid_geom, (2, 2), colspan=2, rowspan=2, \
                              fig=fig)
        self._plot_proposal_count(fig, ax)

        # Output
        out_file = "disper" + str(curve_id).zfill(2) + ".png"
        print("Output: " + out_file)
        fig.savefig(out_file)
        
        
        
    #---------------------------------------------------------------
    def draw_recv_func_set(self, trace_id):
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
        
        # Receiver function
        ax = plt.subplot2grid(grid_geom, (1, 0), colspan=2, fig=fig)
        self._plot_recv_func(fig, ax, trace_id)

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

        # Proposal count
        ax = plt.subplot2grid(grid_geom, (2, 2), colspan=2, rowspan=2, \
                              fig=fig)
        self._plot_proposal_count(fig, ax)

        # Output
        out_file = "recv_func" + str(trace_id).zfill(2) + ".png"
        print("Output: " + out_file)
        fig.savefig(out_file)
               
#-----------------------------------------------------------------------
