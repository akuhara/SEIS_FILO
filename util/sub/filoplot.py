import os 
import sys
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import struct
from matplotlib.ticker import MaxNLocator

class InvResult:
    def __init__(self, param_file, vs_true=None):
        
        # Check if file exists or not
        if os.path.exists(param_file) == False:
            print("ERROR: " + param_file + " is not found", \
                  file=sys.stderr)
            sys.exit(1)
        self._param_file = param_file

        # Get parameters
        self._param = {}
        self._read_param_file()

        if not vs_true is None:
            if os.path.exists(vs_true) == False:
                print("ERROR: " + vs_true + " is not found", \
                      file=sys.stderr)
                sys.exit(1)
            else:
                self._param["vs_true"] = vs_true
        
        return

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
        ylabel = "Probability"
        df = pd.read_csv(file, delim_whitespace=True, header=None, \
                         names=(xlabel, ylabel))
        df.plot(x=xlabel, y=ylabel, ax=ax, kind="area", legend=None)
        ax.set_ylabel(ylabel)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    #---------------------------------------------------------------

    def _plot_sigma(self, fig, ax, mode, trace_id):
        param = self._param
        clabel = "Phase vel. (km/s)"
        c_used = "Phase velocity is used?"
        ulabel = "Group vel. (km/s)"
        u_used = "Group velocity is used?"
        hvlabel = "H/V"
        hv_used = "H/V is used?"
        
        if mode == "group":
            file = "group_sigma" + str(trace_id).zfill(3) + ".ppd"
            used_label = u_used
        elif mode == "phase":
            file = "phase_sigma" + str(trace_id).zfill(3) + ".ppd"
            used_label = c_used
        elif mode == "hv":
            file = "hv_sigma" + str(trace_id).zfill(3) + ".ppd"
            used_label = hv_used
        elif mode == "recv_func":
            file = "rf_sigma" + str(trace_id).zfill(3) + ".ppd"

        
        if mode == "group" or mode == "phase" or mode == "hv":
            df = pd.read_csv(param["obs_disper_file"],\
                             delim_whitespace=True, \
                             header=None, \
                             names=(clabel, c_used, ulabel, u_used, \
                                    hvlabel, hv_used), \
                             comment='#')
            df_obs = df[df[used_label] == "T"]
        else:
            df_obs = (999,)

        xlabel = "STDV of data noise"
        ylabel = "Probability"
        if len(df_obs) > 0:
            df = pd.read_csv(file, delim_whitespace=True, header=None, \
                             names=(xlabel, ylabel))
            df.plot(x=xlabel, y=ylabel, ax=ax, kind="area", legend=None)
        else:
            ax.text(0.5, 0.5, "N/A", size=40, \
                    horizontalalignment="center", \
                    verticalalignment="center")
            ax.set_xlabel(xlabel)
            
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
        plabel = "Probability"
        del_v = (v_max - v_min) / nbin_v
        z_min = 0.0
        z_max = float(param["z_max"])
        nbin_z = int(param["n_bin_z"])
        del_z = (z_max - z_min) / nbin_z
        
        df = pd.read_csv(file, delim_whitespace=True, header=None, \
                         names=(vlabel, zlabel, plabel))
        v_tmp = df[df[vlabel].duplicated()==False][vlabel]
        z_tmp = df[df[zlabel].duplicated()==False][zlabel]
        df_below_sea = df[df[zlabel] > float(param["z_min"])]
        p_max = df_below_sea[plabel].max()
        data = df.pivot(zlabel, vlabel, plabel)
        mappable = ax.pcolormesh(v_tmp, z_tmp, \
                                 data, shading='nearest', \
                                 cmap='hot_r', vmax=p_max)
        cbar = fig.colorbar(mappable, ax=ax)
        cbar.ax.set_ylabel(plabel)
        
        df = pd.read_csv(mean_file, delim_whitespace=True, \
                         header=None, names=(zlabel, vlabel))
        line, = ax.plot(df[vlabel], df[zlabel], color="blue")
        lines.append(line)
        labels.append("Mean model")


        # Reference model (if applicable)
        if "solve_anomaly" in param and \
           (param["solve_anomaly"] == ".true." or \
            param["solve_anomaly"] == "T"):
            print("Anomaly mode")
            
            ref_file = param["ref_vmod_in"]
            
            xlabel2 = "Depth (km)"
            ylabel2 = "P wave velocity (km/s)"
            zlabel2 = "S wave velocity (km/s)"
            df = pd.read_csv(ref_file, delim_whitespace=True, \
                             header=None, \
                             names=(xlabel2, ylabel2, zlabel2))
            if (mode == "vs"):
                ref_data = df[zlabel2]
            elif (mode == "vp"):
                ref_data = df[ylabel2]

            line, = ax.plot(ref_data, df[xlabel2], color="black", \
                            linestyle="dashed")
            
            lines.append(line)
            labels.append("Reference model")
            
        # True model
        if mode == "vs" and "vs_true" in param:
            vp_true, vs_true, rho_true, z_true = self._read_vmod(param["vs_true"])
            line, = ax.plot(vs_true, z_true, color="red")
            lines.append(line)
            labels.append("True model")

        ax.legend(lines, labels, loc="lower left")
        
        ax.set_xlabel(vlabel)
        ax.set_ylabel(zlabel)
        ax.set_ylim([z_max, 0])
        ax.set_xlim([v_min, v_max])
        
    #---------------------------------------------------------------

    def _read_vmod(self, vmod_file):

        vp = []
        vs = []
        rho = []
        z = []
        z_tmp = 0.0
        count = 1
        with open(vmod_file, 'r') as f:
            for line in f:
                tmp_line = line.rstrip()
                # Dealing with comment out
                i = tmp_line.find('#')
                if i >= 0:
                    tmp_line = line[0:i]
                
                # Skip if line is null
                if len(tmp_line) == 0:
                    continue
                
                item = tmp_line.split()
                if count == 1:
                    nlay = item[0]
                else:
                    print(item)
                    print(item[1])
                    vp.append(float(item[0]))
                    vs.append(float(item[1]))
                    rho.append(float(item[2]))
                    z.append(z_tmp)
                    vp.append(float(item[0]))
                    vs.append(float(item[1]))
                    rho.append(float(item[2]))
                    z_tmp += float(item[3])
                    z.append(z_tmp)
                    
                count += 1

        return vp, vs, rho, z

    #---------------------------------------------------------------
    
    def _plot_recv_func(self, fig, ax, trace_id):
        param = self._param
        ppd_file = "syn_rf" + str(trace_id).zfill(3) + ".ppd"
        obs_file = param["recv_func_in"]

        self._read_recv_func_obs(obs_file, trace_id)
        param = self._param
        
        tlabel = "Time (s)"
        alabel = "RF amp."
        plabel = "Probability"
        
        df = pd.read_csv(ppd_file, delim_whitespace=True, \
                         header=None, \
                         names=(alabel, tlabel, plabel))
        if "amp_min" in param:
            amp_min = float(param["amp_min"])
        else:
            amp_min = -0.6

        if "amp_max" in param:
            amp_max = float(param["amp_max"])
        else:
            amp_max = 0.6
            
        t_tmp = df[df[tlabel].duplicated()==False][tlabel]
        a_tmp = df[df[alabel].duplicated()==False][alabel]
        data = df.pivot(tlabel, alabel, plabel)
        mappable = ax.pcolormesh(a_tmp, t_tmp, data, cmap='hot_r', \
                                 shading='nearest')
        cbar = fig.colorbar(mappable, ax=ax)
        cbar.ax.set_ylabel(plabel) 
        ax.set_xlim([float(param["t_start"]), float(param["t_end"])])

        ax.set_xlabel(tlabel)
        ax.set_ylabel(alabel)

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
                        self._param["syn_rf_file"] = tmp_line.replace('\'','')
                        self._param["syn_rf_file"] = self._param["syn_rf_file"].replace('\"','')
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

        clabel = "Phase vel. (km/s)"
        c_used = "Phase velocity is used?"
        ulabel = "Group vel. (km/s)"
        u_used = "Group velocity is used?"
        hvlabel = "H/V"
        hv_used = "H/V is used?"
        if mode == "c":
            ppd_file = "syn_phase" + str(curve_id).zfill(3) + ".ppd"
            vlabel = clabel
            used_label = c_used
        elif mode == "u":
            ppd_file = "syn_group" + str(curve_id).zfill(3) + ".ppd"
            vlabel = ulabel
            used_label = u_used
        elif mode == "hv":
            ppd_file = "syn_hv" + str(curve_id).zfill(3) + ".ppd"
            vlabel = hvlabel
            used_label = hv_used

        if param["freq_or_period"] == "freq":
            xlabel = "Frequency (Hz)"
        elif param["freq_or_period"] == "period":
            xlabel = "Period (s)"
            
        
        plabel = "Probability"
        df = pd.read_csv(ppd_file, delim_whitespace=True, header=None, \
                         names=(xlabel, vlabel, plabel))
        
        x_min = float(param["xmin"])
        del_x = float(param["dx"])
        nx    = int(param["nx"])
        x_max = x_min + (nx -1) * del_x
        v_tmp = df[df[vlabel].duplicated()==False][vlabel]
        x_tmp = df[df[xlabel].duplicated()==False][xlabel]
        data = df.pivot(vlabel, xlabel, plabel)

        # Observation
        df = pd.read_csv(param["obs_disper_file"],\
                         delim_whitespace=True, \
                         header=None, \
                         names=(clabel, c_used, ulabel, u_used, \
                                hvlabel, hv_used), \
                         comment='#')
        df[xlabel] = np.arange(x_min, x_max + 0.5 * del_x, del_x)
        df_obs = df[df[used_label] == "T"]
        
        # Draw
        if len(df_obs) > 0:
            mappable = ax.pcolormesh(x_tmp, v_tmp, data, cmap='hot_r', \
                                     shading='nearest')
            ax.set_xlabel(xlabel)
            ax.set_ylabel(vlabel)
            cbar = fig.colorbar(mappable, ax=ax)
            cbar.ax.set_ylabel(plabel)
            df_obs.plot.scatter(xlabel, vlabel, ax=ax, s=8, \
                                marker=".", c="blue")
        else:
            ax.text(0.5, 0.5, "N/A", size=40, \
                    horizontalalignment="center", \
                    verticalalignment="center")
            ax.set_xlabel(xlabel)
            ax.set_ylabel(vlabel)

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

                if (count - 2) // 7 + 1 == curve_id:
                    iloc = (count - 2) % 7
                    if iloc == 0:
                        self._param["obs_disper_file"] = tmp_line.replace('\'','')
                        self._param["obs_disper_file"] = self._param["obs_disper_file"].replace('\"','')
                    elif iloc == 1:
                        item = tmp_line.split(" ")
                        self._param["freq_or_period"] = item[2]
                    elif iloc == 2:
                        item = tmp_line.split(" ")
                        self._param["nx"]   = item[0]
                        self._param["xmin"] = item[1]
                        self._param["dx"]   = item[2]
                    elif iloc == 3:
                        item = tmp_line.split(" ")
                        self._param["cmin"] = item[0]
                        self._param["cmax"] = item[1]
                        self._param["dc"] = item[2]


    #---------------------------------------------------------------
    
    def _plot_likelihood_history(self, fig, ax):
        param = self._param
        if "diagnostic_mode" in param:
            diagnostic = param["diagnostic_mode"].lower()
        else:
            diagnostic = ".false."

        if diagnostic == ".true." or diagnostic == "t":
            file = "likelihood.history"
            df = pd.read_csv(file, delim_whitespace=True, header=None)
            n_all = len(df.columns)
            if n_all > 5:
                df.plot(ax=ax, legend=None, linewidth=0.5, color="gray")
                n1 = 0
                nn = n_all // 4
                n2 = nn * 1
                n3 = nn * 2
                n4 = nn * 3
                n5 = n_all - 1
                df = df.iloc[:, [n1, n2, n3, n4, n5]]
                df.plot(ax=ax, legend=None, linewidth=2.0)
                ax.set_ylim([-5000,1000])
                ax.set_xlabel("Iteration #")
                ax.set_ylabel("Log-likelihood")
        else:
            ax.text(0.5, 0.5, "N/A (set diagnostic_mode=.true.)", \
                    size=20, \
                    horizontalalignment="center", \
                    verticalalignment="center")
            ax.set_xlabel("Iteration #")
            ax.set_ylabel("Log-likelihood")
            
    #---------------------------------------------------------------

    def _plot_temp_history(self, fig, ax):
        param = self._param
        if "diagnostic_mode" in param:
            diagnostic = param["diagnostic_mode"].lower()
        else:
            diagnostic = ".false."

        if diagnostic == ".true." or diagnostic == "t":
            file = "temp.history"
            df = pd.read_csv(file, delim_whitespace=True, header=None,
                             dtype=np.float64)
            n_all = len(df.columns)
            if n_all > 5:
                df.plot(ax=ax, legend=None, linewidth=0.5, color="gray")
                n1 = 0
                nn = n_all // 4
                n2 = nn * 1
                n3 = nn * 2
                n4 = nn * 3
                n5 = n_all - 1
            df = df.iloc[:, [n1, n2, n3, n4, n5]]
            df.plot(ax=ax, legend=None, linewidth=2.0)
            
            ax.set_xlabel("Iteration #")
            ax.set_ylabel("Temperature")
        else:
            ax.text(0.5, 0.5, "N/A (set diagnostic_mode=.true.)", \
                    size=20, \
                    horizontalalignment="center", \
                    verticalalignment="center")
            ax.set_xlabel("Iteration #")
            ax.set_ylabel("Temperature")

    #---------------------------------------------------------------
    
    def _plot_proposal_count(self, fig, ax):
        param = self._param
        file = "proposal.count"
        df = pd.read_csv(file, delim_whitespace=True, header=None, \
                         index_col=0)
    
        #if param["solve_vp"].lower() == ".true.":
        #    df.index = ['Birth', 'Death', 'Depth', 'Vs', 'Vp']
        #elif param["solve_vp"].lower() == "t":
        #    df.index = ['Birth', 'Death', 'Depth', 'Vs', 'Vp']
        #else:
        #    df.index = ['Birth', 'Death', 'Depth', 'Vs']
        
        df.columns = ['Proposed', 'Accepted']
        df.index.name = 'Proposal type'
        df.plot(kind='bar', ax=ax)


    #---------------------------------------------------------------

    
    def draw_surface_wave_set(self, curve_id):
        grid_geom = (5, 6)
        fig_size = (17, 10)
        param = self._param
        sns.set()
        sns.set_style('ticks')
        fig = plt.figure(figsize=fig_size)
        
        fig.subplots_adjust(wspace=0.7, hspace=0.9)
        
        # Number of layers
        ax = plt.subplot2grid(grid_geom, (0, 0), colspan=3, fig=fig)
        self._plot_n_layers(fig, ax)
        
        # Phase velocity
        ax = plt.subplot2grid(grid_geom, (1, 0), fig=fig)
        self._plot_dispersion(fig, ax, "c", curve_id)

        # Group velocity
        ax = plt.subplot2grid(grid_geom, (1, 1), fig=fig)
        self._plot_dispersion(fig, ax, "u", curve_id)

        # H/V
        ax = plt.subplot2grid(grid_geom, (1, 2), fig=fig)
        self._plot_dispersion(fig, ax, "hv", curve_id)

        # Sigma phase velocity
        ax = plt.subplot2grid(grid_geom, (2, 0), fig=fig)
        self._plot_sigma(fig, ax, "phase", curve_id)

        # Sigma group velocity
        ax = plt.subplot2grid(grid_geom, (2, 1), fig=fig)
        self._plot_sigma(fig, ax, "group", curve_id)

        # Sigma H/V
        ax = plt.subplot2grid(grid_geom, (2, 2), fig=fig)
        self._plot_sigma(fig, ax, "hv", curve_id)

        # Vs-z
        ax = plt.subplot2grid(grid_geom, (3, 0), rowspan=2, fig=fig)
        self._plot_vz(fig, ax, "vs")

        # Vp-z
        if param["solve_vp"].lower() == ".true." or param["solve_vp"].lower() == "t":
            ax = plt.subplot2grid(grid_geom, (3, 1), rowspan=2, fig=fig)
            self._plot_vz(fig, ax, "vp")
            
        # Likelihood history
        ax = plt.subplot2grid(grid_geom, (0, 3), colspan=3, rowspan=1, \
                              fig=fig)
        self._plot_likelihood_history(fig, ax)

        # Temperature hisotry
        ax = plt.subplot2grid(grid_geom, (1, 3), colspan=3, rowspan=1, \
                              fig=fig)
        self._plot_temp_history(fig, ax)

        # Proposal count
        ax = plt.subplot2grid(grid_geom, (2, 3), colspan=3, rowspan=2, \
                              fig=fig)
        self._plot_proposal_count(fig, ax)

        # Output
        out_file = "disper" + str(curve_id).zfill(2) + ".png"
        print("Output: " + out_file)
        fig.savefig(out_file)
        
        
        
    #---------------------------------------------------------------
    def draw_recv_func_set(self, trace_id):
        grid_geom = (5, 4)
        fig_size = (17, 10)
        param = self._param
        sns.set()
        sns.set_style('ticks')
        fig = plt.figure(figsize=fig_size)
        
        fig.subplots_adjust(wspace=0.7, hspace=0.9)
        
        # Number of layers
        ax = plt.subplot2grid(grid_geom, (0, 0), colspan=2, fig=fig)
        self._plot_n_layers(fig, ax)
        
        # Receiver function
        ax = plt.subplot2grid(grid_geom, (1, 0), colspan=1, fig=fig)
        self._plot_recv_func(fig, ax, trace_id)

        # Receiver function sigma
        ax = plt.subplot2grid(grid_geom, (2, 0), colspan=1, fig=fig)
        self._plot_sigma(fig, ax, "recv_func", trace_id)

        # Vs-z
        ax = plt.subplot2grid(grid_geom, (3, 0), rowspan=2, fig=fig)
        self._plot_vz(fig, ax, "vs")

        # Vp-z
        if param["solve_vp"].lower() == ".true." or param["solve_vp"].lower() == "t":
            ax = plt.subplot2grid(grid_geom, (3, 1), rowspan=2, fig=fig)
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
