"""
lamipy project - laminated composites calculations in Python.

plotresults.py - Module containing a class for results plotting.

INPUTS:
lam -- laminate layup characteristics
plot_data -- contains pairs of Load Factor vs Results (stress & strain)
fail_status -- contains data about layer failures
**options -- keywords which set plotting options

OUTPUTS:
.png savefile -- saves figures and plots.
.mp4 animation -- saves video of failure progress.

Joao Paulo Bernhardt - October 2017
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import clt

class PlotResults:
    """ Contains functions for setting options and plotting data.

    Initialization requires:
    lam -- laminate layup characteristics
    plot_data -- contains pairs of Load Factor vs Results (stress & strain)
    fail_status -- contains data about layer failures
    """

    def __init__(self, lam, plot_data, fail_status):
        # Sets obligatory data for the instance.
        self.lam = lam
        self.plot_data = plot_data
        self.fail_status = fail_status

        # Gets number of layers and loadsteps
        self.num_layers = len(lam["ang"])
        self.num_steps = len(plot_data)

        # Isolates stress & strain results from array
        self.res = plot_data[:,1]

        # Initiate default values for Options
        self.save = False
        self.display = True

    def Options(self, **options):
        """ Optional function for determining plot options."""

        # Should lamipy save the plots?
        if isinstance(options["save"], bool):
            self.save = options["save"]

        # Should lamipy display the plots during runtime?
        if isinstance(options["display"], bool):
            self.display = options["display"]


    def ProgAvgStrain(self):
        """ Plots Average Strain (x-axis) vs. Load Factor (y-axis). """

        # Sets plot dimension
        fig = plt.figure()
        #plt.figure(figsize=(10, 8))

        mean_eps = np.zeros((self.num_steps))
        fpf_lf = min(self.fail_status["Load Factor"])
        fpf_strain = 0
        lpf_lf = max(self.fail_status["Load Factor"])
        lpf_strain = 0

        # Iterates the loadsteps
        for step in range(self.num_steps):
            mean_eps[step] = np.mean(np.union1d(
                                    self.res[step]["LCS"]["strain"]["sup"][0], 
                                    self.res[step]["LCS"]["strain"]["inf"][0]))
            # Finds FPF and LPF
            if self.plot_data[step, 0] == fpf_lf and fpf_strain == 0:
                fpf_strain = mean_eps[step]
            elif self.plot_data[step, 0] == lpf_lf and lpf_strain == 0:
                lpf_strain = mean_eps[step]

        # Formats the plot
        plt.plot(mean_eps, self.plot_data[:,0])
        plt.xlabel('Avg. Strain (x axis)')
        plt.ylabel('Load Factor')
        plt.title('Avg. Strain vs Load Factor')
        plt.grid(True)

        # Points LPF and FPF
        plt.annotate('FPF', xy=(fpf_strain, fpf_lf), 
                     xytext=(fpf_strain, fpf_lf),
                     arrowprops=dict(facecolor='green'),
                    )
        plt.annotate('LPF', xy=(lpf_strain, lpf_lf), 
                     xytext=(lpf_strain, lpf_lf),
                     arrowprops=dict(facecolor='red'),
                    )

        # Displays the plot
        if self.display:
            plt.show()

        if self.save:
            fig.savefig("plots/strain_lf.png")
        
        plt.close(fig)
  
    def Profile(self, coord_sys, axis, var, step):
        """ Plots stress/strain in LCS or MCS against Z vector. """
 
        # Sets plot dimension
        fig = plt.figure()
        #plt.figure(figsize=(10, 8))

        Z = clt.assemble_Z(self.lam)
        
        X = { "inf" : np.zeros((self.num_layers)), 
              "sup" : np.zeros((self.num_layers)) }
        Y = { "inf" : np.zeros((self.num_layers)), 
              "sup" : np.zeros((self.num_layers)) }

        P = np.zeros((self.num_layers*2, 2))

        # Sets the proper names
        if coord_sys == "MCS":
            if axis == 0:
                axis_name = "1"
            elif axis == 1:
                axis_name = "2"
            else:
                axis_name = "6"
        else:
            if axis == 0:
                axis_name = "x"
            elif axis == 1:
                axis_name = "y"
            else:
                axis_name = "xy"

        # Formats the plot
        plt.title('Profile ' + coord_sys + '-' + axis_name + ' ' + var)
        plt.xlabel(var + ' (' + coord_sys + '-' + axis_name + ')')
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.ylabel('Z coordinate')        
        plt.grid(True)
        
        # Iterates the layers, add data to plot
        for layer in range(self.num_layers):
            X["inf"][layer] = self.res[step][coord_sys][var]["inf"][axis][layer]
            Y["inf"][layer] = Z[layer]
            X["sup"][layer] = self.res[step][coord_sys][var]["sup"][axis][layer]
            Y["sup"][layer] = Z[layer + 1]
            P[layer*2] = [X["inf"][layer], Y["inf"][layer]]
            P[layer*2 + 1] = [X["sup"][layer], Y["sup"][layer]]
            plt.fill_betweenx([P[layer*2, 1], P[layer*2 + 1, 1]],
                              [P[layer*2, 0], P[layer*2 + 1, 0]],
                              hatch="//", facecolor="none", 
                              edgecolor="r", lw=1.0)

        # Adds main lines    
        plt.plot(P[:,0], P[:,1], color="r", lw=2.5)
        plt.plot([0]*(self.num_layers+1), Z, color="b", lw=2.5)        
        
        # Displays the plot
        if self.display:
            plt.show()

        if self.save:
            fig.savefig("plots/profile.png")
        
        plt.close(fig)

        """

    def ProfileAn(step, coord_sys, axis, var):
    # NOT READY YET.
    # Function for displaying animated plot of profile.

        Z = clt.assemble_Z(lam)
        
        X = { "inf" : np.zeros((num_layers)), "sup" : np.zeros((num_layers)) }
        Y = { "inf" : np.zeros((num_layers)), "sup" : np.zeros((num_layers)) }
        P = np.zeros((num_layers*2, 2))

        if coord_sys == "MCS":
            if axis == 0:
                axis_name = "1"
            elif axis == 1:
                axis_name = "2"
            else:
                axis_name = "6"
        else:
            if axis == 0:
                axis_name = "x"
            elif axis == 1:
                axis_name = "y"
            else:
                axis_name = "xy"

        plt.title('Profile ' + coord_sys + '-' + axis_name + ' ' + var)
        plt.xlabel(var + ' (' + coord_sys + '-' + axis_name + ')')
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.ylabel('Z coordinate')        
        plt.grid(True)
        
        for layer in range(num_layers):
            X["inf"][layer] = res[step][coord_sys][var]["inf"][axis][layer]
            Y["inf"][layer] = Z[layer]
            X["sup"][layer] = res[step][coord_sys][var]["sup"][axis][layer]
            Y["sup"][layer] = Z[layer + 1]
            P[layer*2] = [X["inf"][layer], Y["inf"][layer]]
            P[layer*2 + 1] = [X["sup"][layer], Y["sup"][layer]]
            plt.fill_betweenx([P[layer*2, 1], P[layer*2 + 1, 1]],
                              [P[layer*2, 0], P[layer*2 + 1, 0]],
                              hatch="//", facecolor="none", edgecolor="r", lw=1.0)
        plt.plot(P[:,0], P[:,1], color="r")
        plt.plot([0]*(num_layers+1), Z, color="b")        
        
        return plt.plot(P[:,0], P[:,1], color="r")

        pass

    def AnimatedProfile():
    # NOT READY YET.
    # Function for displaying animated plot of profile.
        fig = plt.figure()

        a = FuncAnimation(fig, ProfileAn, num_steps, fargs=('MCS', 0, 'strain'))
        a.save('temp.mp4')
        

    # Calls plotting functions
    #ProgAvgStrain()
    #Profile("MCS", 1, "strain", 0)
    #Profile("MCS", 2, "strain", 0)
    #Profile("LCS", 0, "strain", 0)
    #Profile("LCS", 1, "strain", 0)

    # Saves specific plot
    sfig = Profile("MCS", 0, "strain", 0)
    sfig.savefig('MCS_1_strain.png', bbox_inches='tight')


"""