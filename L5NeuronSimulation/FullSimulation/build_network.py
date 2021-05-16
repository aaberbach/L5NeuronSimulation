"""
The script used to build the network for the simulation run in run_network.py
"""
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 
sys.path.insert(0, currentdir)

from bmtk.builder import NetworkBuilder
from bmtk.utils.sim_setup import build_env_bionet
import numpy as np
import sys
import synapses   
import h5py
import pandas as pd   
import json   
import scipy.stats as st
from functools import partial

from raster_maker import *
from clustering import *

#np.random.seed(2129)

class SimulationBuilder:
        """Class used to build our BMTK simulation.

        Attributes
        ----------
        params : dict
                contains parameters for the network
        syn : dict
                contains synaptic templates

        n_dend_exc : int
                number of excitatory input cells on the basal dendrites
        n_apic_exc : int
                number of excitatory input cells on the apical dendrites

        n_dend_inh : int
                number of inhibitory (SOM+) input cells on the basal dendrites 
                more than 50 um from the soma.
        n_apic_inh : int
                number of inhibitory (SOM+) input cells on the apical dendrites
        n_prox_dend_inh : int
                number of inhibitory (PV+) input cells on the basal dendrites
                less than 50 um from the soma
        n_soma_inh : int
                number of inhibitory (PV+) input cells on the soma

        clust_per_group : int
                number of clusters per functional group

        net : NetworkBuilder
                the BMTK network for the biophysical cell
        exc_stim : NetworkBuilder
                the BMTK network for excitatory inputs
        prox_inh_stim : NetworkBuilder
                the BMTK network for perisomatic inhibition
        dist_inh_stim : NetworkBuilder
                the BMTK network for dendritic inhibition

        dend_groups : list
                all excitatory functional groups on the basal dendrites
        apic_groups : list
                all excitatory functional groups on the apical dendrites

        Methods
        -------
        build()
                builds the network
        save_groups()
                saves the functional groups to a csv

        _build_exc()
                creates excitatory input nodes and edges
        _build_exc_nodes(segs : pandas.DataFrame, base_name : str, n_cells : int, start=0 : int)
                builds excitatory nodes
        _build_exc_edges(group_list : list)
                builds excitatory edges

        _save_nets()
                builds and saves the BMTK NetworkBuilders

        _build_inh()
                creates inhibitory input nodes and edges
        
        _make_rasters()
                creates the inhibitory and excitatory input rasters
        _gen_exc_spikes(fname : str)
                generates and saves the excitatory spike rasters
        _gen_inh_spikes(n_cells : int, mean_fr : float, std_fr : float, key : str, fname : str)
                creates inhibitory spike rasters, using a noise trace based on averaging excitation and shifting it

        Static Methods
        --------------
        _connector_func(sources : list, targets : list, cells : list)
                sets the number of synapses from the given cells
        _set_location(source : dict, target : dict, cells : list, start_id : int)
                sets the location of the given edge

        _norm_connect(source : dict, target : dict, m : float, s : float, low : int, high : int)
                used to normally distribute connection counts

        _gen_group_spikes(writer : SonataWriter, group : FunctionalGroup, seconds : float, start_time : float, dist : func)
                creates and saves a functional group's spike raster
        _norm_rvs(mean : float, std : float)
                generates a random float from a normal distribution with a near zero minimum
        """

        def __init__(self, params_file):
                """Initializes the simulation builder, 
                setting up attributes but not actually building the BMTK network.

                Parameters
                ----------
                params_file : str
                    path to the JSON file with network parameters
                """     
                #Loads the JSON file with information about the network.           
                with open(params_file) as f:
                        self.params = json.load(f)

                #Loads synapse templates.
                synapses.load()
                self.syn = synapses.syn_params_dicts()

                avg_exc_div = np.mean(list(self.params["divergence"]["exc"].values()))

                self.n_dend_exc = int((self.params["lengths"]["basal_dist"] * self.params["syn_density"]["exc"]) / avg_exc_div)
                self.n_apic_exc = int((self.params["lengths"]["apic"] * self.params["syn_density"]["exc"]) / avg_exc_div)

                self.n_dend_inh = int((self.params["lengths"]["basal_dist"] * self.params["syn_density"]["inh"]) / self.params["divergence"]["basal_inh"]["m"])
                self.n_apic_inh = int((self.params["lengths"]["apic"] * self.params["syn_density"]["inh"]) / self.params["divergence"]["apic_inh"]["m"])

                self.n_prox_dend_inh = int((self.params["lengths"]["basal_prox"] * self.params["syn_density"]["inh"]) / self.params["divergence"]["peri_inh"]["m"])
                self.n_soma_inh = int(self.params["n_soma_syns"] / self.params["divergence"]["peri_inh"]["m"])

                self.clust_per_group = int((self.params["groups"]["cells_per_group"] * avg_exc_div) // (self.params["syn_density"]["exc"] * 10))

        def build(self):
                """Builds the nodes and edges for the network.
                """                
                self.net = NetworkBuilder("biophysical")

                self.net.add_nodes(N=1, pop_name='Pyrc',
                        potental='exc',
                        model_type='biophysical',
                        model_template='hoc:L5PCtemplate',
                        morphology = None)

                self._build_exc()
                self._build_inh()
                self._save_nets()

                self._make_rasters()
                
                #Final build step.
                build_env_bionet(base_dir='./',
                        network_dir='./network',
                        dt = self.params["dt"], 
                        tstop=self.params["time"]["stop"] * 1000.0,
                        report_vars=['v'],
                        dL = self.params["dL"],#target length (um) of segments
                        spikes_threshold=-10,
                        spikes_inputs=[('exc_stim', 'exc_stim_spikes.h5'), ('prox_inh_stim', 'prox_inh_stim_spikes.h5'), ('dist_inh_stim', 'dist_inh_stim_spikes.h5')],
                        components_dir='../biophys_components',
                        compile_mechanisms=True)

        def save_groups(self):
                """saves the apic and dend groups into a csv.
                one row for each node containgin the id of the functional group it is in.
                """        
                all_groups = self.dend_groups + self.apic_groups
                node_ids = []
                func_groups = []

                for func_id, group in enumerate(all_groups):
                        for i in range(group.start_id, group.start_id + group.n_cells):
                                node_ids.append(i)
                                func_groups.append(func_id)

                df = pd.DataFrame()
                df["Node ID"] = node_ids
                df["Functional Group"] = func_groups
                df.to_csv("FunctionalGroups.csv", index=False)  

        def _build_exc(self):
                """Builds the excitatory input cells and their synapses.
                """          

                np.random.seed(1000)      
                # External excitatory inputs
                self.exc_stim = NetworkBuilder('exc_stim')

                #DataFrame of all segments on the cell.
                segs = pd.read_csv("Segments.csv")

                dends = segs[(segs["Type"] == "dend") & (segs["Distance"] >= 50)]
                apics = segs[(segs["Type"] == "apic")]

                apic_start, self.dend_groups = self._build_exc_nodes(dends, "dend", self.n_dend_exc)
        
                _, self.apic_groups = self._build_exc_nodes(apics, "apic", self.n_apic_exc, start=apic_start)

                self._build_exc_edges(self.dend_groups)
                self._build_exc_edges(self.apic_groups)

        #Sets the number of synapses for each input cell.
        def _connector_func(sources, target, cells):
                """Used to set the number of synapses from each excitatory input
                cell in a functional group. Use with "all_to_one" iterator.

                Parameters
                ----------
                sources : list
                        presynaptic nodes (represented as dicts)
                target : dict
                        postsynaptic node
                cells : list
                        list of Cells in the FunctionalGroup

                Returns
                -------
                list
                        list of synapses for each pairing
                """
                return [cell.n_syns for cell in cells]

        #Sets the location of synapses based on the given cell list.
        def _set_location(source, target, cells, start_id):
                """Sets the location of the given synapse.

                Parameters
                ----------
                source : dict
                    source node information
                target : dict
                    target node information
                cells : list
                    Cells in the functional group
                start_id : int
                    start_id for the functional groups the cells come from

                Returns
                -------
                int
                    BMTK section id
                float
                    distance along the section
                """     
                #Gets the proper index within the cell list.           
                index = source.node_id - start_id

                seg = cells[index].get_seg()
                return seg.bmtk_id, seg.x

        #Creates the functional groups and adds the virtual cells to the
        #BMTK NetworkBuilder.
        def _build_exc_nodes(self, segs, base_name, n_cells, start=0):
                """Creates the functional groups and adds the virtual cells to the
                BMTK NetworkBuilder

                Parameters
                ----------
                segs : pandas.DataFrame
                    all the segments available for the functional groups
                base_name : str
                    the string that is appended to to make the group names.
                    groups get 0 - n_groups appended to their names.
                n_cells : int
                    total number of input cells that should be added.
                start : int, optional
                    starting id to be associated with the functional groups, by default 0
                    this is used later to associate cells in functional groups with the correct
                    locations and synapses.

                Returns
                -------
                int
                    what the start parameter should be for the next call to _build_exc_nodes
                list
                    list of functional groups that were created
                """                
                start_id = start

                n_groups = n_cells // self.params["groups"]["cells_per_group"]
                n_extra = n_cells % self.params["groups"]["cells_per_group"] #number of extra cells that don't evenly fit into groups

                group_list = []

                for i in range(n_groups):
                        name = base_name + str(i)

                        #Spreads out the extra cells.
                        N = self.params["groups"]["cells_per_group"]
                        if i < n_extra:
                                N += 1

                        self.exc_stim.add_nodes(N=N,
                                pop_name=name,
                                potential="exc",
                                model_type='virtual')

                        new_group = FunctionalGroup(segs, segs.sample().iloc[0],
                                N, self.clust_per_group, name, start_id, 
                                partial(make_seg_sphere, radius = self.params["groups"]["group_radius"]), 
                                partial(make_seg_sphere, radius = self.params["groups"]["cluster_radius"]))
                        group_list.append(new_group)
                        start_id += N

                return start_id, group_list

        def _build_exc_edges(self, group_list):
                """Creates the connections between each cell in the list of groups
                and the biophysical cell.

                Parameters
                ----------
                group_list : list
                    list of functional groups
                """                
                for i in range(len(group_list)):
                        group = group_list[i]

                        #Creates the edges from each excitatory input cells in the group.
                        conn = self.net.add_edges(source=self.exc_stim.nodes(pop_name=group.name), target=self.net.nodes(),
                                iterator="all_to_one",
                                connection_rule=SimulationBuilder._connector_func,
                                connection_params={'cells': group.cells},
                                syn_weight=1,
                                delay=0.1,
                                dynamics_params='PN2PN.json',
                                model_template=self.syn['PN2PN.json']['level_of_detail'],)

                        #Sets the postsynaptic locations of the connections.
                        conn.add_properties(['sec_id',"sec_x"], 
                                rule=SimulationBuilder._set_location,
                                rule_params={'cells': group.cells, 'start_id': group.start_id},
                                dtypes=[np.int, np.float])


        def _save_nets(self):
                """builds and saves the BMTK NetworkBuilders
                """                
                # Build and save our networks
                self.net.build()
                self.net.save_nodes(output_dir='network')
                self.net.save_edges(output_dir='network')

                self.exc_stim.build()
                self.exc_stim.save_nodes(output_dir='network')

                self.prox_inh_stim.build()
                self.prox_inh_stim.save_nodes(output_dir='network')

                self.dist_inh_stim.build()
                self.dist_inh_stim.save_nodes(output_dir='network')


        def _build_inh(self):
                """Creates inhibitory input nodes and their connections onto the biophysical cell
                """     

                np.random.seed(1000)           
                #####################Perisomatic Inhibition##############################
                self.prox_inh_stim = NetworkBuilder('prox_inh_stim')

                #Nodes that connect to soma.
                self.prox_inh_stim.add_nodes(N=self.n_soma_inh,
                                pop_name='on_soma',
                                potential='exc',
                                model_type='virtual')
                
                #Nodes that connect to proximal dendrites.
                self.prox_inh_stim.add_nodes(N=self.n_prox_dend_inh,
                                pop_name='on_dend',
                                potential='exc',
                                model_type='virtual')

                div_params = self.params["divergence"]["peri_inh"]

                #On soma.
                self.net.add_edges(source=self.prox_inh_stim.nodes(pop_name='on_soma'),
                                target=self.net.nodes(),
                                connection_rule=SimulationBuilder._norm_connect,
                                connection_params={"m":div_params["m"], "s":div_params["s"], "low":div_params["min"], "high":div_params["max"]},
                                syn_weight=1,
                                delay=0.1,
                                dynamics_params='INT2PN.json',
                                model_template=self.syn['INT2PN.json']['level_of_detail'],
                                distance_range=[-2000, 2000.0],
                                target_sections=['somatic'])

                #On dendrites within 50 um
                self.net.add_edges(source=self.prox_inh_stim.nodes(pop_name='on_dend'), 
                                target=self.net.nodes(),
                                connection_rule=SimulationBuilder._norm_connect,
                                connection_params={"m":div_params["m"], "s":div_params["s"], "low":div_params["min"], "high":div_params["max"]},
                                syn_weight=1,
                                delay=0.1,
                                dynamics_params='INT2PN.json',
                                model_template=self.syn['INT2PN.json']['level_of_detail'],
                                distance_range=[0, 50.0],
                                target_sections=['dend'])
                #######################################################################################

                #############################Dendritic Inhibition######################################
                self.dist_inh_stim = NetworkBuilder('dist_inh_stim')

                self.dist_inh_stim.add_nodes(N=self.n_dend_inh,
                                pop_name='dend',
                                potential='exc',
                                model_type='virtual')

                self.dist_inh_stim.add_nodes(N=self.n_apic_inh,
                                pop_name='apic',
                                potential='exc',
                                model_type='virtual')

                div_params = self.params["divergence"]["basal_inh"]

                #Basal edges.
                self.net.add_edges(source=self.dist_inh_stim.nodes(pop_name="dend"),
                                target=self.net.nodes(),
                                connection_rule=SimulationBuilder._norm_connect,
                                connection_params={"m":div_params["m"], "s":div_params["s"], "low":div_params["min"], "high":div_params["max"]},
                                syn_weight=1,
                                delay=0.1,
                                dynamics_params='INT2PN.json',
                                model_template=self.syn['INT2PN.json']['level_of_detail'],
                                distance_range=[50, 2000.0],
                                target_sections=['dend'])
                
                div_params = self.params["divergence"]["apic_inh"]

                #Apic edges.
                self.net.add_edges(source=self.dist_inh_stim.nodes(pop_name="apic"), 
                                target=self.net.nodes(),
                                connection_rule=SimulationBuilder._norm_connect,
                                connection_params={"m":div_params["m"], "s":div_params["s"], "low":div_params["min"], "high":div_params["max"]},
                                syn_weight=1,
                                delay=0.1,
                                dynamics_params='INT2PN.json',
                                model_template=self.syn['INT2PN.json']['level_of_detail'],
                                distance_range=[50, 2000.0],
                                target_sections=['apic'])

        def _norm_connect(source, target, m, s, low, high):
                """Returns a random number of synapses based on
                the given distribution.

                Parameters
                ----------
                source : dict
                    source node
                target : dict
                    target node
                m : float
                    mean number of connections
                s : float
                    standard deviation of number of connections
                low : int
                    minimum number of connections
                high : int
                    maximum number of connections

                Returns
                -------
                int
                    number of connections
                """                
                return int(min(max(np.random.normal(m, s), low), high))

            

        def _make_rasters(self):
                """Generates excitatory and inhibitory input rasters
                """    
                np.random.seed(1000)            
                self._gen_exc_spikes('exc_stim_spikes.h5')

                inh_frs = self.params["inh_frs"]

                #Makes perisomatic inhibitory raster.
                self._gen_inh_spikes(self.n_soma_inh + self.n_prox_dend_inh, 
                        inh_frs["proximal"]["m"], inh_frs["proximal"]["s"], 
                        "prox_inh_stim", 'prox_inh_stim_spikes.h5')
                #Makes dendritic inhibitory raster.
                self._gen_inh_spikes(self.n_apic_inh + self.n_dend_inh, 
                        inh_frs["distal"]["m"], inh_frs["distal"]["s"], 
                        "dist_inh_stim", 'dist_inh_stim_spikes.h5')


        #Generates the spike raster for a given group.
        #The group has the same noise.
        def _gen_group_spikes(writer, group, seconds, start_time, dist):
                """Generates and writes to a h5 file the given functional group's spike trains

                Parameters
                ----------
                writer : SonataWriter
                    how the spike trains are saved
                group : FunctionalGroup
                    the functional group that the spike trains are being made for
                seconds : float
                    length of the spike trains in seconds
                start_time : float
                    what time (ms) the spike trains should start at
                dist : func
                    function for random distribution used for an individual cell's firing rate
                """                
                z = make_noise(num_samples=(int(seconds*1000))-1,num_traces=1)#generates the noise trace common to each cell in the functional group.
                make_save_spikes(writer, True, dist, numUnits=group.n_cells,
                        rateProf=z[0,:],start_id=group.start_id,
                        start_time=start_time)

        #Creates the excitatory input raster from the functional groups.
        def _gen_exc_spikes(self, fname):
                """Generates the excitatory input raster for all of the functional groups

                Parameters
                ----------
                fname : str
                    name of the file to save the rasters in (.h5)
                """    
                #distribution used for generating excitatory firing rates.    
                levy_dist = partial(st.levy_stable.rvs, alpha=1.37, beta=-1.00, loc=0.92, scale=0.44, size=1)

                length = self.params["time"]["stop"] - self.params["time"]["start"]
                buffer = self.params["time"]["start"]

                writer = SonataWriter(fname, ["spikes", "exc_stim"], ["timestamps", "node_ids"], [np.float, np.int])

                for group in (self.dend_groups + self.apic_groups):
                        SimulationBuilder._gen_group_spikes(writer, group, length, buffer*1000, levy_dist)
                

        #Blocks off the bottom of a normal distribution.
        def _norm_rvs(mean, std):
                """Generates a random float from a normal distribution with a near zero minimum

                Parameters
                ----------
                mean : float
                    mean of the distribution
                std : float
                    standard deviation of the distribution

                Returns
                -------
                float
                    random float
                """                
                return max(st.norm.rvs(loc=mean, scale=std, size=1), 0.001)

        # #Makes a spike raster with each cell having its own noise trace.
        # def gen_inh_spikes(n_cells, mean_fr, std_fr, key, file, times):
        #         # node_ids = []
        #         # timestamps = []

        #         length = times[1] - times[0]
        #         buffer = times[0]

        #         writer = SonataWriter(file, ["spikes", key], ["timestamps", "node_ids"], [np.float, np.int])

        #         z = make_noise(num_samples=(int(length*1000))-1,num_traces=1)
        #         make_save_spikes(writer, False, partial(positive_normal, mean=mean_fr, std=std_fr), numUnits=n_cells,rateProf=z[0,:],start_time=buffer*1000)

        #Creates a spike raster with each cell having the same noise coming from the a shifted average of excitation.
        def _gen_inh_spikes(self, n_cells, mean_fr, std_fr, key, fname):
                """Generates a spike raster with each train having the noise trace from
                averaging excitation. Distributes firing rates normally.

                Parameters
                ----------
                n_cells : int
                    number of spike trains
                mean_fr : float
                    mean firing rate
                std_fr : float
                    standard deviation of the firing rate
                key : str
                    name of the second group in the h5 file
                fname : str
                    name of file to save the raster to
                """                
                # node_ids = []
                # timestamps = []

                f = h5py.File("exc_stim_spikes.h5", "r")
                ts = f['spikes']["exc_stim"]['timestamps']
                nid = f['spikes']["exc_stim"]['node_ids']

                #Creates a noise trace based on the excitatory spike raster.
                z = shift_exc_noise(ts, nid, self.params["time"]["stop"], time_shift=self.params["inh_shift"])
                
                writer = SonataWriter(fname, ["spikes", key], ["timestamps", "node_ids"], [np.float, np.int])

                make_save_spikes(writer, False, partial(SimulationBuilder._norm_rvs, mean=mean_fr, std=std_fr), numUnits=n_cells,rateProf=z)


if __name__ == "__main__":
        np.random.seed(2129)
        builder = SimulationBuilder("NetParams.json")
        builder.build()
        builder.save_groups()