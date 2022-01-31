"""
The script used to build the network for the simulation run in run_network.py
"""
import os,sys,inspect
import shutil
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

try:
    np.random.seed(int(sys.argv[-1]))
except:
    np.random.seed(123)

def delete_file(filename):
        try:
                os.remove(filename)
                print(filename + " sucessfully deleted.")
        except:
                print(filename + " not found for deletion.")

def delete_folder(dir_name):
        try:
                shutil.rmtree(dir_name)
                print(dir_name + " sucessfully deleted.")
        except:
                print(dir_name + " not found for deletion.")


def clean_folder():
        delete_file("./node_sets.json")
        delete_file("./run_bionet.py")
        delete_file("./simulation_config.json")
        delete_file("./prox_inh_stim_spikes.h5")
        delete_file("./exc_stim_spikes.h5")
        delete_file("./dist_inh_stim_spikes.h5")
        delete_file("./circuit_config.json")
        delete_folder("./network")
        delete_folder("./output")
        delete_folder("./__pycache__")

class SimulationBuilder:
        """Class used to build our BMTK simulation.

        Attributes
        ----------
        params : dict
                contains parameters for the network
        seed : int
                base random seed for the simulation
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

        _set_prefixed_directory(base_dir_name : str)
                sets up the correct biophy_components structure based on the cell prefix in params for the given directory base

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
        _modify_jsons()
                modifies the various json files however is needed after they are built
        _modify_sim_config()
                modifies the simulation_config.json however is needed
        _update_cellvar_record_locs(sim_config : dict)
                modifies the location of cellvar recordings in the given JSON simulation_config

        Static Methods
        --------------
        _get_directory_prefix(directory : str)
                reads the prefix.txt fil in directory and returns the contents
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

        def __init__(self, params_file, seed=123):
                """Initializes the simulation builder, 
                setting up attributes but not actually building the BMTK network.

                Parameters
                ----------
                params_file : str
                    path to the JSON file with network parameters
                seed : int
                    base random seed for the simulation
                """     
                #Loads the JSON file with information about the network.           
                with open(params_file) as f:
                        self.params = json.load(f)

                self.seed = seed

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
                if self.params["file_current_clamp"]["input_file"]=="None": 
                    self.file_current_clamp = None
                else:
                    self.file_current_clamp = self.params["file_current_clamp"]

        def build(self):
                """Builds the nodes and edges for the network.
                """                
                np.random.seed(self.seed)

                self._set_prefixed_directory("mechanisms")
                self._set_prefixed_directory("templates")

                self.net = NetworkBuilder("biophysical")

                self.net.add_nodes(N=1, pop_name='Pyrc',
                        potental='exc',
                        model_type='biophysical',
                        dynamics_params= self.params["cell"]["dynamic_params"],
                        model_template= self.params["cell"]["model_template"],
                        model_processing = self.params["cell"]["model_processing"],
                        morphology = self.params["cell"]["morphology"])

                self._build_exc()
                self._build_inh()
                self._save_nets()

                self._make_rasters()
                
                #Final build step.
                build_env_bionet(base_dir='./',
                        network_dir='./network',
                        dt = self.params["dt"], 
                        tstop=self.params["time"]["stop"] * 1000.0,
                        report_vars=self.params["record_cellvars"]["vars"],
                        dL = self.params["dL"],#target length (um) of segments
                        spikes_threshold=-10,
                        file_current_clamp=self.file_current_clamp,
                        spikes_inputs=[('exc_stim', 'exc_stim_spikes.h5'), ('prox_inh_stim', 'prox_inh_stim_spikes.h5'), ('dist_inh_stim', 'dist_inh_stim_spikes.h5')],
                        components_dir='../biophys_components',
                        compile_mechanisms=True)

                self._modify_jsons()

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

        def _set_prefixed_directory(self, base_dir_name):
                """Fixes the biophy_components directory. There should be only one directory
                named <base_dir_name> and it should be the one with the prefix.txt file in it
                that has the same prefix as params.

                Parameters
                ----------
                base_dir_name : str
                        base name of the set of directories to be fixed
                """
                #import pdb; pdb.set_trace()
                components_path = "../biophys_components/"
                biophys_subdirs = [ f.name for f in os.scandir(components_path) if f.is_dir() ]
                
                for dir_name in biophys_subdirs:
                        if base_dir_name == dir_name:
                                prefix = SimulationBuilder._get_directory_prefix(components_path + dir_name)
                                if prefix == self.params["cell"]["prefix"]:
                                        return
                                else:
                                        os.rename(components_path + base_dir_name, components_path + prefix + base_dir_name)
                                        


                for dir_name in biophys_subdirs:
                        if base_dir_name in dir_name and self.params["cell"]["prefix"] in dir_name:
                                os.rename(components_path + dir_name, components_path + base_dir_name)

        def _get_directory_prefix(directory):
                """Returns the contents of the prefix.txt file in the given directory.

                Parameters
                ----------
                directory : str
                        directory to look in

                Returns
                -------
                str
                        contents of prefix.txt
                """  
                with open(directory + "/prefix.txt", 'r') as f:
                        return f.read()

        def _build_exc(self):
                """Builds the excitatory input cells and their synapses.
                """          

                # External excitatory inputs
                self.exc_stim = NetworkBuilder('exc_stim')

                #DataFrame of all segments on the cell.
                segs = pd.read_csv(self.params["cell"]["segments_file"])

                dends = segs[(segs["Type"] == "dend") & (segs["Distance"] >= 50)]
                apics = segs[(segs["Type"] == "apic")]

                np.random.seed(self.seed + 1)
                apic_start, self.dend_groups = self._build_exc_nodes(dends, "dend", self.n_dend_exc)
                
                np.random.seed(self.seed + 2)
                _, self.apic_groups = self._build_exc_nodes(apics, "apic", self.n_apic_exc, start=apic_start)

                np.random.seed(self.seed + 3)
                self._build_exc_edges(self.dend_groups)

                np.random.seed(self.seed + 4)
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
                np.random.seed(self.seed + 12)
                self.net.build()
                self.net.save_nodes(output_dir='network')
                np.random.seed(self.seed + 16)
                self.net.save_edges(output_dir='network')

                np.random.seed(self.seed + 13)
                self.exc_stim.build()
                self.exc_stim.save_nodes(output_dir='network')

                np.random.seed(self.seed + 14)
                self.prox_inh_stim.build()
                self.prox_inh_stim.save_nodes(output_dir='network')

                np.random.seed(self.seed + 15)
                self.dist_inh_stim.build()
                self.dist_inh_stim.save_nodes(output_dir='network')


        def _build_inh(self):
                """Creates inhibitory input nodes and their connections onto the biophysical cell
                """     

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
                np.random.seed(self.seed + 5)
                self.net.add_edges(source=self.prox_inh_stim.nodes(pop_name='on_soma'),
                                target=self.net.nodes(),
                                connection_rule=SimulationBuilder._norm_connect,
                                connection_params={"m":div_params["m"], "s":div_params["s"], "low":div_params["min"], "high":div_params["max"]},
                                syn_weight=1,
                                delay=0.1,
                                dynamics_params='PV2PN.json',
                                model_template=self.syn['PV2PN.json']['level_of_detail'],
                                distance_range=[-2000, 2000.0],
                                target_sections=['somatic'])

                #On dendrites within 50 um
                np.random.seed(self.seed + 6)
                self.net.add_edges(source=self.prox_inh_stim.nodes(pop_name='on_dend'), 
                                target=self.net.nodes(),
                                connection_rule=SimulationBuilder._norm_connect,
                                connection_params={"m":div_params["m"], "s":div_params["s"], "low":div_params["min"], "high":div_params["max"]},
                                syn_weight=1,
                                delay=0.1,
                                dynamics_params='PV2PN.json',
                                model_template=self.syn['PV2PN.json']['level_of_detail'],
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
                np.random.seed(self.seed + 7)
                self.net.add_edges(source=self.dist_inh_stim.nodes(pop_name="dend"),
                                target=self.net.nodes(),
                                connection_rule=SimulationBuilder._norm_connect,
                                connection_params={"m":div_params["m"], "s":div_params["s"], "low":div_params["min"], "high":div_params["max"]},
                                syn_weight=1,
                                delay=0.1,
                                dynamics_params='SOM2PN.json',
                                model_template=self.syn['SOM2PN.json']['level_of_detail'],
                                distance_range=[50, 2000.0],
                                target_sections=['dend'])
                
                div_params = self.params["divergence"]["apic_inh"]

                #Apic edges.
                np.random.seed(self.seed + 8)
                self.net.add_edges(source=self.dist_inh_stim.nodes(pop_name="apic"), 
                                target=self.net.nodes(),
                                connection_rule=SimulationBuilder._norm_connect,
                                connection_params={"m":div_params["m"], "s":div_params["s"], "low":div_params["min"], "high":div_params["max"]},
                                syn_weight=1,
                                delay=0.1,
                                dynamics_params='SOM2PN.json',
                                model_template=self.syn['SOM2PN.json']['level_of_detail'],
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
                np.random.seed(self.seed + 9)
                self._gen_exc_spikes('exc_stim_spikes.h5')

                inh_frs = self.params["inh_frs"]

                #Makes perisomatic inhibitory raster.
                np.random.seed(self.seed + 10)
                self._gen_inh_spikes(self.n_soma_inh + self.n_prox_dend_inh, 
                                     inh_frs["proximal"]["m"], 
                                     inh_frs["proximal"]["s"], 
                                     inh_frs["proximal"]["rhythmicity"],
                                     "prox_inh_stim", 
                                     'prox_inh_stim_spikes.h5')
                
                #Makes dendritic inhibitory raster.
                np.random.seed(self.seed + 11)
                self._gen_inh_spikes(self.n_apic_inh + self.n_dend_inh, 
                                     inh_frs["distal"]["m"], 
                                     inh_frs["distal"]["s"],
                                     inh_frs["distal"]["rhythmicity"],
                                     "dist_inh_stim", 
                                     'dist_inh_stim_spikes.h5')


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
                make_save_spikes(writer, True, dist(size=group.n_cells), numUnits=group.n_cells,
                        rateProf=np.tile(z[0,:],(group.n_cells,1)),start_id=group.start_id,
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
        def _gen_inh_spikes(self, n_cells, mean_fr, std_fr, rhythmic_dict, key, fname):
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
                rhythmic_dict : dict
                    dictionary with keys f - frequency, mod - depth of modulation
                key : str
                    name of the second group in the h5 file
                fname : str
                    name of file to save the raster to
                """                
                # node_ids = []
                # timestamps = []
                a, b = (0 - mean_fr) / std_fr, (100 - mean_fr) / std_fr
                d = partial(st.truncnorm.rvs, a=a, b=b, loc=mean_fr, scale=std_fr)
                
                if rhythmic_dict['f'] == "None":
                    f = h5py.File("exc_stim_spikes.h5", "r")
                    ts = f['spikes']["exc_stim"]['timestamps']
                    nid = f['spikes']["exc_stim"]['node_ids']

                    #Creates a noise trace based on the excitatory spike raster.
                    z = shift_exc_noise(ts, nid, self.params["time"]["stop"], time_shift=self.params["inh_shift"])
                    z = np.tile(z,(n_cells,1))
                    
                    writer = SonataWriter(fname, ["spikes", key], ["timestamps", "node_ids"], [np.float, np.int])
                    make_save_spikes(writer, False, d(size=n_cells), numUnits=n_cells,rateProf=z)

                else:
                    # make an array of modulated sin waves
                    # make_save_spikes should be written so that the firing rates are generated
                    #    outside instead of inside the function.
                    frs = d(size=n_cells)
                    
                    t = np.arange(0,self.params["time"]["stop"],0.001)
                    z = np.zeros((n_cells,t.shape[0]))
                    P = 0
                    for i in np.arange(0,n_cells):
                        offset = frs[i]
                        A = offset/((1/rhythmic_dict['mod'])-1)
                        z[i,:] = A*np.sin((2 * np.pi * rhythmic_dict['f'] * t)+P) + offset

                    writer = SonataWriter(fname, ["spikes", key], ["timestamps", "node_ids"], [np.float, np.int])
                    make_save_spikes(writer, False, np.ones((n_cells,1)), numUnits=n_cells,rateProf=z)

        def _modify_jsons(self):
                """modifies the various json files however is needed after they are built"""
                self._modify_sim_config()

        def _modify_sim_config(self):
                """modifies the simulation_config.json however is needed"""
                with open("simulation_config.json", "r") as jsonFile:
                        sim_config = json.load(jsonFile)

                self._update_cellvar_record_locs(sim_config)

                with open("simulation_config.json", "w") as jsonFile:
                        json.dump(sim_config, jsonFile, indent=2)

        def _update_cellvar_record_locs(self, sim_config):
                """modifies the location of cellvar recordings in the given JSON simulation_config
                
                Parameters
                ----------
                sim_config : dict
                    simulation_config to modify
                """
                reports = sim_config["reports"]
                cellvar_reports = [report for report in reports.values() if report["module"] == "membrane_report"]

                for loc, report in zip(self.params["record_cellvars"]["locs"], cellvar_reports):
                        report["sections"] = loc

if __name__ == "__main__":
        try:
                net_params = sys.argv[1]
        except:
                net_params = "L5NetParams.json"

        if ".json" not in net_params:
                net_params = "L5NetParams.json"

        clean_folder()

        builder = SimulationBuilder(net_params)
        builder.build()
        builder.save_groups()
