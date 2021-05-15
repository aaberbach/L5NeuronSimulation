"""Clustering Tools

This script contains classes and functions for implementing synaptic clustering
in the build step of the simulation.
"""

import numpy as np
import pandas as pd
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from functools import partial 

np.random.seed(42)

def calc_dist(p1, p2):
    """Returns the distance between two points.

    Parameters
    ----------
    p1 : list
        [x, y, z] representation of the first coordinate
    p2 : list
        [x, y, z] representation of the second coordinate
    
    Returns
    -------
    float
        The distance between the two given points
    """
    return np.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)

def make_seg_sphere(segments, center, radius=50):
    """Returns segments within a sphere of the 
    given radius centered at the given center.

    Parameters
    ----------
    segments : pandas.DataFrame
        DataFrame with each row representing a segment that can be in the sphere
    center : pandas.Series
        Series representing the segment that is the center of the created sphere
    radius : int, optional
        The radius of sphere to be made in um (default is 50)

    Returns
    -------
    pandas.DataFrame
        DataFrame contiaining the subset of segments that is within the sphere
    """
    o = [center["Coord X"], center["Coord Y"], center["Coord Z"]]
    return segments[calc_dist(o, [segments["Coord X"], segments["Coord Y"], segments["Coord Z"]]) <= radius]


class FunctionalGroup:
    """ A class used to represent a functional group of input cells.
    Each functional group has some number of clusters that it places its cells' synapses into.

    Attributes
    ----------
    center : pd.Series
        the segment that the functional group is centered around
    n_cells : int
        number of cells
    syn_per_cell : list
        list of length two that contains the lower and upper boundary for the
        number of synapses from each cell. [low, high]
    n_clusters : int
        number of clusters
    name : str
        the name associated with the group. Used in the build step of the simulation for proper assignment.
    start_id : int
        the node_id that is associated with the first cell in the cells list
    clustering_func : func
        how segments are chosen for the clusters (formatted like make_seg_sphere above)
    group_segs : pd.DataFrame
        DataFrame of segments that fall within the group
    cells : list
        list of Cell objects that belong to the group
    clusters : list
        list of Cluster objects that belong to the group

    Methods
    -------
    make_clusters()
        Creates the group's clusters
    create_cells()
        Creates the group's cells
    """
    def __init__(self, seg_list, center, num_cells, num_clusters, name, start_id, grouping_func, clustering_func, syn_per_cell=[2, 8]):
        """Initializes the object and creates the clusters and cells.
        
        See above for most parameters

        Parameters
        ----------
        grouping_func : func
            function used to find segments that can be used (formatted like mage_seg_sphere above)
        """
        self.center = center
        self.n_cells = num_cells
        self.syn_per_cell = syn_per_cell
        self.n_clusters = num_clusters
        self.name = name
        self.start_id = start_id
        self.clustering_func = clustering_func
        self.group_segs = grouping_func(seg_list, center)

        self.cells = []
        self.clusters = []

        self.make_clusters()
        self.create_cells()

    #Initializes the clusters of the group.
    #Centers are randomly selected (uniquely) from all possible segs.
    #Might want some purposeful spread?
    def make_clusters(self):
        """Initializes the clusters of the group.
        Cluster centers are randomly selected (uniquely) from group_segs.
        """
        if self.n_clusters > len(self.group_segs):
            raise Exception("The number of segment options should never be bellow n_clusters.")

        center_ids = np.random.choice(range(0, len(self.group_segs)), self.n_clusters, replace=False)
        for id in center_ids:
            self.clusters.append(Cluster(self.group_segs.iloc[id], self.clustering_func, self.group_segs))

    #Creates the group's Cells with randomly assigned number of synapses.
    #Distributes the cell's synapses onto clusters.
    def create_cells(self):
        """Creates the group's cells.
        Gives each cell and random number of synapses within syn_per_cell.
        Distributes each cell's synapses onto the clusters.
        """
        for i in range(self.n_cells):
            cell = Cell(np.random.randint(low=self.syn_per_cell[0], high=self.syn_per_cell[1] + 1))
            cluster_ids = np.random.choice(range(0, self.n_clusters), cell.n_syns, replace=False)
            for id in cluster_ids:
                cell.add_syn(self.clusters[id].random_seg())

            self.cells.append(cell)

            
class Cluster:
    """ A class used to represent a cluster of synapses within a functional group.

    Attributes
    ----------
    center : pd.Series
        the segment that the cluster is centered around
    n_syns : int
        number of synapses attached to the cluster
    cluster_segs : pd.DataFrame
        contains the segments that are within the cluster

    Methods
    -------
    random_seg() -> pd.Series
        Returns a random segment in cluster_segs
    """
    def __init__(self, center, clustering_func, group_segs):
        """
        Parameters
        ----------
        center : pandas.Series
            series representing the segment that is the center of the cluster
        clustering_func : func
            used to select cluster_segs
        group_segs : pandas.DataFrame
            segments that are within the functional group
        """
        self.center = center
        self.n_syns = 0
        self.cluster_segs = clustering_func(group_segs, center)

    #Returns a random segment from cluster segs.
    #Could add count table that spreads synapses.
    def random_seg(self):
        """Returns a random segment in cluster_segs and
        increases the count of synapses in the cluster
        """
        self.n_syns += 1
        return self.cluster_segs.iloc[np.random.choice(len(self.cluster_segs))]
    


#Simple class that contains the number of synapses and the synapse locations for an input cell.
class Cell:
    """ A class used to an input cell within a functional group.

    Attributes
    ----------
    n_syns : int
        number of synapses from the cell
    syn_segs : list
        list of Segments where the cell's synapses should be placed
    num_set : int
        number of synapses that the build step has already created. Functions as a way to iterate through synapses.

    Methods
    -------
    add_syn(seg : pandas.Series)
        Adds the given segment to syn_segs
    get_seg()
        Returns the next synapse to set
        
    """
    def __init__(self, n_syns):
        """Initializes the cell. See above for parameters"""
        self.n_syns = n_syns
        self.syn_segs = [] #List of Segments where the cell's synapses should be placed.
        self.num_set = 0#Number of synapses already created in bmtk.

    def add_syn(self, seg):
        """Adds the given segment to syn_segs

        Parameters
        ----------
        seg : pandas.Series
            represents the segment to added to syn_segs
        """
        if len(self.syn_segs) >= self.n_syns:
            raise Exception("Error: too many synapses added to cell.")

        self.syn_segs.append(Segment(seg["BMTK ID"], seg["X"]))

    #Returns the next syn location.
    def get_seg(self):
        """Returns the next Segment where a synapses should
        be created with this cell as the presynaptic cell."""
        if self.num_set >= self.n_syns:
            raise Exception("Too many synapses created for cell.")

        result = self.syn_segs[self.num_set]
        self.num_set += 1
        return result


#Simple class that contains the BMTK morphology id and x of a segment.
class Segment:    
    """Simple class that represents a single segment in BMTK.

    Attributes
    ----------
    bmtk_id : int
        the id that BMTK associates with the section (not the section id)
    x : float
        distance along the section
    """

    def __init__(self, bmtk_id, x):
        self.bmtk_id = bmtk_id
        self.x = x


def plot_group(segs, group, ax, psoma=[0, 0, 0]):  
    """Plots the given group's center, cells, and clusters on the given plt axes.

    Parameters
    ----------
    segs : pd.DataFrame
        all segments that can possibly be plotted
    group : FunctionalGroup
        the group to be plotted
    ax : Axes3d
        what to plot on
    psoma : list, optional
        used to shift the coordinates of the segments to match the morphology plot.
        formatted as [x, y, z] (default is [0, 0, 0])
    """

    #Plots a black star to represent the center of the functional group.
    ax.scatter3D([group.center["Coord X"] + psoma[0]], [group.center["Coord Y"] + psoma[1]], [group.center["Coord Z"] + psoma[2]], color="Black", marker="*", s=50)
    zdata = []
    xdata = []
    ydata = []
    vdata = []

    #Marks each synapse.
    cell_colors = np.linspace(0, 1, group.n_cells)
    for i, cell in enumerate(group.cells):
        for seg in cell.syn_segs:
            s = segs[(segs["X"] == seg.x) & (segs["BMTK ID"] == seg.bmtk_id)].iloc[0]
            xdata.append(s["Coord X"])
            ydata.append(s["Coord Y"])
            zdata.append(s["Coord Z"])
            vdata.append(cell_colors[i])

    ax.scatter3D(np.array(xdata) + psoma[0], np.array(ydata) + psoma[1], np.array(zdata) + psoma[2], c=vdata, cmap="hsv", alpha = 1, s=10)

    zdata = []
    xdata = []
    ydata = []

    #Marks the center of each synapse.
    for cluster in group.clusters:
        #s = df[(df["X"] == seg.x) & (df["BMTK ID"] == seg.bmtk_id)].iloc[0]
        xdata.append(cluster.center["Coord X"])
        ydata.append(cluster.center["Coord Y"])
        zdata.append(cluster.center["Coord Z"])

    ax.scatter3D(np.array(xdata) + psoma[0], np.array(ydata) + psoma[1], np.array(zdata) + psoma[2], color = "Green", s=30)