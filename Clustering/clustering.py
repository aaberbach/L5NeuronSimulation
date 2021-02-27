import numpy as np
import pandas as pd
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from functools import partial 

np.random.seed(42)

def calc_dist(p1, p2):
    return np.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)

def make_seg_sphere(segments, center, radius=50):
    #o = [center["Coord X"].iloc[0], center["Coord Y"].iloc[0], center["Coord Z"].iloc[0]]
    o = [center["Coord X"], center["Coord Y"], center["Coord Z"]]
    return segments[calc_dist(o, [segments["Coord X"], segments["Coord Y"], segments["Coord Z"]]) <= radius]

# class SegmentList:
#     def __init__(self, file):
#         segs = pd.read_csv(file)

#     def filter(self, part, func, center):
#         return func(self.segs[self.segs["Type"] == part], center)


class FunctionalGroup:
    def __init__(self, seg_list, center, num_cells, num_clusters, name, start_id, grouping_func, clustering_func):
        self.center = center
        self.n_cells = num_cells
        self.n_clusters = num_clusters
        self.name =  name
        self.start_id = start_id
        #self.grouping_func = grouping_func
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
        if self.n_clusters > len(self.group_segs):
            raise Exception("The number of segment options should never be bellow n_clusters.")

        center_ids = np.random.choice(range(0, len(self.group_segs)), self.n_clusters, replace=False)
        for id in center_ids:
            self.clusters.append(Cluster(self.group_segs.iloc[id], self.clustering_func, self.group_segs))

    #Creates the group's Cells with randomly assigned number of synapses.
    #Distributes the cell's synapses onto clusters.
    def create_cells(self):
        for i in range(self.n_cells):
            cell = Cell(np.random.randint(low=2, high=8 + 1))
            cluster_ids = np.random.choice(range(0, self.n_clusters), cell.n_syns, replace=False)
            for id in cluster_ids:
                #import pdb; pdb.set_trace()
                cell.add_syn(self.clusters[id].random_seg())

            self.cells.append(cell)

            

class Cluster:
    def __init__(self, center, clustering_func, group_segs):
        self.center = center
        #import pdb; pdb.set_trace()
        self.cluster_segs = clustering_func(group_segs, center)

    #Returns a random segment from cluster segs.
    #Could add count table that spreads synapses.
    def random_seg(self):
        #import pdb; pdb.set_trace()
        return self.cluster_segs.iloc[np.random.choice(len(self.cluster_segs))]
    


#Simple class that contains the number of synapses and the synapse locations for an input cell.
class Cell:
    def __init__(self, n_syns):
        self.n_syns = n_syns
        self.syn_segs = [] #List of Segments where the cell's synapses should be placed.
        self.num_set = 0#Number of synapses already created in bmtk.

    def add_syn(self, seg):
        if len(self.syn_segs) >= self.n_syns:
            raise Exception("Error: too many synapses added to cell.")

        self.syn_segs.append(Segment(seg["BMTK ID"], seg["X"]))

    #Returns the next syn location.
    def get_seg(self):
        if self.num_set >= self.n_syns:
            raise Exception("Too many synapses created for cell.")

        result = self.syn_segs[self.num_set]
        self.num_set += 1
        return result


#Simple class that contains the BMTK morphology id and x of a segment.
class Segment:
    def __init__(self, bmtk_id, x):
        self.bmtk_id = bmtk_id
        self.x = x

def plot_group(segs, group, ax, psoma=[0, 0, 0]):
    #import pdb; pdb.set_trace()
    ax.scatter3D([group.center["Coord X"] + psoma[0]], [group.center["Coord Y"] + psoma[1]], [group.center["Coord Z"] + psoma[2]], color="Black", marker="*", s=50)
    zdata = []
    xdata = []
    ydata = []
    vdata = []

    cell_colors = np.linspace(0, 1, group.n_cells)
    for i, cell in enumerate(group.cells):
        for seg in cell.syn_segs:
            #import pdb; pdb.set_trace()
            s = segs[(segs["X"] == seg.x) & (segs["BMTK ID"] == seg.bmtk_id)].iloc[0]
            xdata.append(s["Coord X"])
            ydata.append(s["Coord Y"])
            zdata.append(s["Coord Z"])
            vdata.append(cell_colors[i])

    ax.scatter3D(np.array(xdata) + psoma[0], np.array(ydata) + psoma[1], np.array(zdata) + psoma[2], c=vdata, cmap="hsv", alpha = 1, s=10)

    zdata = []
    xdata = []
    ydata = []

    for cluster in group.clusters:
        #import pdb; pdb.set_trace()
        #s = df[(df["X"] == seg.x) & (df["BMTK ID"] == seg.bmtk_id)].iloc[0]
        xdata.append(cluster.center["Coord X"])
        ydata.append(cluster.center["Coord Y"])
        zdata.append(cluster.center["Coord Z"])

    ax.scatter3D(np.array(xdata) + psoma[0], np.array(ydata) + psoma[1], np.array(zdata) + psoma[2], color = "Green", s=30)

# df = pd.read_csv("Segments_small.csv")
# dends = df[df["Type"] == "dend"]
# apics = df[df["Type"] == "apic"]
# #import pdb; pdb.set_trace()
# group = FunctionalGroup(dends, dends.sample().iloc[0], 60, 20, "first", partial(make_seg_sphere, radius = 100), partial(make_seg_sphere, radius = 5))


# ax = plt.axes(projection='3d')

# plot_group(dends, group, ax)
# plot_group(apics, FunctionalGroup(apics, apics.sample().iloc[0], 60, 20, "first", partial(make_seg_sphere, radius = 100), partial(make_seg_sphere, radius = 5)), ax)

# # ax.set_xlim([-200, 300])
# # ax.set_ylim([-200, 1200])
# # ax.set_zlim([-200, 50])
# plt.show()