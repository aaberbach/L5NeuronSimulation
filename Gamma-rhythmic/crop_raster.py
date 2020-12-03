import h5py
import numpy as np

def crop_raster(raster_file, new_file, time, num_cells):
    f = h5py.File(raster_file, 'r')
    raster = f['spikes']#['inh_stim'][]
    key = list(raster.keys())[0]
    raster = f['spikes'][key]
    gids = np.array(raster['node_ids'])
    timestamps = np.array(raster['timestamps'])
    f.close()

    if (num_cells >= np.max(gids)):
        raise Exception("The inputed raster_file has less than {num_cells} cells.")
    if (time > np.max(timestamps)):
        raise Exception("The inputed raster_file does not go for {time} ms.")

    ids = np.where((gids < num_cells) & (timestamps < time))[0]
    gids = gids[ids]
    timestamps = timestamps[ids]
    #import pdb; pdb.set_trace()
    f = h5py.File(new_file, 'w')
    group = f.create_group('spikes')
    group = group.create_group(key)
    group.create_dataset("node_ids", data = gids)
    group.create_dataset("timestamps", data = timestamps)
    f.close()

