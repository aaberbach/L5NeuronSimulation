from bmtk.simulator.bionet.modules.record_cellvars import MembraneReport


def modify_bmtk():
    MembraneReport.step = step
    return


def step(self, sim, tstep):
        # save all necessary cells/variables at the current time-step into memory
        for gid in self._local_gids:
            pop_id = self._gid_map.get_pool_id(gid)
            cell = sim.net.get_cell_gid(gid)
            for var_name in self._variables:
                if '.' in var_name:
                    seg_vals = [getattr(getattr(seg, var_name.split('.')[0], None), var_name.split('.')[1], None) for seg in cell.get_segments()]
                else:
                    seg_vals = [getattr(seg, var_name, None) for seg in cell.get_segments()]
                self._var_recorder.record_cell(pop_id.node_id, population=pop_id.population, vals=seg_vals, tstep=tstep)

            for var_name, fnc in self._transforms.items():
                seg_vals = [fnc(getattr(seg, var_name)) for seg in cell.get_segments()]
                # self._var_recorder.record_cell(gid, var_name, seg_vals, tstep)
                self._var_recorder.record_cell(pop_id.node_id, population=pop_id.population, val=seg_vals, tstep=tstep)

        self._block_step += 1
