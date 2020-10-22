from allensdk.api.queries.biophysical_api import BiophysicalApi

bp = BiophysicalApi()
bp.cache_stimulus = False # change to False to not download the large stimulus NWB file
neuronal_model_id = 497232429    # get this from the web site as above
bp.cache_data(neuronal_model_id, working_directory='neuronal_model')