from .base import VolumeProvider

def get_provider(config):
    if config['volume']['provider']=='stat_table':
        from .stat_table import StatTableProvider
        return StatTableProvider(config)
    if config['volume']['provider']=='raster':
        from .raster import RasterProvider
        return RasterProvider(config)
    raise ValueError('Unknown volume provider')
