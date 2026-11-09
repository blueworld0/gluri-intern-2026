import importlib
import json
import os

results = {}
for name in ['numpy', 'h5py', 'scipy', 'shapely', 'pyproj', 'matplotlib', 'csv', 'shapefile']:
    try:
        module = importlib.import_module(name)
        results[name] = {'ok': True, 'version': getattr(module, '__version__', 'stdlib')}
    except Exception as exc:
        results[name] = {'ok': False, 'error': str(exc)}
print(json.dumps({'imports': results, 'earthdata_env': {
    key: bool(os.environ.get(key))
    for key in ['EARTHDATA_USERNAME', 'EARTHDATA_PASSWORD']
}}, ensure_ascii=False))
