"""Volume suppliers preserve row count, index and order; missing is never zero."""
from abc import ABC, abstractmethod

class VolumeProvider(ABC):
    @abstractmethod
    def provide(self, polygons):
        """Return DataFrame m3_per_ha, volume_source, volume_flag with input index."""
    @abstractmethod
    def provenance(self):
        """Return JSON-serializable input hashes, edition and applied rules."""

def file_record(path, **metadata):
    import hashlib
    from pathlib import Path
    p=Path(path); h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return {'path':str(p.resolve()),'sha256':h.hexdigest(),**metadata}
