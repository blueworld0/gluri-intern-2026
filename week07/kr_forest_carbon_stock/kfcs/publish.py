"""Copy small shareable outputs, with provenance paths relative to the experiment."""
import copy
import json
import os
from pathlib import Path
import re
import shutil

def strip_long_prefix(value):
    return value[4:] if value.startswith('\\\\?\\') else value

def relative_provenance(data,experiment_root):
    root=strip_long_prefix(str(Path(experiment_root).resolve()))
    def visit(value):
        if isinstance(value,dict):return {k:visit(v) for k,v in value.items()}
        if isinstance(value,list):return [visit(v) for v in value]
        if isinstance(value,str) and not value.startswith(('https://','http://')):
            s=strip_long_prefix(value)
            if os.path.isabs(s):return os.path.relpath(s,root).replace('\\','/')
        return value
    out=visit(copy.deepcopy(data))
    out['shared_copy']={'path_base':'../..','path_base_description':'experiment root; paths below are relative to that root',
                        'absolute_paths_removed':True,'full_outputs_location':'Drive local mirror; see README.md'}
    return out

def has_absolute_machine_path(text):
    return bool(re.search(r'(?<![A-Za-z])[A-Za-z]:[\\/]|\\\\\?\\|/(?:Users|home|content|tmp)/',text))

def publish_summaries(run_dir,destination,experiment_root,drive_root,max_bytes=1_000_000):
    src=Path(run_dir);dest=Path(destination);root=Path(experiment_root).resolve();drive=Path(drive_root).resolve()
    dest.mkdir(parents=True,exist_ok=False)
    copied=[];skipped=[]
    patterns=['summary_*.csv','comparison_*.csv','volume_flags.csv','unclassified_codes.csv','report.md','carbon_density_map.png']
    paths=sorted(set(p for pattern in patterns for p in src.glob(pattern)))
    for path in paths:
        if path.stat().st_size>=max_bytes:skipped.append(path.name);continue
        if path.suffix=='.md':
            text=path.read_text(encoding='utf-8')
            for absolute in sorted({root,drive,root.parents[3]},key=lambda p:len(str(p)),reverse=True):
                original=strip_long_prefix(str(absolute));relative=os.path.relpath(original,strip_long_prefix(str(root))).replace('\\','/')
                for spelling in {original,original.replace('\\','/'),original.replace('\\','\\\\')}:
                    text=text.replace(spelling,relative)
            if has_absolute_machine_path(text):raise ValueError('Local absolute path remains in report')
            (dest/path.name).write_text(text,encoding='utf-8')
        else:shutil.copy2(path,dest/path.name)
        copied.append(path.name)
    provenance=json.loads((src/'provenance.json').read_text(encoding='utf-8'))
    shared=json.dumps(relative_provenance(provenance,root),ensure_ascii=False,indent=2,allow_nan=False)
    if has_absolute_machine_path(shared):raise ValueError('Local absolute path remains in shared provenance')
    if len(shared.encode('utf-8'))>=max_bytes:raise ValueError('Shared provenance exceeds size limit')
    (dest/'provenance.json').write_text(shared,encoding='utf-8');copied.append('provenance.json')
    return {'copied':copied,'skipped_over_size':skipped,'bytes':sum(p.stat().st_size for p in dest.iterdir() if p.is_file())}
