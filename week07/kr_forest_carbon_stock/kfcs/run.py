"""python -m kfcs.run --config configs/geumjeongsan_2025.yaml"""
import argparse
import json
from pathlib import Path
from .config import load_config
from .aoi import read_boundary,clip
from .forest_map import read_forest_maps
from .volume import get_provider
from .factors import FactorSet
from .carbon import calculate
from .uncertainty import monte_carlo
from .report import write_outputs

def run(config):
    out=Path(config['output_dir'])
    # Exclusive mkdir is also the concurrency guard. Failed runs remain inspectable.
    out.mkdir(parents=True,exist_ok=False)
    try:
        boundary=read_boundary(config)
        g,audit=read_forest_maps(config,boundary)
        g=clip(g,boundary,config)
        provider=get_provider(config); volume=provider.provide(g)
        if not volume.index.equals(g.index) or len(volume)!=len(g): raise AssertionError('Volume provider changed rows')
        for column in volume: g[column]=volume[column]
        factors=FactorSet(config); g=calculate(factors.attach(g))
        mc=monte_carlo(g,factors,config)
        result=write_outputs(g,config,provider,factors,mc,audit,out)
        print(result.to_string(index=False));print(f'Outputs: {out}')
        return result
    except Exception as exc:
        (out/'FAILED.json').write_text(json.dumps({'error':type(exc).__name__,'message':str(exc)},ensure_ascii=False,indent=2),encoding='utf-8')
        raise

def main():
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--drive-root');p.add_argument('--output-dir')
    args=p.parse_args();run(load_config(args.config,args.drive_root,args.output_dir))

if __name__=='__main__': main()
