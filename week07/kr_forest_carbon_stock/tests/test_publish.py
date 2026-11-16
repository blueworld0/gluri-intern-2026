from pathlib import Path
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from kfcs.publish import publish_summaries,relative_provenance,has_absolute_machine_path

def test_only_small_outputs_are_shared(tmp_path):
    root=tmp_path/'repo'/'experiments'/'experiment';root.mkdir(parents=True)
    source=tmp_path/'drive'/'run';source.mkdir(parents=True)
    (source/'summary_total.csv').write_text('carbon_tC\n1\n',encoding='utf-8')
    (source/'polygons.gpkg').write_bytes(b'not shared')
    (source/'comparison_big.csv').write_bytes(b'x'*1000)
    (source/'report.md').write_text('Done',encoding='utf-8')
    (source/'provenance.json').write_text(json.dumps({'inputs':[{'path':str(source/'polygons.gpkg')}],'config':{'reference_dir':str(root/'reference')}}),encoding='utf-8')
    dest=root/'outputs'/'trial';result=publish_summaries(source,dest,root,tmp_path/'drive',max_bytes=1000)
    assert set(result['copied'])=={'summary_total.csv','report.md','provenance.json'}
    assert result['skipped_over_size']==['comparison_big.csv']
    assert not (dest/'polygons.gpkg').exists()
    content=json.loads((dest/'provenance.json').read_text(encoding='utf-8'))
    assert content['inputs'][0]['path'].startswith('../')
    assert content['config']['reference_dir']=='reference'
    assert not has_absolute_machine_path((dest/'provenance.json').read_text(encoding='utf-8'))

def test_remote_sources_remain_urls(tmp_path):
    out=relative_provenance({'url':'https://example.org/source','path':str(tmp_path/'data')},tmp_path)
    assert out['url']=='https://example.org/source'
    assert out['path']=='data'
