"""Validate the candidate archive with standard library only; no analysis execution."""
from pathlib import Path
import ast, csv, hashlib, json, re, sys

ROOT=Path(__file__).resolve().parents[1]
def load_csv(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
    return h.hexdigest()
def validate():
    failures=[];warnings=[];checked=0;nb_count=0;py_count=0
    mf=ROOT/'manifests/sha256_manifest.csv'
    if not mf.exists():return {'status':'FAIL','failures':['sha256_manifest.csv missing']}
    entries=load_csv(mf)
    for row in entries:
        p=(ROOT/row['Release_file']).resolve()
        if not p.is_relative_to(ROOT):failures.append('Out-of-root manifest path');continue
        if not p.is_file():failures.append('MISSING:'+row['Release_file']);continue
        checked+=1
        if p.stat().st_size!=int(row['Bytes']) or digest(p)!=row['SHA256']:
            failures.append('CHECKSUM_MISMATCH:'+row['Release_file'])
        try:
            if p.suffix=='.py':ast.parse(p.read_text(encoding='utf-8-sig'));py_count+=1
            elif p.suffix=='.ipynb':
                n=json.loads(p.read_text(encoding='utf-8-sig'));nb_count+=1
                if n.get('nbformat')!=4 or not isinstance(n.get('cells'),list):raise ValueError('Notebook structure')
                for i,c in enumerate(n['cells']):
                    if c.get('cell_type') not in ['code','markdown','raw'] or not isinstance(c.get('source'),(list,str)):
                        raise ValueError('Invalid cell structure '+str(i))
                    if c['cell_type']=='code':
                        if not isinstance(c.get('outputs'),list):raise ValueError('Invalid outputs '+str(i))
                        if c['outputs']:failures.append('SAVED_OUTPUT_IN_GITHUB_NOTEBOOK:'+row['Release_file']+':'+str(i))
                        source=''.join(c['source']) if isinstance(c['source'],list) else c['source']
                        if any(re.match(r'^\s*%(?:pip|matplotlib|load_ext|time|timeit)\b',line) for line in source.splitlines()):continue
                        try:ast.parse(source)
                        except SyntaxError:warnings.append(row['Release_file']+':cell'+str(i)+':ordinary_Python_AST_error')
            elif p.suffix=='.json':json.loads(p.read_text(encoding='utf-8-sig'))
            elif p.suffix=='.csv':load_csv(p)
        except Exception as e:failures.append(row['Release_file']+':'+type(e).__name__)
    prov=load_csv(ROOT/'manifests/figure_table_provenance.csv')
    expected={f'Figure{i}' for i in range(1,9)}|{f'FigureS{i}' for i in range(1,12)}|{'Table1','Table2'}|{f'TableS{i}' for i in range(1,26)}
    if {r['Identifier'] for r in prov}!=expected:failures.append('Global Figure/Table identifier coverage')
    for row in prov:
        for rr in row.get('Source_data','').split(';'):
            if rr and not rr.startswith(('EXTERNAL:','WITHHELD_LOCAL:','MISSING:')) and not (ROOT/rr).is_file():
                failures.append('PROVENANCE_FILE_MISSING:'+rr)
    for d in ['main_tables','supplementary_tables']:
        want=2 if d=='main_tables' else 25
        if len(list((ROOT/'source_data'/d).glob('*_DATA.csv')))!=want:failures.append('Table DATA count:'+d)
    if any(x in (ROOT/'LICENSE').read_text(encoding='utf-8') for x in ['DRAFT','SELECTION PENDING']):warnings.append('License pending rights-holder selection')
    if 'authors:' not in (ROOT/'CITATION.cff').read_text(encoding='utf-8'):warnings.append('CITATION draft lacks confirmed authors; not schema-complete')
    warnings.append('Original Step21/25 producers missing; full raw-to-final reproduction is blocked')
    warnings.append('FigureS7A displayed39 versus frozen caption28: confirmed mismatch; scientific assets unchanged')
    return {'status':'FAIL' if failures else ('PASS_WITH_WARNINGS' if warnings else 'PASS'),
      'hashed_files':checked,'python_scripts':py_count,'notebook_copies':nb_count,
      'global_identifiers':len(expected),'failures':failures,'warnings':warnings,
      'scientific_execution':'NOT_PERFORMED','pickle_deserialization':'NOT_PERFORMED'}
if __name__=='__main__':
    result=validate()
    print(json.dumps(result,ensure_ascii=True,indent=2))
    sys.exit(1 if result['status']=='FAIL' else 0)
