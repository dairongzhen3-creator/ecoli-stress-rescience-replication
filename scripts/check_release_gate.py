#!/usr/bin/env python3
"""Fail-closed ReScience prepublication gate. Offline, no API, no upload."""
import argparse, json, pathlib, re, sys, yaml, hashlib
p=argparse.ArgumentParser()
p.add_argument('--scope',choices=['private','public'],default='public')
p.add_argument('--repo',default=None)
p.add_argument('--approval-file',default=None,help='Private author authorization receipt, kept outside public source repository')
a=p.parse_args()
root=pathlib.Path(a.repo).resolve() if a.repo else pathlib.Path(__file__).resolve().parents[1]
issues=[]; warnings=[]
def need(x,msg):
    if not x: issues.append(msg)
try: md=yaml.safe_load((root/'metadata.yaml').read_text(encoding='utf8'))
except Exception as e:
    print(json.dumps({'status':'INVALID_METADATA','reason':str(e)},indent=2));sys.exit(2)
need(isinstance(md,dict),'metadata.yaml is not a YAML mapping')
for k in ('title','authors','affiliations','keywords','code','data','replication','abstract','bibliography','type','domain','language'):
    need(md.get(k),f'metadata missing {k}')
need((root/'paper.md').is_file(),'paper.md missing')
need((root/'paper.bib').is_file(),'paper.bib missing')
need((root/'LICENSE').is_file(),'LICENSE missing')
need(len(list((root/'figures').glob('*.png')))==6,'six manuscript figures not present')
need(md.get('type')=='Replication','wrong article type')
need(md.get('bibliography')=='paper.bib','bibliography must reference paper.bib')
need(md.get('replication',[{}])[0].get('doi')=='10.1038/s44320-026-00226-5','source DOI mismatch')
need(md.get('data',[{}])[0].get('doi')=='10.5281/zenodo.18059846','source data DOI mismatch')
need('ChatGPT' in (root/'AI_DISCLOSURE.md').read_text(encoding='utf8'),'AI disclosure missing')
need(len(list((root/'src').glob('*.py')))>=10,'source scripts missing')
raw=[x for x in root.rglob('*') if x.is_file() and x.stat().st_size>10_000_000]
need(not raw,'unexpected large data file bundled: '+str([x.name for x in raw]))
# External release decisions are checked only for public scope, not smuggled into private CI.
if a.scope=='public':
    need('PRIVATE PRE-SUBMISSION DRAFT' not in (root/'metadata.yaml').read_text(encoding='utf8'),'private draft marker remains in metadata.yaml')
    contact=md.get('authors',[{}])[0].get('email')
    codeurl=md.get('code',[{}])[0].get('url')
    need(bool(contact) and '@' in str(contact),'author-confirmed contact email missing')
    need(str(codeurl).startswith('https://github.com/') and 'template' not in str(codeurl),'public repository URL missing')
    approval=pathlib.Path(a.approval_file).resolve() if a.approval_file else None
    need(approval is not None and approval.is_file(),'explicit private author approval record missing')
    if approval is not None and approval.is_file():
        try:
            decision=json.loads(approval.read_text())
            need(decision.get('author_approved_public_source_release') is True,'author source release approval not confirmed')
            need(decision.get('author_approved_paper_submission') is True,'author paper submission approval not confirmed')
        except Exception as e: issues.append('invalid approval record: '+str(e))
else:
    if not md.get('authors',[{}])[0].get('email'):warnings.append('Contact email unresolved (public release block)')
    if not md.get('code',[{}])[0].get('url'):warnings.append('Public code URL unresolved (public release block)')
result={'scope':a.scope,'status':'PASS' if not issues else 'BLOCKED','issues':issues,'warnings':warnings,'repo_files':sum(x.is_file() for x in root.rglob('*'))}
print(json.dumps(result,indent=2,ensure_ascii=False))
sys.exit(0 if not issues else 2)
