"""Read-only discovery; no result is treated as author code without provenance."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.parse import urlencode
import json,re
queries=['PDROD','"Tu Bowen" peridynamics','peridynamics reinforced concrete']

def query(q):
    url='https://api.github.com/search/repositories?'+urlencode({'q':q,'per_page':20})
    try:
        with urlopen(Request(url,headers={'User-Agent':'Bowen-reproduction-source-audit','Accept':'application/vnd.github+json'}),timeout=20) as r: data=json.load(r)
        return {'query':q,'url':url,'total_count':data.get('total_count'),
                'repositories':[{'name':v['full_name'],'url':v['html_url'],'description':v.get('description')} for v in data.get('items',[])]}
    except Exception as e: return {'query':q,'url':url,'error':str(e)}
with ThreadPoolExecutor(max_workers=3) as pool: results=list(pool.map(query,queries))
q='"基于近场动力学方法的混凝土梁构件冲击损伤研究"'
url='https://www.google.com/search?'+urlencode({'q':q})
try:
    with urlopen(Request(url,headers={'User-Agent':'Mozilla/5.0'}),timeout=20) as r: body=r.read(400000).decode('utf-8','replace')
    Path('outputs/public_search_google.html').write_text(body)
    clean=re.sub('<[^>]+>',' ',body)
    results.append({'query':q,'url':url,'response_text_excerpt':re.sub(r'\s+',' ',clean)[:5000]})
except Exception as e: results.append({'query':q,'url':url,'error':str(e)})
Path('outputs/public_source_search.json').write_text(json.dumps(results,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(results,ensure_ascii=False,indent=2))
