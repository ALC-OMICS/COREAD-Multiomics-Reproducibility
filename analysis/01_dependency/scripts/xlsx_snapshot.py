"""Read frozen XLSX cell values with the Python standard library; do not recalculate formulas."""
from pathlib import Path
from zipfile import ZipFile
import xml.etree.ElementTree as ET
NS={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
def parse_rows(path,sheet):
 with ZipFile(path) as z:
  strings=[]
  if 'xl/sharedStrings.xml' in z.namelist():
   strings=[''.join(t.text or '' for t in el.findall('.//s:t',NS)) for el in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',NS)]
  with z.open(f'xl/worksheets/sheet{sheet}.xml') as f:
   for ev,row in ET.iterparse(f,events=['end']):
    if row.tag!='{'+NS['s']+'}row':continue
    out={}
    for c in row.findall('s:c',NS):
     v=c.find('s:v',NS)
     if c.get('t')=='inlineStr':val=''.join(t.text or '' for t in c.findall('.//s:t',NS))
     elif v is not None:
      val=strings[int(v.text)] if c.get('t')=='s' else v.text
      if c.get('t') not in {'s','str','e'}:
       try:val=float(val);val=int(val) if val.is_integer() else val
       except (ValueError,TypeError):pass
     else:continue
     out[c.get('r')]=val
    if out:yield int(row.get('r')),out
    row.clear()
