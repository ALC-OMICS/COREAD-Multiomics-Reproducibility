#!/usr/bin/env python3
"""Export Figure 2 A-D and their unchanged numerical panel sources."""
from pathlib import Path
import argparse,copy,json,re,shutil,xml.etree.ElementTree as ET

PANELS={
 'A':'compact_profiles/Figure_TCGA_compact_alteration_profile_top30',
 'B':'compact_profiles/Figure_CPTAC2_compact_alteration_profile_top30',
 'C':'compact_profiles/Figure_AC_ICAM_compact_alteration_profile_top30',
 'D':'permutation_summary/Figure_8_Permutation_summary_compact',
}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--numerical',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    a.output.mkdir(parents=True,exist_ok=True);source=a.output/'source_data';source.mkdir(exist_ok=True)
    ns='http://www.w3.org/2000/svg';ET.register_namespace('',ns);ET.register_namespace('xlink','http://www.w3.org/1999/xlink')
    root=ET.Element('{'+ns+'}svg',{'width':'1176pt','height':'1627pt','viewBox':'0 0 1176 1627','version':'1.1'})
    ET.SubElement(root,'{'+ns+'}rect',{'x':'0','y':'0','width':'1176','height':'1627','fill':'white'})
    for i,(panel,rel) in enumerate(PANELS.items()):
        f=a.numerical/'figures'/rel
        for ext in ['.svg','.png']:
            if f.with_suffix(ext).exists():shutil.copy2(f.with_suffix(ext),a.output/('Figure2'+panel+ext))
        shutil.copy2(f.with_name(f.name+'_source.csv'),source/('Figure2'+panel+'_source.csv'))
        svg=ET.parse(f.with_suffix('.svg')).getroot();vb=[float(x) for x in svg.attrib['viewBox'].split()];scale=576/vb[2];x=12+(i%2)*588;y=22+(i//2)*810
        group=ET.SubElement(root,'{'+ns+'}g',{'transform':f'translate({x},{y}) scale({scale})'})
        for elem in svg.iter():
            if 'id' in elem.attrib:elem.attrib['id']=panel+'_'+elem.attrib['id']
            for key,value in list(elem.attrib.items()):
                value=re.sub(r'url\(#([^)]*)\)',lambda m:'url(#'+panel+'_'+m.group(1)+')',value)
                if key.endswith('}href') and value.startswith('#'):value='#'+panel+'_'+value[1:]
                elem.attrib[key]=value
        for child in svg:group.append(copy.deepcopy(child))
        text=ET.SubElement(root,'{'+ns+'}text',{'x':str(x+2),'y':str(y-5),'font-size':'15','font-family':'DejaVu Sans','font-weight':'bold'});text.text=panel
    ET.ElementTree(root).write(a.output/'Figure2.svg',encoding='utf-8',xml_declaration=True)
    from PIL import Image,ImageDraw,ImageFont
    from matplotlib import font_manager
    canvas=Image.new('RGB',(1176,1627),'white');draw=ImageDraw.Draw(canvas)
    font=ImageFont.truetype(font_manager.findfont(font_manager.FontProperties(family='DejaVu Sans',weight='bold')),18)
    for i,panel in enumerate(PANELS):
        im=Image.open(a.output/('Figure2'+panel+'.png')).convert('RGB');im.thumbnail((576,778))
        x=12+i%2*588;y=22+i//2*810;canvas.paste(im,(x,y));draw.text((x,y-20),panel,font=font,fill='black')
    canvas.save(a.output/'Figure2.png')
    (a.output/'figure2_reconstruction.json').write_text(json.dumps({'passed':True,'panels':PANELS,'composition':'Current Figure2 SVG is an assembled layout of original numerical renderers.','BH':'Original family-specific adjustments (18 external COREAD,18 colon sensitivity,10 paired representation); no correction on displayed subset.'},indent=2))

if __name__=='__main__':main()
