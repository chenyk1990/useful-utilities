from rsf.tex import *
import os, sys
os.environ['PSTEXPENOPTS']='color=y fat=3 fatmult=1.5 serifs=n'

env=Environment()
if sys.platform == 'darwin':
    env['ENV']['PATH'] = ':'.join(['/opt/local/bin',env['ENV']['PATH']])

# for a common geophysics paper
Paper('enhemd',lclass='geophysics',use='amsmath,subfigure,listings')

# for revised paper
Paper('topic_revise1',lclass='geophysics',options='manuscript,revised',use='amsmath,subfigure,listings')
Command('topic_revise_plain1.tex','topic_revise1.tex','cp $SOURCE $TARGET')
Paper('topic_revise_plain1',lclass='geophysics',options='manuscript',use='amsmath')

# for reply letter
Paper('reply1',use='hyperref,listings')

# for a common seg abstract
Paper('seg9999_topic',lclass='segabs',use='amsmath,subfigure,listings')

# for seg abstract, without reference page
env.Command('seg2013_emdpf_final.pdf','seg2013_emdpf_final.ltx','pdflatex $SOURCE')

# for a common eage abstract
Paper('topiceage',lclass='mad15',use='listings, amsmath,subfigure') 

## Make sure there is a firstbreak.bst in your directory

# For slides
Paper('slides',lclass='beamer',use='subfigure,amsmath,amsbsy,listings,overpic,color,multicol,amssymb,comment,verbatim', include=r'''
      \mode<presentation>{\usetheme{AnnArbor}}
      \mode<presentation>{\usecolortheme{beaver}}
      \newcommand{\TEXMF}{%s/texmf}
      ''' % os.environ.get('HOME'))

## The following is used to automatically cp IEEE files
## scons ieee & scons paper.pdf 
env.Command(['ieee','IEEEtran.bst','IEEEtran.cls'],[os.getenv('HOME')+'/chenyk.open/temp/ieee/IEEEtran.bst',os.getenv('HOME')+'/chenyk.open/temp/ieee/IEEEtran.cls'],'cp  ${SOURCES[0]}  ${TARGETS[1]} && cp  ${SOURCES[1]}  ${TARGETS[2]} ')
Paper('ieee_filename',lclass='IEEEtran',options='10pt',use='cite,ifthen,seg,color,graphicx,subfigure,amsmath')
      
# For default paper.tex
End(options='manuscript',use='amsmath,subfigure,listing)
