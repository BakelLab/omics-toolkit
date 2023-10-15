#!/usr/bin/env python

import sys
import argparse
import subprocess
import os

def colorstr(rgb): return "#%02x%02x%02x" % (rgb[0],rgb[1],rgb[2])


def hsl_to_rgb(h, s, l):
    c = (1 - abs(2*l - 1)) * s
    x = c * (1 - abs(h *1.0 / 60 % 2 - 1))
    m = l - c/2
    if h < 60:
        r, g, b = c + m, x + m, 0 + m
    elif h < 120:
        r, g, b = x + m, c+ m, 0 + m
    elif h < 180:
        r, g, b = 0 + m, c + m, x + m
    elif h < 240:
        r, g, b, = 0 + m, x + m, c + m
    elif h < 300:
        r, g, b, = x + m, 0 + m, c + m
    else:
        r, g, b, = c + m, 0 + m, x + m
    r = int(r * 255)
    g = int(g * 255)
    b = int(b * 255)
    return (r,g,b)

class scalableVectorGraphics:

    def __init__(self, height, width):
        self.height = height
        self.width = width
        self.out = '''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg
   xmlns:dc="http://purl.org/dc/elements/1.1/"
   xmlns:cc="http://creativecommons.org/ns#"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
   xmlns:svg="http://www.w3.org/2000/svg"
   xmlns="http://www.w3.org/2000/svg"
   xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
   xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
   height="%d"
   width="%d"
   id="svg2"
   version="1.1"
   inkscape:version="0.48.4 r9939"
   sodipodi:docname="easyfig">
  <metadata
     id="metadata122">
    <rdf:RDF>
      <cc:Work
         rdf:about="">
        <dc:format>image/svg+xml</dc:format>
        <dc:type
           rdf:resource="http://purl.org/dc/dcmitype/StillImage" />
        <dc:title>Easyfig</dc:title>
      </cc:Work>
    </rdf:RDF>
  </metadata>
  <defs
     id="defs120" />
  <sodipodi:namedview
     pagecolor="#ffffff"
     bordercolor="#666666"
     borderopacity="1"
     objecttolerance="10"
     gridtolerance="10"
     guidetolerance="10"
     inkscape:pageopacity="0"
     inkscape:pageshadow="2"
     inkscape:window-width="640"
     inkscape:window-height="480"
     id="namedview118"
     showgrid="false"
     inkscape:zoom="0.0584"
     inkscape:cx="2500"
     inkscape:cy="75.5"
     inkscape:window-x="55"
     inkscape:window-y="34"
     inkscape:window-maximized="0"
     inkscape:current-layer="svg2" />
  <title
     id="title4">Easyfig</title>
  <g
     style="fill-opacity:1.0; stroke:black; stroke-width:1;"
     id="g6">''' % (self.height, self.width)

    def drawLine(self, x1, y1, x2, y2, th=1, cl=(0, 0, 0), alpha = 1.0, linecap='round'):
        self.out += '  <line x1="%d" y1="%d" x2="%d" y2="%d"\n        stroke-width="%d" stroke="%s" stroke-opacity="%f" stroke-linecap="%s" />\n' % (x1, y1, x2, y2, th, colorstr(cl), alpha, linecap)

    def drawPath(self, xcoords, ycoords, th=1, cl=(0, 0, 0), alpha=0.9):
        self.out += '  <path d="M%d %d' % (xcoords[0], ycoords[0])
        for i in range(1, len(xcoords)):
            self.out += ' L%d %d' % (xcoords[i], ycoords[i])
        self.out += '"\n        stroke-width="%d" stroke="%s" stroke-opacity="%f" stroke-linecap="butt" fill="none" z="-1" />\n' % (th, colorstr(cl), alpha)


    def writesvg(self, filename):
        outfile = open(filename, 'w')
        outfile.write(self.out + ' </g>\n</svg>')
        outfile.close()

    def drawRightArrow(self, x, y, wid, ht, fc, oc=(0,0,0), lt=1):
        if lt > ht /2:
            lt = ht/2
        x1 = x + wid
        y1 = y + ht/2
        x2 = x + wid - ht / 2
        ht -= 1
        if wid > ht/2:
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d"\n' % (colorstr(fc), colorstr(oc), lt)
            self.out += '           points="%d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d" />\n' % (x, y+ht/4, x2, y+ht/4,
                                                                                                x2, y, x1, y1, x2, y+ht,
                                                                                                x2, y+3*ht/4, x, y+3*ht/4)
        else:
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d"\n' % (colorstr(fc), colorstr(oc), lt)
            self.out += '           points="%d,%d %d,%d %d,%d" />\n' % (x, y, x, y+ht, x + wid, y1)

    def drawLeftArrow(self, x, y, wid, ht, fc, oc=(0,0,0), lt=1):
        if lt > ht /2:
            lt = ht/2
        x1 = x + wid
        y1 = y + ht/2
        x2 = x + ht / 2
        ht -= 1
        if wid > ht/2:
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d"\n' % (colorstr(fc), colorstr(oc), lt)
            self.out += '           points="%d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d" />\n' % (x1, y+ht/4, x2, y+ht/4,
                                                                                                x2, y, x, y1, x2, y+ht,
                                                                                                x2, y+3*ht/4, x1, y+3*ht/4)
        else:
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d"\n' % (colorstr(fc), colorstr(oc), lt)
            self.out += '           points="%d,%d %d,%d %d,%d" />\n' % (x, y1, x1, y+ht, x1, y)

    def drawBlastHit(self, x1, y1, x2, y2, x3, y3, x4, y4, fill=(0, 0, 255), lt=2, alpha=0.1):
        self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d" fill-opacity="%f"\n' % (colorstr(fill), colorstr((0,0,0)), lt, alpha)
        self.out += '           points="%d,%d %d,%d %d,%d %d,%d" />\n' % (x1, y1, x2, y2, x3, y3, x4, y4)

    def drawGradient(self, x1, y1, wid, hei, minc, maxc):
        self.out += '  <defs>\n    <linearGradient id="MyGradient" x1="0%" y1="0%" x2="0%" y2="100%">\n'
        self.out += '      <stop offset="0%%" stop-color="%s" />\n' % colorstr(maxc)
        self.out += '      <stop offset="100%%" stop-color="%s" />\n' % colorstr(minc)
        self.out += '    </linearGradient>\n  </defs>\n'
        self.out += '  <rect fill="url(#MyGradient)" stroke-width="0"\n'
        self.out += '        x="%d" y="%d" width="%d" height="%d"/>\n' % (x1, y1, wid, hei)

    def drawGradient2(self, x1, y1, wid, hei, minc, maxc):
        self.out += '  <defs>\n    <linearGradient id="MyGradient2" x1="0%" y1="0%" x2="0%" y2="100%">\n'
        self.out += '      <stop offset="0%%" stop-color="%s" />\n' % colorstr(maxc)
        self.out += '      <stop offset="100%%" stop-color="%s" />\n' % colorstr(minc)
        self.out += '    </linearGradient>\n</defs>\n'
        self.out += '  <rect fill="url(#MyGradient2)" stroke-width="0"\n'
        self.out += '        x="%d" y="%d" width="%d" height="%d" />\n' % (x1, y1, wid, hei)

    def drawOutRect(self, x1, y1, wid, hei, fill=(255, 255, 255), outfill=(0, 0, 0), lt=1, alpha=1.0, alpha2=1.0):
        self.out += '  <rect stroke="%s" stroke-width="%d" stroke-opacity="%f"\n' % (colorstr(outfill), lt, alpha)
        self.out += '        fill="%s" fill-opacity="%f"\n' % (colorstr(fill), alpha2)
        self.out += '        x="%d" y="%d" width="%d" height="%d" />\n' % (x1, y1, wid, hei)

    def drawAlignment(self, x, y, fill, outfill, lt=1, alpha=1.0, alpha2=1.0):
        self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d" stroke-opacity="%f" fill-opacity="%f"\n' % (colorstr(fill), outfill, lt, alpha, alpha2)
        self.out += '  points="'
        for i, j in zip(x, y):
            self.out += str(i) + ',' + str(j) + ' '
        self.out += '" />\n'
             # print self.out.split('\n')[-2]



    def drawSymbol(self, x, y, size, fill, symbol, alpha=1.0, lt=1):
        x0 = x - size/2
        x1 = size/8 + x - size/2
        x2 = size/4 + x - size/2
        x3 = size*3/8 + x - size/2
        x4 = size/2 + x - size/2
        x5 = size*5/8 + x - size/2
        x6 = size*3/4 + x - size/2
        x7 = size*7/8 + x - size/2
        x8 = size + x - size/2
        y0 = y - size/2
        y1 = size/8 + y - size/2
        y2 = size/4 + y - size/2
        y3 = size*3/8 + y - size/2
        y4 = size/2 + y - size/2
        y5 = size*5/8 + y - size/2
        y6 = size*3/4 + y - size/2
        y7 = size*7/8 + y - size/2
        y8 = size + y - size/2
        if symbol == 'o':
            self.out += '  <circle stroke="%s" stroke-width="%d" stroke-opacity="%f"\n' % (colorstr((0, 0, 0)), lt, alpha)
            self.out += '        fill="%s" fill-opacity="%f"\n' % (colorstr(fill), alpha)
            self.out += '        xc="%d" yc="%d" r="%d" />\n' % (x, y, size/2)
        elif symbol == 'x':
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d" stroke-opacity="%f" fill-opacity="%f"\n' % (colorstr(fill), colorstr((0, 0, 0)), lt, alpha, alpha)
            self.out += '           points="%d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d" />\n' % (x0, y2, x2, y0, x4, y2, x6, y0, x8, y2,
                                                                                                                             x6, y4, x8, y6, x6, y8, x4, y6, x2, y8,
                                                                                                                             x0, y6, x2, y4)
        elif symbol == '+':
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d" stroke-opacity="%f" fill-opacity="%f"\n' % (colorstr(fill), colorstr((0, 0, 0)), lt, alpha, alpha)
            self.out += '           points="%d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d" />\n' % (x2, y0, x6, y0, x6, y2, x8, y2, x8, y6,
                                                                                                                             x6, y6, x6, y8, x2, y8, x2, y6, x0, y6,
                                                                                                                             x0, y2, x2, y2)
        elif symbol == 's':
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d" stroke-opacity="%f" fill-opacity="%f"\n' % (colorstr(fill), colorstr((0, 0, 0)), lt, alpha, alpha)
            self.out += '           points="%d,%d %d,%d %d,%d %d,%d" />\n' % (x0, y0, x0, y8, x8, y8, x8, y0)
        elif symbol == '^':
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d" stroke-opacity="%f" fill-opacity="%f"\n' % (colorstr(fill), colorstr((0, 0, 0)), lt, alpha, alpha)
            self.out += '           points="%d,%d %d,%d %d,%d %d,%d %d,%d %d,%d" />\n' % (x0, y0, x2, y0, x4, y4, x6, y0, x8, y0, x4, y8)
        elif symbol == 'v':
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d" stroke-opacity="%f" fill-opacity="%f"\n' % (colorstr(fill), colorstr((0, 0, 0)), lt, alpha, alpha)
            self.out += '           points="%d,%d %d,%d %d,%d %d,%d %d,%d %d,%d" />\n' % (x0, y8, x2, y8, x4, y4, x6, y8, x8, y8, x4, y0)
        elif symbol == 'u':
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d" stroke-opacity="%f" fill-opacity="%f"\n' % (colorstr(fill), colorstr((0, 0, 0)), lt, alpha, alpha)
            self.out += '           points="%d,%d %d,%d %d,%d" />\n' % (x0, y8, x4, y0, x8, y8)
        elif symbol == 'd':
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d" stroke-opacity="%f" fill-opacity="%f"\n' % (colorstr(fill), colorstr((0, 0, 0)), lt, alpha, alpha)
            self.out += '           points="%d,%d %d,%d %d,%d" />\n' % (x0, y0, x4, y8, x8, y0)
        else:
            sys.stderr.write(symbol + '\n')
            sys.stderr.write('Symbol not found, this should not happen.. exiting')
            sys.exit()








    def drawRightFrame(self, x, y, wid, ht, lt, frame, fill):
        if lt > ht /2:
            lt = ht /2
        if frame == 1:
            y1 = y + ht/2
            y2 = y + ht * 3/8
            y3 = y + ht * 1/4
        elif frame == 2:
            y1 = y + ht * 3/8
            y2 = y + ht * 1/4
            y3 = y + ht * 1/8
        elif frame == 0:
            y1 = y + ht * 1/4
            y2 = y + ht * 1/8
            y3 = y + 1
        x1 = x
        x2 = x + wid - ht/8
        x3 = x + wid
        if wid > ht/8:
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d"\n' % (colorstr(fill), colorstr((0, 0, 0)), lt)
            self.out += '           points="%d,%d %d,%d %d,%d %d,%d %d,%d" />\n' % (x1, y1, x2, y1, x3, y2, x2, y3, x1, y3)
        else:
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d"\n' % (colorstr(fill), colorstr((0, 0, 0)), lt)
            self.out += '           points="%d,%d %d,%d %d,%d" />\n' % (x1, y1, x3, y2, x1, y3)

    def drawRightFrameRect(self, x, y, wid, ht, lt, frame, fill):
        if lt > ht /2:
            lt = ht /2
        if frame == 1:
            y1 = y + ht / 4
        elif frame == 2:
            y1 = y + ht /8
        elif frame == 0:
            y1 = y + 1
        hei = ht /4
        x1 = x
        self.out += '  <rect fill="%s" stroke-width="%d"\n' % (colorstr(fill), lt)
        self.out += '        x="%d" y="%d" width="%d" height="%d" />\n' % (x1, y1, wid, hei)

    def drawLeftFrame(self, x, y, wid, ht, lt, frame, fill):
        if lt > ht /2:
            lt = ht /2
        if frame == 1:
            y1 = y + ht
            y2 = y + ht * 7/8
            y3 = y + ht * 3/4
        elif frame == 2:
            y1 = y + ht * 7/8
            y2 = y + ht * 3/4
            y3 = y + ht * 5/8
        elif frame == 0:
            y1 = y + ht * 3/4
            y2 = y + ht * 5/8
            y3 = y + ht / 2
        x1 = x + wid
        x2 = x + ht/8
        x3 = x
        if wid > ht/8:
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d"\n' % (colorstr(fill), colorstr((0, 0, 0)), lt)
            self.out += '           points="%d,%d %d,%d %d,%d %d,%d %d,%d" />\n' % (x1, y1, x2, y1, x3, y2, x2, y3, x1, y3)
        else:
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d"\n' % (colorstr(fill), colorstr((0, 0, 0)), lt)
            self.out += '           points="%d,%d %d,%d %d,%d" />\n' % (x1, y1, x3, y2, x1, y3)

    def drawLeftFrameRect(self, x, y, wid, ht, lt, frame, fill):
        if lt > ht /2:
            lt = ht /2
        if frame == 1:
            y1 = y + ht * 3/4
        elif frame == 2:
            y1 = y + ht * 5/8
        elif frame == 0:
            y1 = y + ht / 2
        hei = ht /4
        x1 = x
        self.out += '  <rect fill="%s" stroke-width="%d"\n' % (colorstr(fill), lt)
        self.out += '        x="%d" y="%d" width="%d" height="%d" />\n' % (x1, y1, wid, hei)

    def drawPointer(self, x, y, ht, lt, fill):
        x1 = x - int(round(0.577350269 * ht/2))
        x2 = x + int(round(0.577350269 * ht/2))
        y1 = y + ht/2
        y2 = y + 1
        self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d"\n' % (colorstr(fill), colorstr((0, 0, 0)), lt)
        self.out += '           points="%d,%d %d,%d %d,%d" />\n' % (x1, y2, x2, y2, x, y1)

    def drawDash(self, x1, y1, x2, y2, exont):
        self.out += '  <line x1="%d" y1="%d" x2="%d" y2="%d"\n' % (x1, y1, x2, y2)
        self.out += '       style="stroke-dasharray: 5, 3, 9, 3"\n'
        self.out += '       stroke="#000" stroke-width="%d" />\n' % exont

    def drawPolygon(self, x_coords, y_coords, colour=(0,0,255)):
        self.out += '  <polygon points="'
        for i,j in zip(x_coords, y_coords):
            self.out += str(i) + ',' + str(j) + ' '
        self.out += '"\nstyle="fill:%s;stroke=none" />\n'  % colorstr(colour)
    def writeString(self, thestring, x, y, size, ital=False, bold=False, rotate=0, justify='left'):
        if rotate != 0:
            x, y = y, x
        self.out += '  <text\n'
        self.out += '    style="font-size:%dpx;font-style:normal;font-weight:normal;z-index:10\
;line-height:125%%;letter-spacing:0px;word-spacing:0px;fill:#111111;fill-opacity:1;stroke:none;font-family:Sans"\n' % size
        if justify == 'right':
            self.out += '    text-anchor="end"\n'
        elif justify == 'middle':
            self.out += '    text-anchor="middle"\n'
        if rotate == 1:
            self.out += '    x="-%d"\n' % x
        else:
            self.out += '    x="%d"\n' % x
        if rotate == -1:
            self.out += '    y="-%d"\n' % y
        else:
            self.out += '    y="%d"\n' % y
        self.out += '    sodipodi:linespacing="125%"'
        if rotate == -1:
            self.out += '\n    transform="matrix(0,1,-1,0,0,0)"'
        if rotate == 1:
            self.out += '\n    transform="matrix(0,-1,1,0,0,0)"'
        self.out += '><tspan\n      sodipodi:role="line"\n'
        if rotate == 1:
            self.out += '      x="-%d"\n' % x
        else:
            self.out += '      x="%d"\n' % x
        if rotate == -1:
            self.out += '      y="-%d"' % y
        else:
            self.out += '      y="%d"' % y
        if ital and bold:
            self.out += '\nstyle="font-style:italic;font-weight:bold"'
        elif ital:
            self.out += '\nstyle="font-style:italic"'
        elif bold:
            self.out += '\nstyle="font-style:normal;font-weight:bold"'
        self.out += '>' + thestring + '</tspan></text>\n'



def translate_dna(dna):
    code = {     'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
         'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
         'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
         'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
         'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
         'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
         'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
         'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
         'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
         'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
         'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
         'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
         'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
         'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
         'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
         'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'
    }
    protein = ''
    dna = dna.lower()
    for i in range(0, len(dna), 3):
        if dna[i:i+3] in code:
            protein += code[dna[i:i+3]]
        else:
            protein += 'X'
    return protein

def reverse_compliment(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    return "".join(complement.get(base, base) for base in reversed(seq))


def get_genes_gff(gff_file):
    with open(gff_file) as gff:
        gene_dict = {}
        seqDict = {}
        getseq = False
        for line in gff:
            if line.startswith('#'):
                pass
            elif line.startswith('>'):
                getseq = True
                name = line.split()[0][1:]
                seqDict[name] = ''
            elif getseq:
                seqDict[name] += line.rstrip()
            else:
                contig, method, feature, start, stop, score, strand, phase, extra = line.rstrip().split('\t')
                if not contig in gene_dict:
                    gene_dict[contig] = []
                genes = gene_dict[contig]
                gene= 'none'
                uniprot = 'none'
                product = 'none'
                seq = 'none'
                locus_tag = 'none'
                for i in extra.split(';'):
                    if i.startswith('Name='):
                        gene = i.split('=')[1]
                    elif i.startswith('locus_tag'):
                        locus_tag = i.split('=')[1]
                    elif i.startswith('product'):
                        product = i.split('=')[1]
                    elif i.startswith('locus_tag'):
                        locus_tag = i.split('=')[1]
                if feature == 'CDS':
                    genes.append((int(start), int(stop), strand, gene, locus_tag, seq, product, uniprot))
    return gene_dict, seqDict




def get_genes(gbk):
    with open(gbk) as gbk:
        gene_dict = {}
        seqDict = {}
        getseq2 = False
        getseq = False
        getproduct = False
        for line in gbk:
            if line.startswith('LOCUS'):
                contig_name = line.split()[1]
                gene_dict[contig_name] = []
                genes = gene_dict[contig_name]
            elif line.startswith('     CDS             complement('):
                startstop = line.split('(')[1].split(')')[0]
                start, stop = map(int, startstop.split('..'))
                strand = '-'
                gene = 'none'
                uniprot = 'none'
            elif line.startswith('     CDS '):
                startstop = line.split()[1]
                strand = '+'
                start, stop = map(int, startstop.split('..'))
                gene = 'none'
                uniprot = 'none'
            elif line.startswith('                     /locus_tag'):
                locus_tag = line.split('"')[1]
            elif line.startswith('                     /gene='):
                gene = line.split('"')[1]
            elif line.startswith('                     /inference="similar to AA sequence:UniProtKB:'):
                uniprot = line.split(':')[-1].split('"')[0]
            elif line.startswith('                     /product=') or getproduct:
                if line.startswith('                     /product='):
                    getproduct = True
                    product = line.rstrip().split('=')[1]
                else:
                    product += ' ' + line.rstrip()[21:]
                if product.endswith('"'):
                    getproduct = False
            elif line.startswith('                     /translation='):
                seq = line.rstrip().split('"')[1]
                if line.count('"') == 2:
                    genes.append((start, stop, strand, gene, locus_tag, seq, product, uniprot))
                else:
                    getseq = True
            elif getseq:
                seq += line.split()[0]
                if seq.endswith('"'):
                    genes.append((start, stop, strand, gene, locus_tag, seq[:-1], product, uniprot))
                    getseq = False
            elif line.startswith('ORIGIN'):
                getseq2 = True
                seq = ''
            elif line.startswith('//'):
                seqDict[contig_name] = seq
                getseq2 = False
            elif getseq2:
                seq += ''.join(line.split()[1:])
    return gene_dict, seqDict

def read_nucdiff(gffs, query_genbank, ref_genbank, output, working_dir, ref=False, merge=False, get_indel=True, promoter_region=500):
    if ref:
        next = '_ref_'
    else:
        next = '_query_'
    if query_genbank.endswith('.gff'):
        query_genes, query_seq = get_genes_gff(query_genbank)
        ref_genes, ref_seq = get_genes_gff(ref_genbank)
    else:
        query_genes, query_seq = get_genes(query_genbank)
        ref_genes, ref_seq = get_genes(ref_genbank)
    with open(output + '.gff', 'w') as o:
        o.write('# QUERY_GBK=' + query_genbank + '\n')
        o.write('# REF_GBK=' + ref_genbank + '\n')
        for gff in gffs:
            with open(gff + next + 'snps.gff') as snps:
                for line in snps:
                    if line.startswith('#'):
                        o.write(line)
                    if not line.startswith('#'):
                        contig, program, so, query_start, query_stop, score, strand, phase, extra = line.rstrip().split('\t')
                        if not contig in query_genes:
                            sys.exit('You may have switched query and reference genbanks.')
                        query_start = int(query_start)
                        query_stop = int(query_stop)
                        extra_dict = {}
                        for i in extra.split(';'):
                            key, value = i.split('=')
                            extra_dict[key] = value
                        if extra_dict['Name'] in ['deletion', 'insertion'] and not get_indel:
                            continue
                        for i in query_genes[contig]:
                            start, stop, strand, gene, locus, seq, product, uniprot = i
                            if start <= query_start <= stop or start <= query_stop <= stop:
                                gene_seq = query_seq[contig][start-1:stop]
                                gene_seq_altered = query_seq[contig][start-1:query_start-1] + extra_dict['ref_bases'] + query_seq[contig][query_stop:stop]
                                codon_start = query_start - (query_start-start)%3-1
                                codon_seq = query_seq[contig][codon_start:codon_start+3]
                                alt_codon = query_seq[contig][codon_start:query_start-1] + extra_dict['ref_bases'] + query_seq[contig][query_stop:codon_start+3]
                                if strand == '-':
                                    gene_seq = reverse_compliment(gene_seq)
                                    gene_seq_altered = reverse_compliment(gene_seq_altered)
                                    codon_seq = reverse_compliment(codon_seq)
                                    alt_codon = reverse_compliment(alt_codon)
                                aa_seq = translate_dna(gene_seq)
                                aa_seq_altered = translate_dna(gene_seq_altered)
                                if extra_dict['Name'] in ['deletion', 'insertion', 'deletion_promoter', 'insertion_promoter']:
                                    if 'deletion' in extra_dict['Name']:
                                        extra_dict['Name'] = 'frameshift_del'
                                    elif 'insertion' in extra_dict['Name']:
                                        extra_dict['Name'] = 'frameshift_ins'
                                    else:
                                        print('error')
                                elif not '*' in aa_seq_altered:
                                    extra_dict['Name'] = 'stop_gain'
                                elif '*' in aa_seq_altered[:-1]:
                                    extra_dict['Name'] = 'stop_loss'
                                elif aa_seq == aa_seq_altered:
                                    extra_dict['Name'] = 'synonymous'
                                else:
                                    extra_dict['Name'] = 'nonsynonymous'
                                    extra_dict['aa_sub'] = str((query_start-start)/3+1) + '>' + translate_dna(codon_seq) + '>' + translate_dna(alt_codon)
                                extra_dict['in_genes'] = locus
                                extra_dict['in_genes_name'] = gene
                                extra_dict['in_gene_uniprot'] = uniprot
                                break

                            elif start - promoter_region <= query_start < start and strand == '+':
                                if extra_dict['Name'] in ['deletion', 'insertion']:
                                    extra_dict['Name'] = extra_dict['Name'] + '_promoter'
                                else:
                                    extra_dict['Name'] = 'promoter'
                                extra_dict['in_genes'] = locus
                                extra_dict['in_genes_name'] = gene
                                extra_dict['in_gene_uniprot'] = uniprot
                            elif stop + promoter_region >= query_start > stop and strand == '-':
                                if extra_dict['Name'] in ['deletion', 'insertion']:
                                    extra_dict['Name'] = extra_dict['Name'] + '_promoter'
                                else:
                                    extra_dict['Name'] = 'promoter'
                                extra_dict['in_genes'] = locus
                                extra_dict['in_genes_name'] = gene
                                extra_dict['in_gene_uniprot'] = uniprot
                        new_extra = []
                        for i in extra_dict:
                            new_extra.append(i + '=' + extra_dict[i])
                        new_extra.sort()
                        o.write('\t'.join([contig, program, so, str(query_start), str(query_stop), score, strand, phase, ';'.join(new_extra)]) + '\n')
        for gff in gffs:
            with open(gff + next + 'struct.gff') as struct:
                for line in struct:
                    if not line.startswith('#'):
                        contig, program, so, query_start, query_stop, score, strand, phase, extra = line.rstrip().split('\t')
                        if not contig in query_genes:
                            sys.exit('You may have switched query and reference genbanks.')
                        query_start = int(query_start)
                        query_stop = int(query_stop)
                        extra_dict = {}
                        for i in extra.split(';'):
                            key, value = i.split('=')
                            extra_dict[key] = value
                        for i in query_genes[contig]:
                            start, stop, strand, gene, locus, seq, product, uniprot = i
                            if start <= query_start <= query_stop <= stop:
                                if 'in_gene' in extra_dict:
                                    extra_dict['in_genes'] += ',' + locus
                                    extra_dict['in_genes_name'] += ',' + gene
                                    extra_dict['in_gene_uniprot'] = uniprot
                                else:
                                    extra_dict['in_genes'] = locus
                                    extra_dict['in_genes_name'] = gene
                                    extra_dict['in_gene_uniprot'] = uniprot
                            elif query_start <= start <= stop <= query_stop:
                                if query_stop - query_start < 100000:
                                    if 'contains_gene' in extra_dict:
                                        extra_dict['contains_genes'] += ',' + locus
                                        extra_dict['contains_genes_name'] += ',' + gene
                                        extra_dict['contains_gene_uniprot'] += ',' + uniprot
                                    else:
                                        extra_dict['contains_genes'] = locus
                                        extra_dict['contains_genes_name'] = gene
                                        extra_dict['contains_genes_uniprot'] = uniprot
                                # else:
                                    # print line.rstrip()
                            elif start <= query_start <= stop or start <= query_stop <= stop:
                                if 'partial_overlap' in extra_dict:
                                    extra_dict['partial_overlap'] += ',' + locus
                                    extra_dict['partial_overlap_name'] += ',' + gene
                                    extra_dict['partial_overlap_uniprot'] += ',' + uniprot
                                else:
                                    extra_dict['partial_overlap'] = locus
                                    extra_dict['partial_overlap_name'] = gene
                                    extra_dict['partial_overlap_uniprot'] = uniprot
                        ref_contig = extra_dict['ref_sequence']
                        if 'ref_coord' in extra_dict:
                            ref_coord = extra_dict['ref_coord']
                            if '-' in ref_coord:
                                ref_start, ref_stop = map(int, ref_coord.split('-'))
                            else:
                                ref_start = ref_stop = int(ref_coord)
                        elif 'blk_1_ref':
                            ref_coord = extra_dict['blk_1_ref']
                            ref_start, ref_stop = map(int, ref_coord.split('-'))
                        else:
                            print(line.rstrip())
                        for i in ref_genes[ref_contig]:
                            start, stop, strand, gene, locus, seq, product, uniprot = i
                            if start <= ref_start <= ref_stop <= stop:
                                if 'in_gene_ref' in extra_dict:
                                    extra_dict['in_genes_ref'] += ',' + locus
                                    extra_dict['in_genes_ref_name'] += ',' + gene
                                    extra_dict['in_genes_ref_uniprot'] += ',' + uniprot
                                else:
                                    extra_dict['in_genes_ref'] = locus
                                    extra_dict['in_genes_ref_name'] = gene
                                    extra_dict['in_genes_ref_uniprot'] = uniprot
                            elif ref_start <= start <= stop <= ref_stop:
                                if ref_stop - ref_start < 100000:
                                    if 'contains_gene_ref' in extra_dict:
                                        extra_dict['contains_genes_ref'] += ',' + locus
                                        extra_dict['contains_genes_ref_name'] += ',' + gene
                                        extra_dict['contains_genes_ref_uniprot'] += ',' + uniprot
                                    else:
                                        extra_dict['contains_genes_ref'] = locus
                                        extra_dict['contains_genes_ref_name'] = gene
                                        extra_dict['contains_genes_ref_uniprot'] = uniprot

                                # else:
                                    # print line.rstrip()
                            elif start <= ref_start <= stop or start <= ref_stop <= stop:
                                if 'partial_overlap_ref' in extra_dict:
                                    extra_dict['partial_overlap_ref'] += ',' + locus
                                    extra_dict['partial_overlap_ref_name'] += ',' + gene
                                    extra_dict['partial_overlap_ref_uniprot'] += ',' + uniprot
                                else:
                                    extra_dict['partial_overlap_ref'] = locus
                                    extra_dict['partial_overlap_ref_name'] = gene
                                    extra_dict['partial_overlap_ref_uniprot'] = uniprot
                        new_extra = []
                        for i in extra_dict:
                            new_extra.append(i + '=' + extra_dict[i])
                        new_extra.sort()
                        o.write('\t'.join([contig, program, so, str(query_start), str(query_stop), score, strand, phase, ';'.join(new_extra)]) + '\n')
        for i in query_genes:
            plasmid = None
            for j in gffs:
                if os.path.basename(j).startswith(os.path.splitext(os.path.basename(query_genbank))[0] + '.' + i +  'vs'):
                    plasmid = i
            if plasmid is None:
                locus_list = []
                gene_list = []
                uniprot_list = []
                for j in query_genes[i]:
                    start, stop, strand, gene, locus, seq, product, uniprot = j
                    locus_list.append(locus)
                    gene_list.append(gene)
                    uniprot_list.append(uniprot)
                extra = 'Name=Plasmid_loss;contains_gene=' + ','.join(locus_list) + ';contains_gene_name=' +\
                ','.join(gene_list) + ';contains_gene_uniprot=' + ','.join(uniprot_list)
                o.write('\t'.join([i, 'getVar', 'SO:0001059', '1', str(len(query_seq[i])), '.', '.', '.', extra]) + '\n')
        for i in ref_genes:
            plasmid = None
            for j in gffs:
                if os.path.basename(j).endswith('vs' + os.path.splitext(os.path.basename(ref_genbank))[0] + '.' + i):
                    plasmid = i
            if plasmid is None:
                locus_list = []
                gene_list = []
                uniprot_list =[]
                for j in ref_genes[i]:
                    start, stop, strand, gene, locus, seq, product, uniprot = j
                    locus_list.append(locus)
                    gene_list.append(gene)
                    uniprot_list.append(uniprot)
                extra = 'Name=Plasmid_loss;ref_sequence=' + i + ';contains_gene_ref=' + ','.join(locus_list) + ';contains_gene_ref_name=' +\
                ','.join(gene_list) + ';contains_gene_ref_uniprot=' + ','.join(uniprot_list)
                o.write('\t'.join([i, 'getVar', 'SO:0001059', '1', str(len(ref_seq[i])), '.', '.', '.', extra]) + '\n')




def gbk_to_fasta(gbk, out, concat=True):
    length_dict = {}
    getseq = False
    if concat:
        with open(gbk) as f, open(out, 'w') as o:
            for line in f:
                if line.startswith('LOCUS'):
                    name = line.split()[1]
                    length_dict[name] = int(line.split()[2])
                    o.write('>' + name + '\n')
                elif line.startswith('ORIGIN'):
                    getseq = True
                elif line.startswith('//'):
                    getseq = False
                elif getseq:
                    o.write(''.join(line.split()[1:]) + '\n')
    else:
        with open(gbk) as f:
            for line in f:
                if line.startswith('LOCUS'):
                    name = line.split()[1]
                    length_dict[name] = int(line.split()[2])
                    o = open(out + '.' + name + '.fa', 'w')
                    o.write('>' + name + '\n')
                elif line.startswith('ORIGIN'):
                    getseq = True
                elif line.startswith('//'):
                    getseq = False
                    o.close()
                elif getseq:
                    o.write(''.join(line.split()[1:]) + '\n')
    return length_dict


def gff_to_fasta(gff_file, out, concat=True):
    length_dict = {}
    if concat:
        with open(gff_file) as gff, open(out, 'w') as o:
            getseq = False
            for line in gff:
                if line.startswith('>'):
                    name = line.split()[0][1:]
                    o.write(line)
                    getseq = True
                    length_dict[name] = 0
                elif getseq:
                    o.write(line)
                    length_dict[name] += len(line.rstrip())
    else:
        with open(gff_file) as gff:
            getseq = False
            for line in gff:
                if line.startswith('>'):
                    name = line.split()[0][1:]
                    o = open(out + '.' + name + '.fa', 'w')
                    o.write('>' + name + '\n')
                    getseq = True
                    length_dict[name] = 0
                elif getseq:
                    o.write(line)
                    length_dict[name] += len(line.rstrip())
    return(length_dict)




def get_contig_matches(query_gbk, ref_gbk, working_dir):
    if query_gbk.endswith('.gff'):
        query_lengths = gff_to_fasta(query_gbk, working_dir + '/query_all.fa')
        ref_lengths = gff_to_fasta(ref_gbk, working_dir + '/ref_all.fa')
    else:
        query_lengths = gbk_to_fasta(query_gbk, working_dir + '/query_all.fa')
        ref_lengths = gbk_to_fasta(ref_gbk, working_dir + '/ref_all.fa')
    subprocess.Popen('nucmer ' + working_dir + '/query_all.fa ' + working_dir + '/ref_all.fa --prefix ' + working_dir + '/all_v_all', shell=True).wait()
    subprocess.Popen('delta-filter -g ' + working_dir + '/all_v_all.delta > ' +
                      working_dir + '/all_v_all.filter.delta', shell=True, stderr=subprocess.PIPE).wait()
    subprocess.Popen('show-coords ' + working_dir + '/all_v_all.filter.delta > ' + working_dir + '/all_v_all.coords', shell=True).wait()
    get_aligns = False
    matched_bases = {}
    with open(working_dir + '/all_v_all.coords') as f:
        for line in f:
            if line.startswith('=================='):
                get_aligns = True
            elif get_aligns:
                s1, e1, bar, s2, e2, bar, l1, l2, bar, idy, bar, query, ref = line.split()
                s1, e1 = int(s1), int(e1)
                if query in matched_bases:
                    if not ref in matched_bases[query]:
                        matched_bases[query][ref] = set()
                else:
                    matched_bases[query] = {ref:set()}
                for i in range(s1, e1+1):
                    matched_bases[query][ref].add(i)
    query_matches = {}
    ref_matches = {}
    for i in matched_bases:
        for j in matched_bases[i]:
            if i in query_matches and len(matched_bases[i][j]) > query_matches[i][1]:
                query_matches[i] = (j, len(matched_bases[i][j]))
            elif not i in query_matches:# and len(matched_bases[i][j]) > query_lengths[i] /5:
                query_matches[i] = (j, len(matched_bases[i][j]))
            if j in ref_matches and len(matched_bases[i][j]) > ref_matches[j][1]:
                ref_matches[j] = (i, len(matched_bases[i][j]))
            elif not j in ref_matches:# and len(matched_bases[i][j]) > ref_lengths[j] / 5:
                ref_matches[j] = (i, len(matched_bases[i][j]))
    matches = []
    for i in query_matches:
        match = query_matches[i][0]
        if match in ref_matches and ref_matches[match][0] == i:
            matches.append((i, match))
    return matches


def run_nucdiff(matches, working_dir, query_gbk, ref_gbk, nucdiff_path):
    if query_gbk.endswith('.gff'):
        gff_to_fasta(query_gbk, working_dir + '/' + os.path.splitext(os.path.basename(query_gbk))[0], False)
        gff_to_fasta(ref_gbk, working_dir + '/' + os.path.splitext(os.path.basename(ref_gbk))[0], False)
    else:
        gbk_to_fasta(query_gbk, working_dir + '/' + os.path.splitext(os.path.basename(query_gbk))[0], False)
        gbk_to_fasta(ref_gbk, working_dir + '/' + os.path.splitext(os.path.basename(ref_gbk))[0], False)
    gffs = []
    for i in matches:
        print(nucdiff_path + ' ' + working_dir + '/' + os.path.splitext(os.path.basename(ref_gbk))[0]
                         + '.' + i[1] + '.fa ' + working_dir + '/' + os.path.splitext(os.path.basename(query_gbk))[0] + '.' +
                         i[0] + '.fa ' + working_dir + ' ' + os.path.splitext(os.path.basename(query_gbk))[0] + '.' + i[0] +
                         'vs' + os.path.splitext(os.path.basename(ref_gbk))[0] + '.' + i[1])
        subprocess.Popen(nucdiff_path + ' ' + working_dir + '/' + os.path.splitext(os.path.basename(ref_gbk))[0]
                         + '.' + i[1] + '.fa ' + working_dir + '/' + os.path.splitext(os.path.basename(query_gbk))[0] + '.' +
                         i[0] + '.fa ' + working_dir + ' ' + os.path.splitext(os.path.basename(query_gbk))[0] + '.' + i[0] +
                         'vs' + os.path.splitext(os.path.basename(ref_gbk))[0] + '.' + i[1], shell=True).wait()
        gffs.append(os.path.join(working_dir, 'results', os.path.splitext(os.path.basename(query_gbk))[0] + '.' + i[0] +
                         'vs' + os.path.splitext(os.path.basename(ref_gbk))[0] + '.' + i[1]))
    return gffs






parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", help="Will create a gff of changes", required=True)
parser.add_argument("-qg", "--query_genbank", help="Concatenated genbank of genome", metavar="genome.gbk", required=True)
parser.add_argument("-rg", "--ref_genbank", help="Concatenated genbank of genome", metavar="genome.gbk", required=True)
parser.add_argument("-w", '--working_dir', help="Folder to put intermediary files.", required=True)
parser.add_argument("-r", '--reference', action="store_true", default=False, help="Look at changes to reference not query")
parser.add_argument("-n", '--nucdiff', default='nucdiff', help="path to nucdiff.py")
args = parser.parse_args()

if not os.path.exists(args.working_dir):
    os.makedirs(args.working_dir)

matches = get_contig_matches(args.query_genbank, args.ref_genbank, args.working_dir)
gffs = run_nucdiff(matches, args.working_dir, args.query_genbank, args.ref_genbank, args.nucdiff)
read_nucdiff(gffs, args.query_genbank, args.ref_genbank, args.output, args.working_dir)
