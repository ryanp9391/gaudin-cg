# -*- coding: utf-8 -*-
r"""Build Gaudin-sl3-CG-notes.tex and .md from a light markup (src/*.src).
Markup:
  '# Title {#sec:label}', '## ', '### '   headings (numbered); '@@appendix' switches to lettered sections
  paragraphs (blank-line separated); inline $...$ math (LaTeX), **bold**, *em*, `code`, ^[footnote]
  '$$' ... '$$' display math (one number per block; \label{eq:x} allowed)
  '- ' / '1. ' lists (continuation lines indented by 2 spaces)
  '@table label | caption' / '@cols spec' / '@size cmd' / rows '| a | b |' / '@endtable'
  '@verb' ... '@endverb'
  '@include path'
  '@title', '@author', '@date', '@abstract' ... '@endabstract'
"""
import re, sys, os, glob

ROOT = os.path.dirname(os.path.abspath(__file__))

# ------------------------------------------------------------------ unicode math
GREEK = {
 'alpha':'α','beta':'β','gamma':'γ','delta':'δ','epsilon':'ε','varepsilon':'ε','zeta':'ζ','eta':'η','theta':'θ','vartheta':'ϑ',
 'iota':'ι','kappa':'κ','lambda':'λ','mu':'μ','nu':'ν','xi':'ξ','pi':'π','varpi':'ϖ','rho':'ρ','varrho':'ϱ','sigma':'σ','varsigma':'ς',
 'tau':'τ','upsilon':'υ','phi':'φ','varphi':'φ','chi':'χ','psi':'ψ','omega':'ω',
 'Gamma':'Γ','Delta':'Δ','Theta':'Θ','Lambda':'Λ','Xi':'Ξ','Pi':'Π','Sigma':'Σ','Upsilon':'Υ','Phi':'Φ','Psi':'Ψ','Omega':'Ω'}
SYM = {
 'sum':'Σ','prod':'Π','partial':'∂','infty':'∞','cdot':'·','times':'×','otimes':'⊗','oplus':'⊕','to':'→','rightarrow':'→','leftarrow':'←',
 'Rightarrow':'⇒','Leftarrow':'⇐','Leftrightarrow':'⇔','iff':'⇔','mapsto':'↦','leftrightarrow':'↔','le':'≤','leq':'≤','ge':'≥','geq':'≥',
 'ne':'≠','neq':'≠','pm':'±','mp':'∓','approx':'≈','equiv':'≡','sim':'∼','simeq':'≃','cong':'≅','propto':'∝','in':'∈','notin':'∉',
 'subset':'⊂','subseteq':'⊆','cap':'∩','cup':'∪','emptyset':'∅','forall':'∀','exists':'∃','ell':'ℓ','hbar':'ħ','dagger':'†','star':'⋆',
 'circ':'∘','bullet':'•','ldots':'…','dots':'…','cdots':'⋯','vdots':'⋮','ddots':'⋱','langle':'⟨','rangle':'⟩','prime':'′','nabla':'∇',
 'perp':'⊥','parallel':'∥','mid':'∣','vert':'∣','lvert':'∣','rvert':'∣','Vert':'‖','lVert':'‖','rVert':'‖','|':'‖','S':'§',
 'wedge':'∧','vee':'∨','neg':'¬','ll':'≪','gg':'≫','lfloor':'⌊','rfloor':'⌋','lceil':'⌈','rceil':'⌉','setminus':'∖','oplus_':'⊕',
 'sqrt':'√','int':'∫','oint':'∮','aleph':'ℵ','Re':'Re','Im':'Im','angle':'∠','triangle':'△','square':'□','checkmark':'✓',
 'leadsto':'⇝','uparrow':'↑','downarrow':'↓','longrightarrow':'⟶','longmapsto':'⟼','Longrightarrow':'⟹','hookrightarrow':'↪',
 'tfrac':'', 'dfrac':'','displaystyle':'','textstyle':'','nonumber':'','noindent':'',
}
FUNCS = {'det','lim','min','max','log','exp','sin','cos','tan','ln','rank','tr','deg','dim','gcd','sup','inf','arg','ker','Res','mod','id'}
SPACES = {',':' ',';':' ',':':' ','!':'','quad':'  ','qquad':'    ',' ':' ','>':' '}
SUBS = {**{str(i):'₀₁₂₃₄₅₆₇₈₉'[i] for i in range(10)}, '+':'₊','-':'₋','=':'₌','(':'₍',')':'₎',
 'a':'ₐ','e':'ₑ','h':'ₕ','i':'ᵢ','j':'ⱼ','k':'ₖ','l':'ₗ','m':'ₘ','n':'ₙ','o':'ₒ','p':'ₚ','r':'ᵣ','s':'ₛ','t':'ₜ','u':'ᵤ','v':'ᵥ','x':'ₓ',
 'β':'ᵦ','γ':'ᵧ','ρ':'ᵨ','φ':'ᵩ','χ':'ᵪ'}
SUPS = {**{str(i):'⁰¹²³⁴⁵⁶⁷⁸⁹'[i] for i in range(10)}, '+':'⁺','-':'⁻','=':'⁼','(':'⁽',')':'⁾',
 'a':'ᵃ','b':'ᵇ','c':'ᶜ','d':'ᵈ','e':'ᵉ','f':'ᶠ','g':'ᵍ','h':'ʰ','i':'ⁱ','j':'ʲ','k':'ᵏ','l':'ˡ','m':'ᵐ','n':'ⁿ','o':'ᵒ','p':'ᵖ','r':'ʳ','s':'ˢ',
 't':'ᵗ','u':'ᵘ','v':'ᵛ','w':'ʷ','x':'ˣ','y':'ʸ','z':'ᶻ','T':'ᵀ','β':'ᵝ','γ':'ᵞ','φ':'ᶲ','χ':'ᵡ','′':'′','*':'*','∗':'*'}
COMB = {'hat':'\u0302','tilde':'\u0303','bar':'\u0304','overline':'\u0304','widehat':'\u0302','widetilde':'\u0303','dot':'\u0307','ddot':'\u0308','vec':'\u20d7','check':'\u030c'}


def find_group(s, i):
    """s[i] == '{' -> return (content, index after closing)"""
    assert s[i] == '{'
    d = 0
    j = i
    while j < len(s):
        if s[j] == '{' and (j == 0 or s[j-1] != '\\'):
            d += 1
        elif s[j] == '}' and (j == 0 or s[j-1] != '\\'):
            d -= 1
            if d == 0:
                return s[i+1:j], j+1
        j += 1
    raise ValueError('unbalanced braces in: ' + s)


def read_arg(s, i):
    """read one TeX argument starting at i (skipping spaces): braces group or single token"""
    while i < len(s) and s[i] == ' ':
        i += 1
    if i >= len(s):
        return '', i
    if s[i] == '{':
        return find_group(s, i)
    if s[i] == '\\':
        m = re.match(r'\\([A-Za-z]+|.)', s[i:])
        return m.group(0), i + len(m.group(0))
    return s[i], i + 1


def mapscript(txt, table, marker):
    out = []
    for ch in txt:
        if ch in table:
            out.append(table[ch])
        else:
            if len(txt) > 16:
                mark_unclean()
            return marker + ('(' + txt + ')' if len(txt) > 1 else txt)
    return ''.join(out)


UNCLEAN = [False]


def mark_unclean():
    UNCLEAN[0] = True


def tex2uni_checked(s):
    UNCLEAN[0] = False
    r = tex2uni(s)
    return r, UNCLEAN[0]


def tex2uni(s):
    """convert a LaTeX math string to Unicode text"""
    out = []
    i = 0
    n = len(s)
    while i < n:
        c = s[i]
        if c == '\\':
            m = re.match(r'\\([A-Za-z]+)', s[i:])
            if m:
                name = m.group(1)
                i += len(m.group(0))
                if i < n and s[i] == ' ' and name not in ('quad', 'qquad'):
                    i += 1
                if name in GREEK:
                    out.append(GREEK[name])
                elif name in ('frac', 'dfrac', 'tfrac'):
                    a, i = read_arg(s, i); b, i = read_arg(s, i)
                    A, B = tex2uni(a), tex2uni(b)
                    wrapA = not re.fullmatch(r"[\w.′'̂]+", A)
                    wrapB = not re.fullmatch(r"[\w.′'̂]+", B)
                    VULG = {('1', '2'): '½', ('1', '3'): '⅓', ('2', '3'): '⅔', ('1', '4'): '¼', ('3', '4'): '¾', ('1', '5'): '⅕',
                            ('1', '6'): '⅙', ('1', '8'): '⅛', ('1', '9'): '⅑'}
                    if (A, B) in VULG:
                        res = VULG[(A, B)]
                    else:
                        res = ('(' + A + ')' if wrapA else A) + '/' + ('(' + B + ')' if wrapB else B)
                    if i < n and (s[i].isalnum() or s[i] == chr(92)):
                        res += ' '
                    out.append(res)
                elif name == 'binom':
                    a, i = read_arg(s, i); b, i = read_arg(s, i)
                    out.append('C(' + tex2uni(a) + ',' + tex2uni(b) + ')')
                elif name == 'sqrt':
                    if i < n and s[i] == '[':
                        j = s.index(']', i); root = s[i+1:j]; i = j + 1
                        a, i = read_arg(s, i)
                        out.append(tex2uni(root) + '√(' + tex2uni(a) + ')')
                    else:
                        a, i = read_arg(s, i)
                        A = tex2uni(a)
                        out.append('√' + (A if re.fullmatch(r'[\w.]+', A) else '(' + A + ')'))
                elif name in ('mathrm', 'text', 'textrm', 'mathbf', 'boldsymbol', 'mathit', 'operatorname', 'mathsf', 'textbf', 'textit', 'mathcal', 'mathfrak', 'mathscr', 'mathbb', 'emph', 'texttt', 'mbox', 'bm', 'mathop', 'mathtt'):
                    a, i = read_arg(s, i)
                    if name == 'mathbb':
                        mp = {'C': 'ℂ', 'Q': 'ℚ', 'R': 'ℝ', 'Z': 'ℤ', 'N': 'ℕ'}
                        out.append(''.join(mp.get(ch, ch) for ch in a))
                    else:
                        out.append(tex2uni(a) if name not in ('text', 'textrm', 'mbox', 'textbf', 'textit', 'emph', 'texttt') else a)
                elif name in COMB:
                    a, i = read_arg(s, i)
                    A = tex2uni(a)
                    if len(A) >= 1:
                        out.append(A[0] + COMB[name] + A[1:] if len(A) == 1 else '(' + A + ')' + COMB[name])
                elif name == 'begin':
                    env, i = read_arg(s, i)
                    endtag = '\\end{' + env + '}'
                    j = s.index(endtag, i)
                    body = s[i:j]; i = j + len(endtag)
                    if env in ('pmatrix', 'bmatrix', 'matrix', 'vmatrix'):
                        rows = [r for r in re.split(r'\\\\', body) if r.strip()]
                        rr = ['[' + ', '.join(tex2uni(x.strip()) for x in r.split('&')) + ']' for r in rows]
                        out.append('[' + ', '.join(rr) + ']' if len(rr) > 1 else rr[0])
                    elif env in ('cases',):
                        rows = [r for r in re.split(r'\\\\', body) if r.strip()]
                        out.append('{ ' + ' ; '.join(' '.join(tex2uni(x.strip()) for x in r.split('&')) for r in rows) + ' }')
                    else:
                        out.append(tex2uni(body.replace('&', ' ')))
                elif name in ('left', 'right', 'big', 'Big', 'bigg', 'Bigg', 'bigl', 'bigr', 'Bigl', 'Bigr', 'biggl', 'biggr', 'limits', 'nolimits', 'displaystyle', 'textstyle', 'nonumber', 'noindent', 'protect', 'smash', 'phantom', 'vphantom', 'hphantom', 'label', 'rm', 'allowbreak'):
                    if name == 'label':
                        a, i = read_arg(s, i)
                    elif name in ('phantom', 'vphantom', 'hphantom', 'smash'):
                        a, i = read_arg(s, i)
                    # strip; keep following delimiter ('.' dropped)
                    if name in ('left', 'right', 'big', 'Big', 'bigg', 'Bigg', 'bigl', 'bigr', 'Bigl', 'Bigr', 'biggl', 'biggr') and i < n and s[i] == '.':
                        i += 1
                elif name in ('quad', 'qquad'):
                    out.append('  ')
                elif name in SYM:
                    out.append(SYM[name])
                elif name in FUNCS:
                    out.append(name)
                elif name in ('lbrace',):
                    out.append('{')
                elif name in ('rbrace',):
                    out.append('}')
                elif name == 'ref' or name == 'eqref':
                    a, i = read_arg(s, i)
                    out.append(REFS.get(a, '?'))
                elif name == 'underline':
                    a, i = read_arg(s, i); out.append(tex2uni(a))
                elif name in ('colon',):
                    out.append(':')
                else:
                    mark_unclean()
                    out.append(name)  # unknown macro: print name
                continue
            # \ followed by a symbol char
            if i + 1 < n:
                d = s[i+1]
                i += 2
                if d in SPACES:
                    out.append(SPACES[d])
                elif d == '\\':
                    out.append('; ')
                elif d in '{}%&_#$':
                    out.append(d)
                elif d == '|':
                    out.append('‖')
                else:
                    out.append(d)
                continue
            i += 1
            continue
        if c in '_^':
            kind = c
            i += 1
            a, i = read_arg(s, i)
            A = tex2uni(a)
            if re.match(r'\\(mathrm|text|rm)\b', a.strip()):
                out.append(kind + A.strip())
            elif kind == '_':
                out.append(mapscript(A, SUBS, '_'))
            else:
                out.append(mapscript(A, SUPS, '^'))
            continue
        if c == '{':
            a, i = find_group(s, i)
            out.append(tex2uni(a))
            continue
        if c == '}':
            i += 1
            continue
        if c == "'":
            out.append('′'); i += 1; continue
        if c == '~':
            out.append(' '); i += 1; continue
        if c == '&':
            out.append(' '); i += 1; continue
        out.append(c)
        i += 1
    r = ''.join(out)
    r = r.replace('′′', '″')
    return r


# ------------------------------------------------------------------ parsing helpers
REFS = {}


def read_src():
    parts = sorted(glob.glob(os.path.join(ROOT, 'src', '*.src')))
    text = []
    for p in parts:
        text.append(open(p, encoding='utf8').read())
    full = '\n'.join(text)

    def inc(m):
        path = os.path.join(ROOT, m.group(1).strip())
        return open(path, encoding='utf8').read()
    for _ in range(3):
        full = re.sub(r'^@include (.+)$', inc, full, flags=re.M)
    return full


class Block:
    def __init__(self, kind, **kw):
        self.kind = kind
        self.__dict__.update(kw)


def split_cells(line):
    """split a '| a | b |' row on '|' outside $...$ math and not preceded by a backslash"""
    line = line.strip()
    assert line.startswith('|')
    cells = []
    cur = ''
    inmath = False
    i = 1
    while i < len(line):
        c = line[i]
        if c == '\\' and i + 1 < len(line):
            cur += line[i:i+2]
            i += 2
            continue
        if c == '$':
            inmath = not inmath
        if c == '|' and not inmath:
            cells.append(cur.strip())
            cur = ''
        else:
            cur += c
        i += 1
    if cur.strip():
        cells.append(cur.strip())
    return cells


def parse(full):
    lines = full.split('\n')
    blocks = []
    i = 0
    n = len(lines)
    meta = {}
    while i < n:
        ln = lines[i]
        if not ln.strip() or ln.startswith('%%'):
            i += 1; continue
        if ln.startswith('@title '):
            meta['title'] = ln[7:].strip(); i += 1; continue
        if ln.startswith('@author '):
            meta['author'] = ln[8:].strip(); i += 1; continue
        if ln.startswith('@date '):
            meta['date'] = ln[6:].strip(); i += 1; continue
        if ln.startswith('@abstract'):
            j = i + 1; buf = []
            while not lines[j].startswith('@endabstract'):
                buf.append(lines[j]); j += 1
            meta['abstract'] = '\n'.join(buf).strip()
            i = j + 1; continue
        if ln.startswith('@startsec'):
            blocks.append(Block('startsec', n=int(ln.split()[1])))
            i += 1
            continue
        if ln.startswith('@@appendix'):
            blocks.append(Block('appendix')); i += 1; continue
        m = re.match(r'^(#{1,3}) (.*?)(?: \{#([\w:.-]+)\})?\s*$', ln)
        if m:
            blocks.append(Block('head', level=len(m.group(1)), text=m.group(2), label=m.group(3))); i += 1; continue
        if ln.strip() == '$$':
            j = i + 1; buf = []
            while lines[j].strip() != '$$':
                buf.append(lines[j]); j += 1
            blocks.append(Block('eq', body='\n'.join(buf))); i = j + 1; continue
        if ln.startswith('@verb'):
            j = i + 1; buf = []
            while not lines[j].startswith('@endverb'):
                buf.append(lines[j]); j += 1
            blocks.append(Block('verb', body='\n'.join(buf))); i = j + 1; continue
        if ln.startswith('@table'):
            m = re.match(r'@table\s+([\w:.-]+)\s*\|\s*(.*)$', ln)
            label, cap = m.group(1), m.group(2)
            cols = None; size = 'small'; rows = []; j = i + 1
            while not lines[j].startswith('@endtable'):
                l2 = lines[j]
                if l2.startswith('@cols '):
                    cols = l2[6:].strip()
                elif l2.startswith('@size '):
                    size = l2[6:].strip()
                elif l2.strip().startswith('|'):
                    cells = split_cells(l2.strip())
                    if all(re.fullmatch(r':?-{2,}:?', c) for c in cells if c):
                        pass
                    else:
                        rows.append(cells)
                j += 1
            blocks.append(Block('table', label=label, caption=cap, cols=cols, size=size, rows=rows)); i = j + 1; continue
        if re.match(r'^(- |\d+\. )', ln):
            items = []; ordered = bool(re.match(r'^\d+\. ', ln))
            while i < n and re.match(r'^(- |\d+\. )', lines[i]):
                buf = [re.sub(r'^(- |\d+\. )', '', lines[i])]; i += 1
                while i < n and lines[i].startswith('  ') and lines[i].strip():
                    buf.append(lines[i].strip()); i += 1
                items.append(' '.join(buf))
            blocks.append(Block('list', items=items, ordered=ordered)); continue
        # paragraph
        buf = [ln]; i += 1
        while i < n and lines[i].strip() and not re.match(r'^(#{1,3} |\$\$|@|- |\d+\. )', lines[i]):
            buf.append(lines[i]); i += 1
        blocks.append(Block('para', text=' '.join(x.strip() for x in buf)))
    return meta, blocks


# ------------------------------------------------------------------ label pass
def number_labels(blocks):
    sec = [0, 0, 0]; appendix = False
    eq = 0; tab = 0
    for b in blocks:
        if b.kind == 'startsec':
            sec = [b.n - 1, 0, 0]
        elif b.kind == 'appendix':
            appendix = True; sec = [0, 0, 0]
        elif b.kind == 'head':
            L = b.level
            sec[L-1] += 1
            for k in range(L, 3):
                sec[k] = 0
            if appendix:
                num = chr(64 + sec[0]) + ''.join('.' + str(x) for x in sec[1:L])
            else:
                num = '.'.join(str(x) for x in sec[:L])
            b.num = num
            if b.label:
                REFS[b.label] = num
        elif b.kind == 'eq':
            eq += 1; b.num = eq
            for lab in re.findall(r'\\label\{([^}]+)\}', b.body):
                REFS[lab] = str(eq)
        elif b.kind == 'table':
            tab += 1; b.num = tab
            REFS[b.label] = str(tab)


# ------------------------------------------------------------------ TeX emission
def tex_escape_text(t):
    t = re.sub(r'(?<!\\)&', r'\\&', t)
    t = re.sub(r'(?<!\\)%', r'\\%', t)
    t = re.sub(r'(?<!\\)#', r'\\#', t)
    t = re.sub(r'(?<!\\)_', r'\\_', t)
    t = re.sub(r'(?<!\\)>=', r'$\\ge$', t)
    return t


def split_inline(text):
    """split into list of (kind, content): 'm' math, 'c' code, 't' text, 'f' footnote"""
    res = []
    i = 0; n = len(text); buf = ''
    while i < n:
        c = text[i]
        if c == '\\' and i + 1 < n and text[i+1] == '$':
            buf += '\\$'; i += 2; continue
        if c == '$':
            j = i + 1
            while j < n and text[j] != '$':
                if text[j] == '\\':
                    j += 1
                j += 1
            if buf: res.append(('t', buf)); buf = ''
            res.append(('m', text[i+1:j])); i = j + 1; continue
        if c == '`':
            j = text.index('`', i + 1)
            if buf: res.append(('t', buf)); buf = ''
            res.append(('c', text[i+1:j])); i = j + 1; continue
        if c == '^' and i + 1 < n and text[i+1] == '[':
            d = 0; j = i + 1
            while j < n:
                if text[j] == '[': d += 1
                elif text[j] == ']':
                    d -= 1
                    if d == 0: break
                elif text[j] == '$':
                    # skip over math to avoid brackets inside math
                    k = j + 1
                    while text[k] != '$':
                        if text[k] == '\\': k += 1
                        k += 1
                    j = k
                j += 1
            if buf: res.append(('t', buf)); buf = ''
            res.append(('f', text[i+2:j])); i = j + 1; continue
        buf += c; i += 1
    if buf: res.append(('t', buf))
    return res


def style_text_tex(t):
    t = tex_escape_text(t)
    t = re.sub(r'\*\*(.+?)\*\*', r'\\textbf{\1}', t)
    t = re.sub(r'(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])', r'\\emph{\1}', t)
    t = t.replace('\\"', '\x01')
    t = re.sub(r'"(.*?)"', r"``\1''", t)
    t = t.replace('\x01', '\\"')
    t = re.sub(r'(?<=[a-z0-9\)\}])([.!?])\s+(?=[A-Z\\`])', r'\1\n', t)
    return t


def protect(text):
    """replace math/code/footnote segments by placeholders; return (buf, holders)"""
    holders = []
    buf = ''
    for kind, content in split_inline(text):
        if kind == 't':
            buf += content
        else:
            holders.append((kind, content))
            buf += '\x00%d\x00' % (len(holders) - 1)
    return buf, holders


def inline_tex(text):
    buf, holders = protect(text)
    buf = style_text_tex(buf)

    def sub(m):
        kind, content = holders[int(m.group(1))]
        if kind == 'm':
            return '$' + content + '$'
        if kind == 'c':
            esc = re.sub(r'([_&%#$])', r'\\\1', content).replace('{', '\\{').replace('}', '\\}')
            esc = esc.replace('/', '/\\allowbreak{}').replace('\\_', '\\_\\allowbreak{}').replace('-', '-\\allowbreak{}').replace('.', '.\\allowbreak{}')
            return '\\texttt{' + esc + '}'
        return '\\footnote{' + inline_tex(content) + '}'
    return re.sub('\x00(\\d+)\x00', sub, buf)


def heading_tex(text):
    out = ''
    for kind, content in split_inline(text):
        if kind == 'm':
            plain = re.sub(r'[^\x00-\x7f]', '', tex2uni(content)).replace('\\', '')
            plain = re.sub(r'[_^{}]', '', plain)
            out += '\\texorpdfstring{$' + content + '$}{' + plain + '}'
        elif kind == 'c':
            out += '\\texttt{' + content + '}'
        else:
            out += style_text_tex(content)
    return out


def emit_tex(meta, blocks):
    o = []
    o.append(r'''% Generated from Paul/Mathematica/Gaudin/notes-src/ on 2026-10-16; may be edited directly thereafter (the generator is optional).
\documentclass[11pt]{article}
\usepackage{lmodern}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{amsmath,amssymb}
\usepackage{booktabs,longtable,array}
\usepackage[margin=1in]{geometry}
\usepackage{xcolor}
\usepackage[colorlinks=true,linkcolor=blue!50!black,urlcolor=blue!50!black,citecolor=blue!50!black]{hyperref}
\usepackage{microtype}
\setlength{\parskip}{0.4em}
\setlength{\parindent}{0pt}
\setlength{\LTcapwidth}{\textwidth}
\allowdisplaybreaks
\emergencystretch=2em
\sloppy
''')
    o.append('\\title{' + inline_tex(meta['title']) + '}')
    o.append('\\author{' + inline_tex(meta.get('author', '')) + '}')
    o.append('\\date{' + inline_tex(meta.get('date', '')) + '}')
    o.append('\\begin{document}\n\\maketitle')
    if 'abstract' in meta:
        o.append('\\begin{abstract}\n' + '\n\n'.join(inline_tex(p) for p in meta['abstract'].split('\n\n')) + '\n\\end{abstract}')
    o.append('\\tableofcontents\n\\bigskip')
    for b in blocks:
        if b.kind == 'startsec':
            o.append('\\setcounter{section}{%d}' % (b.n - 1))
        elif b.kind == 'appendix':
            o.append('\\appendix')
        elif b.kind == 'head':
            cmd = ['section', 'subsection', 'subsubsection'][b.level-1]
            s = '\\' + cmd + '{' + heading_tex(b.text) + '}'
            if b.label:
                s += '\\label{' + b.label + '}'
            o.append(s)
        elif b.kind == 'para':
            o.append(inline_tex(b.text))
        elif b.kind == 'list':
            env = 'enumerate' if b.ordered else 'itemize'
            o.append('\\begin{' + env + '}\n' + '\n'.join('\\item ' + inline_tex(x) for x in b.items) + '\n\\end{' + env + '}')
        elif b.kind == 'eq':
            body = b.body.strip()
            if re.match(r'\\begin\{(align|gather|multline)\*?\}', body):
                o.append(body)
            else:
                o.append('\\begin{equation}\n' + body + '\n\\end{equation}')
        elif b.kind == 'verb':
            o.append('{\\footnotesize\n\\begin{verbatim}\n' + b.body + '\n\\end{verbatim}}')
        elif b.kind == 'table':
            ncol = max(len(r) for r in b.rows)
            cols = b.cols or ('l' * ncol)
            hdr = b.rows[0]
            s = '{\\' + b.size + '\n\\begin{longtable}{' + cols + '}\n'
            s += '\\caption{' + inline_tex(b.caption) + '}\\label{' + b.label + '}\\\\\n\\toprule\n'
            s += ' & '.join(inline_tex(c) for c in hdr) + ' \\\\\n\\midrule\n\\endfirsthead\n'
            s += '\\toprule\n' + ' & '.join(inline_tex(c) for c in hdr) + ' \\\\\n\\midrule\n\\endhead\n\\bottomrule\n\\endlastfoot\n'
            for r in b.rows[1:]:
                r = r + [''] * (ncol - len(r))
                s += ' & '.join(inline_tex(c) for c in r) + ' \\\\\n'
            s += '\\end{longtable}}\n'
            o.append(s)
    o.append('\\end{document}\n')
    return '\n\n'.join(o)


# ------------------------------------------------------------------ Markdown emission
FOOTNOTES = []


def style_text_md(t):
    t = t.replace('\\"u', 'ü').replace('\\"o', 'ö').replace("\\'e", 'é').replace('\\&', '&').replace('\\%', '%').replace('\\#', '#').replace('\\_', '_').replace('\\$', '$')
    t = re.sub(r'\\(?:emph|textit)\{([^{}]*)\}', r'*\1*', t)
    t = re.sub(r'\\textbf\{([^{}]*)\}', r'**\1**', t)
    t = re.sub(r'\\(?:mbox|text)\{([^{}]*)\}', r'\1', t)
    t = t.replace('---', '—').replace('--', '–').replace('~', ' ')
    t = re.sub(r'\\S\s*\\ref\{([^}]+)\}', lambda m: '§' + REFS.get(m.group(1), '?'), t)
    t = re.sub(r'\\S\s*', '§', t)
    t = re.sub(r'\\(?:eq)?ref\{([^}]+)\}', lambda m: REFS.get(m.group(1), '?'), t)
    t = re.sub(r'\\cite\{[^}]*\}', '', t)
    t = re.sub(r'"(.*?)"', r'“\1”', t)
    return t


def inline_md(text, in_table=False):
    buf, holders = protect(text)
    buf = style_text_md(buf)

    def sub(m):
        kind, content = holders[int(m.group(1))]
        if kind == 'm':
            u, bad = tex2uni_checked(content)
            if bad:
                return '$' + content.replace('|', '\\vert ' if in_table else '|') + '$'
            if in_table:
                u = u.replace('|', '∣')
            return u
        if kind == 'c':
            return '`' + content + '`'
        FOOTNOTES.append(inline_md(content))
        return '[^' + str(len(FOOTNOTES)) + ']'
    return re.sub('\x00(\\d+)\x00', sub, buf)


def split_top_rows(body):
    rows = []
    depth = 0
    cur = ''
    i = 0
    while i < len(body):
        if body.startswith('\\begin{', i):
            depth += 1
        if body.startswith('\\end{', i):
            depth -= 1
        if depth == 0 and body.startswith('\\\\', i):
            rows.append(cur)
            cur = ''
            i += 2
            continue
        cur += body[i]
        i += 1
    rows.append(cur)
    return rows


def emit_md(meta, blocks):
    o = []
    o.append('# ' + inline_md(meta['title']))
    o.append('*' + inline_md(meta.get('author', '')) + (' — ' if meta.get('date') else '') + inline_md(meta.get('date', '')) + '*')
    o.append('*(Markdown twin of `Gaudin-sl3-CG-notes.tex`; same structure and numbering. Math is Unicode; displayed equations are numbered as in the PDF.)*')
    if 'abstract' in meta:
        o.append('## Abstract\n\n' + '\n\n'.join(inline_md(p) for p in meta['abstract'].split('\n\n')))
    toc = ['## Contents', '']
    for b in blocks:
        if b.kind == 'head' and b.level <= 2:
            toc.append(('  ' * (b.level - 1)) + '- ' + b.num + '. ' + inline_md(b.text))
    o.append('\n'.join(toc))
    app = False
    for b in blocks:
        if b.kind == 'appendix':
            app = True
        elif b.kind == 'head':
            hashes = '#' * (b.level + 1)
            prefix = ('Appendix ' + b.num) if (app and b.level == 1) else b.num
            o.append(hashes + ' ' + prefix + '. ' + inline_md(b.text) if not (app and b.level == 1) else hashes + ' ' + prefix + '. ' + inline_md(b.text))
        elif b.kind == 'para':
            o.append(inline_md(b.text))
        elif b.kind == 'list':
            o.append('\n'.join((('%d. ' % (k+1)) if b.ordered else '- ') + inline_md(x) for k, x in enumerate(b.items)))
        elif b.kind == 'eq':
            body = re.sub(r'\\label\{[^}]*\}', '', b.body).strip()
            body = re.sub(r'\\begin\{(aligned|align\*?|gathered|split)\}', '', body)
            body = re.sub(r'\\end\{(aligned|align\*?|gathered|split)\}', '', body)
            rows = [r.strip() for r in split_top_rows(body) if r.strip()]
            conv = [tex2uni_checked(r) for r in rows]
            urows = [c[0].strip() for c in conv]
            bad = any(c[1] for c in conv)
            if bad:
                raw = re.sub(r'\\label\{[^}]*\}', '', b.body).strip()
                o.append('```latex\n' + raw + '\n```\n(' + str(b.num) + ')')
                continue
            urows[0] = urows[0] + '        (' + str(b.num) + ')' if len(urows) == 1 else urows[0]
            if len(urows) > 1:
                urows[-1] = urows[-1] + '        (' + str(b.num) + ')'
            o.append('```text\n' + '\n'.join(urows) + '\n```')
        elif b.kind == 'verb':
            o.append('```\n' + b.body + '\n```')
        elif b.kind == 'table':
            ncol = max(len(r) for r in b.rows)
            o.append('**Table ' + str(b.num) + '.** ' + inline_md(b.caption))
            hdr = [inline_md(c, True) for c in b.rows[0]] + [''] * (ncol - len(b.rows[0]))
            s = '| ' + ' | '.join(hdr) + ' |\n|' + '|'.join(['---'] * ncol) + '|\n'
            for r in b.rows[1:]:
                r = r + [''] * (ncol - len(r))
                s += '| ' + ' | '.join(inline_md(c, True) for c in r) + ' |\n'
            o.append(s.rstrip('\n'))
    txt = '\n\n'.join(o)
    if FOOTNOTES:
        txt += '\n\n---\n\n' + '\n\n'.join('[^%d]: %s' % (k + 1, f) for k, f in enumerate(FOOTNOTES))
    return txt + '\n'


def main():
    full = read_src()
    meta, blocks = parse(full)
    number_labels(blocks)
    outdir = sys.argv[1]
    tex = emit_tex(meta, blocks)
    open(os.path.join(outdir, 'Gaudin-sl3-CG-notes.tex'), 'w', encoding='utf8').write(tex)
    md = emit_md(meta, blocks)
    open(os.path.join(outdir, 'Gaudin-sl3-CG-notes.md'), 'w', encoding='utf8').write(md)
    print('blocks', len(blocks), 'eq', sum(1 for b in blocks if b.kind == 'eq'), 'tables', sum(1 for b in blocks if b.kind == 'table'),
          'footnotes', len(FOOTNOTES))
    # unresolved refs
    for m in re.findall(r'\\(?:eq)?ref\{([^}]+)\}', full):
        if m not in REFS:
            print('UNRESOLVED REF', m)


if __name__ == '__main__':
    main()
