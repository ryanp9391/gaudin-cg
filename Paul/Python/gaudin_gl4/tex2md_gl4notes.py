#!/usr/bin/env python3
"""tex2md_gl4notes.py -- small converter for gl4-slice-engine-notes.tex (restricted LaTeX subset) to its Markdown twin (written by G6, 2026-10-22).
Usage (from Paul/Python/gaudin_gl4):
    python -X utf8 tex2md_gl4notes.py ../../Mathematica/Gaudin/gl4-slice-engine-notes.tex ../../Mathematica/Gaudin/gl4-slice-engine-notes.md
Regenerates the .md byte for byte from the .tex; it understands only the LaTeX constructs used in that file (see below), so run it again after editing the .tex.
Handles: title/author/date, abstract, sections/subsections/appendices with numbering, equation/align/aligned, \\[ \\] displays,
tables (tabular -> pipe tables), itemize/enumerate, footnotes (collected at the end), \\cite/\\ref/\\eqref, the file-name macros, accents, dashes.
Math is kept verbatim as $...$ / $$...$$.
"""
import re
import sys

src = open(sys.argv[1], encoding='utf-8').read()

# ---------------------------------------------------------------- helpers
def strip_comments(s):
    out = []
    for line in s.split('\n'):
        if line.lstrip().startswith('%'):
            continue
        # inline comments: an unescaped % starts a comment
        m = re.search(r'(?<!\\)%', line)
        if m:
            line = line[:m.start()]
        out.append(line)
    return '\n'.join(out)


def match_brace(s, i):
    """s[i] == '{' ; return index of matching '}'"""
    assert s[i] == '{', s[i:i + 20]
    depth = 0
    j = i
    while j < len(s):
        c = s[j]
        if c == '\\':
            j += 2
            continue
        if c == '{':
            depth += 1
        elif c == '}':
            depth -= 1
            if depth == 0:
                return j
        j += 1
    raise ValueError('unbalanced braces near: ' + s[i:i + 60])


src = strip_comments(src)

# ---------------------------------------------------------------- front matter
def grab(cmd, s):
    k = s.find('\\' + cmd + '{')
    if k < 0:
        return ''
    i = s.index('{', k)
    j = match_brace(s, i)
    return s[i + 1:j]

title = grab('title', src)
title = re.sub(r'\\vspace\{[^}]*\}', '', title).strip()
author = grab('author', src)
date = grab('date', src)

body = src[src.index('\\begin{document}') + len('\\begin{document}'):src.index('\\end{document}')]
body = body.replace('\\maketitle', '')

# ---------------------------------------------------------------- bibliography
bib_m = re.search(r'\\begin\{thebibliography\}\{[^}]*\}(.*?)\\end\{thebibliography\}', body, re.S)
bib_text = bib_m.group(1) if bib_m else ''
bib_keys = re.findall(r'\\bibitem\{([^}]*)\}', bib_text)
cite_num = {k: i + 1 for i, k in enumerate(bib_keys)}
# remove the bibliography from body; remember what comes before and after it (after: trailing note)
pre_bib = body[:bib_m.start()] if bib_m else body
post_bib = body[bib_m.end():] if bib_m else ''
# the bibliography may be wrapped in a font-size group "{\footnotesize ... }": strip a dangling group
pre_bib = re.sub(r'\{\\(footnotesize|small|scriptsize)\s*$', '', pre_bib.rstrip())
post_bib = post_bib.strip()
if post_bib.startswith('}'):
    post_bib = post_bib[1:].strip()

body = pre_bib

# ---------------------------------------------------------------- first pass: numbering of sections, tables, equations
labels = {}

sec_counter = 0
sub_counter = 0
app_mode = False
tab_counter = 0
eq_counter = 0

def letter(n):
    return chr(ord('A') + n - 1)

# scan sequentially through body for \appendix, \section, \subsection, tables and equations
token_re = re.compile(r'\\appendix|\\section\{|\\subsection\{|\\begin\{table\}|\\begin\{equation\}|\\begin\{align\}|\\label\{([^}]*)\}')
pos = 0
cur_env = None
last_target = None  # what the next \label attaches to
for m in token_re.finditer(body):
    t = m.group(0)
    if t == '\\appendix':
        app_mode = True
        sec_counter = 0
    elif t == '\\section{':
        sec_counter += 1
        sub_counter = 0
        last_target = ('sec', letter(sec_counter) if app_mode else str(sec_counter))
    elif t == '\\subsection{':
        sub_counter += 1
        last_target = ('sub', (letter(sec_counter) if app_mode else str(sec_counter)) + '.' + str(sub_counter))
    elif t == '\\begin{table}':
        tab_counter += 1
        last_target = ('tab', str(tab_counter))
    elif t == '\\begin{equation}':
        eq_counter += 1
        last_target = ('eq', str(eq_counter))
    elif t == '\\begin{align}':
        last_target = ('align', None)
    elif t.startswith('\\label{'):
        key = m.group(1)
        if last_target is None:
            continue
        kind, val = last_target
        if kind == 'align':
            eq_counter += 1
            labels[key] = str(eq_counter)
        else:
            labels[key] = val
# the number of \label inside an align without a label would be missed; our align lines all carry labels.

# ---------------------------------------------------------------- inline conversion
ACCENTS = {
    "\\v{r}": 'ř', "\\'e": 'é', "\\'{e}": 'é', '\\"o': 'ö', '\\"{o}': 'ö',
}

def conv_code(s):
    s = s.replace('\\_', '_').replace('\\{', '{').replace('\\}', '}').replace('\\%', '%').replace('\\dots', '...').replace('\\&', '&')
    s = s.replace('$\\lambda$', 'λ').replace('~', ' ').replace('\\allowbreak ', '').replace('\\allowbreak', '')
    return s

PREFIX = {'gd': 'GD/', 'gl': 'GL/', 'cpp': 'CPP/', 'py': 'PY3/', 'src': ''}

footnotes = []
MATHS = []

def parse_arg(s, i):
    """s[i] must be '{' ; return (arg, next index)"""
    j = match_brace(s, i)
    return s[i + 1:j], j + 1


def inline(s):
    """convert a text fragment (may contain math) to markdown"""
    # protect math
    def keep(m):
        MATHS.append(m.group(0))
        return '\x00M%d\x00' % (len(MATHS) - 1)
    s = re.sub(r'(?<!\\)\$[^$]*\$', keep, s)
    out = []
    i = 0
    n = len(s)
    while i < n:
        c = s[i]
        if c == '\\':
            # command
            m = re.match(r'\\([A-Za-z]+)', s[i:])
            if m:
                cmd = m.group(1)
                i += m.end()
                # skip one trailing space? (TeX eats spaces after control words; \S{} etc.)
                def next_arg():
                    nonlocal i
                    k = i
                    while k < n and s[k] == ' ':
                        k += 1
                    if k < n and s[k] == '{':
                        a, i = parse_arg(s, k)
                        return a
                    return None
                if cmd in ('textbf',):
                    a = next_arg(); out.append('**' + inline(a) + '**')
                elif cmd == 'emph':
                    a = next_arg(); out.append('*' + inline(a) + '*')
                elif cmd == 'texttt':
                    a = next_arg(); out.append('`' + conv_code(a) + '`')
                elif cmd == 'textsf':
                    a = next_arg(); out.append(inline(a))
                elif cmd in PREFIX:
                    a = next_arg(); out.append('`' + PREFIX[cmd] + conv_code(a) + '`')
                elif cmd == 'path':
                    a = next_arg(); out.append('`' + conv_code(a) + '`')
                elif cmd == 'rg':
                    a = next_arg(); b = next_arg(); out.append('G%s §%s' % (a, b))
                elif cmd == 'cite':
                    a = next_arg(); out.append('[' + ', '.join(str(cite_num[k.strip()]) for k in a.split(',')) + ']')
                elif cmd == 'ref':
                    a = next_arg(); out.append(labels[a.strip()])
                elif cmd == 'eqref':
                    a = next_arg(); out.append('(' + labels[a.strip()] + ')')
                elif cmd == 'label':
                    a = next_arg()
                elif cmd == 'footnote':
                    a = next_arg()
                    footnotes.append(None)
                    idx = len(footnotes)
                    footnotes[idx - 1] = inline(a.strip())
                    out.append('[^%d]' % idx)
                elif cmd == 'texorpdfstring':
                    a = next_arg(); b = next_arg(); out.append(inline(a))
                elif cmd == 'S':
                    out.append('§')
                    # swallow an empty group "{}" after \S
                    if s[i:i + 2] == '{}':
                        i += 2
                elif cmd in ('dots', 'ldots'):
                    out.append('…')
                elif cmd == 'vspace':
                    next_arg()
                elif cmd in ('medskip', 'bigskip', 'smallskip', 'noindent', 'centering', 'small', 'footnotesize', 'scriptsize', 'normalsize', 'par', 'raggedright'):
                    pass
                elif cmd == 'v':
                    a = next_arg(); out.append({'r': 'ř'}.get(a, a))
                else:
                    out.append('\\' + cmd)  # unknown: keep
                continue
            # control symbol
            nxt = s[i + 1] if i + 1 < n else ''
            if nxt in "'\"":
                # accent like \'e or \"o  (possibly braced)
                if s[i + 2:i + 3] == '{':
                    a, k = parse_arg(s, i + 2)
                    ch = a
                    i = k
                else:
                    ch = s[i + 2:i + 3]
                    i += 3
                out.append({("'", 'e'): 'é', ('"', 'o'): 'ö'}.get((nxt, ch), ch))
                continue
            if nxt == '%':
                out.append('%'); i += 2; continue
            if nxt == '&':
                out.append('&'); i += 2; continue
            if nxt == '#':
                out.append('#'); i += 2; continue
            if nxt == ' ':
                out.append(' '); i += 2; continue
            if nxt == ',':
                out.append(' '); i += 2; continue
            if nxt == '\\':
                out.append('\n'); i += 2; continue
            if nxt in '_{}':
                out.append(nxt); i += 2; continue
            out.append(c)
            i += 1
            continue
        if c == '~':
            out.append(' ')
            i += 1
            continue
        if c == '`' and s[i:i + 2] == '``':
            out.append('“'); i += 2; continue
        if c == "'" and s[i:i + 2] == "''":
            out.append('”'); i += 2; continue
        if c == '-' and s[i:i + 3] == '---':
            out.append('—'); i += 3; continue
        if c == '-' and s[i:i + 2] == '--':
            out.append('–'); i += 2; continue
        if c in '{}':
            i += 1  # stray grouping braces
            continue
        out.append(c)
        i += 1
    res = ''.join(out)
    res = re.sub(r'\x00M(\d+)\x00', lambda m: MATHS[int(m.group(1))], res)
    return res

# ---------------------------------------------------------------- block conversion
eq_out_counter = [0]

def fmt_display(math, tag=None):
    math = math.strip()
    if tag:
        return '$$\n%s \\tag{%s}\n$$' % (math, tag)
    return '$$\n%s\n$$' % math


def eq_block(env, text):
    """equation / align environment -> list of display blocks"""
    labs = re.findall(r'\\label\{([^}]*)\}', text)
    if env == 'equation':
        inner = re.sub(r'\\label\{[^}]*\}', '', text).strip()
        lab = labs[0] if labs else None
        return fmt_display(inner, labels[lab] if lab else None)
    # align: split lines on \\ at top level (not inside aligned etc: ours are simple)
    parts = re.split(r'\\\\(?![\[])', text)
    blocks = []
    for p in parts:
        if not p.strip():
            continue
        lab = re.findall(r'\\label\{([^}]*)\}', p)
        q = re.sub(r'\\label\{[^}]*\}', '', p).replace('&', '').strip()
        blocks.append(fmt_display(q, labels[lab[0]] if lab else None))
    return '\n\n'.join(blocks)


def colspec_align(spec):
    out = []
    i = 0
    while i < len(spec):
        c = spec[i]
        if c in '@><!':
            j = match_brace(spec, i + 1)
            i = j + 1
            continue
        if c == 'p':
            j = match_brace(spec, i + 1)
            out.append('l')
            i = j + 1
            continue
        if c in 'lcr':
            out.append(c)
        i += 1
    return out


def tabular_to_md(spec, text):
    aligns = colspec_align(spec)
    rows = []
    header_rows = 1
    seen_mid = False
    # split lines by \\
    text = text.replace('\\toprule', '\n').replace('\\bottomrule', '\n')
    chunks = re.split(r'\\midrule', text)
    head_text = chunks[0]
    rest_text = '\n'.join(chunks[1:])
    def rows_of(t):
        rs = [r for r in re.split(r'\\\\', t) if r.strip()]
        out = []
        for r in rs:
            cells = [c.strip() for c in re.split(r'(?<!\\)&', r)]
            out.append(cells)
        return out
    hr = rows_of(head_text)
    br = rows_of(rest_text)
    ncol = len(aligns)
    def fix(cells):
        cells = cells + [''] * (ncol - len(cells))
        return [inline(c).replace('\n', ' ').replace('|', '\\|') if '$' not in c else inline(c).replace('\n', ' ') for c in cells[:ncol]]
    lines = []
    if not hr:
        hr = [[''] * ncol]
    lines.append('| ' + ' | '.join(fix(hr[0])) + ' |')
    sep = []
    for a in aligns:
        sep.append({'l': ':---', 'r': '---:', 'c': ':---:'}[a])
    lines.append('| ' + ' | '.join(sep) + ' |')
    for r in br:
        lines.append('| ' + ' | '.join(fix(r)) + ' |')
    return '\n'.join(lines)


def table_block(text):
    cap_m = re.search(r'\\caption\{', text)
    cap = ''
    if cap_m:
        i = cap_m.end() - 1
        j = match_brace(text, i)
        cap = text[i + 1:j]
    num = None
    lab = re.search(r'\\label\{([^}]*)\}', text)
    if lab:
        num = labels[lab.group(1)]
    tabs = []
    for m in re.finditer(r'\\begin\{tabular\}\{', text):
        i = m.end() - 1
        j = match_brace(text, i)
        spec = text[i + 1:j]
        end = text.index('\\end{tabular}', j)
        tabs.append(tabular_to_md(spec, text[j + 1:end]))
    out = '**Table %s.** %s\n\n' % (num, inline(cap).strip()) + '\n\n'.join(tabs)
    return out


def list_block(env, text):
    items = re.split(r'\\item\b', text)[1:]
    out = []
    for k, it in enumerate(items, 1):
        it = it.strip()
        # inside an item: possible display equation
        blocks = convert(it, indent=True)
        prefix = ('%d. ' % k) if env == 'enumerate' else '- '
        pad = ' ' * len(prefix)
        lines = blocks.split('\n')
        first = True
        res = []
        for ln in lines:
            if first:
                res.append(prefix + ln)
                first = False
            else:
                res.append((pad + ln) if ln.strip() else '')
        out.append('\n'.join(res))
    return '\n'.join(out)


env_re = re.compile(r'\\begin\{(equation|align|table|itemize|enumerate|abstract)\}')


def convert(text, indent=False):
    out = []
    pos = 0
    while True:
        m = env_re.search(text, pos)
        seg = text[pos:m.start()] if m else text[pos:]
        out.append(paragraphs(seg))
        if not m:
            break
        env = m.group(1)
        endtag = '\\end{%s}' % env
        e = text.index(endtag, m.end())
        inner = text[m.end():e]
        if env in ('equation', 'align'):
            out.append(eq_block(env, inner))
        elif env == 'table':
            inner = re.sub(r'^\[[^\]]*\]', '', inner)
            out.append(table_block(inner))
        elif env in ('itemize', 'enumerate'):
            out.append(list_block(env, inner))
        elif env == 'abstract':
            out.append('## Abstract\n\n' + paragraphs(inner))
        pos = e + len(endtag)
    return '\n\n'.join(b for b in out if b.strip())


def paragraphs(seg):
    """text between environments: handle \\[...\\], \\section, \\subsection, \\appendix, paragraphs"""
    seg = seg.strip('\n')
    if not seg.strip():
        return ''
    # display math \[ ... \]
    parts = re.split(r'(\\\[.*?\\\])', seg, flags=re.S)
    res = []
    for p in parts:
        if p.startswith('\\[') and p.endswith('\\]'):
            res.append(fmt_display(p[2:-2]))
        elif p.strip():
            res.append(para_text(p))
    return '\n\n'.join(res)


sec_state = {'sec': 0, 'sub': 0, 'app': False}


def para_text(p):
    out = []
    # split into paragraphs on blank lines
    for para in re.split(r'\n\s*\n', p.strip('\n')):
        if not para.strip():
            continue
        lines = para.split('\n')
        buf = []
        for ln in lines:
            s = ln.strip()
            if s.startswith('\\appendix'):
                sec_state['app'] = True
                sec_state['sec'] = 0
                continue
            ms = re.match(r'\\(section|subsection)\{', s)
            if ms:
                i = s.index('{')
                j = match_brace(s, i)
                title_ = s[i + 1:j]
                rest = s[j + 1:]
                if ms.group(1) == 'section':
                    sec_state['sec'] += 1
                    sec_state['sub'] = 0
                    num = letter(sec_state['sec']) if sec_state['app'] else str(sec_state['sec'])
                    buf.append('## %s. %s' % (num, inline(title_)))
                else:
                    sec_state['sub'] += 1
                    num = (letter(sec_state['sec']) if sec_state['app'] else str(sec_state['sec'])) + '.' + str(sec_state['sub'])
                    buf.append('### %s. %s' % (num, inline(title_)))
                # the \label follows on the same line; ignore
                continue
            buf.append(s)
        text = '\n'.join(buf)
        # headings must be separated by blank lines
        text = re.sub(r'(^|\n)(#{2,3} [^\n]*)\n', r'\1\n\2\n\n', text)
        # inline conversion line by line is not safe for footnotes spanning lines; convert whole paragraph
        parts = re.split(r'(\n#{2,3} [^\n]*\n)', text)
        for q in parts:
            if re.match(r'\n#{2,3} ', q):
                out.append(q.strip())
            elif q.strip():
                if q.strip().startswith('#'):
                    out.append(q.strip())
                else:
                    out.append(inline(q.strip()))
    return '\n\n'.join(o for o in out if o.strip())


md_body = convert(body)

# bibliography
bib_out = ['## References', '']
for k, m in enumerate(re.finditer(r'\\bibitem\{([^}]*)\}(.*?)(?=\\bibitem\{|$)', bib_text, re.S), 1):
    bib_out.append('[%d] %s' % (k, inline(m.group(2).strip().replace('\n', ' '))))
    bib_out.append('')
post = inline(post_bib.replace('\\end{document}', '').strip()) if post_bib.strip() else ''

header = []
header.append('# ' + inline(title))
header.append('')
auth = inline(author.replace('\\\\', ' ').replace('\\small', '')).strip()
header.append('*%s — %s*' % (re.sub(r'\s+', ' ', auth), inline(date)))
header.append('')
header.append('*(Markdown twin of `gl4-slice-engine-notes.tex`, generated from it by a small script: same sections, same numbering of equations, tables and footnotes. Math is written in `$…$` / `$$…$$`, tables are pipe tables, footnotes are collected at the end.)*')
header.append('')

doc = '\n'.join(header) + '\n' + md_body.strip() + '\n\n' + '\n'.join(bib_out).strip() + '\n'
if post:
    doc += '\n' + post + '\n'
doc += '\n---\n\n'
for k, f in enumerate(footnotes, 1):
    doc += '[^%d]: %s\n\n' % (k, re.sub(r'\s*\n\s*', ' ', f))

open(sys.argv[2], 'w', encoding='utf-8', newline='\n').write(doc)
print('written', sys.argv[2], len(doc), 'chars;', len(footnotes), 'footnotes;', tab_counter, 'tables;', eq_counter, 'equations')
