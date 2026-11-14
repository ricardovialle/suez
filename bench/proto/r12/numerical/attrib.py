# attrib.py <sample.txt> <leaf substring>: attribute self samples of a leaf to the nearest suez-engine ancestor
import re, sys, collections
OURS = ['Profile::operator()', 'grid_values', 'Spectrum::weights', 'Spectrum::wsum', 'Spectrum::mom', 'Block::Block', 'x0_grid', 'rotate_all', 'rotate(', 'make_model', 'newton', 'refine', 'wald', 'fit_gene', 'x0_moments', 'FastLDLT', 'grid_direct']
lines = open(sys.argv[1]).read().split('Sort by top of stack')[0].split('\n')
stack = []; res = collections.Counter()
for ln in lines:
    m = re.match(r'^([ +!:|]*)(\d+) (.*)$', ln)
    if not m: continue
    depth = len(m.group(1)); cnt = int(m.group(2)); name = m.group(3)
    while stack and stack[-1][0] >= depth: stack.pop()
    stack.append((depth, name, cnt))
    if sys.argv[2] in name:
        anc = next((o for d, nm, c in reversed(stack[:-1]) for o in OURS if o in nm), 'other')
        # self count = cnt minus children: approximate by cnt when the next line is not deeper (leaf); keep simple
        res[anc] += cnt
for k, v in res.most_common(12): print(v, k)
