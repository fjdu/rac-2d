def load_latex_tab(fname):
    """
    d: 
    """
    def make_one_line_dic(titles, data):
        return {titles[i]: data[i] for i in xrange(len(titles))}

    d = []
    with open(fname, 'r') as f:
        for line in f:
            a_line = [_.strip() for _ in line.split('&')]
            if not a_line[0].startswith('%'):
                d.append(a_line)
    sources = {}
    for i in xrange(1, len(d)):
        if d[i][0] == '':
            d[i][0] = d[i-1][0]
        else:
            sources.update({d[i][0]: make_one_line_dic(d[0], d[i])})
    #sources     = set([line[0] for line in d[1:]])
    transitions = set([line[1] for line in d[1:]])
    return d, sources, transitions



def load_txt_tab(fname):
    return open(fname, 'r').read().splitlines()
