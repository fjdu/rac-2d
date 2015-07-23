wi = 12
pos = [i*wi for i in xrange(6,2,-1)]

def format_this(s, pos, wi, s_repl, None_if_no_space=True):
    fmtstr = '{0:' + '{0:d}'.format(wi) + 's}'
    for i in xrange(len(pos)-1):
        p = pos[i]
        p_1 = pos[i+1]
        s_this = s[p:(p+wi)].strip()
        s_1    = s[p_1:(p_1+wi)].strip()
        if len(s_this) == 0 and len(s_1) > 0:
            return s[:p] + fmtstr.format(s_repl) + s[p+wi:]
    if None_if_no_space:
        print 'No space: ', s
        return None
    else:
        return s

def format_Grain(fname, grain_name):
    fname_out = 'out_' + fname
    with open(fname, 'r') as f:
        lines = f.readlines()
    with open(fname_out, 'w') as f:
        for s in lines:
            sout = format_this(s, pos, wi, grain_name)
            if sout != None:
                f.write(sout)

format_Grain('tmp1.dat', 'Grain0')
format_Grain('tmp2.dat', 'Grain+')
