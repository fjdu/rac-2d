class contri:
    def __init__(self, fname):
        with open(fname, 'r') as f:
            self.lines = f.readlines()
        self.itime = 0
        self.time = []
        self.keys = []
        s_time_row = 'Time = '
        for i, t in enumerate(self.lines):
            if self.lines[i].startswith(s_time_row):
                self.time.append((float(self.lines[i][len(s_time_row):]),
                                  i))
        for i in xrange(self.time[0][1]+1, self.time[1][1]-1):
            if self.lines[i][0] != ' ':
                self.keys.append(self.lines[i].split()[0])

    def time_match(self, time, i):
        if i == None:
            return False
        if i == (len(self.time)-1):
            return time >= self.time[-1][0]
        return (self.time[i][0] <= time and time < self.time[i+1][0]) or \
               (abs(self.time[i][0]-time) < time * 1e-3)

    def key_match(self, key, i):
        return self.lines[i].startswith(key)

    def locate_time(self, time):
        for i in xrange(len(self.time)):
            if self.time_match(time, i):
                return i

    def get_info(self, key, time):
        if not key in self.keys:
            return None
        if not self.time_match(time, self.itime):
            self.itime = self.locate_time(time)
        if self.itime == None:
            return None
        i_1 = self.time[self.itime][1]
        i0, i1, i2, i3 = None, None, None, None
        for i in xrange(i_1, len(self.lines)):
            if self.key_match(key, i):
                i0 = i
                break
        d = {}
        d['abundance'] = float(self.lines[i0].split()[-1])
        for i in xrange(i0+1, len(self.lines)):
            if self.lines[i].startswith('  Production'):
                i1 = i
                break
        d['prod_rate'] = float(self.lines[i1].split()[-1])
        for i in xrange(i1+1, len(self.lines)):
            if self.lines[i].startswith('  Destruction'):
                i2 = i
                break
        d['productions'] = (i1+1, i2)
        d['destr_rate'] = float(self.lines[i2].split()[-1])
        for i in xrange(i2+1, len(self.lines)):
            if self.lines[i][0] != ' ':
                i3 = i
                break
        d['destructions'] = (i2+1, i3)
        return d

    def parse_one_row():
        return
