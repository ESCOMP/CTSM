import os
from datetime import date
import sys


def get_x(x, lines):
    strips=["'", '"']
    for line in lines:
        if x in line:
            xout = line.split('=')[1]
            for s in strips:
                if s in xout:
                    xout = xout.replace(s, '')
    return xout


def get_tdir(pwd):
    tools = '/tools/'
    if tools in pwd:
        base = pwd.split(tools)[0]
        tdir = base+tools+'tether'
    else:
        tdir = 'TETHER_DIR_NOT_FOUND'
    return tdir


def find_cime(pwd,user):
    ''' search upwards from pwd to find ./cime/scripts '''
    go = True
    thisdir = pwd
    while go:
        ld = os.listdir(thisdir)
        if 'cime' in ld:
            cime_dir = thisdir + '/cime/scripts'
            break
        else:
            last = os.path.basename(thisdir)
            if last == user:
                print('cannot find CIME')
                cime_dir = 'CIME/SCRIPTS_DIR_NOT_FOUND'
                break
            else:
                nx = len(last) + 1
                thisdir = thisdir[:-nx]
            if not os.path.isdir(thisdir):
                print('cannot find CIME')
                cime_dir = 'CIME/SCRIPTS_DIR_NOT_FOUND'
                break
    return cime_dir


def get_defaults(lines):
    user = get_x('USER', lines)
    pwd = get_x('SDIR', lines)
    tag = get_x('TAG',lines).replace('.','')
    tdir = get_tdir(pwd)
    dstr = 'c'+date.today().strftime('%Y%m%d')
    wdir = '/glade/u/home/'+user+'/tethered_sims/'+dstr
    cdir = find_cime(pwd,user)
    defaults = {'PI_COMPSET': 'I1850Clm60BgcCrujra',
                'GRID': 'f19_g17',
                'PROJECT': 'P93300041'}

    for v in ['PI_COMPSET','GRID']:
        x = get_x('PI_COMPSET',lines)
        if x!='UNSET':
            defaults[v]=x
    
    dirs = {'TDIR': tdir, 'WDIR': wdir, 'CDIR': cdir}
    for k in dirs:
        defaults[k] = dirs[k]

    compset = defaults['PI_COMPSET']
    grid = defaults['GRID']
    for suff in ['AD', 'SASU', 'ND']:
        defaults['CASE_' + suff] = '{}.{}.{}.{}'.format(tag, compset, grid, suff)

    return defaults


def write_new_file(lines, defaults):
    with open('config.tmp', 'w') as f:
        for line in lines:
            if '=' in line:
                k, v = line.split('=')
                if v == 'UNSET':
                    if k in defaults:
                        f.write(line.replace(v,defaults[k]) + '\n')
                    else:
                        f.write(line + '\n')
                else:
                    f.write(line + '\n')
            else:
                f.write(line + '\n') 


def main():
    with open('spinup.config', 'r') as f:
        lines = [line.rstrip() for line in f]
    defaults = get_defaults(lines)
    write_new_file(lines, defaults)


if __name__ == '__main__':
    main()
