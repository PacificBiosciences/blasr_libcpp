#!/usr/bin/env python
"""Configure the build.

- Fetch HDF5 headers.
- Create libconfig.h
- Create defines.mk

This is not used by './unittest/'.
"""
import commands
import contextlib
import os
import sys

def log(msg):
    sys.stderr.write(msg)
    sys.stderr.write('\n')

def shell(cmd):
    log(cmd)
    status, output = commands.getstatusoutput(cmd)
    if status:
        raise Exception('%d <- %r' %(status, cmd))
    return output

def update_content(fn, content):
    current_content = open(fn).read() if os.path.exists(fn) else None
    if content != current_content:
        log('writing to %r' %fn)
        log('"""\n' + content + '"""')
        open(fn, 'w').write(content)

def compose_libconfig(pbbam=False):
    if pbbam:
        content = """
#define USE_PBBAM
"""
    else:
        content = """
"""
    return content

def compose_defines_with_hdf_headers(HDF_HEADERS):
    thisdir = os.path.dirname(os.path.abspath(__file__))
    return """
HDF_HEADERS:=%(HDF_HEADERS)s
#HDF5_INCLUDE?=${HDF_HEADERS}/src
CPPFLAGS+=-I${HDF_HEADERS}/src -I${HDF_HEADERS}/c++/src
CPPFLAGS+=-I../pbdata -I../hdf -I../alignment
LIBPBDATA_LIB     ?=../pbdata/libpbdata.so
LIBPBIHDF_LIB     ?=../pbdata/libpbihdf.so
LIBBLASR_LIB      ?=../pbdata/libblasr.so
"""%(dict(thisdir=thisdir, HDF_HEADERS=HDF_HEADERS))

def compose_defines():
    """
    Note that our local 'hdf' subdir will not even build
    in this case.
    """
    thisdir = os.path.dirname(os.path.abspath(__file__))
    return """
LIBPBDATA_INCLUDE ?=../pbdata
LIBPBIHDF_INCLUDE ?=../hdf
LIBBLASR_INCLUDE  ?=../alignment
LIBPBDATA_LIB     ?=%(thisdir)s/pbdata/libpbdata.so
LIBPBIHDF_LIB     ?=%(thisdir)s/pbdata/libpbihdf.so
LIBBLASR_LIB      ?=%(thisdir)s/pbdata/libblasr.so
nohdf             ?=1
"""%(dict(thisdir=thisdir))

def get_OS_STRING():
    G_BUILDOS_CMD = """bash -c 'set -e; set -o pipefail; id=$(lsb_release -si | tr "[:upper:]" "[:lower:]"); rel=$(lsb_release -sr); case $id in ubuntu) printf "$id-%04d\n" ${rel/./};; centos) echo "$id-${rel%%.*}";; *) echo "$id-$rel";; esac' 2>/dev/null"""
    return shell(G_BUILDOS_CMD)
def get_PREBUILT():
    cmd = 'cd ../../../../prebuilt.out 2>/dev/null && pwd || echo -n notfound'
    return shell(cmd)
def get_BOOST_INCLUDE(env):
    key_bi = 'BOOST_INCLUDE'
    key_br = 'BOOST_ROOT'
    if key_bi in env:
        return env[key_bi]
    if key_br in env:
        return env[key_br]
    return '${PREBUILT}/boost/boost_1_55_0'

def get_PBBAM(env, prefix):
    """
    key = 'PBBAM'
    if key in env:
        return env[key]
    cmd = 'cd $(THIRD_PARTY_PREFIX)/../staging/PostPrimary/pbbam 2>/dev/null && pwd || echo -n notfound' %(
            THIRD_PARTY_PREFIX=prefix)
    return shell(cmd)
    """
def get_HTSLIB(env, prefix):
    """
    key = 'HTSLIB'
    if key in env:
        return env[key]
    cmd = 'cd $(THIRD_PARTY_PREFIX)/../staging/PostPrimary/htslib 2>/dev/null && pwd || echo -n notfound' %(
            THIRD_PARTY_PREFIX=prefix)
    return shell(cmd)
    """
def ifenvf(env, key, func):
    if key in env:
        return env[key]
    else:
        return func()
def setifenvf(envout, envin, key, func):
    envout[key] = ifenvf(envin, key, func)
def setifenv(envout, envin, key, val):
    envout[key] = envin.get(key, val)
def setenv(envout, key, val):
    envout[key] = val
def update_env_if(envout, envin, keys):
    for key in keys:
        if key in envin:
            envout[key] = envin[key]
def compose_defs_env(env):
    # We disallow env overrides for anything with a default from GNU make.
    nons = ['CXX', 'CC', 'AR'] # 'SHELL'?
    ovr    = ['%-20s ?= %s' %(k, v) for k,v in env.items() if k not in nons]
    nonovr = ['%-20s := %s' %(k, v) for k,v in env.items() if k in nons]
    return '\n'.join(ovr + nonovr + [''])
def compose_defines_pacbio(envin):
    """
    This is used by mobs via buildcntl.sh.
    """
    env = dict()
    setenv(env, 'SHELL', 'bash')
    setifenvf(env, envin, 'OS_STRING', get_OS_STRING)
    setifenvf(env, envin, 'PREBUILT', get_PREBUILT)
    setifenv(env, envin, 'LIBPBDATA_INCLUDE', '../pbdata')
    setifenv(env, envin, 'LIBPBIHDF_INCLUDE', '../hdf')
    setifenv(env, envin, 'LIBBLASR_INCLUDE', '../alignment')
    setifenv(env, envin, 'LIBPBDATA_LIB', '../pbdata/libpbdata.so')
    setifenv(env, envin, 'LIBPBIHDF_LIB', '../hdf/libpbihdf.so')
    setifenv(env, envin, 'LIBBLASR_LIB', '../alignment/libblasr.so')
    setifenv(env, envin, 'nohdf', '1')
    nondefaults = set([
            'CXX', 'AR',
            'HDF5_INCLUDE', 'HDF5_LIB',
            'PBBAM_INCLUDE', 'PBBAM_LIB',
            'HTSLIB_INCLUDE', 'HTSLIB_LIB',
            'BOOST_INCLUDE',
            'ZLIB_LIB',
    ])
    update_env_if(env, envin, nondefaults)
    return compose_defs_env(env)


@contextlib.contextmanager
def cd(nwd):
    cwd = os.getcwd()
    log('cd %r -> %r' %(cwd, nwd))
    os.chdir(nwd)
    yield
    os.chdir(cwd)
    log('cd %r <- %r' %(cwd, nwd))

def fetch_hdf5_headers():
    """Fetch into ./hdf/HEADERS directory.
    Return actual directory path, relative to subdirs.
    """
    version = 'hdf5-1.8.12-headers'
    if not os.path.isdir(os.path.join('hdf', version)):
        with cd('hdf'):
            cmd = 'curl -k -L https://www.dropbox.com/s/8971bcyy5o42rxb/hdf5-1.8.12-headers.tar.bz2\?dl\=0 | tar xjf -'
            shell(cmd)
    return os.path.join('../hdf', version) # Relative path might help caching.

def update(content_defines_mk, content_libconfig_h):
    """ Write these relative to the same directory as *this* file.
    """
    thisdir = os.path.dirname(os.path.abspath(__file__))
    fn_defines_mk = os.path.join(thisdir, 'defines.mk')
    update_content(fn_defines_mk, content_defines_mk)
    fn_libconfig_h = os.path.join(thisdir, 'pbdata', 'libconfig.h')
    update_content(fn_libconfig_h, content_libconfig_h)

def configure_nopbbam():
    HDF_HEADERS = fetch_hdf5_headers()
    content1 = compose_defines_with_hdf_headers(HDF_HEADERS)
    content2 = compose_libconfig(pbbam=False)
    update(content1, content2)

def configure_nopbbam_nohdf5():
    content1 = compose_defines()
    content2 = compose_libconfig(pbbam=False)
    update(content1, content2)

def configure_pacbio(envin):
    content1 = compose_defines_pacbio(envin)
    content2 = compose_libconfig(pbbam=True)
    update(content1, content2)

def get_make_style_env(envin, args):
    envout = dict()
    for arg in args:
        if '=' in arg:
            k, v = arg.split('=')
            envout[k] = v
    envout.update(envin)
    return envout

def main(prog, *args):
    """We are still deciding what env-vars to use, if any.
    """
    if 'NOPBBAM' in os.environ:
        if 'NOHDF' in os.environ:
            configure_nopbbam_nohdf5()
        else:
            configure_nopbbam()
    else:
        envin = get_make_style_env(os.environ, args)
        configure_pacbio(envin)


if __name__=="__main__":
    main(*sys.argv)
