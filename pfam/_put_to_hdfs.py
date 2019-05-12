#!/usr/bin/env python2.5

import os, sys, subprocess, glob
from copy import copy
from optparse import OptionParser
from download_compaings import get_date_range, get_env_var

__doc__ = """Program 'put's Test and Target Omniture campaings data to hdfs.
You can specify the date range (default for yesterday, right side of the range is not included)
for which dates put data into hdfs."""


class HadoopFS(object):
  def __init__(self, uri=None):
    self.returncode = None
    self.cmd = ['hadoop','fs']
    if uri is not None: self.cmd.extend(['-fs', uri])
    
  def __execute(self, *args, **kwargs):
    cmd = copy(self.cmd)
    cmd.extend(args)
    ## cmd.append(' '.join(['%s=%s' % (k,v) for k,v in kwargs.items()]))
    print ">>>", ' '.join(cmd)
    self.returncode = subprocess.call(cmd)
    ### reverse the value of the return code to make it functional
    return self.returncode == 0 and 1 or 0
  
  def put(self, src, dst):
    if isinstance(src,list):
      src.append(dst)
      return self.__execute('-put', *src)
    return self.__execute('-put', src, dst)
    
  def mkdir(self, dir):
    if not self.test(dir, '-e'):
      return self.__execute('-mkdir', dir)
    else:
      return 1

  def test(self, fname, *args):
    return self.__execute('-test', ' '.join(args), fname)
    

def process(opts, args):
  hadoop_fs = HadoopFS(opts.uri)
  dst_root = get_env_var('TT_HADOOP_DST')
  src_root = get_env_var('TT_DOWNLOAD')

  date_range = get_date_range(args)
  for dt in date_range:
    src = glob.glob(os.path.abspath(os.path.join(src_root, "%.2d" % dt.year, "%.2d" % dt.month, "%.2d" % dt.day, '*')))
    dst = os.path.join(dst_root, "%.2d" % dt.year, "%.2d" % dt.month, "%.2d" % dt.day)
    if len(src):
      if not hadoop_fs.mkdir(dst): sys.exit(hadoop_fs.returncode)
      if not hadoop_fs.put(src, dst): sys.exit(hadoop_fs.returncode)
  return 0

def main (args):
  usage = "%prog [-h| --help] [-u, --uri URI] [yyyy-MM-DD--yyyy-MM-DD]"
  parser = OptionParser(usage=usage, description=__doc__)
  parser.add_option("-u", "--uri", dest="uri", default=None, help="Hadoop URI, default is configured in hadoop conf file.")
  (options, args) = parser.parse_args()
  return process(options, args)

if __name__ == "__main__":
  sys.exit(main(sys.argv))
