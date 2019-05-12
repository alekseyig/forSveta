import subprocess, shlex

class Shell(object):
  def __init__(self, opts=None):
    self.opts = opts
    self._stdout = None
    self._stderr = None
    self._returncode = None

  def run(self, cmd, *args, **kwargs):
    cmd = cmd + ' ' + ' '.join(args)
    if self.opts.test:
      print cmd
      self._stdout, self._stderr = '', ''
      self._returncode = 0
    else:
      cmd_array = shlex.split(cmd)
      ## print cmd
      process = subprocess.Popen(cmd_array, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      self._stdout, self._stderr = process.communicate()
      self._returncode = process.poll()

    ### reverse the value of the return code to make it functional
    return 1 if self._returncode == 0 else 0
