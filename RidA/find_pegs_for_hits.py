#!/usr/bin/env python
from __future__ import with_statement
import sys, os
#from objc._convenience import sort

data_dir = "data"
blast_file = os.path.join("supportive_data", "PubSEED-ridA-Hits.txt")

class QueryHit(object):
  def __init__(self, name, org, contig, start, stop, ln_num, eval):
    # gi|14719469     993047.3        AEYZ01000114    10483   10803
    self.name = name
    self.org = org
    self.contig = contig
    self.start = int(start)
    self.stop = int(stop)
    self.ln_num = ln_num
    self.eval = float(eval)

  def compare(self, peg):
    self.found_hit = False
    if self.contig != peg.contig:
      return
    if (peg.direction == '+' and
      self.start <= peg.start and self.stop >= peg.start or
      self.start <= peg.stop and self.stop >= peg.stop or
      self.start >= peg.start and self.stop <= peg.stop or
      self.start <= peg.start and self.stop >= peg.stop):
      self.found_hit = True
    elif (peg.direction == '-' and
      self.start <= peg.stop and self.stop >= peg.stop or
      self.start <= peg.start and self.stop >= peg.start or
      self.start >= peg.stop and self.stop <= peg.start or
      self.start <= peg.stop and self.stop >= peg.start):
      self.found_hit = True
    elif (peg.direction == '+' and
      self.stop <= peg.start and self.start >= peg.start or
      self.stop <= peg.stop and self.start >= peg.stop or
      self.stop >= peg.start and self.start <= peg.stop or
      self.stop <= peg.start and self.start >= peg.stop):
      self.found_hit = True
    elif (peg.direction == '-' and
      self.start >= peg.stop and self.stop <= peg.stop or
      self.start >= peg.start and self.stop <= peg.start or
      self.start <= peg.start and self.stop >= peg.stop or
      self.start >= peg.start and self.stop <= peg.stop):
      self.found_hit = True

  @property
  def is_matched(self):
    #return  True if len(self.matched_info) else False
    return self.found_hit

  def __str__(self):
    return "QHit:\t%d\t%s\t%s\t%d\t%d\t%.2E" % (self.ln_num, self.name, self.contig, self.start, self.stop, self.eval)
    #return "QHit:\t%d\t%s" % (self.ln_num, self.name)
    #return "%d\t%s" % (self.ln_num, self.name)

  def __cmp__(self, other):
    return cmp(self.eval, other.eval)


class Peg(object):
  def __init__(self, name, org, contig, start, direction, length):
    # fig|993047.3.peg.917    993047.3:AEYZ01000047_36+1302
    self.name = name
    self.org = org
    self.contig = contig
    self.start = int(start)
    self.direction = direction
    self.length = int(length)

    if direction == '-':
      self.stop = self.start - self.length
    else:
      self.stop = self.start + self.length

  def __str__(self):
    #return "Peg: %s %s %d %d %s" % (self.name, self.contig, self.start, self.stop, self.direction)
    #return "Peg:\t%s\t%d" % (self.name, self.length)
    return "%s\t%d" % (self.name, self.length)


def process_pegs_and_qhits(pegs, fcontig, qhits):
  #print fcontig, qhits[fcontig]
  found_qhits = {}
  notfound_qhits = {}
  for_print = {}
  for qhit in qhits[fcontig]:
    for peg in pegs:
      qhit.compare(peg)
      if qhit.is_matched:
        found_qhits[qhit] = None
        for_print.setdefault(peg,[]).append(qhit)
        #sys.stdout.write("%s\t%s\n" % (qhit,peg))
      else:
        notfound_qhits[qhit] = None


  for peg, qhits_list in for_print.iteritems():
    qhits_list.sort()
    sys.stdout.write("%s\t%s\n" % (qhits_list[0],peg))


  for qhit in [qhit for qhit in notfound_qhits if qhit not in found_qhits]:
    sys.stdout.write("%s\t%s\n" % (qhit,"not_found"))


def process(fcontig, qhits):
  pegs = []
  get_pegs_from_fcontig(pegs, fcontig)
  #print len(pegs)
  process_pegs_and_qhits(pegs, fcontig, qhits)


def get_pegs_from_fcontig(pegs, fcontig):
  with open(os.path.join(data_dir, fcontig)) as f:
    for ln in f:
      name, location = ln.strip().split('\t')
      try:
        org, location = location.split(':')
      except ValueError:
        continue

      parts = location.split('_')
      start_length = parts[-1]
      if len(parts[:-1]) > 1:
        contig = '_'.join(parts[:-1])
      else:
        contig = parts[0]

      if '+' in start_length:
        direction = '+'
      else:
        direction = '-'

      start, length = start_length.split(direction)
      pegs.append(Peg(name,org,contig,start,direction,length))


def main(args):
  ln_num = 0
  qhits = {}
  with open(blast_file) as f:
    for ln in f:
      fields = ln.strip().split('\t')
      ln_num += 1
      try:
        if len(fields) == 13 and '.' in fields[1]:
          name, org, contig, hitstart, hitstop, eval = fields[0], fields[1], fields[2], fields[9], fields[10], fields[11]
          q = QueryHit(name,org,contig,hitstart,hitstop,ln_num,eval)
          qhits.setdefault(org,[]).append(q)
        elif len(fields) == 14 and '.' in fields[2]:
          name, org, contig, hitstart, hitstop, eval = fields[1], fields[2], fields[3], fields[10], fields[11], fields[12]
          q = QueryHit(name,org,contig,hitstart,hitstop,ln_num,eval)
          qhits.setdefault(org,[]).append(q)

      except IndexError:
        pass

  for _, _, fcontigs in os.walk(os.path.join(os.path.abspath("."), data_dir)):
    for fcontig in fcontigs:
      process(fcontig, qhits)
      #break

if __name__ == "__main__":
  sys.exit(main(sys.argv))
