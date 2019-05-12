#!/usr/bin/env python2.5

import os, sys, urllib2
from optparse import OptionParser
from datetime import date, timedelta
from xml.dom.minidom import parse, parseString

__doc__ = """Program downloads Test and Target campaings data from 
Omniture. You can specify the campaing id(s) (default all campaings) 
and date range (default for yesterday), right side of the range is not included ."""

BASE_PROTOCOL = 'https://'
BASE_URL = 'testandtarget.omniture.com/admin/api'
EMAIL = 'webanalytics@orbitzworldwide.com'

def get_date_range(args):
  ### determine the date range 
  date_range = []
  if args and '--' in args[-1]:
    start_range, end_range = args[-1].split('--')
    y,m,d = map(lambda x: int(x), start_range.split('-'))
    s = date(y,m,d)
    y,m,d = map(lambda x: int(x), end_range.split('-'))
    e = date(y,m,d)
    if s > e:
      raise AttributeError('The end day in the date range cannot be less than the start date. Are trying to be too smart?')
    one_day = timedelta(days=1)
    while s < e:
      date_range.append(s)
      s += one_day
  else:    
    yesterday = date.today() - timedelta(days=1)
    date_range.append(yesterday)  
  return date_range

def get_env_var(var):
    value = os.environ.get(var, None)
    if value is None:
      raise EnvironmentError('%s environment variable is not set.' % var)
    return value


class OmnitureService(object):
  def __init__(self):
    self.urllib2 = urllib2
    self.omniture_site = None
    self.client_name = get_env_var('TT_CLIENT')
    self.client_password = get_env_var('TT_CLIENT_PASSWORD')
    proxy_handler = self.urllib2.ProxyHandler({'https': get_env_var('http_proxy')})
    opener = self.urllib2.build_opener(proxy_handler)
    self.urllib2.install_opener(opener)

  def __make_base_url(self, site=''):
    url = BASE_PROTOCOL
    if site: url = url + site + '.'
    url = '%s%s?client=%s&email=%s&password=%s' % (url, BASE_URL, self.client_name, EMAIL, self.client_password)
    return url

  def __find_omniture_site(self):
    url = self.__make_base_url()
    for ln in self.urllib2.urlopen(url).readlines():
      if 'form action' in ln:
        self.omniture_site = ln.split('https://')[1].split('.')[0]
        return
    raise RuntimeError("Could not find the 'admin' omniture site to connect to.")

  def __execute_operation(self, **kwargs):
    if self.omniture_site is None:
      self.__find_omniture_site()

    url = self.__make_base_url(self.omniture_site)
    url = url + ''.join(['&%s=%s' % (k,v) for k,v in kwargs.items()])
    #print ">>>>", url
    return self.urllib2.urlopen(url)

  def __get_text(self, nodelist):
    result = []
    for node in nodelist:
      if node.nodeType == node.TEXT_NODE:
        result.append(node.data)
    return ''.join(result)

  def __get_campaings(self, **kwargs):
    operation = 'campaignList'
    return self.__execute_operation(operation = operation, **kwargs)

  def get_campaing_info(self, **kwargs):
    operation = 'viewCampaign'
    return self.__execute_operation(operation = operation, **kwargs)

  def get_campaings_list(self, **kwargs):
    result = self.__get_campaings(**kwargs)
    list = []
    dom = parse(result)
    for camp_element in dom.getElementsByTagName('campaign'):
      id = self.__get_text(camp_element.getElementsByTagName('id')[0].childNodes)
      ##name = self.__get_text(camp_element.getElementsByTagName('name')[0].childNodes)
      #campaign_info = self.get_campaing_info(id=id)
      #campaign_dom = parse(campaign_info)
      #start_date = self.__get_text(campaign_dom.getElementsByTagName('startDate')[0].childNodes)
      #end_date = self.__get_text(campaign_dom.getElementsByTagName('endDate')[0].childNodes)
      #print id, start_date, end_date
      #if kwargs.hasattr('before') and start_date > kwargs.get('before'): continue
      #if kwargs.hasattr('after') and end_date < kwargs.get('after'): continue
      list.append(id)
    return list

  def get_campaing_report(self, **kwargs):
    operation = 'report'
    result = self.__execute_operation(operation = operation, **kwargs)
    return result

  def get_campaing_audit_report(self, **kwargs):
    operation = 'auditReport'
    result = self.__execute_operation(operation = operation, **kwargs)
    return result


class SaveService(object):
  def __init__(self):
    self.storage_root = get_env_var('TT_DOWNLOAD')

  def get_dst_file(self, dt, fname):
    name = os.path.join(self.storage_root, dt)
    if not os.path.exists(name):
      os.makedirs(name)
    dst_file = os.path.abspath(os.path.join(name, str(fname)))
    print ">>>", dst_file
    return dst_file

  def write_dom(self, dom, dt, fname):
    dst = self.get_dst_file(dt, fname)
    dstf = open(dst, 'w')
    dom.writexml(dstf)
    dstf.close()

  def write_data(self, data, dt, fname):
    dst = self.get_dst_file(dt, fname)
    dstf = open(dst, 'w')
    for line in data.readlines():
      dstf.write(line)
    dstf.close()


def process(opts, args):
  save_service = SaveService()
  om_service = OmnitureService()

  campaings_list = []
  if len(args):
    for arg in args:
      cids = [c for c in [c.strip(',') for c in arg.split(',')] if c.isdigit()]
      if len(cids):
        campaings_list.extend(cids)
      else:
        campaings_list = om_service.get_campaings_list()
  else:
    campaings_list = om_service.get_campaings_list()

  date_range = get_date_range(args)

  for d in date_range:
    dt = d.strftime("%Y/%m/%d")
    start = '%sT00:00' % d.strftime("%Y-%m-%d")
    end = '%sT23:59' % d.strftime("%Y-%m-%d")                      

    for campaing_id in campaings_list:
      dom = parse(om_service.get_campaing_report(campaignId=campaing_id, start=start, end=end))

      ### checks if report is not empty and if it is not, fetch and write different reports
      elements = dom.getElementsByTagName('sample')
      if elements:
        save_service.write_dom(dom, dt, '%s_report.xml' % campaing_id)
          
        ### write campaing details
        dom = parse(om_service.get_campaing_info(id=campaing_id))
        save_service.write_dom(dom, dt, '%s_viewCampaign.xml' % campaing_id)
          
        ### write campaing audit report
        audit_report_data = om_service.get_campaing_audit_report(campaignId=campaing_id, start=start, end=end, format='csv')
        save_service.write_data(audit_report_data, dt, '%s_auditReport.csv' % campaing_id)

                    
def main (args):
  usage = "%prog [-h| --help] [cid1,cid2,...] [yyyy-MM-DD--yyyy-MM-DD]"
  parser = OptionParser(usage=usage, description=__doc__)
  #parser.add_option("-a", "--all", dest="all", default=False, action="store_true", help="Download data for all, including deactivated, campaigs.")
  (options, args) = parser.parse_args()
  process(options, args)
  return 0

if __name__ == "__main__":
  sys.exit(main(sys.argv))

