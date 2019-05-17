import os
import sys
import requests
import shutil

uname, password, input_file, output_dir = sys.argv[1:]


login_data = {"page": 'Login', "action": 'perform_login',
              "login": uname,
              "password": password}

url ='http://rast.theseed.org/FIG/rast.cgi'

cookies = requests.post(url, data=login_data).cookies
cookies = {'WebSession': cookies.get('WebSession')}

if not os.path.exists(output_dir):
    os.mkdir(output_dir)

with open(input_file) as f:
    for ln in f:
        _job, _file = ln.strip().split('\t')
        # _file = _job + '_' + _file + '.gbk'
        _file = _file + '.gbk'
        print _job, _file

        download_data = {
            'page': 'DownloadFile',
            'job': _job,
            'file': _file,
            'do_download': 'Download'}

        r = requests.post(data=download_data, url='http://rast.theseed.org/FIG/rast.cgi', cookies=cookies)

        with open(os.path.join(output_dir, _file), 'wb') as f:
            f.write(r.content)
        shutil.move(os.path.join(output_dir, _file), os.path.join(output_dir, _job + '_' + _file))
