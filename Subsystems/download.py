#!/usr/bin/env python

import os
import sys
import requests
import shutil


def main(args):
    uname, password, input_file, output_dir = args

    login_data = {"page": 'Login', "action": 'perform_login',
                  "login": uname,
                  "password": password}

    url = 'http://rast.theseed.org/FIG/rast.cgi'

    cookies = requests.post(url, data=login_data).cookies
    cookies = {'WebSession': cookies.get('WebSession')}

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    with open(input_file) as f:
        for ln in f:
            _job, _file = ln.strip().split('\t')
            _file = _file + '.gbk'
            print _job, _file

            download_data = {
                'page': 'DownloadFile',
                'job': _job,
                'file': _file,
                'do_download': 'Download'}

            r = requests.post(data=download_data, url=url, cookies=cookies)

            with open(os.path.join(output_dir, _file), 'wb') as f2:
                f2.write(r.content)
            shutil.move(os.path.join(output_dir, _file), os.path.join(output_dir, _job + '_' + _file))


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
