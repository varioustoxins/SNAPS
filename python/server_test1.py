# -*- coding: utf-8 -*-
"""
Created on Sat Oct  6 11:01:43 2018

@author: Alex
"""

from http.server import *

def run(server_class=HTTPServer, handler_class=BaseHTTPRequestHandler):
    server_address = ('', 8000)
    httpd = server_class(server_address, handler_class)
    httpd.serve_forever()