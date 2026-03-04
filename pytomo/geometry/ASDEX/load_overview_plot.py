#!/usr/bin/env python 
# -*- coding: utf-8 -*-






def ShotOverview(shotnum):


    url = 'https://www.aug.ipp.mpg.de/cgibin/local_or_sfread/cview.cgi?shot=%d'
    u = 'dG9kc3RyY2k=\n'
    p = 'bWVnYWhlc2xv\n'

    import urllib.request, urllib.error, urllib.parse,base64


    password_mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    password_mgr.add_password(None, "https://www.aug.ipp.mpg.de/cgibin/local_or_sfread/", u.decode('base64'), p.decode('base64'))

    handler = urllib.request.HTTPBasicAuthHandler(password_mgr)
    opener = urllib.request.build_opener(handler)

        
    imgData = opener.open(url%shotnum).read()

    output = open('./overview_plots/overview_%d.gif','wb')
    output.write(imgData)
    output.close()